"""
AUTHOR: Jean-François Nadeau
SCRIPT: GSM_extractor.py
    This script is used to extract pertinent info from a GSM native XML file from the sra database.
    It also search the targeted protein for some type of experiments like ChIP-Seq.
USAGE:
    python3 GSM_extractor.py processes path
    -processes: Integer indicating the number of paralell process to create.
    -path: Path to the directory containing the xml files. If you need to download the XML files, the path should be the word "online"
Examples:
    python3 GSM_extractor 12 xml/
    python3 GSM_extractor 8 online
NOTE: To use the online mode, please fill the config.ini files that contains the information for the Entrez query and the output directory for the xml files.
"""

import os
import sys
from multiprocessing import Pool
import extract_xml
import download_files
import re
from collections import OrderedDict
from functools import partial
from Bio import Entrez
import urllib.request
import tarfile
import xml.etree.ElementTree as ET



""" ###############################################################################
    #                               MAIN FUNCTION                                 #
    ###############################################################################
"""

def main():
    # Read arguments
    lines=[]
    to_process=[]
    try:
        args=read_config()
        cores=int(sys.argv[1])
        path=sys.argv[2]
    except:
        print_help()
    if path=="online":
        download_files.download(args,cores)

    if os.path.isdir(path):
        print("INFO: Extracting samples from the GSE xml files...", file=sys.stderr)
        for root, dirs, files in os.walk(path):
            files=["{}/{}".format(path.rstrip("/"),filename) for filename in files]
            for f in files:
                lines+=extract_xml.get_metadata(f)
            #tsv_pool=Pool(processes=cores)
            #list_of_lists_of_lines=tsv_pool.map(extract_xml.get_metadata,files)
            #lines = list(set([item for sublist in list_of_lists_of_lines for item in sublist]))
            #tsv_pool.close()
        print("INFO: --> Samples extracted succesfully.", file=sys.stderr)
        if "Create_tsv" in args and args["Create_tsv"].lower()=="true":
            print("INFO: Creating tsv file...", file=sys.stderr)
            with open("GSM_extractor_out.tsv","w") as outf:
                outf.write("\n".join(lines))
            print("INFO: --> Tsv file created.", file=sys.stderr)
    elif os.path.isfile(path):
        print("INFO: Reading lines from tsv file", file=sys.stderr)
        with open(path, "r") as inf:
            lines=inf.readlines()
            lines=[x.rstrip("\n") for x in lines]
    
    
    print("INFO: Processing lines to find targeted proteins...", file=sys.stderr)
    process_pool=Pool(processes=cores)
    results=process_pool.map(process_line,lines)
    print(*list(set(results)), sep='\n')
    print ('@INFO: Done!', file=sys.stderr)


###############################################################################
    #                     METADATA PROCESSING MODULES                             #
###############################################################################

def process_line(line):
    
    line=line.rstrip("\n").split("\t")
    #regex_dictio=get_dict(line[3].replace(" ","_").lower())
    regex_dictio=get_common_dict()
    #[corrected_GSM,GSE,GPL,organism,library_strategy,sample_title," | ".join(attributes),sample_source,molecule,platform_technology,library_source,library_selection,series_title,series_sumary,series_design,contributor,organization,release_date,submission_date]
    organism=line[3].lower()
    strategy=line[4].lower()
    sample_title=line[5].lower()
    attributes=line[6].lower()
    source=line[7].lower()
    flags={"FLAG":"(-|::)?[0-9]*x?-?flag","MYC":"(-|::)?[0-9]*x?-?myc|9e10","V5":"(-|::)?[0-9]*x?-?v5","TAP":"(-|::)?[0-9]*x?-?tap","HA":"(-|::)?[0-9]*x?-?ha","GFP":"(-|::)?[0-9]*x?-?e?gfp","T7":"(-|::)?[0-9]*x?-?t7"}
    #flags={"FLAG":r"([^_\s]+-[0-9]*x?flag|[0-9]*x?flag-[^_\s]+|flag)","MYC":r"([^-_\s]+(-c)?-[0-9]*x?myc|(c-)?[0-9]*x?myc-[^_\s]+|(c-)?myc|9e10)","V5":r"([^_\s]+-[0-9]*x?v5|[0-9]*x?v5-[^_\s]+|v5)","TAP":r"([^_\s]+-[0-9]*x?tap|[0-9]*x?tap-[^_\s]+|tap)","HA":r"([^_\s]+-[0-9]*x?ha|[0-9]*x?ha-[^_\s]+|ha)","GFP":r"([^_\s]+-[0-9]*x?gfp|[0-9]*x?gfp-[^_\s]+|gfp)","T7":r"([^_\s]+-[0-9]*x?t7|[0-9]*x?t7-[^_\s]+|t7)"}
    joined=strategy+sample_title+attributes
    if any (x in joined for x in ["chip-seq","chip-exo","mnase-seq","chec-seq","cut-and-run"]) or strategy=="other":
        strategy=check_assay(joined, strategy)
        if any(x in strategy.lower() for x in ["chip-seq","chip-exo","chec-seq","cut-and-run"]):
            if "1st ip:" in attributes:
                target,confidence=custom_case1(regex_dictio,sample_title)
            else:
                target,confidence=check_if_input(sample_title,attributes)
                if confidence==0:
                    target,confidence=search_target(regex_dictio,remove_term(attributes),remove_term(sample_title),source,flags)

        elif strategy.lower() in ["mnase-seq", "other","brdu-seq"]:
            target,confidence=confidence1_only(regex_dictio,remove_term(attributes),remove_term(sample_title),source,flags)
            if confidence==0:
                target=strategy.lower()
        else:
            target,confidence=strategy,0
    else:
        target,confidence=strategy,0
    #print("\t".join([line[0],organism,strategy,target,str(confidence)]))
    return "\t".join([line[0],strategy,organism,str(confidence),target])
    

def custom_case1(regex_dictio,sample_title):
    #custom case for samples with "1st ip" and "2nd ip" attributes
    hits=[]
    info=sample_title.split("-")
    pair=info[1:]
    for n in pair:
        for t in regex_dictio:
            if re.search(regex_dictio[t],n.lower()):
                hits.append(t)
    if "input" in str(pair).lower():
        hits.append("input")
    if len(hits)==1:
        hits.append(hits[0])
    return "/".join(hits),1


#This function checks and correct the library strategy field if needed
def check_assay(info,strategy):
    names={"BREAK-SEQ":"BREAK-?SEQ","5HMC":"5HMC-?SEQ","BRDU-SEQ":"BRDU","ATAC-SEQ":"ATAC-?SEQ","CUT-AND-RUN":"CUT.?AND.?RUN","GRO-SEQ":"GRO-?SEQ","GRO-CAP":"GRO-?CAP","FAIRE-SEQ":"FAIRE","CHIP-EXO":"CHIP-?EXO","CHEC-SEQ":"CHEC-?SEQ","CHIP-ESPAN":"CHIP-?ESPAN"}
    for name in names:
        if re.search(names[name],info.upper())!=None:
            strategy=name
    return strategy.lower()


#This function checks it a control term is found within the sample title and attributes
def check_if_input(title,attributes):
    target="not_found"
    confidence=0
    attributes=re.sub(r'\([^)]*\)', '', attributes.lower().replace("tissue:wce",""))
    for att in attributes.split(" | "):
        if any(x in att for x in ["antibody:","chip:","protein:","flag tagged:","target of ip:"]):
            if any(x in att for x in ["input","mock","wce","control","igg","_in_","wce","whole cell extract", "antibody:none","epitope:none"]):
                target="INPUT"
                confidence=1
        elif "source:" in att:
            if any(x in att for x in ["input","mock","wce","control","igg","_in_","wce","whole cell extract"]):
                target="INPUT"
                confidence=5
    if re.search(r"(input|_in_|wce|whole\scell\sextract|mock|control|untagged|igg-ip)",title.lower())!=None:
        target="INPUT"
        confidence=3
    return target,confidence


# This function recieve the lists generated by the confidence_level function ordered in a confidence level order and return the first non-empty list od size 1. 
# If multiple targets, returns all targets with a confidence of 5.
# If nothing found, return not_found with a confidence of 0


#This function search for flagged terms in a given field and return all hits
def get_flagged(regex_dictio,field,flags):
    targets=[]
    for f in flags:
        if re.search(flags[f],str(field).lower()):
            regex="[a-z0-9-]+(?={})".format(flags[f])
            flagged=re.findall(regex,str(field).lower())
            if flagged:
                for t in regex_dictio:
                    if re.search(regex_dictio[t], str(flagged)):
                        targets.append(t)
    return targets


def get_strain(attributes,title):
    strain=""
    for att in attributes.lower().split(" | "):
        if any(x in att for x in ["strain:","strain name:","strain number:","cell description:"]):
            strain=att.split(":")[1]
                
    return strain

#This function search for targets in a given field
def search_all_targets(regex_dictio,field):
    hits=[]
    if field:
        for t in regex_dictio:
            if re.search(regex_dictio[t], str(field).lower()):
                hits.append(t)
    return list(set(hits))


#This function search for a given target in a given field
def search_this_target(regex_dictio,field,t):
    hits=[]
    if field:
        if re.search(regex_dictio[t], str(field).lower())!=None:
            hits.append(t)
    return hits

def return_best_list(list_of_lists):
    for l in list_of_lists:
        if len(list(set(l)))==1:
            return list(set(l))
    return []

def confidence1_only(regex_dictio,attributes,title,source,flags):
    attributes=" | ".join([attributes,source])
    lvl1=[]
    flags_found=[]
    for att in attributes.lower().split(" | "):
        if any(x in att for x in ["antibody:","chip:"," ","flag tagged:","target of ip:","epitope:","antibody #lot number:", "epitope tags:"]):
            lvl1+=search_all_targets(regex_dictio,att)
            for f in flags.keys():
                if re.search(flags[f],att):
                    flags_found.append(f)
    
    if len(lvl1)==1:
        return lvl1[0],1
    elif flags_found:
        hits=[]
        for f in flags_found:
            regex="[a-z0-9-]+(?={})".format(flags[f])
            if re.search(regex,attributes):
                hits=re.findall(regex,attributes)
            elif re.search(regex,title):
                hits=re.findall(regex,title)
        if hits:
            lvl1=search_all_targets(regex_dictio," ".join(hits))
            
    if len(lvl1)==1:
        return lvl1[0],1
    else:
        return "not_found",0

#This function search for the target and fill different lists according to the confidence level of the research algorithms used"
def search_target(regex_dictio,attributes,title,source,flags):
    #if "Cut-and-Run_H2A" in title:
    attributes=" | ".join([attributes,source])
    lvl1,lvl1_1,lvl1_2,lvl1_3,lvl2,lvl3,lvl4,lvl5=[],[],[],[],[],[],[],[]
    # CONFIDENCE_1
    flags_found=[]
    for att in attributes.lower().split(" | "):
        if any(x in att for x in ["antibody:","chip:","protein:","flag tagged:","target of ip:","epitope:","antibody #lot number:", "epitope tags:"]):
            lvl1+=search_all_targets(regex_dictio,att)
            for f in flags.keys():
                if re.search(flags[f],att):
                    flags_found.append(f)
    if len(lvl1)==1:
        return lvl1[0],1
    elif flags_found:
        hits=[]
        for f in flags_found:
            regex="[a-z0-9-]+(?={})".format(flags[f])
            if re.search(regex,attributes):
                hits=re.findall(regex,attributes)
            elif re.search(regex,title):
                hits=re.findall(regex,title)
        if hits:
            lvl1=search_all_targets(regex_dictio," ".join(hits))
    if len(lvl1)==1:
        return lvl1[0],1
    else:
        for att in attributes.lower().split(" | "):
            if any(x in att for x in ["source name:","source_name:","antibody #lot number:"]):
                source_hits=search_all_targets(regex_dictio,att)
                if len(source_hits)==1:
                    print(title,attributes,source_hits)
                    return source_hits[0],2
        # CONFIDENCE_2_TO_4
        for t in regex_dictio:
            regex=regex_dictio[t].replace(r"(\D+|$)","")+r"(?=(_|\s|-)?(ip|chip|crosslinked|protein\schip|sonication_ip|chromatin\simmuno))"
            if re.search(regex,title.lower()):
                lvl3.append(t)
            if re.search(regex_dictio[t],title.lower()):
                lvl4.append(t)
        lvl3=rm_deltas(regex_dictio,lvl3,title.lower())
        lvl4=rm_deltas(regex_dictio,lvl4,title.lower())
        lvl5=get_flagged(regex_dictio,title,flags)
        lower_conf=[lvl3,lvl4,lvl5]
        for i in range(len(lower_conf)):
            if len(list(set(lower_conf[i])))==1:
                return lower_conf[i][0],i+2
    if len(lvl4)>1 and "strain:" in attributes:
        for att in attributes.lower().split(" | "):
            if "strain:" in att:
                for t in lvl4:
                    if re.search(regex_dictio[t],att):
                        lvl4.remove(t)
        if len(lvl4)==1:
            return lvl4[0],4
    lvl6=" & ".join(lvl4)
    if lvl6:
        return lvl6,6
    else:
        return "not_found",0


# checks if the targets in a given list match a deletion patern. Returns a list of the targets that do not match in any patterns.
def rm_deltas(regex_dictio,target_list,field):
    new=[]
    for t in target_list:
        regex=r'((∆|d|d\(|delta|del)+(-|_)?'+regex_dictio[t]+'|'+regex_dictio[t].replace(r"(\D+|$)","")+r'-?(-as|δ|\stail\sdelete|del|\{delta\}|∆|aa|d|-(\s|_|$)))'
        if not re.search(regex,field.lower()):
            new.append(t)
    return new

def remove_term(field):
    terms=["red1 chip","jhd2_chip_seq_yar","wt_chip_seq_yar","yng2-wa"]
    for term in terms:
        field=field.lower().replace(term,"")
    return field


""" ###############################################################################
    #                              CONFIG MODULES                                 #
    ###############################################################################
"""
# Clear the xml output directory
def clear_outdir(dirPath):
    if not os.path.exists(dirPath):
        os.makedirs(dirPath)
    fileList = os.listdir(dirPath)
    [ os.remove(os.path.abspath(os.path.join(dirPath,fileName))) for fileName in fileList ]


#Build and return the regex_dict
def get_dict(org):
    org_dictio={}
    common_dictio={}
    lines=[]
    f="dictios/{}.dict".format(org)
    with open(f,"r") as f:
        lines=f.readlines()
    for line in lines:
        if line!="\n" and line!="":
            line=line.rstrip("\n").split(",")
            org_dictio[line[0]]=line[1]
    with open("dictios/common.dict","r") as f:
        lines=f.readlines()
    for line in lines:
        if line!="\n" and line!="":
            line=line.rstrip("\n").split(",")
            common_dictio[line[0]]=line[1]
    targets=OrderedDict(org_dictio)
    targets.update(OrderedDict(common_dictio))
    return targets

def get_common_dict():
    common_dictio={}
    with open("dictios/common.dict","r") as f:
        lines=f.readlines()
    for line in lines:
        if line!="\n" and line!="":
            line=line.rstrip("\n").split(",")
            common_dictio[line[0]]=line[1]
    return OrderedDict(common_dictio)

def print_help():
    text="""SCRIPT: GSM_extractor.py
    This script is used to extract pertinent info from a GSM native XML file from the sra database.
    It also search the targeted protein for some type of experiments like ChIP-Seq.

    USAGE:
        python3 GSM_extractor.py processes path
        -processes: Integer indicating the number of paralell process to create. Default is 4.
        -path: Path to the directory containing the xml files. If you need to download the XML files, the path should be the word "online"
    Examples:
        python3 GSM_extractor 12 xml/
        python3 GSM_extractor 8 online

    NOTE: To use the online mode, please fill the config.ini files that contains the information for the Entrez query and the output directory for the xml files."""
    print(text)
    quit()


# Read and return the config file
def read_config():
    params=""
    args={}
    with open("config.ini", "r") as conf:
        params=conf.readlines()
    for line in params:
        if not line.startswith("#"):
            if "=" in line:
                line=line.rstrip("\n").split("=")
                args[line[0]]=line[1]
            else:
                print("ERROR: invalid config file")
                quit()
    return args






if __name__ == "__main__":
    main()
