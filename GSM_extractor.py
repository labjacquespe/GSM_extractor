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
import get_target
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
    regex_dictio=get_common_dict()
    #0:corrected_GSM, 1:GSE, 2:GPL, 3:organism, 4:library_strategy, 5:sample_title, 6:attributes, 7:sample_source, 8:molecule, 9:platform_technology, 10:library_source, 11:library_selection, 12:series_title, 13:series_sumary, 14:series_design, 15:contributor, 16:organization, 17:release_date, 18:submission_date
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
                target,confidence=get_target.custom_case1(regex_dictio,sample_title)
            else:
                target,confidence=check_if_input(sample_title,attributes)
                if confidence==0:
                    target,confidence=get_target.search_target(regex_dictio,remove_term(attributes),remove_term(sample_title),source,flags)

        elif strategy.lower() in ["mnase-seq", "other","brdu-seq"]:
            target,confidence=get_target.confidence1_only(regex_dictio,remove_term(attributes),remove_term(sample_title),source,flags)
            if confidence==0:
                target=strategy.lower()
        else:
            target,confidence=strategy,0
    else:
        target,confidence=strategy,0
    return "\t".join([line[0],strategy,organism,str(confidence),target])
    
def remove_term(field):
    terms=["red1 chip","jhd2_chip_seq_yar","wt_chip_seq_yar","yng2-wa"]
    for term in terms:
        field=field.lower().replace(term,"")
    return field

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




""" ###############################################################################
    #                              CONFIG MODULES                                 #
    ###############################################################################
"""


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
