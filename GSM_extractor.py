"""
AUTHOR: Jean-François Nadeau
SCRIPT: GSM_extractor.py
    This script is used to extract pertinent info from a GSM native XML file from the sra database.
    It also search the targeted protein for some type of experiments like ChIP-Seq.
USAGE:
    python3 GSM_extractor.py processes path
    -processes: Integer indicating the number of paralell process to create. Default is 4.
    -path: Path to the directory containing the xml files. If you need to download the XML files, the path should be the word "online"
Examples:
    python3 GSM_extractor 12 xml/
    python3 GSM_extractor 8 online
NOTE: To use the online mode, please fill the config.ini files that contains the information for the Entrez query and the output directory for the xml files.
"""

import os
import sys
from multiprocessing import Pool
import read_xml
import re
from collections import OrderedDict
from functools import partial
from Bio import Entrez




""" ###############################################################################
    #                               MAIN FUNCTION                                 #
    ###############################################################################
"""
def main():
    # Read arguments
    regex_dictio=get_dict("saccharomyces_cerevisiae")
    args=read_config()
    try:
        cores=sys.argv[1]
        path=sys.argv[2]
    except:
        print_help()
        quit()
    ### ONLINE MODE ###
    if path=="online":
        clear_outdir("xml_directory")
        Entrez.email = args["Entrez_email"]
        xml_out=args["xml_out"]
        for org in args["Organisms"].split(","):
            date=args['Date_range'].split(":")
            for year in range(int(date[0].split("/")[0]),int(date[1].split("/")[0])+1):
                query_list=build_query(args, org, str(year))
                online_mode(regex_dictio,query_list,cores,xml_out) # THIS FUNCTION WILL GET THE HANDLE, AND RUN THE SCRIPT FOR EACH ID IN THE HANDLE IN PARALLEL
    ### LOCAL MODE ###
    else:
        path=path.rstrip("/")
        for root, dirs, files in os.walk(path):
            files=["{}/{}".format(path,filename) for filename in files]
            pool=Pool(processes=int(cores))
            func = partial(get_line,True) 
            lines=pool.map(func,files)
            func = partial(process_line,regex_dictio)
            lines=pool.map(func, lines)




""" ###############################################################################
    #                     METADATA PROCESSING MODULES                             #
    ###############################################################################
"""
#This function checks and correct the library strategy field if needed
def check_assay(STRATEGY,attributes,title,study):
    info=" - ".join([attributes,title,study])
    names={"BREAK-SEQ":"BREAK-?SEQ","BRDU":"BRDU","ATAC-SEQ":"ATAC-?SEQ","CUT-AND-RUN":"CUT.?AND.?RUN","FAIRE-SEQ":"FAIRE-?SEQ","MNASE-SEQ":"MNASE","CHIP-EXO":"CHIP-?EXO","CHEC-SEQ":"CHEC-*SEQ","CHIP-ESPAN":"chip-*espan"}
    for name in names:
        if re.search(names[name],info.upper())!=None:
            STRATEGY=name
    return STRATEGY


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
    if re.search(r"(input|_in_|wce|whole\scell\sextract|mock|control|untagged)",title.lower())!=None:
        target="INPUT"
        confidence=3
    return target,confidence


# This function recieve the lists generated by the confidence_level function ordered in a confidence level order and return the first non-empty list od size 1. 
# If multiple targets, returns all targets with a confidence of 5.
# If nothing found, return not_found with a confidence of 0
def get_best_match(matches_lists):
    for i in range(3):
        target=list(set(matches_lists[i]))
        if len(target)==1:
            return target[0],1
    if len(list(set(matches_lists[3])))==1:
        return matches_lists[3][0],2
    elif len(list(set(matches_lists[4])))==1:
        return matches_lists[4][0],3
    elif len(list(set(matches_lists[5])))==1:
        return matches_lists[5][0],4
    elif(len(list(set(matches_lists[4])))>1):
        return " & ".join(matches_lists[4]),5
    else:
        return "not_found",0


#This function search for flagged terms in a given field and return all hits
def get_flagged(regex_dictio,field,flags):
    targets=[]
    for f in flags:
        if re.search(flags[f],str(field).lower()):
            flagged=re.findall(flags[f],str(field).lower())
            if flagged:
                for t in regex_dictio:
                    if re.search(regex_dictio[t], str(flagged)):
                        targets.append(t)
    return targets


# extract pertinent info from a given xml file
def get_line(local,file):
    return read_xml.fields(file,local)

def get_strain(attributes,title):
    strain=""
    for att in attributes.lower().split(" | "):
        if any(x in att for x in ["strain:","strain name:","strain number:","cell description:"]):
            strain=att.split(":")[1]
                
    return strain

#This function search for targets in a given field
def search_all_targets(regex_dictio,field):
    hits=[]
    for t in regex_dictio:
        if re.search(regex_dictio[t], field):
            hits.append(t)
    return hits


#This function search for a given target in a given field
def search_this_target(regex_dictio,field,t):
    hits=[]
    if field:
        if re.search(regex_dictio[t], str(field)):
            hits.append(t)
    return hits


#This function search for the target and fill different lists according to the confidence level of the research algorithms used"
def search_target(regex_dictio,attributes,title,flags):
    lvl1_1,lvl1_2,lvl1_3,lvl2,lvl3,lvl4=[],[],[],[],[],[]
    # CONFIDENCE_1
    for att in attributes.lower().split(" | "):
        if any(x in att for x in ["antibody:","chip:","protein:","flag tagged:","target of ip:","epitope:","1st ip:","antibody #lot number:", "epitope tags:"]):
            flagged1=flagged2=[]
            for f in flags:
                if re.search(flags[f],att):
                    flagged1+=re.findall(flags[f],attributes.lower())
                if re.search(flags[f],title.lower()):
                    flagged2+=re.findall(flags[f],title.lower())
            for t in regex_dictio:
                lvl1_1+=search_this_target(regex_dictio,att,t)
                lvl1_2+=search_this_target(regex_dictio,flagged1,t)
                lvl1_3+=search_this_target(regex_dictio,flagged2,t)
        elif "source name:" in att or "antibody #lot number:" in att:
            lvl2=rm_deltas(regex_dictio,search_all_targets(regex_dictio,att),att)
    # CONFIDENCE_2_TO_4
    if not lvl1_1+lvl1_2+lvl1_3:
        if len(lvl3)==1:
            return lvl3[0],2
        else:
            for t in regex_dictio:
                regex=regex_dictio[t].replace(r"(\D+|$)","")+r"(?=(_|\s|-)?(ip|chip|crosslinked|protein\schip|sonication_ip|chromatin\simmuno))"
                if re.search(regex,title.lower()):
                    lvl2.append(t)
                if re.search(regex_dictio[t],title.lower()):
                    lvl3.append(t)
            lvl2=rm_deltas(regex_dictio,lvl2,title)
            lvl3=rm_deltas(regex_dictio,lvl3,title)
            lvl4=get_flagged(regex_dictio,title,flags)
    return get_best_match([lvl1_1,lvl1_2,lvl1_3,lvl2,lvl3,lvl4])


# checks if the targets in a given list match a deletion patern. Returns a list of the targets that do not match in any patterns.
def rm_deltas(regex_dictio,target_list,field):
    new=[]
    for t in target_list:
        regex=r'((∆|d|d\(|delta|del)+(-|_)?'+regex_dictio[t]+'|'+regex_dictio[t].replace(r"(\D+|$)","")+r'-?(δ|\stail\sdelete|del|\{delta\}|∆|aa|d|-(\s|_|$)))'
        if re.search(regex,field.lower())==None:
            new.append(t)
    return new


# Process line according to the library strategy
def process_line(regex_dictio,line):
    strain=""
    GSM=list(set(re.findall("GSM[0-9]+",str(line))))
    if GSM:
        if len(GSM)==1:
            flags={"FLAG":"([^_\\s]+-[0-9]*x?flag|[0-9]*x?flag-[^_\\s]+|flag)","MYC":"([^_\\s]+(-c)?-[0-9]*x?myc|(c-)?[0-9]*x?myc-[^_\\s]+|(c-)?myc|9e10)","V5":"([^_\\s]+-[0-9]*x?v5|[0-9]*x?v5-[^_\\s]+|v5)","TAP":"([^_\\s]+-[0-9]*x?tap|[0-9]*x?tap-[^_\\s]+|tap)","HA":"([^_\\s]+-[0-9]*x?ha|[0-9]*x?ha-[^_\\s]+|ha)","GFP":"([^_\\s]+-[0-9]*x?gfp|[0-9]*x?gfp-[^_\\s]+|gfp)","T7":"([^_\\s]+-[0-9]*x?t7|[0-9]*x?t7-[^_\\s]+|t7)"}
            joined=" - ".join([line["SAMPLE_NAME"],line["SAMPLE_TITLE"],line["EXP_TITLE"],line["ATTRIBUTES"],line["STUDY_TITLE"],line["LIB_STRAT"]])
            if any(x in joined.lower() for x in ["chip-seq","chip-exo"]) or line["LIB_STRAT"].lower()=="other":
                strat=check_assay(line["LIB_STRAT"],line["ATTRIBUTES"],line["SAMPLE_TITLE"],line['STUDY_TITLE'])
                if any(x in strat.lower() for x in ["chip-seq","chip-exo","cut-and-run"]):
                    target,confidence=check_if_input(line["SAMPLE_TITLE"],line["ATTRIBUTES"])
                    if confidence==0:
                        target,confidence=search_target(regex_dictio,line["ATTRIBUTES"],line["SAMPLE_TITLE"],flags)
                        strain=get_strain(line["ATTRIBUTES"],line["SAMPLE_TITLE"])
                else:
                    line["LIB_STRAT"]=strat
                    confidence=0
                    target=strat
            else:
                target=line["LIB_STRAT"]
                confidence=0
            line["STRAIN"]=strain
            print("\t".join([GSM[0],str(confidence),target]))
            #print(line["STRAIN"])




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
    f=org+".dict"
    with open(f,"r") as f:
        lines=f.readlines()
    for line in lines:
        if line!="\n" and line!="":
            line=line.rstrip("\n").split(",")
            org_dictio[line[0]]=line[1]
    with open("common.dict","r") as f:
        lines=f.readlines()
    for line in lines:
        if line!="\n" and line!="":
            line=line.rstrip("\n").split(",")
            common_dictio[line[0]]=line[1]
    targets=OrderedDict(org_dictio)
    targets.update(OrderedDict(common_dictio))
    return targets


# Read sets the nomber of cores to use
def parallel_processes(arguments):
    if len(arguments)>1:
        cores=int(arguments[1])
    else:
        cores=4
    return cores


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





""" ###############################################################################
    #                           ONLINE MODE MODULES                               #
    ###############################################################################
"""
#Build and return a query for each month of the given year
def build_query(args, org, year):
    query_list=[]
    months=["01","02","03","04","05","06","07","08","09","10","11","12"]
    #Split the month range if months are specified in the config file
    if "/" in args['Date_range']:
        if year in args['Date_range'].split(":")[0]:
            months=months[months.index(args['Date_range'].split(":")[0].split("/")[1]):]
        if year in args['Date_range'].split(":")[1]:
            months=months[:months.index(args['Date_range'].split(":")[1].split("/")[1])+1]
    for month in months:
        #adding the org
        query=org+"[ORGN] "
        #adding the search terms:
        if "Search_terms" in args:
            query+="AND {}".format(" AND ".join(args['Search_terms'].split(",")))
        #adding the filter out terms:
        if "Filter_out" in args:
            query+=" NOT ( {} )".format(" OR ".join(args["Filter_out"].split(",")))
        #adding the date
        query+=" AND "+year+"/"+month+"[PDAT]"
        if "Custom" in args:
            query+=" "+args['Custom']
        query_list.append(query)
    return query_list


# Dowload and process the given xml file
def efetch(regex_dictio,xml_out,id):
    xml_out=xml_out.rstrip("/")
    sample=Entrez.efetch(db="sra", id=id, format="native").readlines()
    GSM=list(set(re.findall("GSM[0-9]+","".join(sample))))
    if len(GSM)==1:
        filename="{}/{}.xml".format(xml_out,GSM[0])
        if os.path.isfile(filename):
            number=2
            temp=filename
            while os.path.isfile(temp):
                temp="{}_{}.xml".format(filename.rstrip(".xml"),number)
                number+=1
            filename=temp
        with open(filename,"w") as f:
            f.write("".join(sample))
        line=get_line(False,sample)
        process_line(regex_dictio,line)

# Send esearch query and return the resulting dictionary (containing the count, id of the samples etc...)
def get_sra_handle(query, database):
    handle = Entrez.esearch(db=database,retmax=1000000, term=query)
    dic=Entrez.read(handle)
    return dic

# Process queries generated by build_query
def online_mode(regex_dictio,query_list,cores,xml_out):
    #Ajusting the month range if the user specified it
    for query in query_list:
        esearch_dic=get_sra_handle(query, "sra")
        pool=Pool(processes=int(cores))
        if str(esearch_dic["Count"])!="0":
            if esearch_dic['IdList']:
                    func=partial(efetch,regex_dictio,xml_out)
                    pool.map(func,esearch_dic['IdList'])





if __name__ == "__main__":
    main()