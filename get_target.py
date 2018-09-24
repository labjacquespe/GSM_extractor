import re


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
    lvl1,lvl3,lvl4,lvl5=[],[],[],[]
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

