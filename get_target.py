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

#This function search for targets in a given field
def search_all_targets(regex_dictio,field):
    hits=[]
    if field:
        for t in regex_dictio:
            if re.search(regex_dictio[t], str(field).lower()):
                hits.append(t)
    return list(set(hits))


# checks if the targets in a given list match a deletion patern. Returns a list of the targets that do not match in any patterns.
def rm_deltas(regex_dictio,target_list,field):
    new=[]
    for t in target_list:
        regex=r'((∆|d|d\(|delta|del)+(-|_)?'+regex_dictio[t]+'|'+regex_dictio[t].replace(r"(\D+|$)","")+r'-?(-as|δ|\stail\sdelete|del|\{delta\}|∆|aa|d|-(\s|_|$)))'
        if not re.search(regex,field.lower()):
            new.append(t)
    return new

def search_target(regex_dictio,attributes,title,source,flags,series_title):
    attributes=" | ".join([attributes,source])
    c1=confidence1(regex_dictio,attributes,title,flags)
    if len(c1)==1:
        return c1[0],1
    c2=confidence2(regex_dictio,attributes,title)
    if len(c2)==1:
        return c2[0],2
    c3,c4=confidence_3_4(regex_dictio,attributes,title)
    if len(c3)==1:
        return c3[0],3
    if len(c4)==1:
        return c4[0],4
    c5=list(set(search_all_targets(regex_dictio,series_title)))
    if len(c5)==1:
        return c5[0],5
    c6=c1+c2+c3+c4+c5
    if c6:
        return " & ".join(list(set(c6))),6
    else:
        return "not_found",7

def confidence1(regex_dictio,attributes,title,flags):
    #This algorithm returns a target if it directly found in a high confidence field (HCF) or elsewhere associated with a flag found in a HCF
    lvl1=[]
    flags_found=[]
    for att in attributes.lower().split(" | "):
        if any(x in att for x in ["antibody:","chip:"," ","flag tagged:","target of ip:","epitope:","antibody #lot number:", "epitope tags:"]):
            lvl1+=search_all_targets(regex_dictio,att)
            for f in flags.keys():
                if re.search(flags[f],att):
                    flags_found.append(f)
    
    if len(lvl1)==1:
        return lvl1
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
    return list(set(lvl1))

def confidence2(regex_dictio,attributes,title):
    lvl2=[]
    #this algorithm returns a target if a single hit is found in a medium condidence field
    for att in attributes.lower().split(" | "):
            if any(x in att for x in ["source name:","source_name:","antibody #lot number:"]):
                lvl2=search_all_targets(regex_dictio,att)
                if len(lvl2)==1:
                    print(title,attributes,lvl2)
    return list(set(lvl2))

def confidence_3_4(regex_dictio,attributes,title):
    lvl3=[]
    lvl4=[]
    for t in regex_dictio:
        regex=regex_dictio[t].replace(r"(\D+|$)","")+r"(?=(_|\s|-)?(ip|chip|crosslinked|protein\schip|sonication_ip|chromatin\simmuno))"
        if re.search(regex,title.lower()):
            lvl3.append(t)
        if re.search(regex_dictio[t],title.lower()):
            lvl4.append(t)
    lvl3=rm_strain(regex_dictio,attributes,rm_deltas(regex_dictio,lvl3,title.lower()))
    lvl4=rm_strain(regex_dictio,attributes,rm_deltas(regex_dictio,lvl4,title.lower()))
    return list(set(lvl3)),list(set(lvl4))
    

def rm_strain(regex_dictio,attributes,l):
    if len(l)>1 and "strain:" in attributes:
        for att in attributes.lower().split(" | "):
            if "strain:" in att:
                for t in l:
                    if re.search(regex_dictio[t],att):
                        l.remove(t)
        return l
    return l