import xml.etree.ElementTree as ET

def get_samn(exp,xml):
    SAMN='not_available'
    for path in ['EXPERIMENT/DESIGN/SAMPLE_DESCRIPTOR/IDENTIFIERS/EXTERNAL_ID','SAMPLE/IDENTIFIERS/EXTERNAL_ID','Pool/Member/IDENTIFIERS/EXTERNAL_ID','RUN_SET/RUN/Pool/Member/IDENTIFIERS/EXTERNAL_ID']:
        if exp.find(path)!=None and SAMN=='not_available' and 'namespace' in exp.find(path).attrib:
            if exp.find(path).attrib['namespace']=="BioSample":
                SAMN=exp.find(path).text
    return SAMN

def get_srp(exp,xml):
    SRP='not_available'
    for path in ['EXPERIMENT/STUDY_REF', 'STUDY']:
        if exp.find(path)!=None and SRP=='not_available' and 'accession' in exp.find(path).attrib:
            tmp=exp.find(path).attrib['accession']
            if len(tmp)>3 and tmp[1:3]=="RP":
                SRP=tmp
    return SRP

def get_srx(exp,xml):
    SRX='not_available'
    for path in ['EXPERIMENT', 'RUN_SET/RUN/EXPERIMENT_REF']:
        if exp.find(path)!=None and SRX=='not_available' and 'accession' in exp.find(path).attrib:
            tmp=exp.find(path).attrib['accession']
            if len(tmp)>3 and tmp[1:3]=="RX":
                SRX=tmp
    return SRX



def fields(xml,local):
    if local:
        tree = ET.parse(xml)
    else:
        tree = ET.fromstring("".join(xml))
    SRA=SAMN=SRS=SRX=SRP=LIB_NAME=LIB_STRAT=LIB_SOURCE=LIB_SELECT=SRR=PLATFORM=INSTRUMENT_MODEL=LAB_NAME=CENTER_NAME=SUB_ID=ABSTRACT=PDAT=ORGANISM=SAMPLE_NAME=SAMPLE_TITLE=EXP_TITLE=DESCRIPTION="not_available"

    for exp in tree.iter('EXPERIMENT_PACKAGE'):
        #IDENTIFIERS   SRA, SAMN, SRS, SRX, SRP
        if exp.find('SUBMISSION')!=None and 'accession' in exp.find('SUBMISSION').attrib:
            SRA=exp.find('SUBMISSION').attrib['accession']
        SAMN=get_samn(exp,xml)
        if exp.find('SAMPLE')!=None and 'accession' in exp.find('SAMPLE').attrib:
            SRS=exp.find('SAMPLE').attrib['accession']
        SRX=get_srx(exp,xml)
        SRP=get_srp(exp,SAMN)
        if exp.find('RUN_SET/RUN')!=None:
            if 'accession' in exp.find('RUN_SET/RUN').attrib:
                SRR=exp.find('RUN_SET/RUN').attrib['accession']
        #LIBRARY INFO   LIB NAME, LIB STRATEGY, LIB SOURCE, LIB_SELECT
        if exp.find('EXPERIMENT/DESIGN/LIBRARY_DESCRIPTOR/LIBRARY_NAME')!=None:
            LIB_NAME=exp.find('EXPERIMENT/DESIGN/LIBRARY_DESCRIPTOR/LIBRARY_NAME').text
        if exp.find('EXPERIMENT/DESIGN/LIBRARY_DESCRIPTOR/LIBRARY_STRATEGY')!=None:
            LIB_STRAT=exp.find('EXPERIMENT/DESIGN/LIBRARY_DESCRIPTOR/LIBRARY_STRATEGY').text
        if exp.find('EXPERIMENT/DESIGN/LIBRARY_DESCRIPTOR/LIBRARY_SOURCE')!=None:
            LIB_SOURCE=exp.find('EXPERIMENT/DESIGN/LIBRARY_DESCRIPTOR/LIBRARY_SOURCE').text
        if exp.find('EXPERIMENT/DESIGN/LIBRARY_DESCRIPTOR/LIBRARY_SELECTION')!=None:
            LIB_SELECT=exp.find('EXPERIMENT/DESIGN/LIBRARY_DESCRIPTOR/LIBRARY_SELECTION').text
        #PLATFORM INFO   PLATFORM, INSTRUMENT MODEL
        if exp.find('EXPERIMENT/PLATFORM')!=None:
            for platform in exp.find('EXPERIMENT/PLATFORM'):
                PLATFORM=(platform.tag) 
            plat_path='EXPERIMENT/PLATFORM/'+PLATFORM+'/INSTRUMENT_MODEL'
        if exp.find(plat_path)!=None:
            INSTRUMENT_MODEL=exp.find(plat_path).text
        #SUBMISSION INFO    LAB NAME, CENTER NAME, SUBMISSION ID, STUDY ABSTRACT, PUBLICATION DATE
        if exp.find('SUBMISSION')!=None: 
            if 'lab_name' in exp.find('SUBMISSION').attrib:
                LAB_NAME=exp.find('SUBMISSION').attrib['lab_name']
            if 'center_name' in exp.find('SUBMISSION').attrib:
                CENTER_NAME=exp.find('SUBMISSION').attrib['center_name']
            if 'alias' in exp.find('SUBMISSION').attrib:
                SUB_ID=exp.find('SUBMISSION').attrib['alias']
        if exp.find('STUDY/DESCRIPTOR/STUDY_ABSTRACT')!=None:
            ABSTRACT=exp.find('STUDY/DESCRIPTOR/STUDY_ABSTRACT').text
        if exp.find('STUDY/DESCRIPTOR/STUDY_TITLE')!=None:
            STUDY_TITLE=exp.find('STUDY/DESCRIPTOR/STUDY_TITLE').text
        if exp.find('RUN_SET/RUN')!=None and 'published' in exp.find('RUN_SET/RUN').attrib:
            PDAT=exp.find('RUN_SET/RUN').attrib['published'].split()[0]
        
        #SAMPLE INFO   ORGANISM, SAMPLE NAME, SAMPLE_TITLE, EXP_TITLE, ATTRIBUTES
        if exp.find('SAMPLE/SAMPLE_NAME/SCIENTIFIC_NAME')!=None:
            ORGANISM=exp.find('SAMPLE/SAMPLE_NAME/SCIENTIFIC_NAME').text
        if exp.find('SAMPLE')!=None and 'alias' in exp.find('SAMPLE').attrib:
            SAMPLE_NAME=exp.find('SAMPLE').attrib['alias']
        if exp.find('SAMPLE/TITLE')!=None:
            SAMPLE_TITLE=exp.find('SAMPLE/TITLE').text
        if exp.find('EXPERIMENT/TITLE')!=None: 
            EXP_TITLE=exp.find('EXPERIMENT/TITLE').text
        ATTRIBUTES=[]
        if exp.find('SAMPLE/SAMPLE_ATTRIBUTES')!=None:
            for att in exp.find('SAMPLE/SAMPLE_ATTRIBUTES'):
                attribute="-:-"
                if att.find('TAG')!=None and att.find('VALUE')!=None:
                    tag=att.find('TAG').text
                    val=att.find('VALUE').text
                    attribute=str(tag)+":"+str(val)
                    ATTRIBUTES.append(attribute)


        
        line=[str(x) for x in [SRA,SAMN,SRS,SRX,SRP,SRR,LIB_NAME,LIB_STRAT,LIB_SOURCE,LIB_SELECT,PLATFORM,INSTRUMENT_MODEL,LAB_NAME,CENTER_NAME,SUB_ID,ABSTRACT,PDAT,ORGANISM,SAMPLE_NAME,SAMPLE_TITLE,EXP_TITLE," | ".join(ATTRIBUTES),STUDY_TITLE]]
    line_dict={}
    fields=["SRA","SAMN","SRS","SRX","SRP","SRR","LIB_NAME","LIB_STRAT","LIB_SOURCE","LIB_SELECT","PLATFORM","INSTRUMENT_MODEL","LAB_NAME","CENTER_NAME","SUBMISSION_ID","ABSTRACT","PUB_DAT","ORGANISM","SAMPLE_NAME","SAMPLE_TITLE","EXP_TITLE","ATTRIBUTES","STUDY_TITLE"]
    for i in range(len(fields)):
        i=int(i)
        line_dict[fields[i]]=line[i]
    return line_dict
