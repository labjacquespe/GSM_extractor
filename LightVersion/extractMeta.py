import xml.etree.ElementTree as ET
from multiprocessing import Pool
import sys
import os

def get_metadata(xml_path):
    try:
        tree = ET.parse(xml_path)
        root = tree.getroot()
    except:
        return []

    lines=[]
    GSM_list=[]

    contributor="NOT_FOUND"
    organization="NOT_FOUND"
    GPL="NOT_FOUND"
    platform_technology="NOT_FOUND"
    GSE="NOT_FOUND"
    series_title="NOT_FOUND"
    series_sumary="NOT_FOUND"
    series_design="NOT_FOUND"
    GSM="NOT_FOUND"
    submission_date="NOT_FOUND"
    release_date="NOT_FOUND"
    sample_title="NOT_FOUND"
    library_strategy="NOT_FOUND"
    library_source="NOT_FOUND"
    library_selection="NOT_FOUND"
    sample_source="NOT_FOUND"
    organism="NOT_FOUND"
    attributes="NOT_FOUND"
    treatment="NOT_FOUND"
    molecule="NOT_FOUND"


    #Platform common information
    platform=root.find('{http://www.ncbi.nlm.nih.gov/geo/info/MINiML}Platform')
    if platform!=None:
        GPL=platform.attrib["iid"]
        if platform.find("{http://www.ncbi.nlm.nih.gov/geo/info/MINiML}Technology")!=None:
            platform_technology=platform.find("{http://www.ncbi.nlm.nih.gov/geo/info/MINiML}Technology").text

    #Serie common information
    series=root.find('{http://www.ncbi.nlm.nih.gov/geo/info/MINiML}Series')
    if series!=None:
        GSE=series.attrib["iid"]
        if series.find('{http://www.ncbi.nlm.nih.gov/geo/info/MINiML}Title')!=None:
            series_title=series.find('{http://www.ncbi.nlm.nih.gov/geo/info/MINiML}Title').text
        if series.find('{http://www.ncbi.nlm.nih.gov/geo/info/MINiML}Summary')!=None:
            series_sumary=series.find('{http://www.ncbi.nlm.nih.gov/geo/info/MINiML}Summary').text
        if series.find('{http://www.ncbi.nlm.nih.gov/geo/info/MINiML}Overall-Design')!=None:
            series_design=series.find('{http://www.ncbi.nlm.nih.gov/geo/info/MINiML}Overall-Design').text

    #Sample information
    for sample in root.findall('{http://www.ncbi.nlm.nih.gov/geo/info/MINiML}Sample'):
        if sample!=None:
            if sample.find('{http://www.ncbi.nlm.nih.gov/geo/info/MINiML}Type')!=None and sample.find('{http://www.ncbi.nlm.nih.gov/geo/info/MINiML}Type').text=="SRA":
                GSM=sample.attrib["iid"]
                if sample.find('{http://www.ncbi.nlm.nih.gov/geo/info/MINiML}Title')!=None:
                    sample_title=sample.find('{http://www.ncbi.nlm.nih.gov/geo/info/MINiML}Title').text
                if sample.find('{http://www.ncbi.nlm.nih.gov/geo/info/MINiML}Status/{http://www.ncbi.nlm.nih.gov/geo/info/MINiML}Submission-Date')!=None:
                    submission_date=sample.find('{http://www.ncbi.nlm.nih.gov/geo/info/MINiML}Status/{http://www.ncbi.nlm.nih.gov/geo/info/MINiML}Submission-Date').text
                if sample.find('{http://www.ncbi.nlm.nih.gov/geo/info/MINiML}Status/{http://www.ncbi.nlm.nih.gov/geo/info/MINiML}Release-Date')!=None:
                    release_date=sample.find('{http://www.ncbi.nlm.nih.gov/geo/info/MINiML}Status/{http://www.ncbi.nlm.nih.gov/geo/info/MINiML}Release-Date').text
                if sample.find('{http://www.ncbi.nlm.nih.gov/geo/info/MINiML}Library-Strategy')!=None:
                    library_strategy=sample.find('{http://www.ncbi.nlm.nih.gov/geo/info/MINiML}Library-Strategy').text
                if sample.find('{http://www.ncbi.nlm.nih.gov/geo/info/MINiML}Library-Source')!=None:
                    library_source=sample.find('{http://www.ncbi.nlm.nih.gov/geo/info/MINiML}Library-Source').text
                if sample.find('{http://www.ncbi.nlm.nih.gov/geo/info/MINiML}Library-Selection')!=None:
                    library_selection=sample.find('{http://www.ncbi.nlm.nih.gov/geo/info/MINiML}Library-Selection').text
                if sample.find('{http://www.ncbi.nlm.nih.gov/geo/info/MINiML}Contact-Ref')!=None:
                    contact=sample.find('{http://www.ncbi.nlm.nih.gov/geo/info/MINiML}Contact-Ref').attrib["ref"]
                    for author in root.findall('{http://www.ncbi.nlm.nih.gov/geo/info/MINiML}Contributor'):
                        if author.find('{http://www.ncbi.nlm.nih.gov/geo/info/MINiML}Person')!=None and author.attrib["iid"]==contact:
                            contributor="{} {}".format(author.find('{http://www.ncbi.nlm.nih.gov/geo/info/MINiML}Person/{http://www.ncbi.nlm.nih.gov/geo/info/MINiML}First').text,root.findall('{http://www.ncbi.nlm.nih.gov/geo/info/MINiML}Contributor')[0].find('{http://www.ncbi.nlm.nih.gov/geo/info/MINiML}Person/{http://www.ncbi.nlm.nih.gov/geo/info/MINiML}Last').text)
                            if author.find('{http://www.ncbi.nlm.nih.gov/geo/info/MINiML}Organization')!=None:
                                organization=author.find('{http://www.ncbi.nlm.nih.gov/geo/info/MINiML}Organization').text
                channel_count=sample.find('{http://www.ncbi.nlm.nih.gov/geo/info/MINiML}Channel-Count').text
                for channel in sample.findall('{http://www.ncbi.nlm.nih.gov/geo/info/MINiML}Channel'):
                    position=channel.attrib["position"]
                    if channel.find('{http://www.ncbi.nlm.nih.gov/geo/info/MINiML}Source')!=None:
                        sample_source=channel.find('{http://www.ncbi.nlm.nih.gov/geo/info/MINiML}Source').text
                    if channel.find('{http://www.ncbi.nlm.nih.gov/geo/info/MINiML}Organism')!=None:
                        organism=channel.find('{http://www.ncbi.nlm.nih.gov/geo/info/MINiML}Organism').text
                    if channel.find('{http://www.ncbi.nlm.nih.gov/geo/info/MINiML}Treatment-Protocol')!=None:
                        treatment=channel.find('{http://www.ncbi.nlm.nih.gov/geo/info/MINiML}Treatment-Protocol').text
                    if channel.find('{http://www.ncbi.nlm.nih.gov/geo/info/MINiML}Molecule')!=None:
                        molecule=channel.find('{http://www.ncbi.nlm.nih.gov/geo/info/MINiML}Molecule').text
                    attributes=[]
                    if channel.findall('{http://www.ncbi.nlm.nih.gov/geo/info/MINiML}Characteristics')!=None:
                            for attribute in channel.findall('{http://www.ncbi.nlm.nih.gov/geo/info/MINiML}Characteristics'):
                                if "tag" in attribute.attrib:
                                    attributes.append(":".join([attribute.attrib['tag'].rstrip(" "),attribute.text.replace("\n","")]))
                    if channel_count!="1":
                        corrected_GSM="{}-{}".format(GSM,position)
                    else:
                        corrected_GSM=GSM
                    if corrected_GSM not in GSM_list and organism.lower() in ["homo sapiens","mus musculus","drosophila melanogaster","saccharomyces cerevisiae","arabidopsis thaliana","caenorhabditis elegans","schizosaccharomyces pombe","rattus norvegicus","danio rerio","gallus gallus","pan troglodytes"]:
                        GSM_list.append(corrected_GSM)
                        lines.append("\t".join([corrected_GSM,GSE,GPL,organism,library_strategy,sample_title," | ".join(attributes),sample_source,molecule,platform_technology,library_source,library_selection,series_title,series_sumary,series_design,contributor,organization,release_date,submission_date,treatment]).replace("\n","").replace("  ",""))
    return lines


# Usable in standalone, pass the xml path as argument.
def main():
    lines=[]
    if os.path.isdir(sys.argv[1]):
        for root, dirs, files in os.walk(sys.argv[1]):
            files=["{}/{}".format(sys.argv[1].rstrip("/"),filename) for filename in files]
            for file in files:
                lines+=get_metadata(file)
    elif os.path.isfile(sys.argv[1]):
        lines=get_metadata(sys.argv[1])
    print(*list(set(lines)), sep='\n')

if __name__ == "__main__":
    main()
