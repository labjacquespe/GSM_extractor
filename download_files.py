import sys
import os
from Bio import Entrez
from multiprocessing import Pool
from functools import partial
import urllib.request
import tarfile
import GSM_extractor

def download(args, cores):
    print("INFO: Online mode selected, downloading files...", file=sys.stderr)
    try:
        orgs=args["Organisms"].rstrip("\n").split(",")
    except:
        print("ERROR: Organisms not specified in congig.ini file.", file=sys.stderr)
        quit()
    print("INFO: Online mode selected, downloading files...", file=sys.stderr)
    for org in orgs:
        print("INFO: --> Processing {}.".format(org), file=sys.stderr)        
        get_files(org,args,cores)
    print("INFO: --> {} files downloaded succesfully.".format(org), file=sys.stderr)

def get_files(org,args,cores):
    outdir=args["xml_out"]
    print ('@INFO: ---> Accounting for the xml already in the xml_output folder...', file=sys.stderr)
    already_there=[]
    for root, dirs, files in os.walk(outdir):
        for name in files:
            path=outdir+"/"+name
            if os.stat(path).st_size>0:
                already_there.append(name.rstrip(".xml)"))
    if len(already_there)>0:
        print ('@INFO: ---> Fetching the remaining GSE IDs...', file=sys.stderr)
    else:
        print ('@INFO: ---> Fetching GSE IDs ...', file=sys.stderr)
    Entrez.email = args["Entrez_email"]
    query=build_query(args,org)
    GSE_handle=get_handle(query,"gds")
    xml_links={}
    for GSE in GSE_handle['IdList']:
        GSE="GSE{}".format(GSE.lstrip("2").lstrip("0"))
        if GSE not in already_there:
            index=len(GSE)-3
            two_firsts=GSE[3:index]
            xml_links[GSE]="ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE{}{}/{}/miniml/{}_family.xml.tgz".format(two_firsts,"nnn",GSE,GSE)
    print ('@INFO: ---> {} series found.'.format(len(xml_links)), file=sys.stderr)
    print ('@INFO: ---> Downloading xml files ...', file=sys.stderr)
    #for GSE in xml_links.keys():
    download_pool=Pool(processes=cores)
    func=partial(download_GSE,xml_links,outdir)
    download_pool.map(func,list(xml_links.keys()))
    download_pool.close()
    print ('@INFO: ---> Done!\n', file=sys.stderr)

def build_query(args, org):
    #adding the org
    query=org+"[ORGN] "
    #adding the search terms:
    if "Search_terms" in args:
        query+="AND {}".format(" AND ".join(args['Search_terms'].split(",")))
    #adding the filter out terms:
    if "Filter_out" in args:
        query+=" NOT ( {} )".format(" OR ".join(args["Filter_out"].split(",")))
    #adding the date
    query+=" AND {}[PDAT]".format(args['Date_range'])
    if "Custom" in args:
        query+=" "+args['Custom']
    return query

def get_handle(query, database):
    handle = Entrez.esearch(db=database,retmax=1000000, term=query)
    dic=Entrez.read(handle)
    return dic

def download_GSE(xml_links,outdir,GSE):
    path="{}/{}.xml.tar.gz".format(outdir,GSE)
    urllib.request.urlretrieve(xml_links[GSE],path)
    tar=tarfile.open(path)
    content=""
    for member in tar.getmembers():
        if "family.xml" in str(member):
            f=tar.extractfile(member)
            content = f.read()
            with open(path.rstrip(".tar.gz"),"wb") as outf:
                outf.write(content)
    os.remove(path)


def main():
    try:
        args=GSM_extractor.read_config()
        cores=int(sys.argv[1])
    except:
        print("ERROR: problem reading arguments or config file.", file=sys.stderr)
    download(args, cores)

if __name__ == "__main__":
    main()