
import sys
import os
from Bio import Entrez
from multiprocessing import Pool
from functools import partial
import wget
import tarfile

#
##
### Starting function to call when using this module
def Download(args, cores, outdir):
    Entrez.email = args["Entrez_email"]
    existing=[]
    ftpLinks={}

    print(GetMessage("OnlineMode"), file=sys.stderr)
    try:
        orgs=args["Organisms"].rstrip("\n").split(",")
    except:
        print(GetMessage("NoOrgs"), file=sys.stderr)
        sys.exit()
    existing = CheckExisting(outdir)
    query=BuildQuery(args)
    handle=getHandle(query,"gds")
    ftpLinks=GetFTP(handle, existing)

    print (GetMessage("Downloading"), file=sys.stderr)
    downloadPool=Pool(processes=cores)
    func=partial(DownloadFiles,ftpLinks,outdir)
    downloadPool.map(func,list(ftpLinks.keys()))
    downloadPool.close()
    print (GetMessage("DownloadSuccess"), file=sys.stderr)
#
## BUILD QUERY
### Function to create the entrex query, based on the info entered in the config.ini file
def BuildQuery(args):
    query="("+args["Organisms"].replace(",","[ORGN] OR ")+"[ORGN])"
    if "Search_terms" in args:
        query+=" AND {}".format(" AND ".join(args['Search_terms'].split(",")))
    if "Filter_out" in args:
        query+=" NOT ( {} )".format(" OR ".join(args["Filter_out"].split(",")))
    query+=" AND {}[PDAT]".format(args['Date_range'])
    if "Custom" in args:
        query+=" "+args['Custom']
    print(GetMessage("SendingQuery"),query, file=sys.stderr)
    return query
#
## CHECK EXISTING
### Function to check if there is already files in the output directory. The tool will not download the existing xml files
def CheckExisting(outdir):
    existing=[]
    print(GetMessage("CheckExisting"), file=sys.stderr)
    for root, dirs, files in os.walk(outdir):
        for name in files:
            path=outdir+"/"+name
            if os.stat(path).st_size>0:
                existing.append(name.rstrip(".xml"))
    if existing:
        print(str(len(existing)),GetMessage("FilesFound"), file=sys.stderr)
    else:
        print(GetMessage("NoFilesFound"), file=sys.stderr)
    return existing    
#
## DOWNLOAD FILES
### Function to download the GSE xml files
def DownloadFiles(ftpLinks, outdir, GSE):
    path="{}/{}.xml.tar.gz".format(outdir,GSE)
    try:
        wget.download(ftpLinks[GSE], out=path)
        tar=tarfile.open(path)
        content=""
        for member in tar.getmembers():
            if "family.xml" in str(member):
                f=tar.extractfile(member)
                content = f.read()
                with open(path.rstrip(".tar.gz"),"wb") as outf:
                    outf.write(content)
        os.remove(path)
    except:
        print(GetMessage("CannotDownload"),GSE, ftpLinks[GSE], file=sys.stderr)
#
## GET FTP
### Construction of the ftp download links
def GetFTP(handle, existing):
    ftpLinks={}
    for GSE in handle['IdList']:
        GSE="GSE{}".format(GSE.lstrip("2").lstrip("0")) #Removing the database code from the ID (ids in the query starts with the database ID (20 and 200 for GSE) followed by the GEO ID )
        if GSE not in existing:
            short=GSE[3:len(GSE)-3]
            ftpLinks[GSE]="ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE{}{}/{}/miniml/{}_family.xml.tgz".format(short,"nnn",GSE,GSE)
    return ftpLinks
#
## GET HANDLE
### Function to fetch the Entrez query handle (dictio with all the details such as count, ids etc..)
def getHandle(query, database):
    queryResult = Entrez.esearch(db=database,retmax=1000000, term=query)
    handle=Entrez.read(queryResult)
    return handle
#
## GET MESSAGE
### Gunction to regroup messages to make modifications easier
def GetMessage(title):
    messages={
        "OnlineMode":"INFO: Online mode selected, downloading files...",
        "CheckExisting": "INFO: Searching for existing files in the output folder...",
        "FilesFound": " INFO: Files has been found in the output directory. Corresponding files will not be downloaded again.",
        "NoFilesFound": "INFO: No files found in the output directory.",
        "FetchingIDs": "INFO: Fetching GSE IDs matching the query...",
        "SendingQuery": "INFO: Sending query:",
        "Downloading": "INFO: Downloading the files...",
        "DownloadSuccess": "INFO: files downloaded succesfully.",


        "CannotDownload": "ERROR: This GSE cannot be downloaded:",
        "NoOrgs":"ERROR: Organisms not specified in congig.ini file."
        }
    return messages[title]