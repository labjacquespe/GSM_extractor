import sys
import os
import utilities
import downloadFiles
import extractMeta


def main():
    lines=[]
    try:
        args=utilities.readConfig()
        cores=int(sys.argv[1])
        path=sys.argv[2]
    except:
        PrintHelp()
        sys.exit()
    if path=="online":
        path=args["xml_out"]
        downloadFiles.Download(args,cores,path)
    if os.path.isdir(path):
        lines=ExtractMetadata(args, path)
    elif os.path.isfile(path):
        with open(path, "r") as inf:
            lines=[x.rstrip("\n") for x in inf.readlines()]
    lines=RemoveDuplicates(lines)
    #PROCESS LINE TODO

def RemoveDuplicates(lines):
    dictio = {}
    for line in lines:
        line=line.rstrip("\n").split("\t")
        GSM=line[0]
        GSE=line[1]
        if GSM in dictio:
            if GSE not in dictio[GSM][1]:
                dictio[GSM][1]+=","+GSE
        else:
            dictio[GSM]=line
    return dictio.values()

def ExtractMetadata(args, path):
    lines=[]
    for root, dirs, files in os.walk(path):
        files=["{}/{}".format(path.rstrip("/"),filename) for filename in files]
        for f in files:
            lines+=extractMeta.get_metadata(f)
    if args["Create_tsv"].lower()=="true":
        with open("GSM_extractor_out.tsv","w") as outf:
            outf.write("\n".join(lines))
    return lines

def PrintHelp():
    text="""SCRIPT: GSM_extractor.py
    This tool uses GEO metadata to create a simplified and organised TSV database, including the target protein/mark for some of the assay types.

    USAGE:
        python3 GSM_extractor.py processes path
        -processes: Integer indicating the number of paralell process to create.
        -path: Path to the directory containing the xml files. If you need to download the XML files, assign an existing directory to save the files and set this argument to "online".
    Examples:
        python3 GSM_extractor 12 xmlDirectory/
        python3 GSM_extractor 8 online

    NOTE: To use the online mode, please fill the config.ini files that contains the information for the Entrez query and the output directory for the xml files."""
    print(text, file=sys.stderr)
    

if __name__ == "__main__":
    main()