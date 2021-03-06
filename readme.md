# GSM extractor<br>
**AUTHOR:** Jean-François Nadeau<br>
  * This script is used to extract pertinent info in natives XML files from the sra database.<br>
  * For some  type of experiments like ChIP-Seq, it will also search the targeted protein.<br>

**USAGE:**<br>
  * *python3 GSM_extractor.py processes path organism*<br>
  * processes: Integer indicating the number of paralell process to create.<br>
  * path: Path to the directory containing the xml files. If you need to download the XML files, the path should be the word "online". The script will get the xml files from the SRA database by sending Entrez queries using the info specified in the config.ini file.<br>
  * Organism: The alias of the organism you are analysing. This argument is used to load the corresponding gene names dictionary. Currently, you can use saccer, elegans and human.

**EXAMPLES:**<br>
  * python3 GSM_extractor 12 xml/<br>
  * python3 GSM_extractor 8 online<br>

**NOTE:** To use the online mode, fill the config.ini files that contains the information for the Entrez query and the output directory for the xml files.<br>
