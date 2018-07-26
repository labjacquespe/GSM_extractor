# GSM extractor<br>
**AUTHOR:** Jean-François Nadeau<br>
**SCRIPT:** GSM_extractor.py<br>
  * This script is used to extract pertinent info from a GSM native XML file from the sra database.<br>
  *It also search the targeted protein for some type of experiments like ChIP-Seq.<br>

**USAGE:**<br>
  * *python3 GSM_extractor.py processes path*<br>
  * processes: Integer indicating the number of paralell process to create. Default is 4.<br>
  * path: Path to the directory containing the xml files. If you need to download the XML files, the path should be the word "online"<br>

**Examples:**<br>
  * python3 GSM_extractor 12 xml/<br>
  * python3 GSM_extractor 8 online<br>

**NOTE:** To use the online mode, please fill the config.ini files that contains the information for the Entrez query and the output directory for the xml files.<br>
