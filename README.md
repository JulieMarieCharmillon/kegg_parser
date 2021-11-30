# kegg_parser
An HTML parser command line app for "KEGG Mapper – Reconstruct Result" pages.

## Setup

<p>Built with python3.9</p>
<p>kegg_parser.py uses five modules: beautifulsoup4*, pathlib, json, csv and argparse. Be sure to have all of them up-to-date on your environment.</p>
<p> txt_ko_list.py uses pandas*, pathlib and argparse.</p>
<p>* not built-in modules</p> 
<hr>

On your terminal use the following command line:

```
python txt_ko_list.py "initial_tsv_file.tsv"
```

<p>Returns a txt file with the KO list you want to use in the kegg Reconstruct tool (https://www.genome.jp/kegg/mapper/reconstruct.html)</p>
<p>Output file will be find in the folder containing your tsv file under the name "initial_tsv_file_name_KO_LIST.txt"</p>

<hr>

On your terminal use the following command line:

```
python kegg_parser "KEGG Mapper – Reconstruct Result.html" modules
```

<p>Returns a csv file with the list of all the KO from the KEGG Mapper – Reconstruct Result modules pages (under the name: "KEGG_MODULES_Grouped_KO.csv") and a json file with the corresponding dictionary (under the name: "KEGG_MODULES_dictionary.json"). </p>
<p>It works with all "views" : "complete only", "including 1 block missing" and "including any incomplete". But please download the "views" pages separately.</p>

<p>Output files will be found in a "result" subfolder within the folder containing the initial "KEGG Mapper – Reconstruct Result.html" file.</p>

<hr>

On your terminal use the following command line:

```
python kegg_parser "KEGG Mapper – Reconstruct Result.html" pathways
```
Returns a csv file with the list of all the KO from the KEGG Mapper – Reconstruct Result pathways pages (under the name: "KEGG_PATHWAYS_Grouped_KO.csv") and a json file with the corresponding dictionary (under the name: "KEGG_PATHWAYS_dictionary.json").
<p>Output files will be found in a "result" subfolder within the folder containing the initial "KEGG Mapper – Reconstruct Result.html" file.</p>
