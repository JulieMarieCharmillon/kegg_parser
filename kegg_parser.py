from bs4 import BeautifulSoup
from pathlib import Path
import json
import csv
import argparse

parser = argparse.ArgumentParser(description="returns the KEGG Mapper Reconstruction Result as .csv file")
parser.add_argument('path', metavar='html_path', type=Path, help="the 'KEGG Mapper Reconstruction Result.html' path - please don't forget the bracket around the path")
parser.add_argument('results_type', metavar='results_type', type=str, help="enter 'modules' or 'pathways'")

args = parser.parse_args()
 
RECONSTRUCTION_TYPE = args.results_type

ARG_HTML_DIR = Path(args.path)
with open(ARG_HTML_DIR, "r") as f:
    a = f.read()

soup = BeautifulSoup(a, "html.parser")

### KEGG MAPPER RECONSTRUCTION RESULT MODULES PARSER ###

def next_sibling_cat_content(previous_cat): # allows to avoid repetition of "next_sibling" method
    """
        Get to the next relevant sibling you're looking for, ie the forth next.

        Args:
            previous_cat ([bs4.element.Tag]): [the first relevant sibling]

        Returns:
            [bs4.element.Tag]: [next relevant tag]
    """
    cat_contains = previous_cat.next_sibling.next_sibling
    return cat_contains.next_sibling.next_sibling

def find_modules_supercats(): # returns dictionary {KEYS=supercategories: VALUES=soup -> "bs4.element.Tag"}}
    """
        First level of parsing: create the main dictionary.

        Returns:
            [dictionary]: {keys=SUPERCAT:values=soup}
    """
    supercat_list = []
    supercat_soup_list = []
 
    pathways_modules_name = soup.find("div", id='list').find("b")
    pathways_modules_soup = soup.find("div", id='list').find("b").next_sibling.next_sibling
    
    supercat1 = pathways_modules_soup.find("b")
    supercat_list.append(supercat1.string)
        
    for entry in pathways_modules_soup.b.next_siblings:
        
        if str(type(entry)) == "<class 'bs4.element.Tag'>":
            if entry.ul != None:
                supercat_soup_list.append(entry)
            else:
                supercat_list.append(entry.string)
    
    
    
    signatures_modules_name = next_sibling_cat_content(pathways_modules_name)
    signatures_modules_soup = signatures_modules_name.next_sibling.next_sibling

    supercat2 = signatures_modules_soup.find("b")
    supercat_list.append(supercat2.string)

    for entry in signatures_modules_soup.b.next_siblings:
        
        if str(type(entry)) == "<class 'bs4.element.Tag'>":
            if entry.ul != None:
                supercat_soup_list.append(entry)
            else:
                supercat_list.append(entry.string)
    
    d= {supercat_list[i]: supercat_soup_list[i] for i in range(len(supercat_list))}

    # with open("supercat_dict_23Nov21.json", "w") as f:    #checking the dictionary
    #     json.dump(d, f, indent=4)                 # if you want to try: don't forget to add str() around supercat_soup_list[i] in second_dict declaration

    
    return d

def find_modules_cats(): # returns dictionary {key=Categories: values=soup -> "bs4.element.Tag"}
    """
        From supercat_dict's values: create a second dictionary {keys=categories || valeur = soup -> bs4.element.Tag}

        Returns:
            [dict]: [return nested dictionary {keys=categories || valeur = soup -> bs4.element.Tag} inside supercat_dict]
    """
    
    d=find_modules_supercats()
    for supercat in d:
        value = d.get(supercat)
        temp_list_keys=[]
        temp_list_val=[]
        
        for child in value.children:
            if str(type(child))=="<class 'bs4.element.NavigableString'>" and str(child) != "\n":
                temp_list_keys.append(str(child.string).lstrip("\n ").rstrip("\n "))
            else:
                temp_list_val.append(child)
        
        second_dict = {temp_list_keys[i]: temp_list_val[i] for i in range(len(temp_list_keys))}
        d[supercat]=second_dict
        
    return d

def find_modules_sous_cats(): # returns dictionary {keys="MO reaction": values= soup -> "bs4.element.Tag"}
    """
        From find_modules_cats extracts "MO reation name" as keys and soup as values
        This new dictionary is nested inside the the cat_dict

        Returns:
            [dict]: supercat_dict={Supercategories: {Catgories: {MOxxx: soup->bs4.element.Tag}}}
    """
    
    d=find_modules_cats()
    
    for supcat in d:
        second_dict = d.get(supcat)
        for cat in second_dict:
            value = second_dict.get(cat)   
            MOs=value.find("a").string
            
            MOs_name = str(value.a.next_sibling) + str(value.a.find_next("a").string)+ str(value.a.next_sibling.next_sibling.next_sibling)
            
            MOs_full_name = str(str(MOs)+ " " + MOs_name.replace("\xa0", " ").lstrip(" "))
            
            Mo_div_value = value.div

            Mo_div_dict ={}
            Mo_div_dict[MOs_full_name]=Mo_div_value

            for entry in value.li.next_siblings:
                if str(type(entry))=="<class 'bs4.element.NavigableString'>":
                    continue
                elif str(type(entry))== "<class 'bs4.element.Tag'>":
                    MOs = entry.find("a").string
                    MOs_name = str(entry.a.next_sibling) + str(entry.a.find_next("a").string)+ str(entry.a.next_sibling.next_sibling.next_sibling)
                    MOs_full_name = str(str(MOs)+ " " + MOs_name.replace("\xa0", " ").lstrip(" "))
                    Mo_div_value = entry.div
                    Mo_div_dict[MOs_full_name]=Mo_div_value
            
            second_dict[cat]=Mo_div_dict
    
    # with open("find_modules_sous_cats_24Nov21.json", "w") as f:
    #     json.dump(d, f, indent=4)             # if you want to try: don't forget to add str() around Mo_div_value 

    return d

def find_modules_ko_column(): # creates a liste of tuple (KO, gene) for each MO -> returns the full dictionnary
    """
        From find_modules_sous_cats values creates a list of tuples (KO, gene) for each MO

        Returns:
            [dict]: [The FULL dictionary AND a json file containing the full dictionary]
    """
    
    d = find_modules_sous_cats()
    
    for supercat in d:
        second_dict = d.get(supercat)
        for cat in second_dict:
            last_dict = second_dict.get(cat)
            for Mo in last_dict:
                last_soup = last_dict.get(Mo)
                Ko_gene_list=[]
                if last_soup.a != None:
                    KO_name = last_soup.a.string
                else:
                    KO_name = ""
                
                KO_name_list=[KO_name]
                KO_gene_list=[]
                if last_soup.dt != None:
                    for entry in last_soup.dt.next_siblings:
                        if str(type(entry))=="<class 'bs4.element.NavigableString'>":
                            pass
                        elif str(type(entry))== "<class 'bs4.element.Tag'>":
                            
                            if entry.a != None:
                                next_KO_name = entry.a.string
                                KO_name_list.append(next_KO_name)
                            else:
                                next_KO_gene = entry.string
                                KO_gene_list.append(next_KO_gene)
                            
                final_tuple_list = [(a, b) for a,b in zip(KO_name_list, KO_gene_list)]        
                last_dict[Mo]=final_tuple_list
                
    MODULE_BASE_DIR = ARG_HTML_DIR.parent
    RESULTS_DIR = MODULE_BASE_DIR / "RESULTS"
    RESULTS_DIR.mkdir(exist_ok=True, parents=True)

    FILE_DIR = RESULTS_DIR / (ARG_HTML_DIR.stem + "_KEGG_MODULES_dictionary.json")
    
    
    with open(FILE_DIR, "w") as f: # FINAL DICTIONARY
        json.dump(d, f, indent=4)
    return d

def modules_ko_count(): 
    """
        Return the total number of KO.
    """
    
    d=find_modules_ko_column()
    ko_number = 0
    for SUPERCAT in d:
        cat_dict = d.get(SUPERCAT)
        for cat in cat_dict:
            Mo_dict = cat_dict.get(cat)
            for Mo in Mo_dict:
                Ko_list= Mo_dict.get(Mo)
                for ko_tuple in Ko_list:
                    ko_number += 1
    print(f" There are {ko_number} KO.") 
    return ko_number

def write_modules_csv(): # writes the dictionary in a csv file
    """
        Returns the csv file with the KO_list
    """
    
    d=find_modules_ko_column()

    MODULE_BASE_DIR = ARG_HTML_DIR.parent
    RESULTS_DIR = MODULE_BASE_DIR / "RESULTS"
    CSV_FILE_DIR = RESULTS_DIR / (ARG_HTML_DIR.stem + '_KEGG_MODULES_Grouped_KO.csv')
    
    ROWS_list=[]
                                
    for SUPERCAT in d:
        cat_dict = d.get(SUPERCAT)
        for cat in cat_dict:
            Mo_dict = cat_dict.get(cat)
            for Mo in Mo_dict:
                Ko_list= Mo_dict.get(Mo)
                for ko_gene in Ko_list:
                    row_gene = ko_gene[1]
                    row_ko = ko_gene[0] + ";" + row_gene
                    row_Mo = str(Mo) + ";" + row_ko
                    row_cat = str(cat) + ";" + row_Mo
                    row_SUP = str(SUPERCAT) + ";" + row_cat

                    ROWS_list.append(row_SUP)

    ROWS_list_splitted = []
    for row in ROWS_list:
        new_row = row.split(";")
        ROWS_list_splitted.append(new_row)

    with open(CSV_FILE_DIR, "w", newline='') as csvfile:   
        writer = csv.writer(csvfile, delimiter=',')
        
        writer.writerow(['Super Categories', 'Categories', 'MO', 'KO', 'gene'])
        
        writer.writerows(ROWS_list_splitted)

    modules_ko_count()

### END OF KEGG MAPPER RECONSTRUCTION RESULT MODULES PARSER ###



                                        ### KEGG MAPPER RECONSTRUCTION RESULT PATHWAYS PARSER ###

def find_pathways_supercats(): # returns dictionary {KEYS= SUPERCATS : VALUES = soup -> <class 'bs4.element.Tag'>}
    """ 
        Returns: first dictionary with {key=SUPERCAT: values=soup <ul>}
    """
    
    supercat_list = []
    supercat_soup_list = []

    Supercat1 = soup.find("div", id='list').find("b")
    supercat_list.append(Supercat1.string)
    
    for entry in soup.find("div", id='list').find("b").next_siblings:
        
        if str(type(entry)) == "<class 'bs4.element.Tag'>":
            if entry.ul != None:
                supercat_soup_list.append(entry)
            else:
                supercat_list.append(entry.string)
    
    d= {supercat_list[i]: supercat_soup_list[i] for i in range(len(supercat_list))}

    # with open("find_pathways_supercat_24Nov21.json", 'w') as f:       # dictionary checking
    #     json.dump(d, f, indent=4)                         # if you want to try: don't forget to add str() around supercat_soup_list[i] in d declaration

    return d

def find_pathways_cats():    # returns {key=cat: value=soup -> '<class 'bs4.element.Tag'>'}
    """
        From pathway_supercat_dict creates the second dictionary with keys=CATEGORIES and Values = soup <class 'bs4.element.Tag'>

        Returns:
            [dict]: {key = supercategory: value = {key=cat: value=soup -> '<class 'bs4.element.Tag'>'}}
    """
    d = find_pathways_supercats()
    for Cat in d:
        cat_soup = d.get(Cat)

        temp_list_keys=[]
        temp_list_val=[]
        for child in cat_soup.children:
            if str(type(child))=="<class 'bs4.element.NavigableString'>" and str(child) != "\n":
                temp_list_keys.append(str(child.string).lstrip("\n ").rstrip("\n "))
            else:
                temp_list_val.append(child)
        
        second_dict = {temp_list_keys[i]: temp_list_val[i] for i in range(len(temp_list_keys))}
        d[Cat]=second_dict

    # with open("find_pathways_cats_24Nov21.json", 'w') as f:       #dictionary check
    #     json.dump(d, f, indent=4)                             # if you want to try: don't forget to add str() around temp_list_val[i] in second_dict declaration

    return d

def find_pathways(): # returns dictionary {keys="pathways": values= soup -> "bs4.element.Tag"
    """ 
        From pathways_cat_dict extract "pathways name" as keys and soup as values
        This new dictionary is nested inside the the cat_dict

        Returns:
            [dict]: d= {Supercategories: {Catgories: {Pathways: soup->bs4.element.Tag}}}
    """

    d=find_pathways_cats()
    
    for supcat in d:    
        second_dict = d.get(supcat)
        for cat in second_dict:
            value = second_dict.get(cat)
            MOs=value.find("a").string
            
            MOs_name = str(value.a.next_sibling) + str(value.a.find_next("a").string)+ str(value.a.next_sibling.next_sibling.next_sibling)
            
            MOs_full_name = str(str(MOs)+ " " + MOs_name.replace("\xa0", " ").lstrip(" ").rstrip("\n"))
            
            Mo_div_value = value.div

            Mo_div_dict ={}
            Mo_div_dict[MOs_full_name]=Mo_div_value

            for entry in value.li.next_siblings:
                if str(type(entry))=="<class 'bs4.element.NavigableString'>":
                    continue
                elif str(type(entry))== "<class 'bs4.element.Tag'>":
                    MOs = entry.find("a").string
                    MOs_name = str(entry.a.next_sibling) + str(entry.a.find_next("a").string)+ str(entry.a.next_sibling.next_sibling.next_sibling)
                    MOs_full_name = str(str(MOs)+ " " + MOs_name.replace("\xa0", " ").lstrip(" ").rstrip("\n"))
                    Mo_div_value = entry.div
                    Mo_div_dict[MOs_full_name]=Mo_div_value
                    
            second_dict[cat]=Mo_div_dict
            
    # with open("find_pathways_names_24Nov21.json", 'w') as f:       #dictionary check
    #     json.dump(d, f, indent=4)                        # if you want to try: don't forget to add str() around Mo_div_value

    return d

def pathways_ko_column(): # creates a liste of tuple (KO, gene) for each 'pathways' -> returns the full dictionnary
    """
        From pathways_dict values creates a list of tuples (KO, gene) for each pathway.

        Returns:
            [dict]: [The FULL dictionary AND a json file containing the full dictionary]
    """
    d = find_pathways()
    
    for supercat in d:
        second_dict = d.get(supercat)
        for cat in second_dict:
            last_dict = second_dict.get(cat)
            for Mo in last_dict:
                last_soup = last_dict.get(Mo)
                Ko_gene_list=[]
                if last_soup.a != None:
                    KO_name = last_soup.a.string
                else:
                    KO_name = ""
                
                KO_name_list=[KO_name]
                KO_gene_list=[]
                if last_soup.dt != None:
                    for entry in last_soup.dt.next_siblings:
                        if str(type(entry))=="<class 'bs4.element.NavigableString'>":
                            pass
                        elif str(type(entry))== "<class 'bs4.element.Tag'>":
                            
                            if entry.a != None:
                                next_KO_name = entry.a.string
                                KO_name_list.append(next_KO_name)
                            else:
                                next_KO_gene = entry.string
                                KO_gene_list.append(next_KO_gene)
                            
                final_tuple_list = [(a, b) for a,b in zip(KO_name_list, KO_gene_list)]        
                last_dict[Mo]=final_tuple_list
                
    PATHWAYS_BASE_DIR = ARG_HTML_DIR.parent
    RESULTS_DIR = PATHWAYS_BASE_DIR / "RESULTS"
    RESULTS_DIR.mkdir(exist_ok=True, parents=True)
    FILE_DIR = RESULTS_DIR / (ARG_HTML_DIR.stem + "_KEGG_PATHWAYS_dictionary.json")

    with open(FILE_DIR, "w") as f:
        json.dump(d, f, indent=4)

    return d

def pathways_ko_count():
    """
        Returns the total number of KO.
    """
    d=pathways_ko_column()
    ko_number = 0
    for SUPERCAT in d:
        cat_dict = d.get(SUPERCAT)
        for cat in cat_dict:
            Mo_dict = cat_dict.get(cat)
            for Mo in Mo_dict:
                Ko_list= Mo_dict.get(Mo)
                for ko_tuple in Ko_list:
                    ko_number += 1
    print(f" There are {ko_number} KO.") 

def write_pathways_csv():   # writes the dictionary in a csv file
    """
        Returns the csv file with the KO_list
    """
    d=pathways_ko_column()
    ROWS_list=[]
                                
    for SUPERCAT in d:
        cat_dict = d.get(SUPERCAT)
        for cat in cat_dict:
            pathways_dict = cat_dict.get(cat)
            for pathway in pathways_dict:
                Ko_list= pathways_dict.get(pathway)
                for ko_gene in Ko_list:
                    row_gene = ko_gene[1]     
                    row_ko = ko_gene[0] + ";" + row_gene
                    row_pathway = str(pathway) + ";" + row_ko
                    row_cat = str(cat) + ";" + row_pathway
                    row_SUP = str(SUPERCAT) + ";" + row_cat

                    ROWS_list.append(row_SUP)

    ROWS_list_splitted = []
    for row in ROWS_list:
        new_row = row.split(";")
        ROWS_list_splitted.append(new_row)

    PATHWAYS_BASE_DIR = ARG_HTML_DIR.parent
    RESULTS_DIR = PATHWAYS_BASE_DIR / "RESULTS"
    RESULTS_DIR.mkdir(exist_ok=True, parents=True)
    CSV_FILE_DIR = RESULTS_DIR / (ARG_HTML_DIR.stem + '_KEGG_PATHWAYS_Grouped_KO.csv')

    with open(CSV_FILE_DIR, "w", newline='') as csvfile:   
        writer = csv.writer(csvfile, delimiter=',')
        
        writer.writerow(['Super Categories', 'Categories', 'Pathways', 'KO', 'gene'])
        
        writer.writerows(ROWS_list_splitted)

    pathways_ko_count()


# ARGPARSE BEHAVIOR:
if RECONSTRUCTION_TYPE == "modules":       
    write_modules_csv()
elif RECONSTRUCTION_TYPE == "pathways":
    write_pathways_csv()
else:
    print("Please enter a valid option : 'modules' or 'pathways'.")


# if __name__ == "__main__":
#     # find_modules_supercats()
#     # find_modules_cats()
#     # find_modules_sous_cats()
#     # find_modules_ko_column()
#     # write_modules_csv()
#     ### PATHWAYS: ###
#     # find_pathways_supercats()
#     # find_pathways_cats()
#     # find_pathways()
#     write_pathways_csv()