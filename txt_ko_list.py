import pandas as pd
from pathlib import Path
import argparse


# creates the arguments for command line use
parser = argparse.ArgumentParser(description="returns the KO_list from the tsv file in a .txt file.")
parser.add_argument('file_path', metavar='file.tsv path', type=Path, help="The path of the tsv file.")

args = parser.parse_args()

BASE_DIR = Path(args.file_path)


def get_Ko_list():
    """[From file.tsv returns a ko_list_file.txt with comptatible format for KEGG MAPPER Reconstruction]

    Returns:
        [file.txt]: [File containing a list of Ko]
    """
    # creates the output file
    FILE_DIR = BASE_DIR.parent / (BASE_DIR.stem + "_KO_LIST.txt")
    FILE_DIR.touch(exist_ok=True)

    # reads the file
    df = pd.read_csv(BASE_DIR, sep='\t')

    # gets the column containing the KO
    gene_family_column = df["# Gene Family"]    

    # deletes the non-KO rows
    gene_family_column = gene_family_column[~gene_family_column.astype(str).str.startswith('UN')]

    # gets only KO (without the text following)
    gene_family_column = gene_family_column.str.split(":").str.get(0)

    # removes the KO-duplicates
    gene_family_column = gene_family_column.unique()

    Ko_list = gene_family_column.tolist()
    
    with open(FILE_DIR, "a") as f:
        for i,ko in enumerate(Ko_list):
            f.write((f"gene{i}" + "\t" + ko + "\n"))
    
    print(f"There are {len(Ko_list)} KO.")
    return Ko_list
    

get_Ko_list()