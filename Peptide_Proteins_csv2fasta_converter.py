"""
Script for converting excel/csv peptide sequences-protein ids file to .fasta file.


Requires structure on input file:
    - .xlsx or .csv
    - Name of column containing sequences is "Sequence" and is the 2nd column of table.
    - Name of column containing the Description is "# PSMs" and is the 3rd column of the table.
    - Each sequence row is followed by header row, containing "Description"
        column, even if no protein ids have been found for the sequence.
    - Assumes protein id is longer than 20 characters.
    - Each protein id begins with an ">" character in the Description column.
    - There is a recommendation, though no mandatory, to have lines in fasta not longer than 80 characters. For the current use case, the maximal sequence length <50, therefore line clipping to 80 characters is not implemented.


Requires the file to be located in the same directory as the script. A single filename
    ending in .csv or .xlsx has to be found in the directory to avoid ambiguities.


Output files are written to the same directory. These are the following:
    - summary.txt: summarizes the conversion
    - peptide_protein.fasta: converted fasta file.

Note: There can (and will) be duplicated entries in the sequences' headers in the fasta file, since a protein can have multiple sequences in the input that belong only to itself. Unique entries is not a requirement of the fasta format, but is important to be acknowledged for further processing of the fasta file.
"""


from pathlib import Path
import pandas as pd


FILENAME_FASTA = "peptide_protein.fasta"
FILENAME_SUMMARY = "summary.txt"
COL_SEQUENCE = "Sequence"
COL_DESCRIPTION = "# PSMs"


def check_directory():
    """
    Checks directory for .xlsx and .csv Path().glob("*.csv").name.
    If more than one is identified, raises OSError.
    """
    
    num_csv = len(list(Path().glob("*.csv")))
    num_xlsx = len(list(Path().glob("*.xlsx")))
    
    if (num_csv + num_xlsx) == 1:
        return (num_csv, num_xlsx)
    else:
        raise OSError(f"Expected a single .xlsx"
                       " OR .csv file in the directory."
                      f" Found {num_csv} .csv and "
                      f"{num_xlsx} .xlsx files.")
        

def get_filename():
    num_csv, num_xlsx = check_directory()
    if num_csv == 1:
        filename = str(list(Path().glob("*.csv"))[0])
    else:
        filename = str(list(Path().glob("*.xlsx"))[0])
    return filename


def open_file(path: str) -> pd.DataFrame:
    """
    Reads file into pd.DataFrame.
    Requires input file to be excel (xlsx) or csv (csv).
    """
    
    extension = path.split(".")[-1]
    if extension == "xlsx":
        reader = pd.read_excel
    elif extension == "csv":
        reader = pd.read_csv
    else:
        raise OSError("Invalid extension or filetype."
                      " Please make sure the file format"
                      " is excel (ending with .xlsx) or"
                      " csv (ending with .csv).")
    
    df = reader(Path(path))
    return df


def parse_dataframe(df: pd.DataFrame) -> dict:
    """
    Parses the dataframe into dictionary with
     keys as peptides and values as protein names.
    Used subsequently to write to fasta file.
    """
    
    nrows = df.shape[0]
    pept_prot = {}
    for row_tuple in df.iterrows():
        
        i, row_data = row_tuple
        
        if row_data[COL_DESCRIPTION] == "Description":
            peptide_seq = df.loc[i-1, COL_SEQUENCE]
            protein_ids = []
            
            for j in range(i+1, nrows):
                value = df.loc[j, COL_DESCRIPTION]
                
                if (type(value) != str) or (type(value) == str and value[0] != ">"):
                    break
                else:
                    # Remove an unexpected whitespace at the end;
                    value = value.split(" ")[0]
                    assert len(value) > 20, "Protein name assumed to have > 20 chars."
                    protein_ids.append(value)
            
            pept_prot[peptide_seq] = protein_ids
    
    return pept_prot
            

def summarize_dict(pept_prot_dict: dict, filename: str):
    """
    Writes to a text file a summary. For validation and debugging.
    """
    
    n_seq = len(pept_prot_dict)
    n_prot_ids = sum([len(val) for val in pept_prot_dict.values()])
    seqs_with_no_ids = [ seq for seq in pept_prot_dict.keys()
                         if pept_prot_dict[seq]==[] ]
    seq_lengths = [len(seq) for seq in pept_prot_dict.keys()]
    
    line1 = f"Number of sequences: {n_seq}\n"
    line2 = f"Number of protein ids extracted: {n_prot_ids}\n"
    line3 = f"Max sequence length: {max(seq_lengths)}\n"
    line4 = f"Number of sequences with no proteins ids identified: {len(seqs_with_no_ids)}\n"
    line5 = f"Sequences with no protein ids: {seqs_with_no_ids}\n"
    
    lines = [line1, line2, line3, line4, line5]
    with open(Path(filename), "w+") as f:
        f.writelines(lines)
        
    
def write_to_fasta(pept_prot_dict: dict, filename: str):
    """
    Writes the peptide-protein dict to fasta.
    """
    path = Path(filename)
    with open(path, "w+") as file:
        
        for seq, ids in pept_prot_dict.items():

            ids_edited = []
            for id_ in ids:
                if id_[0] == ">":
                    id_str = id_[1:]
                else:
                    id_str = id_
                ids_edited.append(id_str)
            seq_header = "; ".join(ids_edited)
            
            file.write("> " + seq_header + "\n")
            file.write(seq + "\n")
     
            
if __name__ == "__main__":
    
    print("Start.")
    
    filename = get_filename()
    print(f"Identified file: {filename}\nConverting.")
    
    df = open_file(filename)
    peptide_protein_dict = parse_dataframe(df)
    
    print(f"Writing summary to: {FILENAME_SUMMARY}")
    summarize_dict(peptide_protein_dict, FILENAME_SUMMARY)
    
    print(f"Writing fast to: {FILENAME_FASTA}")
    write_to_fasta(peptide_protein_dict, FILENAME_FASTA)
    
    print("Finished.")