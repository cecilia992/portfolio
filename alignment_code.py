#modules to load 

import pandas as pd
from io import StringIO
import urllib.parse
import urllib.request
import Bio.PDB
import numpy as np

#Here you can input a UNIPROT
# print('Please input a valid uniprot ID:')
# query = input()
query="P24941"


#download the uniprot mapping file (this is not required and can be uncommented)
urllib.request.urlretrieve("https://ftp.ebi.ac.uk/pub/databases/msd/sifts/flatfiles/csv/pdb_chain_uniprot.csv.gz", "pdb_chain_uniprot_ebi.csv.gz")

#getting the mapping from uniprot to pdb containing the query uniprot
url = 'https://www.uniprot.org/uploadlists/'

params = {
'from': 'ACC+ID',
'to': 'PDB_ID',
'format': 'tab',
'query': str(query)
}
#parse according to the parameters established
data = urllib.parse.urlencode(params)
data = data.encode('utf-8')
req = urllib.request.Request(url, data)
with urllib.request.urlopen(req) as f:
    response = f.read()
result = response.decode('utf-8')
#get the data into a dataframe for easy access
df = pd.read_csv(StringIO(result), sep='\t', header=[0])
df["To"]=df["To"].str.lower()
pdb_list = df["To"].to_list()

#download the pdb files and make a make a list of structures 
structures=[]
for i in range(0, len(df)):
    pdb_code=df.loc[i, "To"].lower()
    pdb_filename=pdb_code+".pdb"
    urllib.request.urlretrieve('http://files.rcsb.org/download/'+str(PDB)+".pdb", str(PDB.lower())+".pdb")
    structure = Bio.PDB.PDBParser().get_structure(pdb_code, pdb_filename)
    structure = structure[0] #select the first model - some pdbs have only one model, some have multiple models (eg NMR structures)
    structures.append(structure)

#open the pdb to uniprot mapping file into a dataframe to get initial position and final position of the protein of interest (query) in the pdb
df2 = pd.read_csv('pdb_chain_uniprot_ebi.csv.gz', compression='gzip', header=[0],low_memory=False)
df2=df2.reset_index()
df2.columns = df2.iloc[0]
df2=df2.drop(0, axis=0)

#select only the rows with the PDBs corresponding to the uniprot of interest 
pdb_uniprot_mapping = df2[df2["PDB"].isin(pdb_list)]
pdb_uniprot_mapping = pdb_uniprot_mapping.reset_index()

#getting the protein sequence from uniprot
urllib.request.urlretrieve("http://www.uniprot.org/uniprot/"+str(query)+".fasta", str(query)+".fasta")
with open(str(query)+".fasta") as file:
    list_lines=[]
    for line in file: 
        line=line.strip()
        list_lines.append(line)
protein_seq = ""
for i in range(1, len(list_lines)):
    protein_seq+=list_lines[i]
print(protein_seq)

#establish the reference PDB
#print("please input a reference PDB ID:")
#reference_PDB=input()
reference_PDB="2vta"
reference_file = reference_PDB+".pdb"
ref_structure = Bio.PDB.PDBParser().get_structure(reference_PDB, reference_file)
ref_model = ref_structure[0]

#iterating through the structures I obtained before 
for i in range (0, len(pdb_uniprot_mapping)):
    sample_pdb = pdb_uniprot_mapping.loc[i, 'PDB']
    sample_model = structures[i]
    #get initial and final position of the uniprot in the pdb file, using the mapping file opened before 
    start_id1 = int(pdb_uniprot_mapping.loc[(pdb_uniprot_mapping['PDB'] == reference_PDB)&(pdb_uniprot_mapping['SP_PRIMARY'] == query), 'PDB_BEG'])
    end_id1 = int(pdb_uniprot_mapping.loc[(pdb_uniprot_mapping['PDB'] == reference_PDB)&(pdb_uniprot_mapping['SP_PRIMARY'] == query), 'PDB_END'])
    atoms_to_be_aligned1 = range(start_id1, end_id1 + 1)
    start_id2 = int(pdb_uniprot_mapping.loc[(pdb_uniprot_mapping['PDB'] == sample_pdb)&(pdb_uniprot_mapping['SP_PRIMARY'] == query), 'PDB_BEG'])
    end_id2 = int(pdb_uniprot_mapping.loc[(pdb_uniprot_mapping['PDB'] == sample_pdb)&(pdb_uniprot_mapping['SP_PRIMARY'] == query), 'PDB_END'])
    atoms_to_be_aligned2 = range(start_id2, end_id2 + 1)
    #Build paired lists of c-alpha atoms:
    ref_atoms = []
    sample_atoms = []
    # Iterate of all chains in the model in order to find all residues
    for ref_chain in ref_model:
        # Iterate of all residues in each model in order to find proper atoms
        for ref_res in ref_chain:
        # Check if residue number ( .get_id() ) is in the list
        if ref_res.get_id()[1] in atoms_to_be_aligned1:
            # Append CA atom to list
            ref_atoms.append(ref_res['CA'])

    # Do the same for the sample structure
    for sample_chain in sample_model:
        for sample_res in sample_chain:
            if sample_res.get_id()[1] in atoms_to_be_aligned2:
                sample_atoms.append(sample_res['CA'])

    #We need to make sure that the two lists of atoms are of the same length; the Superimposer does not work otherwise
    ref_atoms = ref_atoms[:len(sample_atoms)]

    # Now we initiate the superimposer:

    super_imposer = Bio.PDB.Superimposer()
    super_imposer.set_atoms(ref_atoms, sample_atoms)
    super_imposer.apply(sample_model.get_atoms())

    # Print RMSD:
    print(pdb)
    print (super_imposer.rms)

    # Save the aligned version of the sample pdb
    io = Bio.PDB.PDBIO()
    io.set_structure(sample_model) 
    io.save(sample_pdb+"aligned.pdb")


