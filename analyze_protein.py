from Bio.PDB import PDBParser, is_aa
from Bio.SeqUtils import seq1
from collections import Counter
import numpy as np

pdb_file = "/home/rajat/Documents/Bio_python_Exp/Protein_Exp/protein_sample/AF-C0NA21-F1-model_v6.pdb"

parser = PDBParser(QUIET=True)
structure = parser.get_structure("Protein", pdb_file)

amino_acids_3 = []
sequence_1 = ""
atom_coords = []
chain_ids = set()

print("\n🔬 AUTOMATED PROTEIN ANALYSIS\n")

for model in structure:
    for chain in model:
        chain_ids.add(chain.id)

        for residue in chain:
            # Automatically check if it's an amino acid
            if is_aa(residue, standard=True):
                resname = residue.get_resname()
                amino_acids_3.append(resname)

                # Convert automatically to 1-letter
                sequence_1 += seq1(resname)

                for atom in residue:
                    atom_coords.append(atom.coord)

# ---- Stats ----
residue_count = len(amino_acids_3)
unique_aas = set(amino_acids_3)
freq = Counter(amino_acids_3)

print("Chains:", len(chain_ids))
print("Total Residues:", residue_count)
print("Unique Amino Acids:", len(unique_aas))
print("Amino Acids:", unique_aas)

# ---- Sequence ----
print("\n🧬 Sequence (first 100):")
print(sequence_1[:100])

# ---- Frequency ----
print("\n📊 Amino Acid Frequency:")
for aa, count in freq.most_common():
    print(f"{aa}: {count}")

# ---- Geometry ----
coords = np.array(atom_coords)

print("\n📐 Geometry:")
print("Center of Mass:", coords.mean(axis=0))
print("Min coords:", coords.min(axis=0))
print("Max coords:", coords.max(axis=0))
print("Size (Å):", coords.max(axis=0) - coords.min(axis=0))