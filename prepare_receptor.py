"""Instructor step (needs internet): download PDB 6OIM, split protein and sotorasib (ligand MOV),
write data/receptor.pdbqt, data/crystal_ligand.pdbqt and fill data/box.yaml."""
import subprocess, sys, urllib.request, yaml
import numpy as np

urllib.request.urlretrieve("https://files.rcsb.org/download/6OIM.pdb", "data/6OIM.pdb")
protein, ligand, coords = [], [], []
for line in open("data/6OIM.pdb"):
    if line.startswith("ATOM") and line[21] == "A":
        protein.append(line)
    elif line.startswith("HETATM") and line[17:20].strip() == "MOV" and line[21] == "A":
        ligand.append(line)
        coords.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
open("data/protein.pdb", "w").writelines(protein + ["END\n"])
open("data/crystal_ligand.pdb", "w").writelines(ligand + ["END\n"])
c = np.array(coords)
yaml.safe_dump({"center": [round(float(x), 2) for x in c.mean(0)],
                "size": [round(float(x), 1) for x in (c.max(0) - c.min(0) + 8)]},
               open("data/box.yaml", "w"))
# receptor: Meeko's receptor preparation (protonates and writes PDBQT)
subprocess.run([sys.executable, "-m", "meeko.cli.mk_prepare_receptor", "--read_pdb", "data/protein.pdb",
                "-o", "data/receptor", "-p", "-a"], check=True)
# crystal ligand as PDBQT for RMSD (template from SMILES so bond orders are right)
subprocess.run(["mk_prepare_ligand.py", "-i", "data/crystal_ligand.pdb", "-o", "data/crystal_ligand.pdbqt"], check=True)
print("done:", open("data/box.yaml").read())
