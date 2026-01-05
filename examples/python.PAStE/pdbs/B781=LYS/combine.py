from pathlib import Path

pdb_files = sorted(Path('.').glob("*.pdb"))  # 当前目录下的所有 PDB 文件
with open("multi_model.pdb", "w") as outfile:
    for i, pdb_file in enumerate(pdb_files, 1):
        outfile.write(f"MODEL     {i:>4}\n")
        with open(pdb_file) as infile:
            for line in infile:
                if line.startswith(("ATOM", "HETATM", "TER")):
                    outfile.write(line)
        outfile.write("ENDMDL\n")

