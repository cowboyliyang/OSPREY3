from resources.lookup_dicts import *


class VALRotamer:
    def __init__(self, name, H, HA, CB, HB, CG1, HG11, HG12, HG13, CG2, HG21, HG22, HG23, CA_back, N_back, C_back):
        self.name = name
        self.H = H
        self.HA = HA
        self.CB = CB
        self.HB = HB
        self.CG1 = CG1
        self.HG11 = HG11
        self.HG12 = HG12
        self.HG13 = HG13
        self.CG2 = CG2
        self.HG21 = HG21
        self.HG22 = HG22
        self.HG23 = HG23
        self.CA_back = CA_back
        self.N_back = N_back
        self.C_back = C_back

    def print_pdb(self, atom_num, chain, res_num, file):
        pdb_info = [""] * 9
        pdb_info[0] = "ATOM".ljust(6)
        pdb_info[1] = str(atom_num).rjust(5)
        pdb_info[3] = "VAL".ljust(3)
        pdb_info[4] = chain.rjust(1)
        pdb_info[5] = str(res_num).rjust(4)

        for atom in vars(self):
            xyz = getattr(self, atom)
            if isinstance(xyz, list):
                if atom == "H":
                    continue
                if atom in ["CA_back", "N_back", "C_back"]:
                    small_name = atom[:-5]
                    pdb_info[2] = small_name.center(4)
                else:
                    pdb_info[2] = atom.center(4)
                pdb_info[6] = str('%8.3f' % (xyz[0])).rjust(8)
                pdb_info[7] = str('%8.3f' % (xyz[1])).rjust(8)
                pdb_info[8] = str('%8.3f' % (xyz[2])).rjust(8)
                file.write("%s%s %s %s %s%s    %s%s%s\n" % (pdb_info[0], pdb_info[1], pdb_info[2], pdb_info[3],
                                                            pdb_info[4], pdb_info[5], pdb_info[6], pdb_info[7],
                                                            pdb_info[8]))
                atom_num += 1
                pdb_info[1] = str(atom_num).rjust(5)
        file.write("TER\n")

class CYSRotamer:
    def __init__(self, name, H, HA, CB, HB2, HB3, SG, HG, CA_back, N_back, C_back):
        self.name = name
        self.H = H
        self.HA = HA
        self.CB = CB
        self.HB2 = HB2
        self.HB3 = HB3
        self.SG = SG
        self.HG = HG
        self.CA_back = CA_back
        self.N_back = N_back
        self.C_back = C_back

    def print_pdb(self, atom_num, chain, res_num, file):
        pdb_info = [""] * 9
        pdb_info[0] = "ATOM".ljust(6)
        pdb_info[1] = str(atom_num).rjust(5)
        pdb_info[3] = "CYS".ljust(3)
        pdb_info[4] = chain.rjust(1)
        pdb_info[5] = str(res_num).rjust(4)

        for atom in vars(self):
            xyz = getattr(self, atom)
            if isinstance(xyz, list):
                if atom == "H":
                    continue
                if atom in ["CA_back", "N_back", "C_back"]:
                    small_name = atom[:-5]
                    pdb_info[2] = small_name.center(4)
                else:
                    pdb_info[2] = atom.center(4)
                pdb_info[6] = str('%8.3f' % (xyz[0])).rjust(8)
                pdb_info[7] = str('%8.3f' % (xyz[1])).rjust(8)
                pdb_info[8] = str('%8.3f' % (xyz[2])).rjust(8)
                file.write("%s%s %s %s %s%s    %s%s%s\n" % (pdb_info[0], pdb_info[1], pdb_info[2], pdb_info[3],
                                                            pdb_info[4], pdb_info[5], pdb_info[6], pdb_info[7],
                                                            pdb_info[8]))
                atom_num += 1
                pdb_info[1] = str(atom_num).rjust(5)
        file.write("TER\n")

class PRORotamer:
    def __init__(self, name, CD, HD2, HD3, CG, HG2, HG3, CB, HB2, HB3, HA, CA_back, N_back, C_back):
        self.name = name
        self.CD = CD
        self.HD2 = HD2
        self.HD3 = HD3
        self.CG = CG
        self.HG2 = HG2
        self.HG3 = HG3
        self.CB = CB
        self.HB2 = HB2
        self.HB3 = HB3
        self.HA = HA
        self.CA_back = CA_back
        self.N_back = N_back
        self.C_back = C_back

    def print_pdb(self, atom_num, chain, res_num, file):
        pdb_info = [""] * 9
        pdb_info[0] = "ATOM".ljust(6)
        pdb_info[1] = str(atom_num).rjust(5)
        pdb_info[3] = "PRO".ljust(3)
        pdb_info[4] = chain.rjust(1)
        pdb_info[5] = str(res_num).rjust(4)

        for atom in vars(self):
            xyz = getattr(self, atom)
            if isinstance(xyz, list):
                if atom == "H":
                    continue
                if atom in ["CA_back", "N_back", "C_back"]:
                    small_name = atom[:-5]
                    pdb_info[2] = small_name.center(4)
                else:
                    pdb_info[2] = atom.center(4)
                pdb_info[6] = str('%8.3f' % (xyz[0])).rjust(8)
                pdb_info[7] = str('%8.3f' % (xyz[1])).rjust(8)
                pdb_info[8] = str('%8.3f' % (xyz[2])).rjust(8)
                file.write("%s%s %s %s %s%s    %s%s%s\n" % (pdb_info[0], pdb_info[1], pdb_info[2], pdb_info[3],
                                                            pdb_info[4], pdb_info[5], pdb_info[6], pdb_info[7],
                                                            pdb_info[8]))
                atom_num += 1
                pdb_info[1] = str(atom_num).rjust(5)
        file.write("TER\n")

class LEURotamer:
    def __init__(self, name, H, HA, CB, HB2, HB3, CG, HG, CD1, HD11, HD12, HD13, CD2, HD21, HD22, HD23, CA_back,
                 N_back, C_back):
        self.name = name
        self.H = H
        self.HA = HA
        self.CB = CB
        self.HB2 = HB2
        self.HB3 = HB3
        self.CG = CG
        self.HG = HG
        self.CD1 = CD1
        self.HD11 = HD11
        self.HD12 = HD12
        self.HD13 = HD13
        self.CD2 = CD2
        self.HD21 = HD21
        self.HD22 = HD22
        self.HD23 = HD23
        self.CA_back = CA_back
        self.N_back = N_back
        self.C_back = C_back

    def print_pdb(self, atom_num, chain, res_num, file):
        pdb_info = [""] * 9
        pdb_info[0] = "ATOM".ljust(6)
        pdb_info[1] = str(atom_num).rjust(5)
        pdb_info[3] = "LEU".ljust(3)
        pdb_info[4] = chain.rjust(1)
        pdb_info[5] = str(res_num).rjust(4)

        for atom in vars(self):
            xyz = getattr(self, atom)
            if isinstance(xyz, list):
                if atom == "H":
                    continue
                if atom in ["CA_back", "N_back", "C_back"]:
                    small_name = atom[:-5]
                    pdb_info[2] = small_name.center(4)
                else:
                    pdb_info[2] = atom.center(4)
                pdb_info[6] = str('%8.3f' % (xyz[0])).rjust(8)
                pdb_info[7] = str('%8.3f' % (xyz[1])).rjust(8)
                pdb_info[8] = str('%8.3f' % (xyz[2])).rjust(8)
                file.write("%s%s %s %s %s%s    %s%s%s\n" % (pdb_info[0], pdb_info[1], pdb_info[2], pdb_info[3],
                                                            pdb_info[4], pdb_info[5], pdb_info[6], pdb_info[7],
                                                            pdb_info[8]))
                atom_num += 1
                pdb_info[1] = str(atom_num).rjust(5)
        file.write("TER\n")

class ILERotamer:
    def __init__(self, name, H, HA, CB, HB, CG2, HG21, HG22, HG23, CG1, HG12, HG13, CD1, HD11, HD12, HD13,
                 CA_back, N_back, C_back):
        self.name = name
        self.H = H
        self.HA = HA
        self.CB = CB
        self.HB = HB
        self.CG2 = CG2
        self.HG21 = HG21
        self.HG22 = HG22
        self.HG23 = HG23
        self.CG1 = CG1
        self.HG12 = HG12
        self.HG13 = HG13
        self.CD1 = CD1
        self.HD11 = HD11
        self.HD12 = HD12
        self.HD13 = HD13
        self.CA_back = CA_back
        self.N_back = N_back
        self.C_back = C_back

    def print_pdb(self, atom_num, chain, res_num, file):
        pdb_info = [""] * 9
        pdb_info[0] = "ATOM".ljust(6)
        pdb_info[1] = str(atom_num).rjust(5)
        pdb_info[3] = "ILE".ljust(3)
        pdb_info[4] = chain.rjust(1)
        pdb_info[5] = str(res_num).rjust(4)

        for atom in vars(self):
            xyz = getattr(self, atom)
            if isinstance(xyz, list):
                if atom == "H":
                    continue
                if atom in ["CA_back", "N_back", "C_back"]:
                    small_name = atom[:-5]
                    pdb_info[2] = small_name.center(4)
                else:
                    pdb_info[2] = atom.center(4)
                pdb_info[6] = str('%8.3f' % (xyz[0])).rjust(8)
                pdb_info[7] = str('%8.3f' % (xyz[1])).rjust(8)
                pdb_info[8] = str('%8.3f' % (xyz[2])).rjust(8)
                file.write("%s%s %s %s %s%s    %s%s%s\n" % (pdb_info[0], pdb_info[1], pdb_info[2], pdb_info[3],
                                                            pdb_info[4], pdb_info[5], pdb_info[6], pdb_info[7],
                                                            pdb_info[8]))
                atom_num += 1
                pdb_info[1] = str(atom_num).rjust(5)
        file.write("TER\n")

class METRotamer:
    def __init__(self, name, H, HA, CB, HB2, HB3, CG, HG2, HG3, SD, CE, HE1, HE2, HE3, CA_back, N_back, C_back):
        self.name = name
        self.H = H
        self.HA = HA
        self.CB = CB
        self.HB2 = HB2
        self.HB3 = HB3
        self.CG = CG
        self.HG2 = HG2
        self.HG3 = HG3
        self.SD = SD
        self.CE = CE
        self.HE1 = HE1
        self.HE2 = HE2
        self.HE3 = HE3
        self.CA_back = CA_back
        self.N_back = N_back
        self.C_back = C_back

    def print_pdb(self, atom_num, chain, res_num, file):
        pdb_info = [""] * 9
        pdb_info[0] = "ATOM".ljust(6)
        pdb_info[1] = str(atom_num).rjust(5)
        pdb_info[3] = "MET".ljust(3)
        pdb_info[4] = chain.rjust(1)
        pdb_info[5] = str(res_num).rjust(4)

        for atom in vars(self):
            xyz = getattr(self, atom)
            if isinstance(xyz, list):
                if atom == "H":
                    continue
                if atom in ["CA_back", "N_back", "C_back"]:
                    small_name = atom[:-5]
                    pdb_info[2] = small_name.center(4)
                else:
                    pdb_info[2] = atom.center(4)
                pdb_info[6] = str('%8.3f' % (xyz[0])).rjust(8)
                pdb_info[7] = str('%8.3f' % (xyz[1])).rjust(8)
                pdb_info[8] = str('%8.3f' % (xyz[2])).rjust(8)
                file.write("%s%s %s %s %s%s    %s%s%s\n" % (pdb_info[0], pdb_info[1], pdb_info[2], pdb_info[3],
                                                            pdb_info[4], pdb_info[5], pdb_info[6], pdb_info[7],
                                                            pdb_info[8]))
                atom_num += 1
                pdb_info[1] = str(atom_num).rjust(5)
        file.write("TER\n")

class TRPRotamer:
    def __init__(self, name, H, HA, CB, HB2, HB3, CG, CD1, HD1, NE1, HE1, CE2, CZ2, HZ2, CH2, HH2, CZ3, HZ3, CE3,
                 HE3, CD2, CA_back, N_back, C_back):
        self.name = name
        self.H = H
        self.HA = HA
        self.CB = CB
        self.HB2 = HB2
        self.HB3 = HB3
        self.CG = CG
        self.CD1 = CD1
        self.HD1 = HD1
        self.NE1 = NE1
        self.HE1 = HE1
        self.CE2 = CE2
        self.CZ2 = CZ2
        self.HZ2 = HZ2
        self.CH2 = CH2
        self.HH2 = HH2
        self.CZ3 = CZ3
        self.HZ3 = HZ3
        self.CE3 = CE3
        self.HE3 = HE3
        self.CD2 = CD2
        self.CA_back = CA_back
        self.N_back = N_back
        self.C_back = C_back

    def print_pdb(self, atom_num, chain, res_num, file):
        pdb_info = [""] * 9
        pdb_info[0] = "ATOM".ljust(6)
        pdb_info[1] = str(atom_num).rjust(5)
        pdb_info[3] = "TRP".ljust(3)
        pdb_info[4] = chain.rjust(1)
        pdb_info[5] = str(res_num).rjust(4)

        for atom in vars(self):
            xyz = getattr(self, atom)
            if isinstance(xyz, list):
                if atom == "H":
                    continue
                if atom in ["CA_back", "N_back", "C_back"]:
                    small_name = atom[:-5]
                    pdb_info[2] = small_name.center(4)
                else:
                    pdb_info[2] = atom.center(4)
                pdb_info[6] = str('%8.3f' % (xyz[0])).rjust(8)
                pdb_info[7] = str('%8.3f' % (xyz[1])).rjust(8)
                pdb_info[8] = str('%8.3f' % (xyz[2])).rjust(8)
                file.write("%s%s %s %s %s%s    %s%s%s\n" % (pdb_info[0], pdb_info[1], pdb_info[2], pdb_info[3],
                                                            pdb_info[4], pdb_info[5], pdb_info[6], pdb_info[7],
                                                            pdb_info[8]))
                atom_num += 1
                pdb_info[1] = str(atom_num).rjust(5)
        file.write("TER\n")

class PHERotamer:
    def __init__(self, name, H, HA, CB, HB2, HB3, CG, CD1, HD1, CE1, HE1, CZ, HZ, CE2, HE2, CD2, HD2,
                 CA_back, N_back, C_back):
        self.name = name
        self.H = H
        self.HA = HA
        self.CB = CB
        self.HB2 = HB2
        self.HB3 = HB3
        self.CG = CG
        self.CD1 = CD1
        self.HD1 = HD1
        self.CE1 = CE1
        self.HE1 = HE1
        self.CZ = CZ
        self.HZ = HZ
        self.CE2 = CE2
        self.HE2 = HE2
        self.CD2 = CD2
        self.HD2 = HD2
        self.CA_back = CA_back
        self.N_back = N_back
        self.C_back = C_back

    def print_pdb(self, atom_num, chain, res_num, file):
        pdb_info = [""] * 9
        pdb_info[0] = "ATOM".ljust(6)
        pdb_info[1] = str(atom_num).rjust(5)
        pdb_info[3] = "PHE".ljust(3)
        pdb_info[4] = chain.rjust(1)
        pdb_info[5] = str(res_num).rjust(4)

        for atom in vars(self):
            xyz = getattr(self, atom)
            if isinstance(xyz, list):
                if atom == "H":
                    continue
                if atom in ["CA_back", "N_back", "C_back"]:
                    small_name = atom[:-5]
                    pdb_info[2] = small_name.center(4)
                else:
                    pdb_info[2] = atom.center(4)
                pdb_info[6] = str('%8.3f' % (xyz[0])).rjust(8)
                pdb_info[7] = str('%8.3f' % (xyz[1])).rjust(8)
                pdb_info[8] = str('%8.3f' % (xyz[2])).rjust(8)
                file.write("%s%s %s %s %s%s    %s%s%s\n" % (pdb_info[0], pdb_info[1], pdb_info[2], pdb_info[3],
                                                            pdb_info[4], pdb_info[5], pdb_info[6], pdb_info[7],
                                                            pdb_info[8]))
                atom_num += 1
                pdb_info[1] = str(atom_num).rjust(5)
        file.write("TER\n")

class LYSRotamer:
    def __init__(self, name, H, HA, CB, HB2, HB3, CG, HG2, HG3, CD, HD2, HD3, CE, HE2, HE3, NZ, HZ1, HZ2, HZ3,
                 CA_back, N_back, C_back):
        self.name = name
        self.H = H
        self.HA = HA
        self.CB = CB
        self.HB2 = HB2
        self.HB3 = HB3
        self.CG = CG
        self.HG2 = HG2
        self.HG3 = HG3
        self.CD = CD
        self.HD2 = HD2
        self.HD3 = HD3
        self.CE = CE
        self.HE2 = HE2
        self.HE3 = HE3
        self.NZ = NZ
        self.HZ1 = HZ1
        self.HZ2 = HZ2
        self.HZ3 = HZ3
        self.CA_back = CA_back
        self.N_back = N_back
        self.C_back = C_back

    def print_pdb(self, atom_num, chain, res_num, file):
        pdb_info = [""] * 9
        pdb_info[0] = "ATOM".ljust(6)
        pdb_info[1] = str(atom_num).rjust(5)
        pdb_info[3] = "LYS".ljust(3)
        pdb_info[4] = chain.rjust(1)
        pdb_info[5] = str(res_num).rjust(4)

        for atom in vars(self):
            xyz = getattr(self, atom)
            if isinstance(xyz, list):
                if atom == "H":
                    continue
                if atom in ["CA_back", "N_back", "C_back"]:
                    small_name = atom[:-5]
                    pdb_info[2] = small_name.center(4)
                else:
                    pdb_info[2] = atom.center(4)
                pdb_info[6] = str('%8.3f' % (xyz[0])).rjust(8)
                pdb_info[7] = str('%8.3f' % (xyz[1])).rjust(8)
                pdb_info[8] = str('%8.3f' % (xyz[2])).rjust(8)
                file.write("%s%s %s %s %s%s    %s%s%s\n" % (pdb_info[0], pdb_info[1], pdb_info[2], pdb_info[3],
                                                            pdb_info[4], pdb_info[5], pdb_info[6], pdb_info[7],
                                                            pdb_info[8]))
                atom_num += 1
                pdb_info[1] = str(atom_num).rjust(5)
        file.write("TER\n")

class ARGRotamer:
    def __init__(self, name, H, HA, CB, HB2, HB3, CG, HG2, HG3, CD, HD2, HD3, NE, HE, CZ, NH1, HH11, HH12, NH2,
                 HH21, HH22, CA_back, N_back, C_back):
        self.name = name
        self.H = H
        self.HA = HA
        self.CB = CB
        self.HB2 = HB2
        self.HB3 = HB3
        self.CG = CG
        self.HG2 = HG2
        self.HG3 = HG3
        self.CD = CD
        self.HD2 = HD2
        self.HD3 = HD3
        self.NE = NE
        self.HE = HE
        self.CZ = CZ
        self.NH1 = NH1
        self.HH11 = HH11
        self.HH12 = HH12
        self.NH2 = NH2
        self.HH21 = HH21
        self.HH22 = HH22
        self.CA_back = CA_back
        self.N_back = N_back
        self.C_back = C_back

    def print_pdb(self, atom_num, chain, res_num, file):
        pdb_info = [""] * 9
        pdb_info[0] = "ATOM".ljust(6)
        pdb_info[1] = str(atom_num).rjust(5)
        pdb_info[3] = "ARG".ljust(3)
        pdb_info[4] = chain.rjust(1)
        pdb_info[5] = str(res_num).rjust(4)

        for atom in vars(self):
            xyz = getattr(self, atom)
            if isinstance(xyz, list):
                if atom == "H":
                    continue
                if atom in ["CA_back", "N_back", "C_back"]:
                    small_name = atom[:-5]
                    pdb_info[2] = small_name.center(4)
                else:
                    pdb_info[2] = atom.center(4)
                pdb_info[6] = str('%8.3f' % (xyz[0])).rjust(8)
                pdb_info[7] = str('%8.3f' % (xyz[1])).rjust(8)
                pdb_info[8] = str('%8.3f' % (xyz[2])).rjust(8)
                file.write("%s%s %s %s %s%s    %s%s%s\n" % (pdb_info[0], pdb_info[1], pdb_info[2], pdb_info[3],
                                                            pdb_info[4], pdb_info[5], pdb_info[6], pdb_info[7],
                                                            pdb_info[8]))
                atom_num += 1
                pdb_info[1] = str(atom_num).rjust(5)
        file.write("TER\n")

class HIDRotamer:
    def __init__(self, name, H, HA, CB, HB2, HB3, CG, ND1, HD1, CE1, HE1, NE2, CD2, HD2, CA_back, N_back, C_back):
        self.name = name
        self.H = H
        self.HA = HA
        self.CB = CB
        self.HB2 = HB2
        self.HB3 = HB3
        self.CG = CG
        self.ND1 = ND1
        self.HD1 = HD1
        self.CE1 = CE1
        self.HE1 = HE1
        self.NE2 = NE2
        self.CD2 = CD2
        self.HD2 = HD2
        self.CA_back = CA_back
        self.N_back = N_back
        self.C_back = C_back

    def print_pdb(self, atom_num, chain, res_num, file):
        pdb_info = [""] * 9
        pdb_info[0] = "ATOM".ljust(6)
        pdb_info[1] = str(atom_num).rjust(5)
        pdb_info[3] = "HID".ljust(3)
        pdb_info[4] = chain.rjust(1)
        pdb_info[5] = str(res_num).rjust(4)

        for atom in vars(self):
            xyz = getattr(self, atom)
            if isinstance(xyz, list):
                if atom == "H":
                    continue
                if atom in ["CA_back", "N_back", "C_back"]:
                    small_name = atom[:-5]
                    pdb_info[2] = small_name.center(4)
                else:
                    pdb_info[2] = atom.center(4)
                pdb_info[6] = str('%8.3f' % (xyz[0])).rjust(8)
                pdb_info[7] = str('%8.3f' % (xyz[1])).rjust(8)
                pdb_info[8] = str('%8.3f' % (xyz[2])).rjust(8)
                file.write("%s%s %s %s %s%s    %s%s%s\n" % (pdb_info[0], pdb_info[1], pdb_info[2], pdb_info[3],
                                                            pdb_info[4], pdb_info[5], pdb_info[6], pdb_info[7],
                                                            pdb_info[8]))
                atom_num += 1
                pdb_info[1] = str(atom_num).rjust(5)
        file.write("TER\n")

class HIERotamer:
    def __init__(self, name, H, HA, CB, HB2, HB3, CG, ND1, CE1, HE1, NE2, HE2, CD2, HD2, CA_back, N_back, C_back):
        self.name = name
        self.H = H
        self.HA = HA
        self.CB = CB
        self.HB2 = HB2
        self.HB3 = HB3
        self.CG = CG
        self.ND1 = ND1
        self.CE1 = CE1
        self.HE1 = HE1
        self.NE2 = NE2
        self.HE2 = HE2
        self.CD2 = CD2
        self.HD2 = HD2
        self.CA_back = CA_back
        self.N_back = N_back
        self.C_back = C_back

    def print_pdb(self, atom_num, chain, res_num, file):
        pdb_info = [""] * 9
        pdb_info[0] = "ATOM".ljust(6)
        pdb_info[1] = str(atom_num).rjust(5)
        pdb_info[3] = "HIE".ljust(3)
        pdb_info[4] = chain.rjust(1)
        pdb_info[5] = str(res_num).rjust(4)

        for atom in vars(self):
            xyz = getattr(self, atom)
            if isinstance(xyz, list):
                if atom == "H":
                    continue
                if atom in ["CA_back", "N_back", "C_back"]:
                    small_name = atom[:-5]
                    pdb_info[2] = small_name.center(4)
                else:
                    pdb_info[2] = atom.center(4)
                pdb_info[6] = str('%8.3f' % (xyz[0])).rjust(8)
                pdb_info[7] = str('%8.3f' % (xyz[1])).rjust(8)
                pdb_info[8] = str('%8.3f' % (xyz[2])).rjust(8)
                file.write("%s%s %s %s %s%s    %s%s%s\n" % (pdb_info[0], pdb_info[1], pdb_info[2], pdb_info[3],
                                                            pdb_info[4], pdb_info[5], pdb_info[6], pdb_info[7],
                                                            pdb_info[8]))
                atom_num += 1
                pdb_info[1] = str(atom_num).rjust(5)
        file.write("TER\n")

class HIPRotamer:
    def __init__(self, name, H, HA, CB, HB2, HB3, CG, ND1, HD1, CE1, HE1, NE2, HE2, CD2, HD2, CA_back, N_back, C_back):
        self.name = name
        self.H = H
        self.HA = HA
        self.CB = CB
        self.HB2 = HB2
        self.HB3 = HB3
        self.CG = CG
        self.ND1 = ND1
        self.HD1 = HD1
        self.CE1 = CE1
        self.HE1 = HE1
        self.NE2 = NE2
        self.HE2 = HE2
        self.CD2 = CD2
        self.HD2 = HD2
        self.CA_back = CA_back
        self.N_back = N_back
        self.C_back = C_back

    def print_pdb(self, atom_num, chain, res_num, file):
        pdb_info = [""] * 9
        pdb_info[0] = "ATOM".ljust(6)
        pdb_info[1] = str(atom_num).rjust(5)
        pdb_info[3] = "HIP".ljust(3)
        pdb_info[4] = chain.rjust(1)
        pdb_info[5] = str(res_num).rjust(4)

        for atom in vars(self):
            xyz = getattr(self, atom)
            if isinstance(xyz, list):
                if atom == "H":
                    continue
                if atom in ["CA_back", "N_back", "C_back"]:
                    small_name = atom[:-5]
                    pdb_info[2] = small_name.center(4)
                else:
                    pdb_info[2] = atom.center(4)
                pdb_info[6] = str('%8.3f' % (xyz[0])).rjust(8)
                pdb_info[7] = str('%8.3f' % (xyz[1])).rjust(8)
                pdb_info[8] = str('%8.3f' % (xyz[2])).rjust(8)
                file.write("%s%s %s %s %s%s    %s%s%s\n" % (pdb_info[0], pdb_info[1], pdb_info[2], pdb_info[3],
                                                            pdb_info[4], pdb_info[5], pdb_info[6], pdb_info[7],
                                                            pdb_info[8]))
                atom_num += 1
                pdb_info[1] = str(atom_num).rjust(5)
        file.write("TER\n")

class SERRotamer:
    def __init__(self, name, H, HA, CB, HB2, HB3, OG, HG, CA_back, N_back, C_back):
        self.name = name
        self.H = H
        self.HA = HA
        self.CB = CB
        self.HB2 = HB2
        self.HB3 = HB3
        self.OG = OG
        self.HG = HG
        self.CA_back = CA_back
        self.N_back = N_back
        self.C_back = C_back

    def print_pdb(self, atom_num, chain, res_num, file):
        pdb_info = [""] * 9
        pdb_info[0] = "ATOM".ljust(6)
        pdb_info[1] = str(atom_num).rjust(5)
        pdb_info[3] = "SER".ljust(3)
        pdb_info[4] = chain.rjust(1)
        pdb_info[5] = str(res_num).rjust(4)

        for atom in vars(self):
            xyz = getattr(self, atom)
            if isinstance(xyz, list):
                if atom == "H":
                    continue
                if atom in ["CA_back", "N_back", "C_back"]:
                    small_name = atom[:-5]
                    pdb_info[2] = small_name.center(4)
                else:
                    pdb_info[2] = atom.center(4)
                pdb_info[6] = str('%8.3f' % (xyz[0])).rjust(8)
                pdb_info[7] = str('%8.3f' % (xyz[1])).rjust(8)
                pdb_info[8] = str('%8.3f' % (xyz[2])).rjust(8)
                file.write("%s%s %s %s %s%s    %s%s%s\n" % (pdb_info[0], pdb_info[1], pdb_info[2], pdb_info[3],
                                                            pdb_info[4], pdb_info[5], pdb_info[6], pdb_info[7],
                                                            pdb_info[8]))
                atom_num += 1
                pdb_info[1] = str(atom_num).rjust(5)
        file.write("TER\n")

class THRRotamer:
    def __init__(self, name, H, HA, CB, HB, CG2, HG21, HG22, HG23, OG1, HG1, CA_back, N_back, C_back):
        self.name = name
        self.H = H
        self.HA = HA
        self.CB = CB
        self.HB = HB
        self.CG2 = CG2
        self.HG21 = HG21
        self.HG22 = HG22
        self.HG23 = HG23
        self.OG1 = OG1
        self.HG1 = HG1
        self.CA_back = CA_back
        self.N_back = N_back
        self.C_back = C_back

    def print_pdb(self, atom_num, chain, res_num, file):
        pdb_info = [""] * 9
        pdb_info[0] = "ATOM".ljust(6)
        pdb_info[1] = str(atom_num).rjust(5)
        pdb_info[3] = "THR".ljust(3)
        pdb_info[4] = chain.rjust(1)
        pdb_info[5] = str(res_num).rjust(4)

        for atom in vars(self):
            xyz = getattr(self, atom)
            if isinstance(xyz, list):
                if atom == "H":
                    continue
                if atom in ["CA_back", "N_back", "C_back"]:
                    small_name = atom[:-5]
                    pdb_info[2] = small_name.center(4)
                else:
                    pdb_info[2] = atom.center(4)
                pdb_info[6] = str('%8.3f' % (xyz[0])).rjust(8)
                pdb_info[7] = str('%8.3f' % (xyz[1])).rjust(8)
                pdb_info[8] = str('%8.3f' % (xyz[2])).rjust(8)
                file.write("%s%s %s %s %s%s    %s%s%s\n" % (pdb_info[0], pdb_info[1], pdb_info[2], pdb_info[3],
                                                            pdb_info[4], pdb_info[5], pdb_info[6], pdb_info[7],
                                                            pdb_info[8]))
                atom_num += 1
                pdb_info[1] = str(atom_num).rjust(5)
        file.write("TER\n")

class TYRRotamer:
    def __init__(self, name, H, HA, CB, HB2, HB3, CG, CD1, HD1, CE1, HE1, CZ, OH, HH, CE2, HE2, CD2, HD2,
                 CA_back, N_back, C_back):
        self.name = name
        self.H = H
        self.HA = HA
        self.CB = CB
        self.HB2 = HB2
        self.HB3 = HB3
        self.CG = CG
        self.CD1 = CD1
        self.HD1 = HD1
        self.CE1 = CE1
        self.HE1 = HE1
        self.CZ = CZ
        self.OH = OH
        self.HH = HH
        self.CE2 = CE2
        self.HE2 = HE2
        self.CD2 = CD2
        self.HD2 = HD2
        self.CA_back = CA_back
        self.N_back = N_back
        self.C_back = C_back

    def print_pdb(self, atom_num, chain, res_num, file):
        pdb_info = [""] * 9
        pdb_info[0] = "ATOM".ljust(6)
        pdb_info[1] = str(atom_num).rjust(5)
        pdb_info[3] = "TYR".ljust(3)
        pdb_info[4] = chain.rjust(1)
        pdb_info[5] = str(res_num).rjust(4)

        for atom in vars(self):
            xyz = getattr(self, atom)
            if isinstance(xyz, list):
                if atom == "H":
                    continue
                if atom in ["CA_back", "N_back", "C_back"]:
                    small_name = atom[:-5]
                    pdb_info[2] = small_name.center(4)
                else:
                    pdb_info[2] = atom.center(4)
                pdb_info[6] = str('%8.3f' % (xyz[0])).rjust(8)
                pdb_info[7] = str('%8.3f' % (xyz[1])).rjust(8)
                pdb_info[8] = str('%8.3f' % (xyz[2])).rjust(8)
                file.write("%s%s %s %s %s%s    %s%s%s\n" % (pdb_info[0], pdb_info[1], pdb_info[2], pdb_info[3],
                                                            pdb_info[4], pdb_info[5], pdb_info[6], pdb_info[7],
                                                            pdb_info[8]))
                atom_num += 1
                pdb_info[1] = str(atom_num).rjust(5)
        file.write("TER\n")

class ASNRotamer:
    def __init__(self, name, H, HA, CB, HB2, HB3, CG, OD1, ND2, HD21, HD22, CA_back, N_back, C_back):
        self.name = name
        self.H = H
        self.HA = HA
        self.CB = CB
        self.HB2 = HB2
        self.HB3 = HB3
        self.CG = CG
        self.OD1 = OD1
        self.ND2 = ND2
        self.HD21 = HD21
        self.HD22 = HD22
        self.CA_back = CA_back
        self.N_back = N_back
        self.C_back = C_back

    def print_pdb(self, atom_num, chain, res_num, file):
        pdb_info = [""] * 9
        pdb_info[0] = "ATOM".ljust(6)
        pdb_info[1] = str(atom_num).rjust(5)
        pdb_info[3] = "ASN".ljust(3)
        pdb_info[4] = chain.rjust(1)
        pdb_info[5] = str(res_num).rjust(4)

        for atom in vars(self):
            xyz = getattr(self, atom)
            if isinstance(xyz, list):
                if atom == "H":
                    continue
                if atom in ["CA_back", "N_back", "C_back"]:
                    small_name = atom[:-5]
                    pdb_info[2] = small_name.center(4)
                else:
                    pdb_info[2] = atom.center(4)
                pdb_info[6] = str('%8.3f' % (xyz[0])).rjust(8)
                pdb_info[7] = str('%8.3f' % (xyz[1])).rjust(8)
                pdb_info[8] = str('%8.3f' % (xyz[2])).rjust(8)
                file.write("%s%s %s %s %s%s    %s%s%s\n" % (pdb_info[0], pdb_info[1], pdb_info[2], pdb_info[3],
                                                            pdb_info[4], pdb_info[5], pdb_info[6], pdb_info[7],
                                                            pdb_info[8]))
                atom_num += 1
                pdb_info[1] = str(atom_num).rjust(5)
        file.write("TER\n")

class GLNRotamer:
    def __init__(self, name, H, HA, CB, HB2, HB3, CG, HG2, HG3, CD, OE1, NE2, HE21, HE22, CA_back, N_back, C_back):
        self.name = name
        self.H = H
        self.HA = HA
        self.CB = CB
        self.HB2 = HB2
        self.HB3 = HB3
        self.CG = CG
        self.HG2 = HG2
        self.HG3 = HG3
        self.CD = CD
        self.OE1 = OE1
        self.NE2 = NE2
        self.HE21 = HE21
        self.HE22 = HE22
        self.CA_back = CA_back
        self.N_back = N_back
        self.C_back = C_back

    def print_pdb(self, atom_num, chain, res_num, file):
        pdb_info = [""] * 9
        pdb_info[0] = "ATOM".ljust(6)
        pdb_info[1] = str(atom_num).rjust(5)
        pdb_info[3] = "GLN".ljust(3)
        pdb_info[4] = chain.rjust(1)
        pdb_info[5] = str(res_num).rjust(4)

        for atom in vars(self):
            xyz = getattr(self, atom)
            if isinstance(xyz, list):
                if atom == "H":
                    continue
                if atom in ["CA_back", "N_back", "C_back"]:
                    small_name = atom[:-5]
                    pdb_info[2] = small_name.center(4)
                else:
                    pdb_info[2] = atom.center(4)
                pdb_info[6] = str('%8.3f' % (xyz[0])).rjust(8)
                pdb_info[7] = str('%8.3f' % (xyz[1])).rjust(8)
                pdb_info[8] = str('%8.3f' % (xyz[2])).rjust(8)
                file.write("%s%s %s %s %s%s    %s%s%s\n" % (pdb_info[0], pdb_info[1], pdb_info[2], pdb_info[3],
                                                            pdb_info[4], pdb_info[5], pdb_info[6], pdb_info[7],
                                                            pdb_info[8]))
                atom_num += 1
                pdb_info[1] = str(atom_num).rjust(5)
        file.write("TER\n")

class ASPRotamer:
    def __init__(self, name, H, HA, CB, HB2, HB3, CG, OD1, OD2, CA_back, N_back, C_back):
        self.name = name
        self.H = H
        self.HA = HA
        self.CB = CB
        self.HB2 = HB2
        self.HB3 = HB3
        self.CG = CG
        self.OD1 = OD1
        self.OD2 = OD2
        self.CA_back = CA_back
        self.N_back = N_back
        self.C_back = C_back

    def print_pdb(self, atom_num, chain, res_num, file):
        pdb_info = [""] * 9
        pdb_info[0] = "ATOM".ljust(6)
        pdb_info[1] = str(atom_num).rjust(5)
        pdb_info[3] = "ASP".ljust(3)
        pdb_info[4] = chain.rjust(1)
        pdb_info[5] = str(res_num).rjust(4)

        for atom in vars(self):
            xyz = getattr(self, atom)
            if isinstance(xyz, list):
                if atom == "H":
                    continue
                if atom in ["CA_back", "N_back", "C_back"]:
                    small_name = atom[:-5]
                    pdb_info[2] = small_name.center(4)
                else:
                    pdb_info[2] = atom.center(4)
                pdb_info[6] = str('%8.3f' % (xyz[0])).rjust(8)
                pdb_info[7] = str('%8.3f' % (xyz[1])).rjust(8)
                pdb_info[8] = str('%8.3f' % (xyz[2])).rjust(8)
                file.write("%s%s %s %s %s%s    %s%s%s\n" % (pdb_info[0], pdb_info[1], pdb_info[2], pdb_info[3],
                                                            pdb_info[4], pdb_info[5], pdb_info[6], pdb_info[7],
                                                            pdb_info[8]))
                atom_num += 1
                pdb_info[1] = str(atom_num).rjust(5)
        file.write("TER\n")

class GLURotamer:
    def __init__(self, name, H, HA, CB, HB2, HB3, CG, HG2, HG3, CD, OE1, OE2, CA_back, N_back, C_back):
        self.name = name
        self.H = H
        self.HA = HA
        self.CB = CB
        self.HB2 = HB2
        self.HB3 = HB3
        self.CG = CG
        self.HG2 = HG2
        self.HG3 = HG3
        self.CD = CD
        self.OE1 = OE1
        self.OE2 = OE2
        self.CA_back = CA_back
        self.N_back = N_back
        self.C_back = C_back

    def print_pdb(self, atom_num, chain, res_num, file):
        pdb_info = [""] * 9
        pdb_info[0] = "ATOM".ljust(6)
        pdb_info[1] = str(atom_num).rjust(5)
        pdb_info[3] = "GLU".ljust(3)
        pdb_info[4] = chain.rjust(1)
        pdb_info[5] = str(res_num).rjust(4)

        for atom in vars(self):
            xyz = getattr(self, atom)
            if isinstance(xyz, list):
                if atom == "H":
                    continue
                if atom in ["CA_back", "N_back", "C_back"]:
                    small_name = atom[:-5]
                    pdb_info[2] = small_name.center(4)
                else:
                    pdb_info[2] = atom.center(4)
                pdb_info[6] = str('%8.3f' % (xyz[0])).rjust(8)
                pdb_info[7] = str('%8.3f' % (xyz[1])).rjust(8)
                pdb_info[8] = str('%8.3f' % (xyz[2])).rjust(8)
                file.write("%s%s %s %s %s%s    %s%s%s\n" % (pdb_info[0], pdb_info[1], pdb_info[2], pdb_info[3],
                                                            pdb_info[4], pdb_info[5], pdb_info[6], pdb_info[7],
                                                            pdb_info[8]))
                atom_num += 1
                pdb_info[1] = str(atom_num).rjust(5)
        file.write("TER\n")
