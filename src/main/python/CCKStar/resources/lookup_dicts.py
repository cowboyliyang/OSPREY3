# all identities (ignores Ala and Gly)
AllIdentities = ['VAL', 'CYS', 'PRO', 'LEU', 'ILE', 'MET', 'TRP', 'PHE', 'LYS', 'ARG', 'HID', 'HIE', 'HIP',
                 'SER', 'THR', 'TYR', 'ASN', 'GLN', 'ASP', 'GLU']

# conversion between 3 and 1 letter abbreviations
SingleAAtoTripleAA = {'GLY': 'G', 'ALA': 'A', 'ala': 'A', 'VAL': 'V', 'CYS': 'C', 'PRO': 'P', 'LEU': 'L', 'ILE': 'I',
                      'MET': 'M',
                      'TRP': 'W', 'PHE': 'F', 'LYS': 'K', 'ARG': 'R', 'HIS': 'H', 'SER': 'S', 'THR': 'T', 'TYR': 'Y',
                      'ASN': 'N', 'GLN': 'Q', 'ASP': 'D', 'GLU': 'E'}

# atoms of each identity
ResAtoms = {'VAL': '[H, HA, CB, HB, CG1, HG11, HG12, HG13, CG2, HG21, HG22, HG23]',
            'CYS': '[H, HA, CB, HB2, HB3, SG, HG]',
            'PRO': '[CD, HD2, HD3, CG, HG2, HG3, CB, HB2, HB3, HA]',
            'LEU': '[H, HA, CB, HB2, HB3, CG, HG, CD1, HD11, HD12, HD13, CD2, HD21, HD22, HD23]',
            'ILE': '[H, HA, CB, HB, CG2, HG21, HG22, HG23, CG1, HG12, HG13, CD1, HD11, HD12, HD13]',
            'MET': '[H, HA, CB, HB2, HB3, CG, HG2, HG3, SD, CE, HE1, HE2, HE3]',
            'TRP': '[H, HA, CB, HB2, HB3, CG, CD1, HD1, NE1, HE1, CE2, CZ2, HZ2, CH2, HH2, CZ3, HZ3, CE3, HE3, CD2]',
            'PHE': '[H, HA, CB, HB2, HB3, CG, CD1, HD1, CE1, HE1, CZ, HZ, CE2, HE2, CD2, HD2]',
            'LYS': '[H, HA, CB, HB2, HB3, CG, HG2, HG3, CD, HD2, HD3, CE, HE2, HE3, NZ, HZ1, HZ2, HZ3]',
            'ARG': '[H, HA, CB, HB2, HB3, CG, HG2, HG3, CD, HD2, HD3, NE, HE, CZ, NH1, HH11, HH12, NH2, HH21, HH22]',
            'HID': '[H, HA, CB, HB2, HB3, CG, ND1, HD1, CE1, HE1, NE2, CD2, HD2]',
            'HIE': '[H, HA, CB, HB2, HB3, CG, ND1, CE1, HE1, NE2, HE2, CD2, HD2]',
            'HIP': '[H, HA, CB, HB2, HB3, CG, ND1, HD1, CE1, HE1, NE2, HE2, CD2, HD2]',
            'SER': '[H, HA, CB, HB2, HB3, OG, HG]',
            'THR': '[H, HA, CB, HB, CG2, HG21, HG22, HG23, OG1, HG1]',
            'TYR': '[H, HA, CB, HB2, HB3, CG, CD1, HD1, CE1, HE1, CZ, OH, HH, CE2, HE2, CD2, HD2]',
            'ASN': '[H, HA, CB, HB2, HB3, CG, OD1, ND2, HD21, HD22]',
            'GLN': '[H, HA, CB, HB2, HB3, CG, HG2, HG3, CD, OE1, NE2, HE21, HE22]',
            'ASP': '[H, HA, CB, HB2, HB3, CG, OD1, OD2]',
            'GLU': '[H, HA, CB, HB2, HB3, CG, HG2, HG3, CD, OE1, OE2]'}

# element to rotamer name
NameToElement = {'H': 'H', 'HA': 'H', 'CB': 'C', 'HB': 'H', 'CG1': 'C', 'HG11': 'H', 'HG12': 'H', 'HG13': 'H',
                 'CG2': 'C', 'HG21': 'H', 'HG22': 'H', 'HG23': 'H', 'HB2': 'H', 'HB3': 'H', 'OG': 'O', 'HG': 'H',
                 'CG': 'C', 'ND1': 'N', 'HD1': 'H', 'CE1': 'C', 'HE1': 'H', 'NE2': 'N', 'CD2': 'C', 'HD2': 'H',
                 'CD1': 'C', 'HD11': 'H', 'HD12': 'H', 'HD13': 'H', 'HE2': 'H', 'HG2': 'H', 'HG3': 'H', 'CD': 'C',
                 'HD3': 'H', 'CE': 'C', 'HE3': 'H', 'NZ': 'N', 'HZ1': 'H', 'HZ2': 'H', 'HZ3': 'H', 'OE1': 'O',
                 'HE21': 'H', 'HE22': 'H', 'CZ': 'C', 'HZ': 'H', 'CE2': 'C', 'OH': 'O', 'HH': 'H', 'OE2': 'O',
                 'NE1': 'N', 'CZ2': 'C', 'CH2': 'C', 'HH2': 'H', 'CZ3': 'C', 'CE3': 'C', 'HA2': 'H', 'HA3': 'H',
                 'NE': 'N', 'HE': 'H', 'NH1': 'N', 'HH11': 'H', 'HH12': 'H', 'NH2': 'N', 'HH21': 'H', 'HH22': 'H',
                 'HB1': 'H', 'SG': 'S', 'OD1': 'O', 'ND2': 'N', 'HD21': 'H', 'HD22': 'H', 'HD23': 'H', 'SD': 'S',
                 'OD2': 'O', 'OG1': 'O', 'HG1': 'H'}

# chain alphabet. PDB differentiates between upper and lowercase for chain ID.
ChainAlphabet = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T',
                 'U', 'V', 'W', 'X', 'Y', 'Z', 'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n',
                 'o', 'p', 'q', 'r', 's', 't', 'u', 'v', 'w', 'x', 'y', 'z']

# python vars can't have "-", so relation btw lovell name (with -) and class name (no -)
ClassToLovellRotamer = {'p': 'p', 't': 't', 'm': 'm', 'p__60': 'p_-60', 'p_60': 'p_60', 'p_180': 'p_180', 'p_0': 'p_0',
                        'p_120': 'p_120', 'p__120': 'p_-120', 't__60': 't_-60', 't_60': 't_60', 't_180': 't_180',
                        't_0': 't_0', 't_120': 't_120', 't__120': 't_-120', 'm__60': 'm_-60', 'm_60': 'm_60',
                        'm_180': 'm_180', 'm_0': 'm_0', 'm_120': 'm_120', 'm__120': 'm_-120', 'p_80': 'p-80',
                        'p80': 'p80', 't_160': 't-160', 't_80': 't-80', 't60': 't60', 'm_70': 'm-70', 'm170': 'm170',
                        'm80': 'm80', 'pp': 'pp', 'pt': 'pt', 'tp': 'tp', 'tt': 'tt', 'mp': 'mp', 'mt': 'mt',
                        'mm': 'mm', 'ptpt': 'ptpt', 'pttp': 'pttp', 'pttt': 'pttt', 'pttm': 'pttm', 'ptmt': 'ptmt',
                        'tptp': 'tptp', 'tptt': 'tptt', 'tptm': 'tptm', 'ttpp': 'ttpp', 'ttpt': 'ttpt', 'tttp': 'tttp',
                        'tttt': 'tttt', 'tttm': 'tttm', 'ttmt': 'ttmt', 'ttmm': 'ttmm', 'mptt': 'mptt', 'mtpp': 'mtpp',
                        'mtpt': 'mtpt', 'mttp': 'mttp', 'mttt': 'mttt', 'mttm': 'mttm', 'mtmt': 'mtmt', 'mtmm': 'mtmm',
                        'mmtp': 'mmtp', 'mmtt': 'mmtt', 'mmtm': 'mmtm', 'mmmt': 'mmmt', 'pt20': 'pt20', 'pm0': 'pm0',
                        'tp_100': 'tp-100', 'tp60': 'tp60', 'tt0': 'tt0', 'mp0': 'mp0', 'mt_30': 'mt-30',
                        'mm_40': 'mm-40', 'mm100': 'mm100', 'down': 'down', 'up': 'up', 'p90': 'p90', 't80': 't80',
                        'm_85': 'm-85', 'm_30': 'm-30', 'p90_0': 'p90_0', 'p90_180': 'p90_180', 't80_0': 't80_0',
                        't80_180': 't80_180', 'm_85_0': 'm-85_0', 'm_85_180': 'm-85_180', 'm_30_0': 'm-30_0',
                        'm_30_180': 'm-30_180', 'pt_20': 'pt-20', 'tp10': 'tp10', 'tm_20': 'tm-20', 'mt_10': 'mt-10',
                        'p_90': 'p-90', 't_105': 't-105', 't90': 't90', 'm_90': 'm-90', 'm0': 'm0', 'm95': 'm95',
                        'ptp85': 'ptp85', 'ptp180': 'ptp180', 'ptt85': 'ptt85', 'ptt180': 'ptt180', 'ptt_85': 'ptt-85',
                        'ptm180': 'ptm180', 'ptm_85': 'ptm-85', 'tpp85': 'tpp85', 'tpp180': 'tpp180', 'tpt85': 'tpt85',
                        'tpt180': 'tpt180', 'ttp85': 'ttp85', 'ttp180': 'ttp180', 'ttp_105': 'ttp-105',
                        'ttt85': 'ttt85', 'ttt180': 'ttt180', 'ttt_85': 'ttt-85', 'ttm105': 'ttm105',
                        'ttm180': 'ttm180', 'ttm_85': 'ttm-85', 'mtp85': 'mtp85', 'mtp180': 'mtp180',
                        'mtp_105': 'mtp-105', 'mtt85': 'mtt85', 'mtt180': 'mtt180', 'mtt_85': 'mtt-85',
                        'mtm105': 'mtm105', 'mtm180': 'mtm180', 'mtm_85': 'mtm-85', 'mmt85': 'mmt85',
                        'mmt180': 'mmt180', 'mmt_85': 'mmt-85', 'mmm180': 'mmm180', 'mmm_85': 'mmm-85',
                        'p_10': 'p-10', 'p30': 'p30', 't_20': 't-20', 't30': 't30', 'm_20': 'm-20', 'm_80': 'm-80',
                        'm120': 'm120', 'ptp': 'ptp', 'ptm': 'ptm', 'tpp': 'tpp', 'tpt': 'tpt', 'ttp': 'ttp',
                        'ttt': 'ttt', 'ttm': 'ttm', 'mtp': 'mtp', 'mtt': 'mtt', 'mtm': 'mtm', 'mmp': 'mmp',
                        'mmt': 'mmt', 'mmm': 'mmm', 't0': 't0', 't70': 't70'}


# VALp coordinates when translated to origin (useful for hull alignments)
VALp_D_origin = {'CA': [0, 0, 0],
                   'HA': [0.8059999999999974, -0.5640000000000001, 0.4930000000000003],
                   'N': [0.605000000000004, 0.8560000000000016, -1.016],
                   'C': [-0.7420000000000044, 0.8540000000000028, 1.0280000000000005]}

VALp_L_origin = {'CA': [0, 0, 0],
                   'HA': [0.8059999999999974, -0.5640000000000001, -0.4930000000000003],
                   'N': [0.605000000000004, 0.8560000000000016, 1.016],
                   'C': [-0.7420000000000044, 0.8540000000000028, -1.0280000000000005]}
