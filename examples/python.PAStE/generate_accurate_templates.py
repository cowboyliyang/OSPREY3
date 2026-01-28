#!/usr/bin/env python3
"""
基于PDB CIF文件中的真实几何参数生成准确的OSPREY模板文件
"""

# 从CIF文件提取的真实几何参数
ALY_GEOMETRY = {
    'NZ-CE': 1.465,
    'NZ-CH': 1.347,  # 酰胺键 - 比我之前用的1.522短！
    'CH-OH': 1.211,  # C=O双键
    'CH-CH3': 1.507,
    'CE-NZ-CH': 120.0,  # sp2杂化角度，不是109.5°
    'NZ-CH-OH': 120.0,
    'NZ-CH-CH3': 120.0,
}

SLL_GEOMETRY = {
    'NZ-CE': 1.465,
    'NZ-CQ': 1.347,  # 酰胺键（我的命名是CQ，CIF里是CX）
    'CQ-OQ1': 1.213,  # C=O双键（我的命名是OQ1，CIF里是OX）
    'CQ-CS': 1.507,   # 我的命名是CS，CIF里是CK
    'CS-CT': 1.530,
    'CT-CU': 1.507,
    'CU-OU1': 1.208,  # C=O双键
    'CU-OU2': 1.342,  # C-OH单键
    'CE-NZ-CQ': 120.0,
    'NZ-CQ-OQ1': 120.0,
    'NZ-CQ-CS': 120.0,
}

def generate_aly_template():
    """生成更准确的ALY template.in文件"""
    content = """
Custom template for acetylated lysine - ACCURATE VERSION
Based on PDB Chemical Component Dictionary (ALY.cif)
ACETYL-LYSINE
 ALY  INT     1
 CORR OMIT DU   BEG
   0.00000
   1  DUMM  DU    M    0  -1  -2     0.000     0.000     0.000   0.00000
   2  DUMM  DU    M    1   0  -1     1.449     0.000     0.000   0.00000
   3  DUMM  DU    M    2   1   0     1.522   111.100     0.000   0.00000
   4  N     N     M    3   2   1     1.335   116.600   180.000  -0.3479
   5  H     H     E    4   3   2     1.010   119.800     0.000   0.2747
   6  CA    CT    M    4   3   2     1.449   121.900   180.000  -0.2400
   7  HA    H1    E    6   4   3     1.090   109.500   300.000   0.1426
   8  CB    CT    3    6   4   3     1.525   111.100    60.000  -0.0094
   9  HB2   HC    E    8   6   4     1.090   109.500   300.000   0.0362
  10  HB3   HC    E    8   6   4     1.090   109.500    60.000   0.0362
  11  CG    CT    3    8   6   4     1.525   109.470   180.000   0.0187
  12  HG2   HC    E   11   8   6     1.090   109.500   300.000   0.0103
  13  HG3   HC    E   11   8   6     1.090   109.500    60.000   0.0103
  14  CD    CT    3   11   8   6     1.525   109.470   180.000  -0.0479
  15  HD2   HC    E   14  11   8     1.090   109.500   300.000   0.0621
  16  HD3   HC    E   14  11   8     1.090   109.500    60.000   0.0621
  17  CE    CT    3   14  11   8     1.525   109.470   180.000  -0.0143
  18  HE2   HP    E   17  14  11     1.090   109.500   300.000   0.1135
  19  HE3   HP    E   17  14  11     1.090   109.500    60.000   0.1135
  20  NZ    N     3   17  14  11     1.465   109.470   180.000  -0.4157
  21  HZ    H     E   20  17  14     1.010   120.000    60.000   0.2719
  22  CH    C     3   20  17  14     1.347   120.000   180.000   0.5973
  23  OH    O2    E   22  20  17     1.211   120.000     0.000  -0.5679
  24  CH3   CT    3   22  20  17     1.507   120.000   180.000  -0.1825
  25  HH31  HC    E   24  22  20     1.090   109.500    60.000   0.0603
  26  HH32  HC    E   24  22  20     1.090   109.500   180.000   0.0603
  27  HH33  HC    E   24  22  20     1.090   109.500   300.000   0.0603
  28  C     C     M    6   4   3     1.522   111.100   180.000   0.7341
  29  O     O     E   28   6   4     1.229   120.500     0.000  -0.5894

IMPROPER
 -M   CA   N    H
 CA   +M   C    O

DONE
"""
    return content

def generate_sll_template():
    """生成更准确的SLL template.in文件"""
    content = """
Custom template for succinyl-lysine - ACCURATE VERSION
Based on PDB Chemical Component Dictionary (SLL.cif)
SUCCINYL-LYSINE
 SLL  INT     1
 CORR OMIT DU   BEG
   0.00000
   1  DUMM  DU    M    0  -1  -2     0.000     0.000     0.000   0.00000
   2  DUMM  DU    M    1   0  -1     1.449     0.000     0.000   0.00000
   3  DUMM  DU    M    2   1   0     1.522   111.100     0.000   0.00000
   4  N     N     M    3   2   1     1.335   116.600   180.000  -0.3479
   5  H     H     E    4   3   2     1.010   119.800     0.000   0.2747
   6  CA    CT    M    4   3   2     1.449   121.900   180.000  -0.2400
   7  HA    H1    E    6   4   3     1.090   109.500   300.000   0.1426
   8  CB    CT    3    6   4   3     1.525   111.100    60.000  -0.0094
   9  HB2   HC    E    8   6   4     1.090   109.500   300.000   0.0362
  10  HB3   HC    E    8   6   4     1.090   109.500    60.000   0.0362
  11  CG    CT    3    8   6   4     1.525   109.470   180.000   0.0187
  12  HG2   HC    E   11   8   6     1.090   109.500   300.000   0.0103
  13  HG3   HC    E   11   8   6     1.090   109.500    60.000   0.0103
  14  CD    CT    3   11   8   6     1.525   109.470   180.000  -0.0479
  15  HD2   HC    E   14  11   8     1.090   109.500   300.000   0.0621
  16  HD3   HC    E   14  11   8     1.090   109.500    60.000   0.0621
  17  CE    CT    3   14  11   8     1.525   109.470   180.000  -0.0143
  18  HE2   HP    E   17  14  11     1.090   109.500   300.000   0.1135
  19  HE3   HP    E   17  14  11     1.090   109.500    60.000   0.1135
  20  NZ    N     3   17  14  11     1.465   109.470   180.000  -0.4157
  21  HZ    H     E   20  17  14     1.010   120.000    60.000   0.2719
  22  CQ    C     3   20  17  14     1.347   120.000   180.000   0.5973
  23  OQ1   O2    E   22  20  17     1.213   120.000     0.000  -0.5679
  24  CS    CT    3   22  20  17     1.507   120.000   180.000  -0.0094
  25  HS2   HC    E   24  22  20     1.090   109.500    60.000   0.0362
  26  HS3   HC    E   24  22  20     1.090   109.500   180.000   0.0362
  27  CT    CT    3   24  22  20     1.530   111.100   300.000  -0.0094
  28  HT2   HC    E   27  24  22     1.090   109.500    60.000   0.0362
  29  HT3   HC    E   27  24  22     1.090   109.500   180.000   0.0362
  30  CU    C     3   27  24  22     1.507   111.100   300.000   0.7341
  31  OU1   O2    E   30  27  24     1.208   120.500     0.000  -0.5894
  32  OU2   O2    E   30  27  24     1.342   120.500   180.000  -0.5894
  33  C     C     M    6   4   3     1.522   111.100   180.000   0.7341
  34  O     O     E   33   6   4     1.229   120.500     0.000  -0.5894

IMPROPER
 -M   CA   N    H
 CA   +M   C    O

DONE
"""
    return content

# 生成新的模板文件
print("Generating accurate ALY template...")
with open('K03_template_accurate.in', 'w') as f:
    f.write(generate_aly_template())
print("Generated K03_template_accurate.in")

print("\nGenerating accurate SLL template...")
with open('K04_template_accurate.in', 'w') as f:
    f.write(generate_sll_template())
print("Generated K04_template_accurate.in")

print("\n=== Key improvements in accurate templates ===")
print("\nALY (Acetyl-lysine):")
print("  NZ-CH bond: 1.347 Å (was 1.522 Å) - correct amide bond length")
print("  CH-OH bond: 1.211 Å (was ~1.229 Å) - correct C=O double bond")
print("  CE-NZ-CH angle: 120.0° (was 109.5°) - correct sp2 geometry")
print("  NZ-CH-OH angle: 120.0° (was ~120.5°) - correct planar geometry")

print("\nSLL (Succinyl-lysine):")
print("  NZ-CQ bond: 1.347 Å (was 1.522 Å) - correct amide bond length")
print("  CQ-OQ1 bond: 1.213 Å (was ~1.229 Å) - correct C=O double bond")
print("  CE-NZ-CQ angle: 120.0° (was 109.5°) - correct sp2 geometry")
print("  NZ-CQ-CS angle: 120.0° - correct planar geometry")
print("  CU-OU1 bond: 1.208 Å - correct carboxyl C=O")
print("  CU-OU2 bond: 1.342 Å - correct carboxyl C-OH")

print("\nThese templates now use geometry from the PDB Chemical Component Dictionary!")
