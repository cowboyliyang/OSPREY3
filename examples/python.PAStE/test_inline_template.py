import osprey
osprey.start(heapSizeMiB=8000)

# Test inline template like the example
customTemplates = """
Custom template for N6-acetyl-lysine
Modified from LYS with acetyl group on NZ
ACETYL-LYSINE
 ALY  INT     1
 CORR OMIT DU   BEG
   0.00000
   1  DUMM  DU    M    0  -1  -2     0.000     0.000     0.000   0.00000
   2  DUMM  DU    M    1   0  -1     1.449     0.000     0.000   0.00000
   3  DUMM  DU    M    2   1   0     1.522   111.100     0.000   0.00000
   4  N     N     M    3   2   1     1.335   116.600   180.000  -0.4157
   5  H     H     E    4   3   2     1.010   119.800     0.000   0.2719
   6  CA    CT    M    4   3   2     1.449   121.900   180.000  -0.0252
   7  HA    H1    E    6   4   3     1.090   109.500   300.000   0.0843
   8  CB    CT    3    6   4   3     1.525   111.100    60.000  -0.0377
   9  HB2   H1    E    8   6   4     1.090   109.500   300.000   0.0352
  10  HB3   H1    E    8   6   4     1.090   109.500    60.000   0.0352
  11  CG    CT    3    8   6   4     1.525   109.470   180.000   0.0358
  12  HG2   HC    E   11   8   6     1.090   109.500   300.000   0.0103
  13  HG3   HC    E   11   8   6     1.090   109.500    60.000   0.0103
  14  CD    CT    3   11   8   6     1.525   109.470   180.000  -0.0461
  15  HD2   HC    E   14  11   8     1.090   109.500   300.000   0.0621
  16  HD3   HC    E   14  11   8     1.090   109.500    60.000   0.0621
  17  CE    CT    3   14  11   8     1.525   109.470   180.000   0.0275
  18  HE2   HP    E   17  14  11     1.090   109.500   300.000   0.1135
  19  HE3   HP    E   17  14  11     1.090   109.500    60.000   0.1135
  20  NZ    N     3   17  14  11     1.470   109.470   180.000  -0.509
  21  HZ    H     E   20  17  14     1.010   109.470    60.000   0.231
  22  CH    C     3   20  17  14     1.335   120.000   180.000   0.597
  23  OH    O     E   22  20  17     1.229   120.500     0.000  -0.557
  24  CH3   CT    3   22  20  17     1.522   116.600   180.000  -0.182
  25  HH31  HC    E   24  22  20     1.090   109.500    60.000   0.060
  26  HH32  HC    E   24  22  20     1.090   109.500   180.000   0.060
  27  HH33  HC    E   24  22  20     1.090   109.500   300.000   0.060
  28  C     C     M    6   4   3     1.522   111.100   180.000   0.5973
  29  O     O     E   28   6   4     1.229   120.500     0.000  -0.5679

IMPROPER
 -M   CA   N    H
 CA   +M   C    O

DONE
"""

customTemplateCoords = """
ALY 26
N  1.755f  0.276f  -2.802f  -0.4157f  N
H  1.583f  1.269f  -2.842f  0.2719f  H
CA  0.442f  -0.381f  -2.823f  -0.0252f  CT
HA  0.575f  -1.462f  -2.781f  0.0843f  H1
CB  -0.377f  0.078f  -1.615f  -0.0377f  CT
HB2  -0.511f  1.159f  -1.657f  0.0352f  H1
HB3  -1.352f  -0.409f  -1.631f  0.0352f  H1
CG  0.358f  -0.294f  -0.327f  0.0358f  CT
HG2  0.492f  -1.375f  -0.285f  0.0103f  HC
HG3  1.333f  0.192f  -0.312f  0.0103f  HC
CD  -0.461f  0.164f  0.879f  -0.0461f  CT
HD2  -0.594f  1.245f  0.837f  0.0621f  HC
HD3  -1.435f  -0.322f  0.864f  0.0621f  HC
CE  0.275f  -0.208f  2.167f  0.0275f  CT
HE2  0.409f  -1.289f  2.209f  0.1135f  HP
HE3  1.250f  0.278f  2.183f  0.1135f  HP
NZ  -0.509f  0.231f  3.324f  -0.509f  N
HZ  -1.357f  0.683f  3.191f  0.231f  H
CH  -0.053f  0.007f  4.572f  0.597f  C
OH  1.006f  -0.557f  4.738f  -0.557f  O
CH3  -0.861f  0.460f  5.761f  -0.182f  CT
HH31  -0.336f  0.194f  6.679f  0.060f  HC
HH32  -0.994f  1.541f  5.719f  0.060f  HC
HH33  -1.835f  -0.027f  5.746f  0.060f  HC
C  -0.283f  -0.013f  -4.092f  0.5973f  C
O  -0.079f  1.056f  -4.616f  -0.5679f  O
ENDRES
"""

customRotamers = """
1
ALY 5 27
N CA CB CG
CA CB CG CD
CB CG CD CE
CG CD CE NZ
CE NZ CH CH3
62 180 68 180 180
62 180 180 65 180
62 180 180 180 0
62 180 180 -65 180
62 180 -68 180 180
-177 68 180 65 180
-177 68 180 180 0
-177 68 180 -65 180
-177 180 68 65 180
-177 180 68 180 0
-177 180 180 65 180
-177 180 180 180 0
-177 180 180 -65 180
-177 180 -68 180 0
-177 180 -68 -65 180
-90 68 180 180 0
-90 68 180 -65 180
-67 180 -68 -65 180
-62 180 68 65 180
-62 180 68 180 0
-62 180 180 65 180
-62 180 180 180 0
-62 180 180 -65 180
"""

print("Testing inline template...")
try:
    lib = osprey.TemplateLibrary(
        extraTemplates=[customTemplates],
        extraTemplateCoords=[customTemplateCoords],
        extraRotamers=[customRotamers]
    )
    print("✓ SUCCESS! Inline template loaded!")
except Exception as e:
    print(f"❌ Error: {e}")
    import traceback
    traceback.print_exc()
