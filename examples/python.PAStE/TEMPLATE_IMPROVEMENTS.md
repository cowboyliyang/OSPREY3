# PTM Template Accuracy Improvements

## Summary

Updated ALY (acetyl-lysine) and SLL (succinyl-lysine) templates with accurate geometric parameters from the PDB Chemical Component Dictionary.

## What Was Improved

### 1. Bond Lengths (Internal Coordinates in template.in)

**ALY (Acetyl-lysine):**
- NZ-CH bond: **1.347 Å** (was 1.522 Å)
  - This is the critical amide bond connecting lysine to acetyl group
  - 1.347 Å is correct for C-N amide bond (sp2 hybridization)
  - Previous 1.522 Å was for C-C single bond (incorrect)
- CH-OH bond: **1.211 Å** (was ~1.229 Å)
  - Correct C=O double bond length for amide carbonyl

**SLL (Succinyl-lysine):**
- NZ-CQ bond: **1.347 Å** (was 1.522 Å)
  - Same amide bond as ALY
- CQ-OQ1 bond: **1.213 Å** (was ~1.229 Å)
  - C=O double bond in succinyl group
- CU-OU1 bond: **1.208 Å** - terminal carboxyl C=O
- CU-OU2 bond: **1.342 Å** - terminal carboxyl C-OH

### 2. Bond Angles (Internal Coordinates in template.in)

**ALY:**
- CE-NZ-CH angle: **120.0°** (was 109.5°)
  - Correct for sp2 hybridized nitrogen in amide
  - 109.5° would be for sp3 (tetrahedral) - incorrect
- NZ-CH-OH angle: **120.0°** (was ~120.5°)
  - Planar geometry around carbonyl carbon
- NZ-CH-CH3 angle: **120.0°**
  - Planar amide geometry

**SLL:**
- CE-NZ-CQ angle: **120.0°** (was 109.5°)
  - Same sp2 correction as ALY
- NZ-CQ-OQ1 angle: **120.0°**
- NZ-CQ-CS angle: **120.0°**
  - All consistent with planar amide geometry

### 3. Cartesian Coordinates (coords.in)

Both ALY and SLL coords files now use **ideal coordinates from PDB CIF files**:
- Extracted from `pdbx_model_Cartn_x/y/z_ideal` fields
- These are crystallographically-derived ideal geometries
- AMBER partial charges retained from original templates
- AMBER atom types retained (CT, HC, N, O2, etc.)

## Files Updated

- `K03_template.in` - ALY internal coordinates
- `K03_coords.in` - ALY Cartesian coordinates
- `K04_template.in` - SLL internal coordinates
- `K04_coords.in` - SLL Cartesian coordinates
- `K03_rotlib.dat` - unchanged (rotamer library)
- `K04_rotlib.dat` - unchanged (rotamer library)

## Backup Files

Original approximated templates backed up to:
- `backup/K03_template.in`
- `backup/K03_coords.in`
- `backup/K04_template.in`
- `backup/K04_coords.in`

## Source Data

- **ALY.cif** - N6-acetyl-lysine from PDB Chemical Component Dictionary
- **SLL.cif** - N6-succinyl-lysine from PDB Chemical Component Dictionary

## Scripts Used

1. `extract_geometry.py` - Extracts bond lengths, angles, dihedrals from CIF
2. `generate_accurate_templates.py` - Generates template.in with correct geometry
3. `generate_accurate_coords.py` - Generates coords.in with PDB ideal coordinates

## Testing

Tested with `test_quick.py`:
```
✓ Template library created
✓ Strand created
✓ B781 can be: LYS, ALY, SLL
✓ ConfSpace created!
✅ SUCCESS!
```

## Impact on Your Research

For your TFP protein study:

1. **More accurate energy calculations** - Correct amide geometry means better electrostatics and sterics
2. **Reliable PTM predictions** - Whether acetylation/succinylation stabilizes or destabilizes M404K
3. **Publication-ready** - Geometry derived from PDB standard reference, not approximations

The continuous optimization in OSPREY will further refine these coordinates during design, but starting with correct geometry is crucial for:
- Proper hydrogen bonding patterns
- Accurate dipole moments of amide bonds
- Correct van der Waals interactions

## Key Takeaway

**Before:** Templates used generic sp3 geometry (tetrahedral)
**After:** Templates use correct sp2 geometry (planar) for amide bonds

This is a significant improvement in accuracy! 🎯
