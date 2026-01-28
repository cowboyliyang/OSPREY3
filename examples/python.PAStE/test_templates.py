#!/usr/bin/env python3
"""测试acetyl和succinyl lysine模板是否能正确加载"""

import osprey
osprey.start(heapSizeMiB=8000)

print("Testing template library loading...")

try:
    # 测试只加载K03
    print("\n1. Testing K03 (acetyl-lysine) only...")
    lib1 = osprey.TemplateLibrary(
        extraTemplates=['K03_template.in'],
        extraTemplateCoords=['K03_coords.in'],
        extraRotamers=['K03_rotlib.dat']
    )
    print("✓ K03 loaded successfully!")
    
    # 测试只加载K04
    print("\n2. Testing K04 (succinyl-lysine) only...")
    lib2 = osprey.TemplateLibrary(
        extraTemplates=['K04_template.in'],
        extraTemplateCoords=['K04_coords.in'],
        extraRotamers=['K04_rotlib.dat']
    )
    print("✓ K04 loaded successfully!")
    
    # 测试同时加载K03和K04
    print("\n3. Testing K03 + K04 together...")
    lib3 = osprey.TemplateLibrary(
        extraTemplates=['K03_template.in', 'K04_template.in'],
        extraTemplateCoords=['K03_coords.in', 'K04_coords.in'],
        extraRotamers=['K03_rotlib.dat', 'K04_rotlib.dat']
    )
    print("✓ K03 + K04 loaded successfully!")
    
    print("\n✅ All template tests passed!")
    print("\nYou can now run PASTE_6dv2.py with modified lysine residues.")
    
except Exception as e:
    print(f"\n❌ Error: {e}")
    import traceback
    traceback.print_exc()
