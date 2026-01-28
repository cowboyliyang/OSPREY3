import osprey
osprey.start(heapSizeMiB=8000)

print("Step 1: Creating template library...")
customizedTemplateLib = osprey.TemplateLibrary(
    extraTemplates=['K03_template.in', 'K04_template.in'],
    extraTemplateCoords=['K03_coords.in', 'K04_coords.in'],
    extraRotamers=['K03_rotlib.dat', 'K04_rotlib.dat']
)
print("✓ Template library created")

print("\nStep 2: Creating Strand...")
protein = osprey.Strand('6dv2_strip_reduce_prep_rc_add_rc.pdb', 
                       templateLib=customizedTemplateLib, 
                       residues=['B734', 'B782'])
print("✓ Strand created")

print("\nStep 3: Setting flexibility...")
protein.flexibility['B781'].setLibraryRotamers(osprey.WILD_TYPE, 'ALY', 'SLL').addWildTypeRotamers().setContinuous()
print("✓ B781 can be: LYS, ALY, SLL")

print("\nStep 4: Creating ConfSpace...")
confSpace = osprey.ConfSpace([protein])
print(f"✓ ConfSpace created!")
print(f"✅ SUCCESS! Script works correctly.")
