from pymol import cmd
import numpy as np

def rodrigues_rotate(v, k, theta):
    """使用Rodrigues公式将向量v绕单位向量k旋转theta弧度"""
    v = np.array(v)
    k = np.array(k) / np.linalg.norm(k)
    return (v * np.cos(theta) +
            np.cross(k, v) * np.sin(theta) +
            k * np.dot(k, v) * (1 - np.cos(theta)))

def angle_between(v1, v2):
    v1 = np.array(v1) / np.linalg.norm(v1)
    v2 = np.array(v2) / np.linalg.norm(v2)
    return np.degrees(np.arccos(np.clip(np.dot(v1, v2), -1.0, 1.0)))

def get_ring_centroid(resi_sel):
    """Phe/Tyr/Trp"""
    ring_atoms = ["CG","CD1","CD2","CE1","CE2","CZ","CH2","CZ2","CZ3"]
    coords = []
    for atom in ring_atoms:
        sel = f"({resi_sel}) and name {atom}"
        if cmd.count_atoms(sel) > 0:
            coords.append(cmd.get_atom_coords(sel))
    if len(coords) < 4:
        return None
    return np.mean(coords, axis=0)

def detect_met_aromatic(met_sel="chain A and resi 45", cutoff=6.0):
    """检测单个Met残基与附近芳香残基的Met-aromatic桥"""
    # Get Met
    coords = {}
    for atom in ["CG","SD","CE"]:
        coords[atom] = np.array(cmd.get_atom_coords(f"({met_sel}) and name {atom}"))

    v1 = -(coords["SD"] - coords["CG"])
    v2 = -(coords["SD"] - coords["CE"])
    axis = np.cross(v1, v2)
    a = rodrigues_rotate(v1, axis, np.pi/2)
    g = rodrigues_rotate(v2, axis, np.pi/2)

    nearby = cmd.get_model(f"(byres ({met_sel} around {cutoff})) and resn PHE+TYR+TRP")
    checked = set()
    for atom in nearby.atom:
        aro_id = (atom.chain, atom.resi)
        if aro_id in checked:
            continue
        checked.add(aro_id)
        aro_sel = f"chain {atom.chain} and resi {atom.resi}"

 
        centroid = get_ring_centroid(aro_sel)
        if centroid is None:
            continue
        
        dist = np.linalg.norm(coords["SD"] - centroid)
        if dist > cutoff:
            continue

        aro_vec = centroid - coords["SD"]
        theta = angle_between(a, aro_vec)
        phi   = angle_between(g, aro_vec)

        status = "not satisfy" if (theta >= 109.5 and phi >= 109.5) else "satisfy"
        print(f"Met {met_sel} → {aro_sel} | distance {dist:.2f} Å | θ={theta:.1f}° φ={phi:.1f}° | {status}")


cmd.extend("detect_met_aromatic", detect_met_aromatic)

