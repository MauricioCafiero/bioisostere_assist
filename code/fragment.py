"""RDKit-only fragmentation for bioisostere search: single-bond cuts
(rdMMPA -- the same algorithm mmpdb wraps) producing dummy-tagged
substituent-like pieces, then capping the open valence with a methyl so the
fragment is a complete, embeddable molecule.

No OpenEye/BROOD involved. Unlike a BROOD-native fragment database, the
methyl cap here is permanent -- there's no proprietary format requiring the
dummy atom's identity to survive, so the capped fragment IS the shape
representation used downstream (see CafChemShape.py).
"""
from rdkit import Chem
from rdkit.Chem import rdMMPA

MAX_FRAG_HEAVY_ATOMS = 8  # keep "substituent-like" pieces, not the big remainder


def cap_dummy_with_methyl(mol):
    """Replace every dummy atom (atomic num 0) with carbon. Returns a new,
    sanitized RDKit Mol, or None on failure.
    """
    rw = Chem.RWMol(mol)
    for atom in rw.GetAtoms():
        if atom.GetAtomicNum() == 0:
            atom.SetAtomicNum(6)
            atom.SetNoImplicit(False)
            atom.SetIsotope(0)
    capped = rw.GetMol()
    try:
        Chem.SanitizeMol(capped)
    except Exception:
        return None
    return capped


def fragments_from_smiles(smiles, source_label, max_frag_heavy_atoms=MAX_FRAG_HEAVY_ATOMS):
    """Single-bond-cut fragmentation (rdMMPA, maxCuts=1). Returns a list of
    (fragment_smiles_with_dummy_tag, source_label) for pieces small enough to
    be substituent-like (<= max_frag_heavy_atoms heavy atoms, including the
    dummy) and containing exactly one attachment point.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return []
    cuts = rdMMPA.FragmentMol(mol, maxCuts=1, resultsAsMols=False)

    seen = set()
    out = []
    for core, chains in cuts:
        for piece in chains.split(".") + ([core] if core else []):
            piece = piece.strip()
            if not piece or piece in seen:
                continue
            seen.add(piece)
            pmol = Chem.MolFromSmiles(piece)
            if pmol is None:
                continue
            n_dummies = sum(1 for a in pmol.GetAtoms() if a.GetAtomicNum() == 0)
            n_heavy = pmol.GetNumAtoms()
            if n_dummies == 1 and n_heavy <= max_frag_heavy_atoms:
                out.append((piece, source_label))
    return out


def build_fragment_pool(seed_smiles_with_labels, max_frag_heavy_atoms=MAX_FRAG_HEAVY_ATOMS):
    """seed_smiles_with_labels: iterable of (smiles, label) pairs.
    Returns a deduped list of (fragment_smiles, source_label) -- first
    occurrence's label wins on duplicates.
    """
    pool = []
    for smi, label in seed_smiles_with_labels:
        pool.extend(fragments_from_smiles(smi, label, max_frag_heavy_atoms))
    seen = {}
    for smi, label in pool:
        seen.setdefault(smi, label)
    return list(seen.items())
