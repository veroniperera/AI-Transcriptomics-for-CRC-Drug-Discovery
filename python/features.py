import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, rdFingerprintGenerator, rdMolDescriptors

from common import SEED

DESC_2D = (
    ["MolWt", "MolLogP", "MolMR", "TPSA", "NumHDonors", "NumHAcceptors",
     "NumRotatableBonds", "HeavyAtomCount", "RingCount", "FractionCSP3"]
    + ["Chi0", "Chi1", "Chi0v", "Chi1v", "Chi2v", "Chi3v"]
    + ["Kappa1", "Kappa2", "Kappa3", "HallKierAlpha"]
    + ["BalabanJ", "BertzCT", "LabuteASA", "NumAromaticRings", "NumAliphaticRings", "NumHeteroatoms"]
    + [f"EState_VSA{i}" for i in range(1, 10)]
    + [f"PEOE_VSA{i}" for i in range(1, 10)]
)
DESC_3D = ["PMI1", "PMI2", "PMI3", "NPR1", "NPR2", "RadiusOfGyration",
           "InertialShapeFactor", "Eccentricity", "Asphericity", "SpherocityIndex"]
assert len(DESC_2D) == 44 and len(DESC_3D) == 10

_FN = {n: getattr(Descriptors, n) for n in DESC_2D}
_FN3 = {n: getattr(rdMolDescriptors, "Calc" + n) for n in DESC_3D}

_ecfp4 = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
_ecfp6 = rdFingerprintGenerator.GetMorganGenerator(radius=3, fpSize=2048)
_rdk = rdFingerprintGenerator.GetRDKitFPGenerator(maxPath=7, fpSize=2048)


def embed(mol, seed=SEED):
    mh = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = seed
    if AllChem.EmbedMolecule(mh, params) != 0:
        params.useRandomCoords = True
        if AllChem.EmbedMolecule(mh, params) != 0:
            return None
    try:
        AllChem.MMFFOptimizeMolecule(mh, mmffVariant="MMFF94", maxIters=2000)
    except Exception:
        pass
    return mh


def descriptors(mol):
    row = {n: f(mol) for n, f in _FN.items()}
    mh = embed(mol)
    for n, f in _FN3.items():
        try:
            row[n] = f(mh) if mh is not None else np.nan
        except Exception:
            row[n] = np.nan
    return row


def featurize(smiles):
    rows = []
    for s in smiles:
        m = Chem.MolFromSmiles(s)
        rows.append(descriptors(m) if m is not None else {n: np.nan for n in DESC_2D + DESC_3D})
    return pd.DataFrame(rows, columns=DESC_2D + DESC_3D)


def bits(smiles, kind="ecfp4"):
    gen = {"ecfp4": _ecfp4, "ecfp6": _ecfp6, "rdkit": _rdk}[kind]
    out = np.zeros((len(smiles), 2048), dtype=np.uint8)
    for i, s in enumerate(smiles):
        m = Chem.MolFromSmiles(s)
        if m is not None:
            out[i] = gen.GetFingerprintAsNumPy(m)
    return out


def tanimoto(a, b):
    a = a.astype(np.float32)
    b = b.astype(np.float32)
    inter = a @ b.T
    union = a.sum(1)[:, None] + b.sum(1)[None, :] - inter
    return np.divide(inter, union, out=np.zeros_like(inter), where=union > 0)


def ad_score(query_bits, train_bits, k=5, exclude_self=False):
    sim = tanimoto(query_bits, train_bits)
    if exclude_self:
        np.fill_diagonal(sim, -1.0)
    return np.sort(sim, axis=1)[:, -k:].mean(1)
