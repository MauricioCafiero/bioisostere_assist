# bioisostere_assist

A CafChem project.

Finds shape-bioisosteric substituent fragments — pieces of a molecule that
could be swapped for something else while preserving 3D shape (and
optionally pharmacophore/electrostatic character) at the attachment point.
Built entirely on open-source tools: **RDKit for fragmentation and
conformer generation, ODDT for the 3D shape descriptors.**

## How it works

1. **Fragment** (`code/fragment.py`) — cuts each input molecule at single
   bonds with RDKit's `rdMMPA` (the same algorithm the `mmpdb` package
   wraps, used directly here with no extra dependency), keeping
   substituent-like pieces (≤ 8 heavy atoms by default, configurable).
2. **Cap** — each fragment's open valence (a dummy `*` atom) is permanently
   capped with a methyl carbon, turning it into a complete, embeddable
   molecule.
3. **Embed** (`code/CafChemShape.py`) — up to 10 conformers per fragment
   (RDKit ETKDGv3), MMFF94-optimized, near-duplicates pruned by RMSD during
   embedding, and any conformer more than 10 kcal/mol above the lowest-energy
   one discarded as too strained. The *rest* are kept — not just the single
   best — since the conformer that overlays best against a query usually
   isn't the global energy minimum.
4. **Shape descriptor** — each surviving conformer gets an ODDT shape
   descriptor. **Swappable via `--method`:**
   - `usr` — Ultrafast Shape Recognition: pure 3D shape, fastest.
   - `usr_cat` — USR extended with pharmacophoric atom-type constraints
     (hydrophobic/aromatic/acceptor/donor atoms compared separately).
   - `electroshape` — shape + chirality + partial-charge (electrostatic)
     information.
5. **Compare** — similarity between two fragments is the **max over every
   (conformer of A, conformer of B) pair**, i.e. "do these two fragments
   have *any* pair of conformations that overlay well."

## Library storage and the "new lead" workflow

`build_library_cli.py` is the expensive step (run once, or whenever you want
to rebuild): fragmenting the whole 100K-compound pool and embedding+shaping
every resulting fragment took **4 minutes 16 seconds** end to end (4,994
unique fragments, 4,981 successfully shaped). It saves **two** files from
`--out outputs/library_100k_usr`:

- `library_100k_usr.csv` — a human-readable manifest (fragment SMILES,
  source compound, conformer count).
- `library_100k_usr.pkl` — the actual library: every surviving conformer's
  shape-descriptor vector, keyed by fragment SMILES, plus the build
  parameters used (method, conformer count, energy window, RMSD threshold).
  **This is what makes it a reusable library rather than a one-off report**
  — the 4-minute computation is captured here and never needs to be redone.

When a new lead molecule comes in, `query_lead_cli.py` loads that pickle
(instant), fragments *only the lead* (fast — one molecule), embeds+shapes
*only the lead's own fragments* using the same parameters the library was
built with (read back out of the pickle automatically, overridable via
flags), and ranks every lead fragment against the full library. Querying
rosuvastatin (10 substituent fragments) against the full 100K-derived
library (4,981 fragments) took **under 2 seconds** total — the whole point
of separating "build" from "query" is that the second step never re-touches
the first step's work.

```bash
# once (or whenever the pool changes):
python3 build_library_cli.py --method usr --out ../outputs/library_100k_usr

# per new lead, as many times as you like:
python3 query_lead_cli.py --library ../outputs/library_100k_usr.pkl \
    --lead "CC(C)c1nc(N(C)S(C)(=O)=O)nc(-c2ccc(F)cc2)c1/C=C/[C@H](O)C[C@H](O)CC(=O)O" \
    --top 8
```

A couple of its 10 fragments, and the top shape-similar hits found in the
library (real output from the run above):

```
Lead fragment: CC(C)[*:1]                  (the isopropyl group on the pyrimidine ring)
  1.000  CC(C)[*:1]           (from ZINC000100783622)
  0.975  CN(C)[*:1]           (from ZINC000070666542)
  0.962  N[C@H](CO)[*:1]      (from ZINC000005273845)
  0.954  CC(C)(Cl)[*:1]       (from ZINC000014490904)
  0.951  C1CC1[*:1]           (from ZINC000065259628)

Lead fragment: Fc1ccc([*:1])cc1            (the 4-fluorophenyl ring)
  1.000  Fc1ccc([*:1])cc1     (from ZINC000021944692)
  0.943  Fc1ccc([*:1])nc1     (from ZINC000049155361)
  0.932  Fc1cccc([*:1])c1     (from ZINC000514289382)

Lead fragment: CN(S(C)(=O)=O)[*:1]         (the N-methyl-N-mesyl group)
  1.000  CN(C)S(=O)(=O)[*:1]  (from ZINC000263585619)
```

## Data

`data/human_druglike_100k.*` — 100,000 unique compounds pulled from ZINC's
`special/current` subsets (`in-man`, `world`, `investigational` — ~79.7K,
the full drug/clinical-status set — topped up with a ~20.2K random sample
of `biogenic` natural products to reach 100K without drowning the drug
compounds in natural-product diversity). All SMILES-valid.

- `human_druglike_100k.smi` — SMILES + ZINC ID, one per line.
- `human_druglike_100k.info.tsv` — same compounds with their full ZINC tag list.
- `human_druglike_100k.scaffolds.tsv` / `.scaffold_freq.tsv` — Bemis-Murcko
  scaffold decomposition of the same pool.

## Installation

`oddt` is unmaintained (last released 2019) and needs an older-NumPy-era
build path, so this uses Python 3.11 specifically and a specific install
order:

```bash
uv venv --python 3.11 bio-env
source bio-env/bin/activate

uv pip install -r requirements-base.txt
uv pip install --no-build-isolation -r requirements-oddt.txt
uv pip install -r requirements.txt
```

`pip install` works the same way in place of `uv pip install` if you're not
using `uv`. Install in this order: `oddt` must be built without isolation
after the base build tools are present, or it can pull an incompatible NumPy
version and fail to compile.

**`oddt` is vendored, not pulled from PyPI** (`vendor/oddt-0.7.tar.gz`, a
verified byte-identical copy of the official PyPI sdist) so installation
keeps working even if the PyPI release or upstream GitHub repo disappears.
A fork of `oddt/oddt` is also maintained at
[github.com/MauricioCafiero/oddt](https://github.com/MauricioCafiero/oddt)
as a second, patchable copy.

`CafChemShape.py` also shims `np.in1d` onto `np.isin` before importing
oddt: oddt calls `np.in1d`, which NumPy 2.x removed.

## Usage

### `build_library_cli.py` — build (or rebuild) the library

```bash
cd code
python3 build_library_cli.py --csv ../data/human_druglike_100k.smi \
    --method usr --out ../outputs/library_100k_usr
```

| Flag | Default | Description |
|---|---|---|
| `--csv` | `data/human_druglike_100k.smi` | Whitespace-separated `SMILES [ID]` file, one molecule per line |
| `--n-mols` | all | Limit to the first N seed molecules |
| `--method` | `usr` | `usr`, `usr_cat`, or `electroshape` |
| `--max-frag-atoms` | 8 | Max heavy atoms (incl. dummy) for a fragment to be kept |
| `--n-confs` | 10 | Conformers embedded per fragment before pruning |
| `--energy-window` | 10.0 | kcal/mol above the lowest-energy conformer before a conformer is discarded as strained |
| `--prune-rms-thresh` | 0.5 | Å; conformers this close to one already kept are treated as duplicates |
| `--out` | `outputs/library` | Output path stem — writes `<out>.csv` and `<out>.pkl` |

### `query_lead_cli.py` — rank a new lead's fragments against a built library

```bash
python3 query_lead_cli.py --library ../outputs/library_100k_usr.pkl \
    --lead "CC(C)c1nc(N(C)S(C)(=O)=O)nc(-c2ccc(F)cc2)c1/C=C/[C@H](O)C[C@H](O)CC(=O)O" \
    --top 10 --out ../outputs/lead_matches.csv
```

| Flag | Default | Description |
|---|---|---|
| `--library` | — (required) | Path to a `.pkl` from `build_library_cli.py` |
| `--lead` | — (required) | SMILES of the new lead molecule |
| `--max-frag-atoms`, `--n-confs`, `--energy-window`, `--prune-rms-thresh` | match library | Override the library's own build parameters if needed |
| `--top` | 15 | Top matches per lead fragment to report |
| `--out` | `outputs/lead_matches.csv` | Output CSV path |
