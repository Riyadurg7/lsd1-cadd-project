# LSD1 CADD Mini Project

Computer-aided drug design study of LSD1 (KDM1A, PDB 6TUY): ChEMBL IC50 data
curation, descriptor/fragment analysis, similarity screening, and docking with
AutoDock Vina and AutoDock4.

## Layout

| Folder | Contents |
| --- | --- |
| `notebooks/` | Colab notebooks (`LSD1_FULL_Project.ipynb` is the complete pipeline) and teammates' notebooks |
| `docking/` | Colab working directory: curated datasets, figures/tables, receptor prep, Vina ligands and results |
| `enamine/` | Enamine similarity search: reference compounds, hits, selected compounds |
| `ad4/` | AutoDock4 run: grid/docking parameter files, `.dlg` logs, complexes, PLIP interactions |
| `dockingflow/` | Course AutoDock4 example (2XKN); source of the `autodock4`/`autogrid4` binaries |
| `scripts/` | Local helper scripts (paths resolved relative to the project root) |
| `report/` | Project report and the report template |

## Pipeline order

1. `notebooks/LSD1_FULL_Project.ipynb` writes the dataset, figures and Vina docking files into `docking/`
2. `scripts/enamine_similarity.py`, then `scripts/select_and_prepare_hits.py`, `make_minimal_set.py`, `fix_macrocycle_ligand.py`
3. `docking/run_docking.py` runs Vina on the minimal set
4. `scripts/setup_autodock4.py`, `fix_dpf_about.py`, then `ad4/run_ad4.sh` (in WSL)
5. `scripts/parse_ad4_results.py`, `extract_poses.py`, `run_plip.py`, `check_fad_contacts.py`

## Not tracked in git

These are kept locally but listed in `.gitignore`:

- `enamine/Enamine_SMILES.csv` (209 MB), `enamine/Enamine_similarity_hits.csv`, `docking/screening_library.csv`
- `*.map` AutoGrid maps (regenerate with `autogrid4 -p grid.gpf`)
- the `vina.exe`, `autodock4` and `autogrid4` binaries
- `archive/` backup zips
