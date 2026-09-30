    # Group 5 - Drug docking: does Vina recover the KRAS G12C binder?

    **Research question:** Does AutoDock Vina rank sotorasib above five property-matched decoys when docked into the KRAS G12C switch-II pocket (PDB 6OIM)?

    **Data:** `data/ligands/*.pdbqt` - sotorasib and five decoy kinase inhibitors, already prepared with Meeko. `data/receptor.pdbqt` and `data/box.yaml` are produced once by the instructor with `prepare_receptor.py` (needs internet).

    Everything happens in **one file, `analysis.py`**. Each function is a stub that raises
    `NotImplementedError`. The file header says who implements what.

    | Owner | Functions |
    |---|---|
    | Student A | `dock_ligand`, `pose_rmsd` |
    | Student B | `dock_all`, `plot_ranking` |
    | **Both** (this is where the merge conflict happens) | `load_box`, `main`, `summary_sentence` |

    Shared parameters live in `config.yaml`. Both students must set the value marked
    `# BOTH` — you will disagree, and Git cannot decide for you.

    Run with `python analysis.py`. It must run without error after your merge.

**Before class (instructor):** `pip install -r requirements.txt && python prepare_receptor.py` downloads 6OIM
from the RCSB, writes `data/receptor.pdbqt`, `data/crystal_ligand.pdbqt` (sotorasib as bound) and fills `data/box.yaml`.
Docking six ligands at exhaustiveness 4 takes about a minute per ligand on a laptop.

    ## Results (fill in after the merge)

    <!-- one sentence answering the research question, one figure -->

    ## Reflection

    1. What caused each merge conflict?
    2. How could branching strategy or file layout have avoided it?
    3. What is the difference between the history produced by `git pull` and `git pull --rebase`?
