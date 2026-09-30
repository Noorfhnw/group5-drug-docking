"""
Shared analysis script. Both students edit THIS file on their own branches.

Student A implements : dock_ligand, pose_rmsd
Student B implements : dock_all, plot_ranking
BOTH implement       : load_box, main, summary_sentence   <- expect a merge conflict here

Run: python analysis.py
"""
import yaml
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

CONFIG = yaml.safe_load(open("config.yaml"))

import glob, os
import numpy as np
from vina import Vina

# ---------- BOTH ----------
def load_box(config):
    """Return (center [x,y,z], size [x,y,z]) from data/box.yaml, using CONFIG for overrides."""
    box_path = config.get("box", "data/box.yaml")
    box = yaml.safe_load(open(box_path)) or {}

    # config.yaml wins over box.yaml when it carries an explicit center/size
    center = config.get("center") or box.get("center")
    size = config.get("size") or box.get("size")

    if center is None or size is None:
        raise ValueError(
            f"{box_path} has no center/size - run prepare_receptor.py first"
        )

    center = [float(x) for x in center]
    size = [float(x) for x in size]
    if len(center) != 3 or len(size) != 3:
        raise ValueError("center and size must each have three values")
    if any(s <= 0 for s in size):
        raise ValueError(f"box size must be positive, got {size}")

    return center, size


# ---------- Student A ----------
def dock_ligand(receptor, ligand, center, size, exhaustiveness):
    """Dock one ligand with Vina. Return (best_score, pose_pdbqt_string).
    Hint: v = Vina(sf_name='vina'); v.set_receptor(...); v.set_ligand_from_file(...);
    v.compute_vina_maps(center=..., box_size=...); v.dock(exhaustiveness=..., n_poses=5)"""
    v = Vina(sf_name="vina", verbosity=0)
    v.set_receptor(receptor)
    v.set_ligand_from_file(ligand)
    v.compute_vina_maps(center=list(center), box_size=list(size))
    v.dock(exhaustiveness=exhaustiveness, n_poses=5)

    # energies(): one row per pose, first column = total score, best pose first
    best_score = float(np.atleast_2d(v.energies())[0, 0])
    pose = v.poses(n_poses=1)
    return best_score, pose


def _pdbqt_heavy_atoms(pdbqt):
    """Parse a PDBQT (file path or text) -> (elements list, Nx3 coords) of heavy atoms."""
    if os.path.exists(pdbqt):
        with open(pdbqt) as fh:
            text = fh.read()
    else:
        text = pdbqt

    elements, coords = [], []
    for line in text.splitlines():
        if not line.startswith(("ATOM", "HETATM")):
            continue
        # PDBQT: the autodock type sits in the trailing columns (77-79)
        atype = line[77:79].strip() or line[12:16].strip()
        if atype.upper() in ("CL", "BR"):
            element = atype.capitalize()
        else:
            element = atype[0].upper()
        if element == "H":
            continue  # heavy atoms only
        # AutoDock types: A = aromatic carbon, NA/OA/SA = H-bond acceptor flavours
        element = {"A": "C"}.get(element, element)
        elements.append(element)
        coords.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
    return elements, np.asarray(coords, dtype=float)


def pose_rmsd(pose_pdbqt, crystal_pdbqt):
    """Heavy-atom RMSD between the docked pose and the crystal ligand (same atom order
    is NOT guaranteed - match by element and nearest neighbour, or use the symmetry-
    corrected RMSD from meeko/rdkit)."""
    pose_el, pose_xyz = _pdbqt_heavy_atoms(pose_pdbqt)
    xtal_el, xtal_xyz = _pdbqt_heavy_atoms(crystal_pdbqt)

    if len(pose_el) == 0 or len(xtal_el) == 0:
        raise ValueError("no heavy atoms parsed from one of the PDBQT inputs")
    if len(pose_el) != len(xtal_el):
        raise ValueError(
            f"atom count mismatch: pose has {len(pose_el)}, crystal has {len(xtal_el)}"
        )

    # Optimal one-to-one assignment, but only between atoms of the same element:
    # forbidden pairs get an infinite cost so the solver never picks them.
    from scipy.optimize import linear_sum_assignment

    d2 = ((pose_xyz[:, None, :] - xtal_xyz[None, :, :]) ** 2).sum(axis=2)
    same = np.equal(np.array(pose_el)[:, None], np.array(xtal_el)[None, :])
    cost = np.where(same, d2, 1e9)
    rows, cols = linear_sum_assignment(cost)
    if not same[rows, cols].all():
        raise ValueError("pose and crystal ligand have different element compositions")

    return float(np.sqrt(d2[rows, cols].mean()))


# ---------- Student B ----------
def dock_all(receptor, ligand_dir, center, size, exhaustiveness):
    """Dock every *.pdbqt in ligand_dir. Return a DataFrame: ligand, score."""
    ligands = sorted(glob.glob(os.path.join(ligand_dir, "*.pdbqt")))
    if not ligands:
        raise FileNotFoundError(f"no *.pdbqt ligands found in {ligand_dir}")

    rows = []
    for path in ligands:
        name = os.path.splitext(os.path.basename(path))[0]
        score, pose = dock_ligand(receptor, path, center, size, exhaustiveness)
        rows.append({"ligand": name, "score": float(score)})

        os.makedirs("results/poses", exist_ok=True)
        with open(f"results/poses/{name}.pdbqt", "w") as fh:
            fh.write(pose)

    scores = pd.DataFrame(rows).sort_values("score").reset_index(drop=True)
    scores["rank"] = np.arange(1, len(scores) + 1)
    os.makedirs("results", exist_ok=True)
    scores.to_csv("results/scores.csv", index=False)
    return scores


def plot_ranking(scores, out="results/ranking.png"):
    """Bar plot of Vina scores sorted best-first, sotorasib highlighted."""
    ranked = scores.sort_values("score").reset_index(drop=True)
    colours = ["#c0392b" if "sotorasib" in n.lower() else "#7f8c8d"
               for n in ranked["ligand"]]

    fig, ax = plt.subplots(figsize=(7, 4))
    ax.bar(ranked["ligand"], ranked["score"], color=colours)
    ax.set_ylabel("Vina score (kcal/mol)")
    ax.set_title("Docking ranking, KRAS G12C switch-II pocket (6OIM)")
    ax.tick_params(axis="x", rotation=45)
    for label in ax.get_xticklabels():
        label.set_ha("right")
    ax.invert_yaxis()  # more negative = better, so best bar sits highest
    fig.tight_layout()

    os.makedirs(os.path.dirname(out) or ".", exist_ok=True)
    fig.savefig(out, dpi=150)
    plt.close(fig)
    return out


# ---------- BOTH ----------
def summary_sentence(rmsd, scores):
    """One sentence: redocking RMSD, and the rank of sotorasib among the six ligands."""
    ranked = scores.sort_values("score").reset_index(drop=True)
    hit = ranked[ranked["ligand"].str.lower().str.contains("sotorasib")]
    if hit.empty:
        raise ValueError("sotorasib not found among the docked ligands")

    row = hit.iloc[0]
    rank = int(row.get("rank", hit.index[0] + 1))
    best_decoy = ranked.drop(index=hit.index[0])["score"].min()
    return (
        f"Redocking sotorasib into the KRAS G12C switch-II pocket reproduced the "
        f"crystal pose to {rmsd:.2f} A heavy-atom RMSD, and Vina ranked it "
        f"{rank} of {len(ranked)} ligands at {row['score']:.1f} kcal/mol "
        f"(best decoy: {best_decoy:.1f} kcal/mol)."
    )


def main():
    center, size = load_box(CONFIG)
    receptor = CONFIG["receptor"]
    exhaustiveness = int(CONFIG["exhaustiveness"])
    os.makedirs("results", exist_ok=True)

    # Student A: redock the crystal ligand and measure how well the pose is recovered
    sotorasib = os.path.join(CONFIG["ligand_dir"], "sotorasib.pdbqt")
    redock_score, redock_pose = dock_ligand(
        receptor, sotorasib, center, size, exhaustiveness
    )
    rmsd = pose_rmsd(redock_pose, CONFIG["crystal_ligand"])
    pd.DataFrame([{"ligand": "sotorasib", "score": redock_score, "rmsd": rmsd}]).to_csv(
        "results/redock.csv", index=False
    )

    # Student B: dock the whole set and plot the ranking
    scores = dock_all(receptor, CONFIG["ligand_dir"], center, size, exhaustiveness)
    plot_ranking(scores)

    print(summary_sentence(rmsd, scores))


if __name__ == "__main__":
    main()
