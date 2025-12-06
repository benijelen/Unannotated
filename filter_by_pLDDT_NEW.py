#!/usr/bin/env python3

import os
import shutil

INPUT_DIR = "/home1/10930/bjelen/esmfold_out"
OUTPUT_DIR = os.path.join(INPUT_DIR, "high_conf")
SUMMARY_TSV = os.path.join(INPUT_DIR, "pLDDT_summary.tsv")

# IMPORTANT:
# ESMFold (HuggingFace infer_pdb) is giving pLDDT in 0-1 range in *your* setup.
# But we’ll auto-detect 0–100 just in case.

GLOBAL_MEAN_CUTOFF = 0.75      # keep whole protein if mean pLDDT >= 0.75
DOMAIN_MIN_LEN = 33            # min contiguous confident region length (aa)
DOMAIN_PLDDT_CUTOFF = 0.70     # every residue in that region must be >=0.70

def parse_plddt_from_pdb(pdb_path):
    """
    Extract per-residue pLDDT for CA atoms from an ESMFold PDB.
    pLDDT is stored in the B-factor column.
    We return a list of floats (in whatever range the file uses).
    """
    plDDT_vals = []

    with open(pdb_path) as fh:
        for line in fh:
            if not line.startswith("ATOM"):
                continue

            atom_name = line[12:16].strip()
            if atom_name != "CA":
                continue

            bfac_str = line[60:66].strip()
            try:
                plddt_val = float(bfac_str)
            except ValueError:
                continue

            plDDT_vals.append(plddt_val)

    return plDDT_vals

def normalize_plddt(plddt_list):
    """
    Return pLDDT scaled to 0–1.
    If max > 1.5, we assume the input is 0–100 and divide by 100.
    """
    if not plddt_list:
        return []

    max_val = max(plddt_list)
    if max_val > 1.5:
        # Looks like 0–100
        return [v / 100.0 for v in plddt_list]
    else:
        # Already 0–1
        return plddt_list

def has_confident_domain(plddt_list, min_len=DOMAIN_MIN_LEN, cutoff=DOMAIN_PLDDT_CUTOFF):
    """
    Return True if there's at least one contiguous run of length >= min_len
    where *every* residue pLDDT >= cutoff.
    """
    run_len = 0
    for val in plddt_list:
        if val >= cutoff:
            run_len += 1
            if run_len >= min_len:
                return True
        else:
            run_len = 0
    return False

def evaluate_structure(plddt_list):
    """
    Decide if this structure is worth keeping.
    keep if:
      global mean pLDDT >= GLOBAL_MEAN_CUTOFF
      OR
      there exists a contiguous >=33 aa domain with per-residue pLDDT >=0.70
    Returns (keep_bool, mean_plddt, domain_bool)
    """
    if not plddt_list:
        return (False, 0.0, False)

    p = normalize_plddt(plddt_list)

    mean_val = sum(p) / len(p)
    domain_ok = has_confident_domain(p)

    keep = (mean_val >= GLOBAL_MEAN_CUTOFF) or domain_ok
    return keep, mean_val, domain_ok

def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    results = []
    pdb_files = [f for f in os.listdir(INPUT_DIR) if f.endswith(".pdb")]

    for pdb_name in sorted(pdb_files):
        pdb_path = os.path.join(INPUT_DIR, pdb_name)

        raw_plddt_list = parse_plddt_from_pdb(pdb_path)
        keep, mean_plddt, domain_ok = evaluate_structure(raw_plddt_list)

        seq_len = len(raw_plddt_list)

        if keep:
            status = "KEEP"
            dest_path = os.path.join(OUTPUT_DIR, pdb_name)
            shutil.copy2(pdb_path, dest_path)
        else:
            status = "DROP"

        results.append({
            "pdb": pdb_name,
            "length": seq_len,
            "mean_pLDDT": round(mean_plddt, 3),
            "has_conf_domain": domain_ok,
            "status": status
        })

    with open(SUMMARY_TSV, "w") as out:
        header = [
            "pdb",
            "length",
            "mean_pLDDT_0to1",
            f"has_conf_domain>={DOMAIN_MIN_LEN}aa@>={DOMAIN_PLDDT_CUTOFF}",
            "status"
        ]
        out.write("\t".join(header) + "\n")
        for r in results:
            out.write(
                f"{r['pdb']}\t"
                f"{r['length']}\t"
                f"{r['mean_pLDDT']}\t"
                f"{r['has_conf_domain']}\t"
                f"{r['status']}\n"
            )

    print(f"Summary written to {SUMMARY_TSV}")
    print(f"High-confidence PDBs copied to {OUTPUT_DIR}")

if __name__ == "__main__":
    main()
