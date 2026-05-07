#!/usr/bin/env python3
"""
Find possible duplicate entries in the odorants database by comparing canonical SMILES.
Requires: obabel
"""

import os
import sys
import re
import subprocess

# Enforce correct path resolution for imports
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from data.odorutils import load_odors


def num_heavy_atoms(smiles: str) -> int:
    """Count non-bracket, non-digit, non-parenthesis characters in SMILES."""
    if not smiles:
        return 0
    return len(re.sub(r"[^A-Za-z]", "", smiles))


def get_canonical(odorant: dict) -> str:
    """Get canonical SMILES via obabel subprocess execution."""
    raw_smiles = odorant.get('smiles', '')
    if not raw_smiles:
        return None
        
    # --- SANITIZATION ---
    # Strip project-specific arguments (e.g., |rflp:C4) before passing to OpenBabel
    base_smiles = raw_smiles.split('|')[0].strip()
        
    try:
        # Secure subprocess execution with capture and timeout
        result = subprocess.run(
            ["obabel", f"-:{base_smiles}", "-ocan"],
            capture_output=True,
            text=True,
            timeout=5,
        )

        if result.returncode == 0 and result.stdout.strip():
            # obabel output contains the canonical SMILES as the first token
            return result.stdout.strip().split()[0]
        else:
            print(f"ERROR: obabel failed to generate output for {smiles}", file=sys.stderr)
            return None
    except subprocess.TimeoutExpired:
        print(f"ERROR: obabel timed out processing {smiles}", file=sys.stderr)
        return None
    except Exception as e:
        print(f"FATAL: Subprocess execution failed for {smiles} - {e}", file=sys.stderr)
        return None


def main():
    """Main execution loop for duplicate detection."""
    odors = load_odors()

    # -----------------------------------------
    # ARCHITECTURAL ENFORCEMENT: VALIDATE PAYLOAD
    # -----------------------------------------
    if odors is None:
        print("FATAL ERROR: load_odors() returned None.", file=sys.stderr)
        print("The system failed to extract the dataset. Check the file paths and I/O logic in data.odorutils.", file=sys.stderr)
        sys.exit(1)

    if not isinstance(odors, list) or len(odors) == 0:
        print("WARNING: load_odors() returned an empty list. No data to process.", file=sys.stderr)
        sys.exit(0)
    # -----------------------------------------

    # Execute the comparison matrix
    canonical_map = {}
    
    for odor in odors:
        # 1. Calculate canonical SMILES exactly once per molecule
        canonical = get_canonical(odor)
        if not canonical:
            continue
            
        # 2. Hash the molecule into the map
        if canonical not in canonical_map:
            canonical_map[canonical] = []
        canonical_map[canonical].append(odor)
        
    print("Pre-computation complete. Identifying duplicates...", file=sys.stderr)

    # 3. Extract duplicates directly from the map
    duplicate_count = 0
    for canonical, molecules in canonical_map.items():
        if len(molecules) > 1:
            duplicate_count += 1
            names = [m.get('full_name', 'Unknown') for m in molecules]
            
            # Format: Name1 and Name2 and Name3 may be the same...
            joined_names = " and ".join(names)
            print(f"{joined_names} may be the same molecule. The canonical SMILES is {canonical}.")
            
    if duplicate_count == 0:
        print("No duplicates detected in the dataset.")

if __name__ == "__main__":
    main()