# REEC — Restrained Ensemble Energy Conformer

Build with `make bin/reec`; run `make check-reec` for numerical and command-line regression checks.
Run the executable from the repository root so the shared amino acid and force data resolve correctly.

```sh
bin/reec --alignment hm/OR1A1.hm.ali --target OR1A1 \
  --template OR1A1_tpl=hm/OR1A1_tpl.pdb \
  --contacts tmp/receptor.contacts --output tmp/receptor.pdb \
  --iterations 6000 --models 3
```

A target-only PIR `.ali` file also works: REEC aligns each template sequence to the target.
Explicit template rows are preferable for GPCR alignments. Rows must have equal column counts,
match the actual template residue sequence, and end in `*`. Supply repeated `--template` options
for multiple templates. Subsequent templates are rigidly fitted to the first by shared CA anchors.
The local backbone and matching sidechain with the highest sequence identity initialize each residue.
The current input contract covers one protein chain, the twenty standard amino acids, and template
residue numbers without insertion codes. MSE template atoms map to methionine. Ambiguous inputs fail
explicitly. Protein output uses target sequence positions starting at 1 and chain A.

Contact files are whitespace separated, accept `#` comments, and use target residue numbers:

```text
disulfide 97 179
contact hard 258 OH 180 OD1|OD2 2.5 3.2 100
contact soft 220 OH 290 OH 3.0 4.6 10
```

`|` lists chemically equivalent atom alternatives; satisfying one pair satisfies the contact.
A mandatory contact uses `hard`, an important preference uses `soft`. Bounds are in Angstroms;
weights are relative restraint strengths. Declared template SSBOND and sulfur CONECT records map
through the sequence alignment and become mandatory disulfides. REEC never pairs all cysteines
indiscriminately. A sulfur can participate in only one disulfide.

REEC uses the project's `Protein`, `AminoAcid`, `Atom`, and `Bond` classes for residue construction,
topology, and reference geometry. The optimizer combines bond and angle distances, rigid-bond
planarity, signed chiral volumes, trans peptide geometry, template positional restraints, contact
bounds, and a repulsive steric term. A spatial neighbor list accelerates clashes; an analytical-gradient
L-BFGS optimizer relaxes coordinates in three stages with increasing mandatory-contact strength.
Untemplated loops start as arched bridges between anchors; missing sidechains start from the residue
library. Deterministic alternative starts are available with `--models`.

The score is a restraint objective, **not** a calibrated enthalpy, free energy, DOPE score, or docking
score. REEC is a local homology modeler, with no claim of MODELLER-level sampling or accuracy.
Long insertions and conflicting contacts may fail to converge. Output contains all standard heavy
atoms; the existing HYDRO workflow can subsequently add hydrogens. Solvent, membrane, protonation,
Ramachandran distributions, and experimental uncertainty are not modeled explicitly.

Every run writes a JSON report (default `OUTPUT.pdb.json`, or `--report PATH`). Validation covers
mandatory contact violations, covalent length deviations, geometry-distance deviations, chiral
inversions, and severe nonbonded overlaps. Exit 0 publishes the PDB by atomic rename; exit 2 writes
a failed report and retains any previous PDB; exit 1 indicates input, library, or I/O errors.
The publication limits are 0.08 A maximum bond deviation, 0.18 A maximum geometry-distance deviation,
zero inverted reference chiral centers, zero overlaps above 1.0 A, and mandatory bounds within 0.05 A.
The overlap measure uses project van der Waals radii and is not the MolProbity clashscore.

External validation is still valuable. The [wwPDB validation guide](https://www.wwpdb.org/validation/2016/XrayValidationReportHelp)
describes independent bond, angle, chirality, torsion, and clash assessments; REEC's report does not
replace those assessments.

## Python workflows

MODELLER remains the default. Select REEC explicitly, or set `AROMA_HM_BACKEND=reec`:

```sh
python3 hm/dohm.py OR1A1 --backend reec --output tmp/OR1A1.reec.pdb
python3 hm/fixfail.py OR5K1 \
  out/OR5/OR5K1/OR5K1~hazelnut_pyrazine.active.model1.pdb \
  --backend reec --skip-reshape --output tmp/OR5K1.reec.pdb
```

The dohm adapter uses the existing alignment builder and composite template creator. Explicit
`--template FILE` arguments bypass composite construction; `--alignment FILE` optionally supplies
the corresponding PIR rows. The fixfail adapter runs the existing receptor-class IC reshape by
default, then rebuilds and relaxes the complete target sequence. `--skip-reshape` repairs the input
conformation directly. Both adapters accept `--iterations`, `--models`, and repeated `--contacts`.

The adapters preserve declared source disulfides, receptor-specific cysteine topology,
Y/H6.55–D/E45.51, and the Class I R/K4.52–D/E5.50 salt bridge. The Y5.58–Y7.53 interaction is a soft
preference. Output carries receptor BW and helix metadata. Publication occurs after final REEC
validation; no subsequent unconstrained IC relaxation alters the result.

Fixed ligand heavy atoms from a repair input participate in steric scoring and retain their exact
coordinates. The executable also accepts `--environment FILE.pdb` for that purpose. Ligand contact
pairs in fixfail use the existing `LIGAND_ATOM RES:ATOM` syntax. `reset`, `fit`, `predock`, `dock`, and
`loop` actions are available; docking actions publish to the normal receptor path. For reviewable
local experiments, use `--output` with no docking action. Failed numerical reports are retained as
`tmp/RECEPTOR.reec-failure.json`.
