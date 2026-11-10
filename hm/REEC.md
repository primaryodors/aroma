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
With `--nofail` (`-N`), exit 2 also exports the least-bad candidate while its report retains
`valid: false`. The default native behavior still preserves the previous PDB on rejection.
The publication limits are 0.08 A maximum bond deviation, 0.18 A maximum geometry-distance deviation,
zero inverted reference chiral centers, zero overlaps above 1.0 A, and mandatory bounds within 0.05 A.
The overlap measure uses project van der Waals radii and is not the MolProbity clashscore.

External validation is still valuable. The [wwPDB validation guide](https://www.wwpdb.org/validation/2016/XrayValidationReportHelp)
describes independent bond, angle, chirality, torsion, and clash assessments; REEC's report does not
replace those assessments.

## Python workflows

BW lookups account for receptor insertion/deletion annotations. Missing BW positions resolve to
zero rather than an adjacent residue. Python lookup strings and native residue placeholders
address inserted residues as `7.36+1`, `7.36+2`, and so on. PDB metadata stores event counts as
`REMARK 800 SITE BW_INSERTION 7.36 2` or
`REMARK 800 SITE BW_DELETION 5.49 1`; native lookup, reverse lookup, and residue placeholders
use these records. Docking and modeling entry points refresh legacy PDB numbering annotations
without changing coordinates. Source fingerprints are refreshed only for metadata-only edits.

Run the numbering checks with `make test/bw_indel_test`, followed by
`python3 -m unittest discover -s test -p 'bw_numbering_test.py'`.

`dohm.py` defaults to REEC. Select MODELLER explicitly with `--backend modeller`,
or set `AROMA_HM_BACKEND=modeller`. The explicit command-line option takes precedence
over the environment. `fixfail.py` retains its existing backend default.

OR5K receptors use `pdbs/_validated/OR5K1.active.pdb` as their default template,
with explicit GPCR-aligned template and target rows. This also applies when selecting
MODELLER. In the REEC workflow, an explicit `--template` overrides the default.
Other receptor subfamilies retain composite template selection.

Batch refreshes through `dohm.py` try three generation attempts by default. Exhausted geometry
retries emit a protein model labeled `DOHM FALLBACK`, with `status: fallback` and
`refresh_accepted: false` in its JSON sidecar. A rejected model retains `valid: false`.
`--max-attempts 0` explicitly selects unbounded retries. `fixfail.py` retains its rejection behavior.

Validated sources receive a whole-protein geometry check. A failed source or postprocessing result
triggers REEC repair with source-conformation priors. The validated source is replaced only after
independent geometry review passes. Geometry-repaired sources carry a content hash and can be
copied on later calls without repeating radial repointing or unconstrained IC minimization.

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

## Sampling and geometry details

Every supplied template contributes soft CA distance bounds for nearby residues; these supplement
local initialization and positional anchors. Alkyl bonds have a threefold torsion preference.
Alternate starts rotate selected sidechains around CA–CB in 120-degree increments before relaxation.
Proline rings are excluded from this move. The final terminal carboxyl OXT is constructed explicitly.
Disulfides also restrain C–S–S angles and the central torsion through cross-link distances, allowing
either handedness. The approximate 104-degree angle and 90-degree torsion are consistent with
[published stereochemical modeling criteria](https://repository.ias.ac.in/4123/1/312.pdf).

The report includes selected atom alternatives, fixed ligand counts, and template CA RMS displacement.
`--keep-work` retains alignment, template, contact, and result files; the existing `nodel` option does
the same in dohm. Failed runs retain their workspace for diagnosis. `--fit-environment` rigidly fits
a ligand's source protein to the first template before adding fixed ligand atoms. The dohm `riglig`
path applies that fit automatically. Ligand heavy-atom CONECT records are carried into the output.
Hydrogen bonds use heavy-atom distance bounds; hydrogen placement and directional scoring are outside
the current optimizer. Inactive fixfail repair is supported; the `reset` action constructs active
models and is rejected for inactive repair.

## Independent publication audit

REEC normalizes backbone bond and angle references and enforces named L alpha centers and the
standard isoleucine/threonine beta centers. Peptide planarity has an explicit signed-volume term;
publication also limits its maximum volume deviation to 0.25 cubic Angstroms. Disulfides have an
explicit 2.05 A covalent bond target included in the maximum bond-error metric.

The Python adapters independently audit the written PDB before replacing a destination. This
checks standard heavy-atom inventories, serialized backbone lengths, peptide angles and twist,
named stereochemistry, disulfide geometry, and mandatory distances. The dohm path additionally
applies a rigid membrane orientation and BW centering using the same anchor definitions as the
project's UPRIGHT/BWCENTER workflow. It transforms both receptor and ligand together and audits
those final coordinates. Fixfail retains the supplied coordinate frame.

```sh
python3 hm/reec_review.py tmp/receptor.pdb --output tmp/receptor.review.json
```

Add `--plot tmp/receptor.png` for a validation plot; Matplotlib is optional. Regression checks also
apply this independent audit to every successful CLI fixture. The [worked examples](reec_examples/README.md)
include two receptor classes, a composite template, an IC reshape, and a ligand-bearing failed dock.

Explicit ligand hydrogens can participate as fixed steric obstacles with
`--environment-hydrogens`. Their supplied coordinates and connectivity are retained; the report's
heavy-atom count excludes them. This supports a ligand conformer prepared with explicit hydrogen
geometry before receptor relaxation. Connection records precede `TER`, matching the native
protein reader's chain termination behavior, so disulfides reload as covalent bonds and retain
correct sulfur valence during subsequent hydrogen preparation.
