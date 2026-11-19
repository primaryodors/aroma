# Methodology for Aroma's public OR receptor tree

Analysis performed 17–18 November 2026; catalog integration performed 19 November 2026.

## Objective and interpretation

This document explains how the publicly shared `btree` values in `data/receptor.json` were derived. The analysis estimated relationships among the olfactory-receptor (OR) proteins in Aroma and tested whether established family and subfamily names correspond to clades. Family labels were evaluated against sequence-derived trees rather than imposed during the principal inference. The result is a **sampled human paralog gene tree**, not a mammalian species tree or a time-calibrated reconstruction of OR evolution.

The inferred master phylogram contains the full maximum-likelihood (ML) topology, including weakly supported internal relationships. Its fitted branch lengths represent estimated amino-acid substitutions per alignment position. Support and sensitivity analyses distinguished robust relationships from uncertain ones without removing uncertain branches from the topology used for the public btree.

**The OR receptor tree presented in Aroma is the master phylogram**, represented publicly by its binary ancestry paths and node labels. More precisely, the btree preserves the full inferred topology but does not encode the fitted branch lengths or support values. The distinction matters when interpreting the public tree: a binary branching pattern does not imply that every divergence is well resolved.

The analysis was motivated partly by published mammalian OR phylogenies that combine some traditional numeric families into broader groups [1]. Those groupings were hypotheses for testing, not constraints on the master tree.

## Sequence sampling and provenance

Input sequences came from the raw `sequence` fields of `data/receptor.json`. The catalog snapshot contained 434 entries: 406 ORs, six trace amine-associated receptors (TAARs), five vomeronasal type-1 receptors (VN1Rs), and 17 MS4A proteins. The ORs represented 17 numeric families and 180 named subfamilies. Ninety-five subfamilies had only one sampled receptor and therefore supplied no internal subfamily divergence to test.

The principal analysis used all 406 ORs. The 11 TAAR/VN1R sequences were used for external-rooting sensitivity analyses. MS4As were excluded from the GPCR alignment because their membrane architecture does not provide an appropriate homologous seven-transmembrane receptor core for this comparison.

The original catalog snapshot had SHA-256:

```text
a0575d84ec458243eaf9b1778b2c8522227efcb937ac237d3d1e1597e1a05129
```

The digest identifies the input before later btree edits; it is not the digest of the current catalog. Whitespace was removed from sequence strings in the analytical snapshot, and amino-acid alphabets were checked. Existing catalog alignments were not reused: some stored rows differed in length or did not preserve the entire source sequence.

Sequence versions were independently checked against official UniProt reference sequences [9]. Of 406 OR sequences, 386 matched exactly. Twenty differed somewhere in the protein, and ten had differences affecting positions represented in the primary alignment. OR4C45 was a conspicuous case, with 26 affected primary-alignment positions. Project sequences were retained in the principal analysis; a separate canonical-reference analysis tested the consequences of these differences.

## Homology alignment and positional quality control

A receptor core was extracted from five residues before the annotated start of transmembrane region 1 through 12 residues after the annotated end of transmembrane region 7, inclusive and clipped to the protein boundaries. The aim was to retain the homologous transmembrane scaffold and adjacent conserved regions while reducing dependence on poorly alignable terminal extensions.

Multiple sequence alignment used **ClustalW 2.1**, with protein input and FASTA output [2]. The principal core alignment was generated with the following option structure; other alignment settings used program defaults:

```bash
clustalw -INFILE=or_core.fasta -TYPE=PROTEIN -ALIGN -OUTPUT=FASTA -OUTFILE=or_core.aligned.fasta -QUIET
```

ClustalW guide trees were intermediate alignment products and were not treated as the final evolutionary inference. All subsequent ML analyses used the resulting amino-acid alignments.

Quality control required an unchanged taxon inventory, equal aligned row lengths, and exact recovery of every input sequence after removing alignment gaps. Raw residue coordinates were mapped onto alignment columns. Conserved Ballesteros–Weinstein (BW) `.50` positions supplied an independent check of structural register: positions 3.50, 4.50, and 7.50 agreed across all 406 ORs, and most other anchors agreed across nearly all receptors. Annotation disagreements, including an OR4C45 offset, were recorded rather than silently corrected. These checks assess consistency with project annotations; they do not establish perfect homology at every variable loop position.

The OR core alignment contained 334 columns. The primary mask retained columns occupied by an amino acid in at least 80% of OR sequences, giving **287 analyzed positions**. A stricter 95% occupancy mask retained 273 positions. A transmembrane mask retained 219 positions at which at least 80% of ORs contributed a residue annotated within a transmembrane region. This is a transmembrane-occupancy criterion, not merely a majority vote among nongap residues. Gaps were not modeled as an additional evolutionary character state.

The joint OR/outgroup core alignment contained 417 sequences and 347 columns. Occupancy thresholds were calculated from the OR rows, preventing divergent outgroups from determining which OR sites were retained. After masking, the 406 OR sequences were identical, position for position, to the primary OR-only alignment.

## Substitution-model selection and ML inference

Tree inference used **IQ-TREE 3.1.4** [3]. ModelFinder selected among 16 candidate models by the Bayesian information criterion (BIC) [4]: four empirical amino-acid exchangeability matrices, LG, WAG, JTT, and VT; matrix-associated or observed amino-acid frequencies; and four-category discrete-gamma or FreeRate among-site rate heterogeneity. This was a specified candidate set, not an exhaustive search of all protein models.

The selected model was **JTT+F+R4**. JTT specifies the exchangeability matrix, `+F` uses observed amino-acid frequencies, and `+R4` estimates four rate categories and their proportions. The model was fitted to the entire masked receptor-core alignment without structural partitioning.

The principal invocation had the following option structure [10]:

```bash
iqtree3 -s or_core.occupancy80.fasta -m MFP --mset LG,WAG,JTT,VT --mrate G,R --cmin 4 --cmax 4 -B 1000 --bnni --alrt 1000 --prefix primary
```

IQ-TREE performed a stochastic ML topology search with branch-length and model-parameter optimization. Interrupted calculations resumed from checkpoints. Previously saved candidates and constrained-search results were checked against the selected unconstrained fit; the workflow includes a refinement step if a tested candidate improves its likelihood. This guards against an obviously inferior selected tree but does not prove that the global likelihood maximum was located. The completed workflow did not constitute a large, independent multistart convergence experiment.

Support was assessed with **1,000 ultrafast bootstrap replicates**, including bootstrap-tree nearest-neighbor-interchange optimization (`--bnni`), and **1,000 SH-aLRT replicates** [5,10]. These are distinct support statistics. Neither a bootstrap percentage nor an SH-aLRT value is a posterior probability that a clade is biologically correct.

## Alignment, sequence-version, and model sensitivity

Five alternative OR analyses evaluated dependence on alignment coverage, sequence versions, and evolutionary model assumptions:

| Analysis | Data or model change | Analyzed positions | Support calculation |
| --- | --- | ---: | --- |
| Strict occupancy | Core alignment, at least 95% occupancy | 273 | 1,000 SH-aLRT and 1,000 NNI-corrected ultrafast bootstrap replicates |
| Transmembrane | At least 80% transmembrane occupancy | 219 | 1,000 SH-aLRT replicates |
| Full length | Independent full-length ClustalW alignment, at least 80% occupancy | 298 | 1,000 SH-aLRT replicates |
| Canonical reference | UniProt-reference sequences, independent ClustalW alignment, at least 80% occupancy | 286 | 1,000 SH-aLRT replicates |
| Profile heterogeneity | Primary alignment; LG+C20+F+G4-derived posterior mean site frequencies | 287 | 1,000 SH-aLRT replicates |

Except for the profile-model analysis, these runs used JTT+F+R4. The primary tree supplied a starting topology, not a fixed topology: subsequent searches could rearrange it. The canonical-reference core boundaries were transferred by global pairwise sequence mapping using BLOSUM62, gap-open score −10, and gap-extension score −0.5; its multiple alignment was still generated by ClustalW.

The posterior mean site frequency (PMSF) analysis used the primary tree as its guide through `-ft`, with `LG+C20+F+G4` as the generating profile mixture [6,11]. This tests sensitivity to heterogeneous amino-acid preferences across sites. It is a guide-dependent approximation, not a full unrestricted C20 topology search or an independent Bayesian analysis.

For conservative branch retention, the principal split had to achieve SH-aLRT ≥80 and ultrafast bootstrap ≥95, with no strongly supported incompatible split in any of the five sensitivity analyses. In alternatives with SH-aLRT support alone, ≥80 was the operational conflict threshold. Mere absence of a branch in an alternative tree was distinguished from a strongly supported contradictory branch. Split comparisons were unrooted and keyed by taxon membership, so they were independent of tip order and display orientation.

## Rooting and direction of divergence

External-rooting analyses used the joint 417-sequence alignment, followed by separate TAAR-only and VN1R-only outgroup analyses. The combined-outgroup analysis included both support statistics; the separate-outgroup analyses used SH-aLRT support. A further joint analysis used the PMSF profile model to assess model dependence of the apparent root.

Root assessment distinguished support for an outgroup cluster from support for its attachment to the OR tree. A well-supported outgroup cluster alone does not identify the earliest divergence within ORs. The workflow compared the flanking OR splits and agreement among distinct analyses, rather than counting multiple extractions from one joint tree as independent corroboration.

These analyses did **not** establish a reliable external root. The displayed Class I/Class II orientation is therefore conditional. Reversible amino-acid models estimate unrooted relationships; they do not establish direction or absolute divergence time. Outgroup divergence, compositional heterogeneity, and long-branch attraction remain relevant concerns.

## Family and subfamily monophyly tests

Numeric families, named subfamilies, and the proposed higher groups OR1/3/7, OR2/13, OR4/12, OR5/8/9, OR6/10/11, OR5/14, and Class I were evaluated by their actual sampled membership. The assessment distinguished an unrooted separating edge from rooted monophyly under the display orientation. Maximal pure components were recorded when a named group did not form a single clade.

For numeric families and proposed higher groups absent from the principal topology, additional ML searches imposed only the requested monophyly constraint while allowing other relationships to vary. These constrained fits and the unconstrained tree formed the candidate set for approximately unbiased (AU) topology tests using **10,000 resampling replicates** [7]. The workflow reports raw AU values and Holm-adjusted values across the tested constrained candidates. Multi-member subfamilies were evaluated by topology, support, and sensitivity; they did not all receive separate constrained AU tests. Single-member subfamilies were recorded as untestable for internal monophyly.

Rejection applies to the tested constrained history under the chosen data and model. Failure to reject is not evidence that a group is monophyletic, and rejection is not a proof over every possible constrained topology. The finite candidate set, stochastic search, limited sequence length, and substitution-model adequacy constrain the interpretation.

## Ancillary sequence analyses

After tree inference, IQ-TREE reconstructed marginal ancestral amino-acid states on the fixed master topology under JTT+F+R4 [10]. These calculations were conditional on the fitted alignment, topology, and model; they did not provide an independently established common OR ancestor or determine the root. Ancestral sequences are not encoded in the public btree.

Descriptive residue and seven-residue motif profiles were calculated for named groups and retained clades. Raw frequencies and Henikoff-style position-weighted frequencies were reported to reduce redundancy effects [8]. Candidate diagnostic residues required high within-group frequency and low outside-group frequency; candidate derived characters additionally required different, high-confidence parent and child ancestral states. Motifs were not used to force the ML topology. Shared states were not automatically interpreted as synapomorphies, and per-sequence enrichment tests assuming phylogenetic independence were avoided.

## Encoding and interpreting the public btree

The November 19 catalog update encoded the **full master phylogram topology** into the `btree` values of all 406 OR entries. A binary root was placed on the conditional Class I/Class II edge, resolving the three-way display root without adding an internal unrooted split. The `00` OR prefix was retained. Verification showed exact agreement of all **403 internal unrooted splits** between the encoded paths and the master tree, with all non-btree receptor data and non-OR records preserved. Binary paths encode topology only: they do not carry branch lengths, support, or the degree of confidence in a divergence.

Each `btree` string is a path from the application root to a receptor tip. After the shared `00` OR prefix, each successive `0` or `1` selects one of the two descendant branches. A prefix shared by several receptors identifies their enclosing subtree; their longest common prefix identifies their most recent shared node under the chosen rooting convention. The choice of digit for either sister branch is an ordering convention without biochemical or evolutionary significance.

Path length counts represented branching events, not substitutions or elapsed time. Similarly, the number of shared leading digits is not a calibrated measure of sequence divergence, functional similarity, or ligand specificity. Receptors in a small common subtree can have distinct ligand profiles, and related functions can evolve independently.

Major-node labels in `data/tree_nodes.json` were independently rebuilt against the new prefix memberships and retained-split evidence. Nonmonophyletic named taxa were represented by supported component clades or explicit mixed-group labels rather than assigned to misleading common prefixes. Labels such as “core” or “clade 1” identify the represented subset rather than asserting that an entire traditional family is monophyletic. These provisional clade labels are application annotations, not proposed replacements for official receptor nomenclature.

## Scope of the public record and limitations

The public result described here is the btree and its associated receptor-tree labels. This document records the sequence source, alignment rules, model-selection scope, inference settings, sensitivity analyses, and encoding conventions. It does not assume public access to intermediate analytical files. The binary paths alone cannot recover fitted branch lengths, support statistics, alignment masks, or ancestral states, and this methods description is not a complete computational replication package.

Support was evaluated for splits identified by receptor membership, rather than inferred from an internal node's position after rerooting. The same unrooted split can have different rooted interpretations, especially at the conditional Class I/Class II root. The public binary topology retains relationships below the adopted support thresholds; the presence of an internal node should therefore not be read as a claim of uniformly strong evidence.

The study did not include broad mammalian ortholog sampling, a codon model, molecular-clock calibration, gene-tree/species-tree reconciliation, or a formal recombination/gene-conversion analysis. An early-versus-late transmembrane neighbor screen was descriptive and cannot establish mosaic ancestry. Likewise, amino-acid composition diagnostics flag possible model violations but do not identify corrupted sequences by themselves. Dense paralog sampling cannot compensate fully for a short alignment or uncertain ancient homology. The resulting tree is an evidence-based working organization of the sampled receptor repertoire, with explicit uncertainty at poorly resolved divergences.

## References

1. Hayden S, Bekaert M, Crider TA, Mariani S, Murphy WJ, Teeling EC. 2010. Ecological adaptation determines functional mammalian olfactory subgenomes. *Genome Research* 20:1–9. [doi:10.1101/gr.099416.109](https://doi.org/10.1101/gr.099416.109).
2. Thompson JD, Higgins DG, Gibson TJ. 1994. CLUSTAL W: improving the sensitivity of progressive multiple sequence alignment through sequence weighting, position-specific gap penalties and weight matrix choice. *Nucleic Acids Research* 22:4673–4680. [doi:10.1093/nar/22.22.4673](https://doi.org/10.1093/nar/22.22.4673).
3. Wong TKF et al. 2026. IQ-TREE 3: phylogenomic inference software using complex evolutionary models. *Molecular Biology and Evolution* 43:msag117. [doi:10.1093/molbev/msag117](https://doi.org/10.1093/molbev/msag117).
4. Kalyaanamoorthy S, Minh BQ, Wong TKF, von Haeseler A, Jermiin LS. 2017. ModelFinder: fast model selection for accurate phylogenetic estimates. *Nature Methods* 14:587–589. [doi:10.1038/nmeth.4285](https://doi.org/10.1038/nmeth.4285).
5. Hoang DT, Chernomor O, von Haeseler A, Minh BQ, Vinh LS. 2018. UFBoot2: improving the ultrafast bootstrap approximation. *Molecular Biology and Evolution* 35:518–522. [doi:10.1093/molbev/msx281](https://doi.org/10.1093/molbev/msx281).
6. Wang H-C, Minh BQ, Susko E, Roger AJ. 2018. Modeling site heterogeneity with posterior mean site frequency profiles accelerates accurate phylogenomic estimation. *Systematic Biology* 67:216–235. [Author-hosted article](https://www.mathstat.dal.ca/~hcwang/sysbio2018.pmsf.pdf).
7. Shimodaira H. 2002. An approximately unbiased test of phylogenetic tree selection. *Systematic Biology* 51:492–508. [doi:10.1080/10635150290069913](https://doi.org/10.1080/10635150290069913).
8. Henikoff S, Henikoff JG. 1994. Position-based sequence weights. *Journal of Molecular Biology* 243:574–578. [doi:10.1016/0022-2836(94)90032-9](https://doi.org/10.1016/0022-2836(94)90032-9).
9. UniProt Consortium. UniProt protein sequence database. Canonical reference sequences consulted November 2026. [UniProt](https://www.uniprot.org/).
10. IQ-TREE developers. Command reference: support estimation, topology tests, and ancestral reconstruction. Documentation consulted November 2026. [Command reference](https://iqtree.github.io/doc/Command-Reference).
11. IQ-TREE developers. Complex models: profile mixtures and posterior mean site frequency inference. Documentation consulted November 2026. [Complex models](https://iqtree.github.io/doc/Complex-Models).
