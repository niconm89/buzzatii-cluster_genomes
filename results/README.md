# Results

Main results reported in Moreyra et al. (2023), *Molecular Phylogenetics and Evolution*. Species are abbreviated following FlyBase (see Table 1 of the paper).

## 1. Genome assembly and annotation
| File | Description | Paper |
|---|---|---|
| [Assembly contiguity statistics](01_genome_assembly_annotation/Table%201.xlsx) | Contiguity statistics for the 13 genome assemblies | Table 1 |
| [Genomic data by species](01_genome_assembly_annotation/Table%20S1%20-%20Species_genomic_data.xlsx) | Assembly accessions and versions used for each species | Table S1 |
| [Sequencing yields](01_genome_assembly_annotation/Table%20S2%20-%20Genome_sequencing_yields_and_statistics.xlsx) | Illumina and PacBio sequencing statistics | Table S2 |
| [Scaffold vs contig statistics](01_genome_assembly_annotation/Table%20S3%20-%20Scaffold%20vs%20Contigs_.xlsx) | Scaffold and contig N50/L50 and BUSCO completeness | Table S3 |
| [Repeat content](01_genome_assembly_annotation/Table%20S4%20-%20Genome_repeat_content.xlsx) | Repetitive elements per genome | Table S4 |
| [Annotation statistics](01_genome_assembly_annotation/Table%20S5%20-%20Genome_annotations_stats.xlsx) | Gene number and length per annotation | Table S5 |
| [AED distribution](01_genome_assembly_annotation/AED-distribution_annotations.txt) | Annotation edit distance distribution per genome | Figure S3 |
| BUSCO summaries: [Dald](01_genome_assembly_annotation/BUSCO_short_summaries/short_summary_Dald.txt), [Dari](01_genome_assembly_annotation/BUSCO_short_summaries/short_summary_Dari.txt), [Dato](01_genome_assembly_annotation/BUSCO_short_summaries/short_summary_Dato.txt), [Dbrb](01_genome_assembly_annotation/BUSCO_short_summaries/short_summary_Dbrb.txt), [Dbuz](01_genome_assembly_annotation/BUSCO_short_summaries/short_summary_Dbuz.txt), [Dhyd](01_genome_assembly_annotation/BUSCO_short_summaries/short_summary_Dhyd.txt), [DkoeA](01_genome_assembly_annotation/BUSCO_short_summaries/short_summary_DkoeA.txt), [DkoeB](01_genome_assembly_annotation/BUSCO_short_summaries/short_summary_DkoeB.txt), [Dmel](01_genome_assembly_annotation/BUSCO_short_summaries/short_summary_Dmel.txt), [Dmoj](01_genome_assembly_annotation/BUSCO_short_summaries/short_summary_Dmoj.txt), [Dnav](01_genome_assembly_annotation/BUSCO_short_summaries/short_summary_Dnav.txt), [Drep](01_genome_assembly_annotation/BUSCO_short_summaries/short_summary_Drep.txt), [Dvir](01_genome_assembly_annotation/BUSCO_short_summaries/short_summary_Dvir.txt) | BUSCO completeness (3,285 dipteran genes) | Figure S2 |

## 2. Phylogeny and concordance factors
| File | Description | Paper |
|---|---|---|
| [Species tree (ML)](02_phylogeny/Drosophila.treefile) | IQ-TREE maximum likelihood tree, 1,866 BUSCO groups, bootstrap support | Figure 2 |
| [Species tree with gCF/sCF](02_phylogeny/concordance_factors/concord.cf.tree) | Node labels: bootstrap / gCF / sCF | Figure 2 |
| [Concordance factor statistics](02_phylogeny/concordance_factors/concordance_gCF-sCF.stat) | IQ-TREE output: gCF, gDF1-2, gDFP, sCF, sDF1-2 per branch | Table S6 |
| [Concordance factors and ILS tests](02_phylogeny/concordance_factors/output.csv) | Concordance statistics per named clade with ILS tests (gEF_p, sEF_p) | Table S6 |
| [Concordance factors table](02_phylogeny/concordance_factors/Table%20S6%20-%20Concordance_factors.xlsx) | Table as published | Table S6 |

## 3. Divergence times
| File | Description | Paper |
|---|---|---|
| [Divergence times table](03_divergence_times/Table%20S6%20-%20Mean_age_%2B_CI_95__HPD.xlsx) | Node ages and 95% HPD for both approaches | Table S7 |
| [Node age-calibrated tree (Mya)](03_divergence_times/mcmctree_node-age/FigTree.MYA.tre) | MCMCTree dated tree, posterior means and 95% HPD | Figure 3 |
| [MCMCTree control file](03_divergence_times/mcmctree_node-age/mcmctree.ctl) | Run settings | Text S6 |
| [Input tree with calibrations](03_divergence_times/mcmctree_node-age/Drosophila.nwk.treefile) | Species tree with node age constraints | Text S6 |
| [Replicate run tree](03_divergence_times/mcmctree_node-age/replicate_run/FigTree.tre) | Independent replicate used to check convergence | Text S6 |
| [Mutation rate-calibrated tree](03_divergence_times/beast_mutation-rate/4FDS_strict_BD_dated.tre) | BEAST maximum clade credibility tree, posterior medians and 95% HPD | Figure 3 |
| [BEAST configuration](03_divergence_times/beast_mutation-rate/4FDS_strict_BD.xml) | Strict clock, Birth-Death prior | Text S6 |
| [Four-fold degenerate sites matrix](03_divergence_times/beast_mutation-rate/supermatrix.phylip) | Concatenated 4FDS alignment | Text S6 |
| [Partitions](03_divergence_times/beast_mutation-rate/partitions.tsv) · [PartitionFinder best scheme](03_divergence_times/beast_mutation-rate/best_scheme.txt) | Partitioning scheme and substitution models | Text S6 |

## 4. Orthologs and taxonomically restricted genes (TRGs)
| File | Description | Paper |
|---|---|---|
| [Orthogroups](04_orthologs_TRGs/Table%20S8%20-%20Orthogroups.tsv) | OrthoMCL orthogroups across the 13 proteomes | Table S8 |
| [TRG validation](04_orthologs_TRGs/Table%20S9%20-%20validation_TRGs.tsv) | Candidate TRGs per branch classified as validated or divergent | Figure 4, Table S9 |

## 5. Molecular evolution of TRGs
| File | Description | Paper |
|---|---|---|
| [Selection tests on TRGs](05_molecular_evolution/Table%20S10%20-%20TRGs_SelectionAnalysis.tsv) | codeml ω estimates and likelihood ratio tests | Figure 5, Table S10 |

## 6. Functional characterization of TRGs
| File | Description | Paper |
|---|---|---|
| [Enriched GO terms](06_functional_characterization/Table%20S11%20-%20enriched_GOs.tsv) | GO enrichment of TRGs per lineage | Table S11 |
| [Functional summary](06_functional_characterization/Table%20S12%20-%20Functional_summary.xlsx) | Summary of functional annotation of TRGs | Table S12 |
| [GOseq, subgenus *Drosophila*](06_functional_characterization/GOseq_per_lineage/sDrosophila.GOseq.enriched.tsv)<br>[GOseq, *mojavensis* cluster](06_functional_characterization/GOseq_per_lineage/cmoj.GOseq.enriched.tsv)<br>[GOseq, *mulleri* subgroup](06_functional_characterization/GOseq_per_lineage/smulleri.GOseq.enriched.tsv)<br>[GOseq, *buzzatii* cluster](06_functional_characterization/GOseq_per_lineage/cbuzz.GOseq.enriched.tsv)<br>[GOseq, *serido* sibling set](06_functional_characterization/GOseq_per_lineage/sibling.GOseq.enriched.tsv)<br>[GOseq, generalists (Dhyd + Dmex)](06_functional_characterization/GOseq_per_lineage/generalist.GOseq.enriched.tsv) | GOseq enrichment results per lineage | Table S11 |
| subgenus *Drosophila*: [BP](06_functional_characterization/Revigo_per_lineage/sDrosophila_Revigo_Dmel_BP.csv) · [CC](06_functional_characterization/Revigo_per_lineage/sDrosophila_Revigo_Dmel_CC.csv) · [MF](06_functional_characterization/Revigo_per_lineage/sDrosophila_Revigo_Dmel_MF.csv)<br>*buzzatii* cluster: [BP](06_functional_characterization/Revigo_per_lineage/cbuzz_Revigo_Dmel_BP.csv) · [CC](06_functional_characterization/Revigo_per_lineage/cbuzz_Revigo_Dmel_CC.csv) · [MF](06_functional_characterization/Revigo_per_lineage/cbuzz_Revigo_Dmel_MF.csv)<br>*serido* sibling set: [BP](06_functional_characterization/Revigo_per_lineage/sibling_Revigo_Dmel_BP.csv) · [CC](06_functional_characterization/Revigo_per_lineage/sibling_Revigo_Dmel_CC.csv) · [MF](06_functional_characterization/Revigo_per_lineage/sibling_Revigo_Dmel_MF.csv)<br>generalists (Dhyd + Dmex): [BP](06_functional_characterization/Revigo_per_lineage/generalist_Revigo_Dmel_BP.csv) · [CC](06_functional_characterization/Revigo_per_lineage/generalist_Revigo_Dmel_CC.csv) · [MF](06_functional_characterization/Revigo_per_lineage/generalist_Revigo_Dmel_MF.csv) | Revigo-reduced GO terms (BP, CC, MF) | Table S11 |

## Figures
| File | Description |
|---|---|
| [Figure 1](../figures/main/Figure_1.pdf) | Scaffold vs contig N50 |
| [Figure 2](../figures/main/Figure_2.pdf) | Species tree with concordance factors |
| [Figure 3](../figures/main/Figure_3.pdf) | Divergence times with both calibration approaches |
| [Figure 4](../figures/main/Figure_4.pdf) | Candidate TRGs per branch and incomplete TRGs |
| [Figure 5](../figures/main/Figure_5.pdf) | ω values for candidate TRGs |
| [Figure S1](../figures/supplementary/Figure%20S1%20-%20Assembly-protocol.png) | Genome assembly protocol |
| [Figure S2](../figures/supplementary/Figure%20S2%20-%20BUSCO%20completeness.pdf) | BUSCO completeness |
| [Figure S3](../figures/supplementary/Figure%20S3%20%20-%20AED_plot.pdf) | AED distribution |
| [Figure S4](../figures/supplementary/Figure%20S4%20-%20AED_comparison.png) | AED of TRGs vs toolkit genes |
| [Figure S5](../figures/supplementary/Figure%20S5%20-%20Protein-length-comparison.png) | Protein length of TRGs vs toolkit genes |

Genome assemblies and annotations are available at NCBI (see main README).
