# Table output folders

Top-level folders in this directory begin with the numbered workflow step that
generates or explains their contents. Filenames retain analysis and manuscript
identifiers where those labels help connect an output to the tutorial.

| Workflow step | Folder | Contents |
| --- | --- | --- |
| 04 | `04_QC/` | Spot-, sample-, domain-, and gene-level QC summaries |
| 04 | `04_SCT_Harmony_integration/` | PCA variance and recomputed-versus-published Harmony cluster comparison |
| 05 | `05_tissue_supergroups_Figure_12/` | Cluster markers, domain counts, and tissue-supergroup annotation records |
| 06 | `06_pseudobulk_structural_domains_Figure_13/` | Replicate-aware pseudobulk matrices, metadata, and PCA coordinates |
| 07 | `07_developmental_trends_GO/` | Developmental expression trends and GO-enrichment inputs and results |
| 08 | `08_RNA_velocity/` | RNA-velocity metadata, QC summaries, run records, and package versions |
| 08 | `08_RNA_velocity_pilot/` | Retained pilot RNA-velocity run record |
| 09 | `09_monocle3_pseudotime/` | Monocle 3 root-selection, barcode-matching, and spot-pseudotime tables |
| 10 | `10_scRNA_reference_integration/` | scRNA-seq reference QC, library, cluster, and plotting summaries |
| 10 | `10_scRNA_SCINA_annotation/` | SCINA signatures, assignments, summaries, and run record |
| 10 | `10_scRNA_Visium_mapping/` | Seurat-transfer and SPOTlight checkpoints, diagnostics, proportions, and agreement tables |
| 11 | `11_TO_GCN/` | TO-GCN inputs, original program outputs, post-processing results, and GO gene lists |

The workflow-step prefix identifies where an output is generated. It does not
replace manuscript figure numbers such as Figure 12 or Figure 13.
