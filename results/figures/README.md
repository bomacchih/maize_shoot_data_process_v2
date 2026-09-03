# Figure output folders

Top-level folders in this directory begin with the numbered workflow step that
generates or explains their contents. Figure filenames retain the manuscript
figure and panel identifiers so they can be matched directly to the paper and
tutorial text.

| Workflow step | Folder | Contents |
| --- | --- | --- |
| 05 | `05_tissue_supergroups_Figure_12/` | Current tissue-supergroup and anatomical annotation panels for Figure 12 |
| 05 | `05_tissue_supergroups_legacy_Figure_7_8_4/` | Earlier copies retained for provenance; current documentation uses the Figure 12 folder |
| 06 | `06_pseudobulk_structural_domains_Figure_13/` | Structural-domain pseudobulk diagnostics and Figure 13 |
| 07 | `07_developmental_trends_GO/` | Developmental expression trends and GO enrichment |
| 08 | `08_RNA_velocity/` | Stochastic and dynamical RNA-velocity outputs |
| 09 | `09_monocle3_pseudotime/` | Monocle 3 root selection, trajectories, and pseudotime |
| 10 | `10_scRNA_reference_integration/` | scRNA-seq reference QC and Harmony integration |
| 10 | `10_scRNA_SCINA_annotation/` | SCINA cell-type annotation of the scRNA-seq reference |
| 10 | `10_scRNA_Visium_mapping/` | Seurat transfer and SPOTlight mapping to Visium |
| 11 | `11_TO_GCN/` | TO-GCN expression, network, and profile panels |

The workflow-step prefix identifies where an output is generated. It is not a
replacement for manuscript figure numbers such as Figure 12 or Figure 13.
