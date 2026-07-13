# Multi-tissue chromatin activity and human-specific TF-binding gains in human accelerated regions

Static website and reproducibility bundle for the Claude Science computational analysis manuscript.

## 🌐 Website

Published via GitHub Pages: **https://hypnolynx.github.io/har-multitissue-manuscript/**

- **[index.html](index.html)** — main manuscript (rendered)
- **[supplementary.html](supplementary.html)** — lncRNA supplementary analysis (rendered)

## Repository structure

```
.
├── index.html                 # rendered main manuscript
├── supplementary.html         # rendered supplementary
├── manuscript/                # source Markdown
│   ├── HAR_manuscript.md
│   └── HAR_LncRNA_supplementary.md
├── figures/                   # all 7 manuscript & supplementary figures (PNG)
│   ├── matrix_heatmap.png                 (Fig 1)
│   ├── nonbrain_ranking.png               (Fig 2)
│   ├── fig7_zoohar189.png                 (Fig 3)
│   ├── ZOOHAR167_flagship_mechanism.png   (Fig 4)
│   ├── evidence_tier_classification.png   (Supp)
│   ├── zooHAR_lncRNA_expression_by_tier.png (Supp)
│   └── tierA_exon_intron_summary.png      (Supp)
└── data/                      # 11 per-locus source tables (CSV)
```

## Summary

Human accelerated regions (HARs) are conserved sequences carrying a dense set of human-lineage
substitutions. Building a 312 × 16 single-cell chromatin-accessibility matrix for the Zoonomia
HARs, this work finds 45 zooHARs dominantly accessible outside the brain and screens them for
human-specific substitutions inside experimentally bound transcription-factor motifs, controlling
each affinity change against a per-locus GC-biased gene-conversion null. Eight loci show an
affinity gain that exceeds the null (P < 0.05), led by a gained CTCF site at ZOOHAR.189 near *MYLK2*.

## Data availability

All values were re-derived from primary per-locus outputs. UniBind bound TFBS are the hg38 Robust
collection; substitutions were recovered from the Ensembl EPO primate alignment and human alleles
verified against hg38. See `data/` for the full source tables.
