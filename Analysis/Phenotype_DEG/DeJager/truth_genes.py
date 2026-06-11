"""
Curated truth-set and negative-control gene lists for evaluating batch correction
quality via known-phenotype DEG effect sizes.

A good batch correction should AMPLIFY effects on truth-set genes (high signal)
while NOT inflating effects on housekeeping genes (low noise).
SNR = mean(|log2FC|_truth) / mean(|log2FC|_negative_control)
"""

# ----------------------------------------------------------------------------
# SEX-DIFFERENTIAL TRUTH SET
# ----------------------------------------------------------------------------
# Genes with extreme sex-biased expression. XIST is X-inactivation-specific
# (female-only); Y-chromosome genes are male-only. These should produce
# |log2FC| > 5 in any reasonably-corrected dataset.
SEX_TRUTH_GENES = [
    "XIST",     # X-inactive specific transcript (female-biased)
    "RPS4Y1",   # Y-linked ribosomal protein
    "DDX3Y",    # Y-linked DEAD-box helicase
    "KDM5D",    # Y-linked histone demethylase
    "UTY",      # Y-linked ubiquitously transcribed
    "EIF1AY",   # Y-linked translation initiation
    "NLGN4Y",   # Y-linked neuroligin (brain-relevant)
    "USP9Y",    # Y-linked ubiquitin protease
    "ZFY",      # Y-linked zinc finger
    "TXLNGY",   # Y-linked taxilin gamma
    "TMSB4Y",   # Y-linked thymosin beta 4
    "PRKY",     # Y-linked protein kinase
]

# ----------------------------------------------------------------------------
# ALZHEIMER'S DISEASE TRUTH SET
# ----------------------------------------------------------------------------
# Top-tier AD GWAS hits + canonical expression markers from
# Mathys et al. 2019 (Nature), Bellenguez et al. 2022 (Nat Genet),
# Wightman et al. 2021 (Nat Genet). Most concentrated in microglia.
AD_TRUTH_GENES = [
    "APOE",      # strongest AD risk gene; astrocyte/microglia
    "BIN1",      # 2nd strongest AD GWAS hit
    "MS4A6A",    # microglia-specific AD GWAS
    "MS4A4A",    # microglia AD GWAS
    "CR1",       # complement receptor 1, AD GWAS
    "PICALM",    # AD GWAS, endocytosis
    "CLU",       # clusterin (apoJ), AD GWAS
    "ABCA7",     # AD GWAS, lipid transport
    "TREM2",     # microglia AD risk gene
    "CD33",      # microglia AD GWAS
    "EPHA1",     # AD GWAS
    "PLCG2",     # microglia AD GWAS
    "ABI3",      # microglia AD GWAS
    "INPP5D",    # microglia AD GWAS
    "SPI1",      # microglia master TF, AD GWAS
    "MEF2C",     # AD GWAS, neuronal/microglia
    "HLA-DRB1",  # AD GWAS, antigen presentation
    "SORL1",     # AD risk, APP processing
    "SLC24A4",   # AD GWAS
    "ZCWPW1",    # AD GWAS
    "CASS4",     # AD GWAS
    "ADAM10",    # APP processing
    "PSEN1",     # familial AD
    "PSEN2",     # familial AD
    "APP",       # amyloid precursor
    "MAPT",      # tau
    "GFAP",      # reactive astrocyte marker (elevated in AD)
    "C1QA",      # complement, microglial activation in AD
    "C1QB",      # complement
    "C1QC",      # complement
]

# ----------------------------------------------------------------------------
# NEGATIVE CONTROL: HOUSEKEEPING GENES
# ----------------------------------------------------------------------------
# Standard housekeeping genes that should NOT show differential expression
# by sex or AD status. If a batch correction inflates their |log2FC|, it's
# adding spurious noise rather than amplifying biology.
HOUSEKEEPING_GENES = [
    "ACTB",      # beta-actin
    "GAPDH",     # glyceraldehyde-3-phosphate dehydrogenase
    "PPIA",      # peptidylprolyl isomerase A
    "HPRT1",     # hypoxanthine phosphoribosyltransferase
    "B2M",       # beta-2 microglobulin
    "RPL13A",    # ribosomal protein L13a
    "RPS18",     # ribosomal protein S18
    "RPL27",     # ribosomal protein L27
    "RPL30",     # ribosomal protein L30
    "RPS27",     # ribosomal protein S27
    "RPL10",     # ribosomal protein L10
    "RPL11",     # ribosomal protein L11
    "RPLP0",     # ribosomal protein lateral stalk subunit P0
    "TBP",       # TATA-box binding protein
    "TFRC",      # transferrin receptor
    "PGK1",      # phosphoglycerate kinase
    "YWHAZ",     # 14-3-3 zeta
    "SDHA",      # succinate dehydrogenase complex flavoprotein
    "UBC",       # ubiquitin C
    "PUM1",      # pumilio homolog 1
    "POLR2A",    # RNA polymerase II
    "EEF1A1",    # eukaryotic elongation factor 1 alpha
    "EEF2",      # eukaryotic elongation factor 2
    "TUBA1B",    # tubulin alpha 1B
    "TUBB",      # beta-tubulin
    "HSP90AB1",  # heat shock protein 90
    "PSMB2",     # proteasome subunit beta 2
    "PSMB4",     # proteasome subunit beta 4
    "ATP5F1B",   # ATP synthase F1 beta
    "VPS29",     # vacuolar protein sorting 29
]


# Phenotype → truth gene set mapping
PHENOTYPE_TRUTH = {
    "msex": SEX_TRUTH_GENES,
    "cogdx_binary": AD_TRUTH_GENES,
}

# Single negative-control set used for both phenotypes
NEGATIVE_CONTROL = HOUSEKEEPING_GENES
