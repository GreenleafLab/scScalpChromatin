# Labels for clusters

BroadClust <- list(
    "Tc" = "Lymphoid",
    "My" = "Myeloid",
    "Ma" = "Mast",
    "Kc" = "Keratinocytes",
    "Fb" = "Fibroblasts",
    "Ve" = "Vascular",
    "Le" = "Lymphatic",
    "Mu" = "Muscle", 
    "Me" = "Melanocytes",
    "Bc" = "Plasma",
    "Other" = "Other"
)

rna.NamedClust <- list(
    "rFb1" = "D.Fib", # Papillary/Reticular dermal fibroblasts
    "rFb2" = "D.Sheath",
    "rTc1" = "Tregs",
    "rTc2" = "CD4.Tc",
    "rTc3" = "CD8.Tc",
    "rMy1" = "DCs_1",
    "rMy2" = "Macs_1",
    "rMy3" = "CLEC9a.DC",
    "rMy4" = "M1.macs",
    "rMa1" = "Mast", # Mast cells
    "rKc1" = "Spinous.Kc_1",
    "rKc2" = "Spinous.Kc_2",
    "rKc3" = "HF.Kc_1",
    "rKc4" = "Basal.Kc_1",
    "rKc5" = "Basal.Kc_2",
    "rMu1" = "Muscle",
    "rMu2" = "Pericytes",
    "rVe1" = "Vas.Endo",
    "rLe1" = "Lymph.Endo",
    "rMe1" = "Melanocytes_1",
    "rMe2" = "Melanocytes_2",
    "rBc1" = "Plasma",
)

rna.FineClust <- list(
    # Lymphoid / T-cells
    "rTc1" = "CD4.Tc_1", # T-helper: NR3C1, RORA, IL7R, CREM, etc. 
    "rTc2" = "CD4.Tc_2", # T-helper: JUN, FOS, HSP, CD69 etc.
    "rTc3" = "CD8_Tc_1",  # Cytotoxic T-cells: CCL4, CCL5, CD8A, GZMK, IFNG  
    "rTc4" = "Treg", # Regulatory T cells: IKZF2, IL2RA, CTLA4, FOXP3 
    "rTc5" = "NK", # Natural Killer cells: XCL1, XCL2, GNLY, NKG7, KLRD1
    "rTc6" = "Cyc.Tc", # Cycling T cells: MKI67, TOP2A, etc.
    # Myeloid
    "rMy1" = "cDC2", # CD1c, CLEC10a (conventional DCs - type 2)
    "rMy2" = "M2.macs_1", # C1Qa/b/c, FOLR2, CD14, CD163, (CCL13)
    "rMy3" = "M2.macs_2", # CXCL2, CXCL3, (CCL20, S100A8/9) 
    "rMy4" = "CLEC9a.DC", # CLEC9a, CLEC4C, XCR1 https://www.frontiersin.org/articles/10.3389/fimmu.2014.00239/full
    "rMy5" = "M1.macs", # IL15, IL32, CCR7 (CCL19, CCL17)
    "rMy6" = "TCR.macs", # CD3, TCR gene positive macrophages
    "rMy7" = "TREM2.macs", # TREM2
    "rMy8" = "Plasma.contam", # (MS4A1, IGHM, CD79A; likely contaminating plasma cells / doublets)
    # Keratinocytes
    "rKc1" = "Basal.Kc_1",
    "rKc2" = "Spinous.Kc_1",
    "rKc3" = "Spinous.Kc_2", 
    "rKc4" = "Spinous.Kc_3", # Borderline HF (SOX9 / KRT75 weak)
    "rKc5" = "Infundibulum", # SOX9, DKK3 (infundibulum)
    "rKc6" = "Cyc.Kc_1", # Cycling keratinocytes
    "rKc7" = "Inf.Segment", # Lhx2, LGR5 high (Inferior segment)
    "rKc8" = "Sebaceous", # KRT17, KRT18, KRT23, CDH11, S100P, RUNX3
    "rKc9" = "Isthmus", # CD200, AR, PPARG, KRT7 (Isthmus/Sebaceous)
    "rKc10" = "Eccrine", # AQP5, KRT7
    # Fibroblasts
    "rFb1" = "D.Fib_1", # CXCL1,2,3
    "rFb2" = "D.Sheath", # COL11A1, EDNRA
    "rFb3" = "D.Fib_2", # CCL19, CXCL12
    "rFb4" = "D.Fib_3", # APCDD1, COL18A1, F13A1
    "rFb5" = "D.Fib_4", # WISP2, AOX1, ARFGEF3
    "rFb6" = "D.Fib_5", # NECAB1, SCN7A
    "rFb7" = "D.Papilla", # HHIP, PTCH1, etc.
    # Endothelial
    "rVe1" = "Vas.Endo_1",
    "rVe2" = "Vas.Endo_2",
    "rVe3" = "Vas.Endo_3",
    "rLe1" = "Lymph.Endo_1",
    "rVe4" = "Vas.Endo_4",
    "rVe5" = "Unknown",
    # Non-subclustered
    "rMa1" = "Mast",
    "rMu1" = "Muscle",
    "rMu2" = "Pericytes",
    "rMe1" = "Melanocytes_1",
    "rMe2" = "Melanocytes_2",
    "rBc1" = "Plasma",
    "Other" = "Other",
    # Hair Folicle Subclustered
    "rHF1" = "Sheath_1",
    "rHF2" = "Sheath_2",
    "rHF3" = "HFSCs",
    "rHF4" = "HG",
    "rHF5" = "Matrix"
)

atac.NamedClust <- list(
    "aFb1" = "D.Fib", # Papillary/Reticular dermal fibroblasts
    "aFb2" = "D.Sheath",
    "aTc1" = "CD4.Tc",
    "aTc2" = "CD8.Tc",
    "aTc3" = "Tregs",
    "aMy1" = "DCs_1",
    "aMy2" = "Macs_1",
    "aMy3" = "CLEC9a.DC",
    "aKc1" = "Basal.Kc_1",
    "aKc2" = "Spinous.Kc_1",
    "aKc3" = "Spinous.Kc_2",
    "aKc4" = "HF.Kc_1",
    "aKc5" = "HF.Kc_2",
    "aKc6" = "HF.Kc_3",
    "aKc7" = "HF.Kc_4",
    "aMu1" = "Muscle",
    "aMu2" = "Pericytes",
    "aVe1" = "Vas.Endo_1",
    "aVe2" = "Vas.Endo_2",
    "aLe1" = "Lymph.Endo",
    "aMe1" = "Melanocytes",
    "aBc1" = "Plasma"
)

atac.FineClust <- list(
    # Lymphoid / T-cells
    "aTc1" = "CD4.Tc_1", # T-helper: NR3C1, RORA, IL7R, CREM, etc. 
    "aTc2" = "CD4.Tc_2", # T-helper: JUN, FOS, HSP, CD69 etc.
    "aTc3" = "CD8.Tc",  # Cytotoxic T-cells: CCL4, CCL5, CD8A, GZMK, IFNG  (Also NK cells absorbed here)
    "aTc4" = "Treg", # Regulatory T cells: IKZF2, IL2RA, CTLA4, FOXP3 
    "aTc5" = "CD4.Tc_3",  
    # Myeloid
    "aMy1" = "M2.macs_1",
    "aMy2" = "cDC2_1", # CD1c, CLEC10a (conventional DCs - type 2)
    "aMy3" = "M2.macs_2", # (between M2 macs and DCs)
    "aMy4" = "M2.macs_3", # CXCL8
    "aMy5" = "CLEC9a.DC", # CLEC9a, CLEC4C, XCR1 https://www.frontiersin.org/articles/10.3389/fimmu.2014.00239/full
    "aMy6" = "M1.macs", # IL15, IL32, CCR7 (CCL19, CCL17)
    "aMy7" = "cDC2_2",
    # Keratinocytes
    "aKc1" = "Basal.Kc_1",
    "aKc2" = "Spinous.Kc_2",
    "aKc3" = "Spinous.Kc_1",
    "aKc4" = "Infundibulum", # SOX9, DKK3
    "aKc5" = "Inf.Segment_1", # Lhx2, LGR5 high
    "aKc6" = "Sebaceous", # RUNX3, KRT23, KRT18
    "aKc7" = "Inf.Segment_2", # Lhx2, LGR5 high
    "aKc8" = "Isthmus", # CD200, AR high
    "aKc9" = "Matrix", 
    "aKc10" = "Eccrine", # AQP5
    "aKc11" = "Unknown_1", # (Suspected doublet)
    # Fibroblasts
    "aFb1" = "D.Fib_1", 
    "aFb2" = "D.Fib_2", 
    "aFb3" = "D.Sheath", # COL11A1
    "aFb4" = "D.Fib_3", # NECAB1, SCN7A (rFb5)
    "aFb5" = "D.Fib_4",
    "aFb6" = "D.Papilla", # HHIP, PTCH1, etc.
    # Endothelial
    "aVe1" = "Vas.Endo_1",
    "aVe2" = "Vas.Endo_2",
    "aVe3" = "Vas.Endo_3",
    "aVe4" = "Unknown_2",
    "aLe1" = "Lymph.Endo",
    # Non-subclustered
    "aMu1" = "Muscle",
    "aMu2" = "Pericytes",
    "aMe1" = "Melanocytes",
    "aBc1" = "Plasma",
    "Other" = "Other",
    # Hair Folicle Subclustered
    "aHF1" = "Sheath_1",
    "aHF2" = "Sheath_2",
    "aHF3" = "Migrating",
    "aHF4" = "HG",
    "aHF5" = "Matrix",
    "aHF6" = "HFSCs"
)

