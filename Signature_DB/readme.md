Curate all important signatures in human cancers

# Immune escape
```
# human
MHC_I_core <- c(
  "HLA-A",
  "HLA-B",
  "HLA-C",
  "B2M",
  "TAP1",
  "TAP2",
  "TAPBP",
  "NLRC5",
  "PSMB8",
  "PSMB9",
  "IFNGR1",
  "JAK1",
  "STAT1",
  "IRF1"
)

HALLMARK_INTERFERON_GAMMA_RESPONSE <- c("ADAR","APOL6","ARID5B","ARL4A","AUTS2","B2M","BANK1","BATF2","BPGM","BST2","BTG1","C1R","C1S","CASP1","CASP3","CASP4","CASP7","CASP8","CCL2","CCL5","CCL7","CD274","CD38","CD40","CD69","CD74","CD86","CDKN1A","CFB","CFH","CIITA","CMKLR1","CMPK2","CMTR1","CSF2RB","CXCL10","CXCL11","CXCL9","DDX60","DHX58","EIF2AK2","EIF4E3","EPSTI1","FAS","FCGR1A","FGL2","FPR1","GBP4","GBP6","GCH1","GPR18","GZMA","HELZ2","HERC6","HIF1A","HLA-A","HLA-B","HLA-DMA","HLA-DQA1","HLA-DRB1","HLA-G","ICAM1","IDO1","IFI27","IFI30","IFI35","IFI44","IFI44L","IFIH1","IFIT1","IFIT2","IFIT3","IFITM2","IFITM3","IFNAR2","IL10RA","IL15","IL15RA","IL18BP","IL2RB","IL4R","IL6","IL7","IRF1","IRF2","IRF4","IRF5","IRF7","IRF8","IRF9","ISG15","ISG20","ISOC1","ITGB7","JAK2","KLRK1","LAP3","LATS2","LCP2","LGALS3BP","LY6E","LYSMD2","MARCHF1","MT2A","MTHFD2","MVP","MX1","MX2","MYD88","NAMPT","NCOA3","NFKB1","NFKBIA","NLRC5","NMI","NOD1","NUP93","OAS2","OAS3","OASL","OGFR","P2RY14","PARP12","PARP14","PDE4B","PELI1","PFKP","PIM1","PLA2G4A","PLSCR1","PML","PNP","PNPT1","PSMA2","PSMA3","PSMB10","PSMB2","PSMB8","PSMB9","PSME1","PSME2","PTGS2","PTPN1","PTPN2","PTPN6","RAPGEF6","RBCK1","RIGI","RIPK1","RIPK2","RNF213","RNF31","RSAD2","RTP4","SAMD9L","SAMHD1","SECTM1","SELP","SERPING1","SLAMF7","SLC25A28","SOCS1","SOCS3","SOD2","SP110","SPPL2A","SRI","SSPN","ST3GAL5","ST8SIA4","STAT1","STAT2","STAT3","STAT4","TAP1","TAPBP","TDRD7","TMT1B","TNFAIP2","TNFAIP3","TNFAIP6","TNFSF10","TOR1B","TRAFD1","TRIM14","TRIM21","TRIM25","TRIM26","TXNIP","UBE2L6","UPP1","USP18","VAMP5","VAMP8","VCAM1","WARS1","XAF1","XCL1","ZBP1","ZNFX1")

HALLMARK_INTERFERON_ALPHA_RESPONSE <- c("ADAR","B2M","BATF2","BST2","C1S","CASP1","CASP8","CCRL2","CD47","CD74","CMPK2","CMTR1","CNP","CSF1","CXCL10","CXCL11","DDX60","DHX58","EIF2AK2","ELF1","EPSTI1","GBP2","GBP4","GMPR","HELZ2","HERC6","HLA-C","IFI27","IFI30","IFI35","IFI44","IFI44L","IFIH1","IFIT2","IFIT3","IFITM1","IFITM2","IFITM3","IL15","IL4R","IL7","IRF1","IRF2","IRF7","IRF9","ISG15","ISG20","LAMP3","LAP3","LGALS3BP","LPAR6","LY6E","MOV10","MVB12A","MX1","NCOA7","NMI","NUB1","OAS1","OASL","OGFR","PARP12","PARP14","PARP9","PLSCR1","PNPT1","PROCR","PSMA3","PSMB8","PSMB9","PSME1","PSME2","RIPK2","RNF31","RSAD2","RTP4","SAMD9","SAMD9L","SELL","SLC25A28","SP110","STAT2","TAP1","TDRD7","TENT5A","TMEM140","TRAFD1","TRIM14","TRIM21","TRIM25","TRIM26","TRIM5","TXNIP","UBA7","UBE2L6","USP18","WARS1")

sig_antigen_processing_loading_nonoverlap <- c(
  "PSMB10",   # immunoproteasome
  "PSME1",    # PA28 alpha
  "PSME2",    # PA28 beta
  "PSME3",    # PA28 gamma
  
  "ERAP1",    # peptide trimming
  "ERAP2",
  
  "CALR",     # calreticulin
  "CANX",     # calnexin
  "PDIA3",    # ERp57
  
  "HSPA5",    # BiP / GRP78 (ER folding)
  "HSP90B1",  # gp96
  
  "SEC61A1",  # peptide translocation into ER
  "SEC61B",
  "SEC61G",
  
  "TAPBPL",   # TAPBPR (TAPBP-like, editing)
  
  "BAG6",     # peptide loading complex-associated
  "UBE2L6",   # IFN-induced ubiquitin conjugation
  
  "FBXO6",    # ER-associated degradation linkage
  "DERL1",    # retrotranslocation (ERAD)
  
  "SYVN1"     # HRD1 E3 ligase (ERAD)
)
```

```
MHC_I_core_mouse <- c(
  "H2-K1",
  "H2-D1",
  "B2m",
  "Tap1",
  "Tap2",
  "Tapbp",
  "Nlrc5",
  "Psmb8",
  "Psmb9",
  "Ifngr1",
  "Jak1",
  "Stat1",
  "Irf1"
)
HALLMARK_INTERFERON_GAMMA_RESPONSE_mouse <- unique(c(
  "Adar","Arid5b","Arl4a","Auts2","B2m","Bank1","Batf2","Bpgm","Bst2","Btg1",
  "C1ra","C1s1","Casp1","Casp3","Casp4","Casp7","Casp8",
  "Ccl2","Ccl5","Ccl7",
  "Cd274","Cd38","Cd40","Cd69","Cd74","Cd86","Cdkn1a",
  "Cfb","Cfh","Ciita","Cmklr1","Cmpk2","Cmtr1","Csf2rb",
  "Cxcl9","Cxcl10","Cxcl11",
  "Ddx60","Dhx58","Eif2ak2","Eif4e3","Epsti1","Fas","Fcgr1","Fgl2","Fpr1",
  "Gbp4","Gbp6","Gch1","Gpr18","Gzma","Helz2","Herc6","Hif1a",
  "H2-K1","H2-D1","H2-DMa","H2-Aa","H2-Ab1",
  "Icam1","Ido1",
  "Ifi27","Ifi30","Ifi35","Ifi44","Ifi44l","Ifih1",
  "Ifit1","Ifit2","Ifit3","Ifitm2","Ifitm3",
  "Ifnar2","Il10ra","Il15","Il15ra","Il18bp","Il2rb","Il4r","Il6","Il7",
  "Irf1","Irf2","Irf4","Irf5","Irf7","Irf8","Irf9",
  "Isg15","Isg20","Isoc1","Itgb7",
  "Jak2","Klrk1","Lap3","Lats2","Lcp2","Lgals3bp","Ly6e","Lysmd2",
  "Marchf1","Mt2a","Mthfd2","Mvp","Mx1","Mx2","Myd88",
  "Nampt","Ncoa3","Nfkb1","Nfkbia","Nlrc5","Nmi","Nod1","Nup93",
  "Oas2","Oas3","Oasl1","Oasl2","Ogfr",
  "P2ry14","Parp12","Parp14","Pde4b","Peli1","Pfkp","Pim1","Pla2g4a","Plscr1","Pml","Pnp","Pnpt1",
  "Psma2","Psma3","Psmb2","Psmb8","Psmb9","Psmb10","Psme1","Psme2",
  "Ptgs2","Ptpn1","Ptpn2","Ptpn6","Rapgef6","Ripk1","Ripk2","Rnf213","Rnf31","Rsad2","Rtp4",
  "Samd9l","Samhd1","Sectm1a","Selp","Serping1","Slamf7","Slc25a28","Socs1","Socs3","Sod2",
  "Sp110","Sppl2a","Sri","Sspn","St3gal5","St8sia4",
  "Stat1","Stat2","Stat3","Stat4",
  "Tap1","Tapbp","Tdrd7","Tmt1b",
  "Tnfaip2","Tnfaip3","Tnfaip6","Tnfsf10","Tor1b","Trafd1",
  "Trim14","Trim21","Trim25","Trim26","Txnip","Ube2l6","Upp1","Usp18",
  "Vamp5","Vamp8","Vcam1","Wars1","Xaf1","Xcl1","Zbp1","Znfx1"
))
HALLMARK_INTERFERON_ALPHA_RESPONSE_mouse <- unique(c(
  "Adar","B2m","Batf2","Bst2","C1s1","Casp1","Casp8","Ccrl2","Cd47","Cd74",
  "Cmpk2","Cmtr1","Cnp","Csf1","Cxcl10","Cxcl11","Ddx60","Dhx58","Eif2ak2","Elf1","Epsti1",
  "Gbp2","Gbp4","Gmpr","Helz2","Herc6",
  "H2-D1",
  "Ifi27","Ifi30","Ifi35","Ifi44","Ifi44l","Ifih1",
  "Ifit2","Ifit3","Ifitm1","Ifitm2","Ifitm3",
  "Il15","Il4r","Il7",
  "Irf1","Irf2","Irf7","Irf9",
  "Isg15","Isg20","Lamp3","Lap3","Lgals3bp","Lpar6","Ly6e",
  "Mov10","Mvb12a","Mx1","Ncoa7","Nmi","Nub1",
  "Oas1a","Oasl1","Oasl2","Ogfr",
  "Parp12","Parp14","Parp9","Plscr1","Pnpt1","Procr",
  "Psma3","Psmb8","Psmb9","Psme1","Psme2",
  "Ripk2","Rnf31","Rsad2","Rtp4","Samd9l","Sell","Slc25a28","Sp110",
  "Stat2","Tap1","Tdrd7","Tent5a","Tmem140","Trafd1",
  "Trim14","Trim21","Trim25","Trim26","Txnip","Uba7","Ube2l6","Usp18","Wars1"
))
sig_antigen_processing_loading_nonoverlap_mouse <- c(
  "Psmb10",
  "Psme1",
  "Psme2",
  "Psme3",
  "Erap1",
  "Erap2",
  "Calr",
  "Canx",
  "Pdia3",
  "Hspa5",
  "Hsp90b1",
  "Sec61a1",
  "Sec61b",
  "Sec61g",
  "Tapbpl",
  "Bag6",
  "Ube2l6",
  "Fbxo6",
  "Derl1",
  "Syvn1"
)
```

