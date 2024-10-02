
# http://localhost:17449/lab/tree/projects/public_resources/2023_SOX9_Sci_Adv_PB/SOX9_Sci_Adv.ipynb
getAllMarkerList.Seurat <- function(tmp.seuset, tmp.control, tmp.disease, tmp.Stem, tmp.TA, tmp.Enterocyte) {
    marker.list <- list()
    # # must have group and celltype
    tmp.seuset@active.ident <- factor(tmp.seuset$celltype)
    tmp.control.seuset <- subset(tmp.seuset, group==tmp.control)
    tmp.disease.seuset <- subset(tmp.seuset, group==tmp.disease)
    #
    normal.stem.markers <- FindMarkers(tmp.control.seuset, ident.1 = tmp.Stem, logfc.threshold = 0, min.pct = 0)
    normal.TA.markers <- FindMarkers(tmp.control.seuset, ident.1 = tmp.TA, logfc.threshold = 0, min.pct = 0)
    normal.Enterocyte.markers <- FindMarkers(tmp.control.seuset, ident.1 = tmp.Enterocyte, logfc.threshold = 0, min.pct = 0)
    #
    crc.stem.markers <- FindMarkers(tmp.disease.seuset, ident.1 = tmp.Stem, logfc.threshold = 0, min.pct = 0)
    crc.TA.markers <- FindMarkers(tmp.disease.seuset, ident.1 = tmp.TA, logfc.threshold = 0, min.pct = 0)
    crc.Enterocyte.markers <- FindMarkers(tmp.disease.seuset, ident.1 = tmp.Enterocyte, logfc.threshold = 0, min.pct = 0)
    #
    marker.list[["normal.stem"]] <- normal.stem.markers
    marker.list[["normal.TA"]] <- normal.TA.markers
    marker.list[["normal.Enterocyte"]] <- normal.Enterocyte.markers
    marker.list[["crc.stem"]] <- crc.stem.markers
    marker.list[["crc.TA"]] <- crc.TA.markers
    marker.list[["crc.Enterocyte"]] <- crc.Enterocyte.markers
    #
    tmp.seuset@active.ident <- factor(tmp.seuset$group)
    tmp.Stem.seuset <- subset(tmp.seuset, celltype==tmp.Stem)
    tmp.TA.seuset <- subset(tmp.seuset, celltype==tmp.TA)
    tmp.Enterocyte.seuset <- subset(tmp.seuset, celltype==tmp.Enterocyte)
    crc.vs.normal.stem.markers <- FindMarkers(tmp.Stem.seuset, ident.1 = tmp.disease, logfc.threshold = 0, min.pct = 0)
    crc.vs.normal.TA.markers <- FindMarkers(tmp.TA.seuset, ident.1 = tmp.disease, logfc.threshold = 0, min.pct = 0)
    crc.vs.normal.Enterocyte.markers <- FindMarkers(tmp.Enterocyte.seuset, ident.1 = tmp.disease, logfc.threshold = 0, min.pct = 0)
    #
    marker.list[["crc.vs.normal.stem"]] <- crc.vs.normal.stem.markers
    marker.list[["crc.vs.normal.TA"]] <- crc.vs.normal.TA.markers
    marker.list[["crc.vs.normal.Enterocyte"]] <- crc.vs.normal.Enterocyte.markers
    marker.list
}

