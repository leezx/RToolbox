#' description of the example function
#' @param param1 param1
#'
#' @return the return object
#' @export
#' @examples
#' example code
#'
example.function <- function(param1) {
	NULL
}

#' Get the genes of interested GO terms by some key words
#' @param organism organism (hs and mm)
#' @param keyWords keyWords (for grey function, e.g. cell cycle|DNA replication)
#'
#' @return a list with the GO terms and genes
#' @export
#' @examples
#' tmpList <- toC_get_genes_of_GO_by_keyWords(organism = "hs", keyWords = "axon", returnAll = F)
#' tmpList <- toC_get_genes_of_GO_by_keyWords(organism = "hs", 
#'                                            keyWords = "cell cycle|DNA replication|cell division|segregation")
#'
toC_get_genes_of_GO_by_keyWords <- function(organism="hs", keyWords="axon", returnAll=F) {
    ######################
    # functions from clusterProfiler
    get_GO_data <- function(OrgDb, ont, keytype) {
    GO_Env <- get_GO_Env()
    use_cached <- FALSE

    if (exists("organism", envir=GO_Env, inherits=FALSE) &&
        exists("keytype", envir=GO_Env, inherits=FALSE)) {

        org <- get("organism", envir=GO_Env)
        kt <- get("keytype", envir=GO_Env)

        if (org == DOSE:::get_organism(OrgDb) &&
            keytype == kt &&
            exists("goAnno", envir=GO_Env, inherits=FALSE)) {
            ## https://github.com/GuangchuangYu/clusterProfiler/issues/182
            ## && exists("GO2TERM", envir=GO_Env, inherits=FALSE)){

            use_cached <- TRUE
        }
    }

    if (use_cached) {
        goAnno <- get("goAnno", envir=GO_Env)
    } else {
        OrgDb <- GOSemSim:::load_OrgDb(OrgDb)
        kt <- keytypes(OrgDb)
        if (! keytype %in% kt) {
            stop("keytype is not supported...")
        }

        kk <- keys(OrgDb, keytype=keytype)
        goAnno <- suppressMessages(
            select(OrgDb, keys=kk, keytype=keytype,
                   columns=c("GOALL", "ONTOLOGYALL")))

        goAnno <- unique(goAnno[!is.na(goAnno$GOALL), ])

        assign("goAnno", goAnno, envir=GO_Env)
        assign("keytype", keytype, envir=GO_Env)
        assign("organism", DOSE:::get_organism(OrgDb), envir=GO_Env)
    }

    if (ont == "ALL") {
        GO2GENE <- unique(goAnno[, c(2,1)])
    } else {
        GO2GENE <- unique(goAnno[goAnno$ONTOLOGYALL == ont, c(2,1)])
    }

    GO_DATA <- DOSE:::build_Anno(GO2GENE, get_GO2TERM_table())

    goOnt.df <- goAnno[, c("GOALL", "ONTOLOGYALL")] %>% unique
    goOnt <- goOnt.df[,2]
    names(goOnt) <- goOnt.df[,1]
    assign("GO2ONT", goOnt, envir=GO_DATA)
    return(GO_DATA)
    }

    get_GO_Env <- function () {
        if (!exists(".GO_clusterProfiler_Env", envir = .GlobalEnv)) {
            pos <- 1
            envir <- as.environment(pos)
            assign(".GO_clusterProfiler_Env", new.env(), envir=envir)
        }
        get(".GO_clusterProfiler_Env", envir = .GlobalEnv)
    }

    get_GO2TERM_table <- function() {
        library(dplyr)
        GOTERM.df <- get_GOTERM()
        GOTERM.df[, c("go_id", "Term")] %>% unique
    }

    get_GOTERM <- function() {
        library(GO.db) # GO.db::GOTERM
        pos <- 1
        envir <- as.environment(pos)
        if (!exists(".GOTERM_Env", envir=envir)) {
            assign(".GOTERM_Env", new.env(), envir)
        }
        GOTERM_Env <- get(".GOTERM_Env", envir = envir)
        if (exists("GOTERM.df", envir = GOTERM_Env)) {
            GOTERM.df <- get("GOTERM.df", envir=GOTERM_Env)
        } else {
            GOTERM.df <- toTable(GOTERM)
            assign("GOTERM.df", GOTERM.df, envir = GOTERM_Env)
        }
        return(GOTERM.df)
    }
    #######################
    # get_GO_data <- rvcheck::get_fun_from_pkg("clusterProfiler", "get_GO_data") ## for easy installation
    if (organism=="hs") {
        GO_DATA <- get_GO_data("org.Hs.eg.db", "BP", "SYMBOL") # org.Mm.eg.db
    } else if (organism=="mm") {
        GO_DATA <- get_GO_data("org.Mm.eg.db", "BP", "SYMBOL") # org.Mm.eg.db
    } else {
        stop("currently only support hs and mm!")
    }
    returnList <- list()
    if (returnAll==T) {
        GO_DATA
    } else {
        # keyWords <- "cell cycle|DNA replication|cell division|segregation"
        tmpGO <- GO_DATA$PATHID2NAME[grep(keyWords, GO_DATA$PATHID2NAME)]
        returnList[["PATHID2NAME"]] <- tmpGO[names(tmpGO) %in% names(GO_DATA$PATHID2EXTID)]
        # the length not matched
        returnList[["PATHID2EXTID"]] <- GO_DATA$PATHID2EXTID[names(returnList[["PATHID2NAME"]])]
        returnList
    }
}
