
# RNA-seq DEG Explorer with GSEA — Gang Lab

library(shiny)
library(plotly)
library(DT)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(DOSE)
library(enrichplot)
library(edgeR)
library(pheatmap)
library(grid)
library(tools)
library(clusterProfiler)
library(org.Mm.eg.db) # Mouse
library(org.Dr.eg.db) # Zebrafish
library(org.Hs.eg.db) # Human
library(org.Rn.eg.db) # Rat
library(org.Ss.eg.db) # Pig
# library(org.Oc.eg.db) # Rabbit

options(
  shiny.maxRequestSize          = 500 * 1024^2,
  shiny.launch.browser          = TRUE,
  shiny.maxWebsocketMessageSize = 500 * 1024^2
)


# Helper functions


# Groups always come from comparison_info in RData — no fallback needed
extract_groups_from_comparison <- function(cmp, comparison_info = NULL) {
  if (!is.null(comparison_info) && cmp %in% names(comparison_info))
    return(trimws(comparison_info[[cmp]]))
  return(NULL)
}

# Detect all annotation columns in DE_output dynamically
# Returns id, symbol, gene_name, entrez, and all_annot (all annotation cols found)
detect_gene_columns <- function(df) {
  cols          <- colnames(df)
  stat_patterns <- "_logFC$|_FDR$|_PValue$|_FC$|_logCPM$|_log2CPM$"
  
  char_cols  <- cols[sapply(df, is.character)]
  annot_cols <- char_cols[!grepl(stat_patterns, char_cols)]
  
  pick_col <- function(options) {
    m <- options[options %in% annot_cols]
    if (length(m) > 0) m[1] else NULL
  }
  
  gene_id_col     <- pick_col(c("ENSEMBL"))
  gene_symbol_col <- pick_col(c("GeneSymbol", "geneName"))
  
  # only use geneName separately if not already used as gene_symbol_col
  gene_name_col <- if (!is.null(gene_symbol_col) && gene_symbol_col == "geneName") NULL
  else pick_col(c("geneName"))
  
  entrez_col <- pick_col(c("ENTREZID", "ENTREZ"))
  
  list(
    id        = gene_id_col,
    symbol    = gene_symbol_col,
    gene_name = gene_name_col,
    entrez    = entrez_col,
    all_annot = annot_cols   # all annotation cols — used to populate MA/Volcano/DEG table dynamically
  )
}


# Species / OrgDb lookup

SPECIES_MAP <- list(
  "M. musculus"   = list(OrgDb = org.Mm.eg.db, KEGGOrg = "mmu"),
  "D. rerio"      = list(OrgDb = org.Dr.eg.db, KEGGOrg = "dre"),
  "H. sapiens"    = list(OrgDb = org.Hs.eg.db, KEGGOrg = "hsa"),
  "R. norvegicus" = list(OrgDb = org.Rn.eg.db, KEGGOrg = "rno"),
  "S. scrofa"     = list(OrgDb = org.Ss.eg.db, KEGGOrg = "ssc")
)


# Run GSEA for one comparison (GO BP + KEGG)

run_gsea_for_cmp <- function(cmp, df, gene_symbol_col, OrgDb, KEGGOrg, gsea_seed) {
  logFC_col <- paste0(cmp, "_logFC")
  if (!logFC_col %in% colnames(df)) return(list(GO_BP = NULL, KEGG = NULL))
  
  gene_list <- df[[logFC_col]]
  names(gene_list) <- df[[gene_symbol_col]]
  gene_list <- na.omit(gene_list)
  gene_list <- sort(gene_list, decreasing = TRUE)
  gene_list <- gene_list[!duplicated(names(gene_list)) &
                           !is.na(names(gene_list)) & names(gene_list) != ""]
  
  go_obj <- tryCatch({
    set.seed(gsea_seed)
    gseGO(geneList      = gene_list,
          ont           = "BP",
          keyType       = "ALIAS",
          exponent      = 1,
          minGSSize     = 10,
          maxGSSize     = 500,
          pvalueCutoff  = 1.0,
          pAdjustMethod = "fdr",
          OrgDb         = OrgDb,
          seed          = TRUE,
          verbose       = FALSE)
  }, error = function(e) { cat("gseGO error for", cmp, ":", e$message, "\n"); NULL })
  
  kegg_obj <- tryCatch({
    ids       <- bitr(names(gene_list), fromType = "ALIAS", toType = "ENTREZID",
                      OrgDb = OrgDb, drop = TRUE)
    dedup_ids <- ids[!duplicated(ids$ALIAS), ]
    kegg_list <- gene_list[names(gene_list) %in% dedup_ids$ALIAS]
    names(kegg_list) <- dedup_ids$ENTREZID[match(names(kegg_list), dedup_ids$ALIAS)]
    kegg_list <- sort(na.omit(kegg_list), decreasing = TRUE)
    kegg_list <- kegg_list[!duplicated(names(kegg_list))]
    set.seed(gsea_seed)
    gseKEGG(geneList      = kegg_list,
            organism      = KEGGOrg,
            exponent      = 1,
            minGSSize     = 10,
            maxGSSize     = 500,
            pvalueCutoff  = 1.0,
            pAdjustMethod = "fdr",
            keyType       = "ncbi-geneid",
            seed          = TRUE,
            verbose       = FALSE)
  }, error = function(e) { cat("gseKEGG error for", cmp, ":", e$message, "\n"); NULL })
  
  make_gsea_entry <- function(gsea_obj) {
    if (is.null(gsea_obj) || nrow(gsea_obj@result) == 0) return(NULL)
    gsea_obj@result$Description <- gsub(" - .*$", "", gsea_obj@result$Description)
    list(gsea_object = gsea_obj, result = gsea_obj@result, gene_list = gene_list)
  }
  
  list(GO_BP = make_gsea_entry(go_obj), KEGG = make_gsea_entry(kegg_obj))
}


# Timepoint helpers

is_timepoint_cmp <- function(cmp) grepl("_Timepoint$", cmp)

# Gets ordered groups for a timepoint comparison from comparison_info
# Filters to only groups that actually have samples in sample_group_map
get_timepoint_groups <- function(cmp, sample_group_map, comparison_info) {
  groups <- comparison_info[[cmp]]
  valid  <- groups[groups %in% unique(unname(sample_group_map))]
  return(valid)
}


# process_rdata

process_rdata <- function(DE_output, sinfo, comparison_info, species = NULL, gsea_seed = NULL) {
  
  df <- DE_output %>% dplyr::mutate(dplyr::across(dplyr::where(is.character), trimws))
  
  # detect annotation columns dynamically
  gene_cols       <- detect_gene_columns(df)
  gene_id_col     <- gene_cols$id
  gene_symbol_col <- gene_cols$symbol
  all_annot_cols  <- gene_cols$all_annot
  cat("Gene ID col:", gene_id_col, "| Symbol col:", gene_symbol_col, "\n")
  
  # comparisons from logFC columns e.g. "NOR_KO_Day0_vs_NOR_WT_Day0_logFC" → "NOR_KO_Day0_vs_NOR_WT_Day0"
  comparisons  <- gsub("_logFC$", "", grep("_logFC$", colnames(df), value = TRUE))
  cat("Comparisons:", paste(comparisons, collapse = ", "), "\n")
  
  # bare_col_names = log2CPM columns with _log2CPM stripped
  # e.g. "1-MCKO-HYP-Day-12_log2CPM" → "1-MCKO-HYP-Day-12"
  log2CPM_cols   <- grep("_log2CPM$", colnames(df), value = TRUE)
  bare_col_names <- gsub("_log2CPM$", "", log2CPM_cols, ignore.case = TRUE)
  
  # detect Group column — support both "Group" and "group"
  group_col <- if ("Group" %in% colnames(sinfo)) "Group"
  else if ("group" %in% colnames(sinfo)) "group"
  else stop("sinfo must contain a 'Group' or 'group' column")
  
  # map samples to groups using sinfo$Samples
  # strip _S1, _S2 etc suffix from sinfo$Samples if present (pig datasets)
  sinfo_clean  <- gsub("_S[0-9]+$", "", sinfo$Samples)
  sinfo_group  <- sinfo[[group_col]]
  group_labels <- sinfo_group[match(bare_col_names, sinfo_clean)]
  
  if (any(is.na(group_labels)))
    stop("Some samples in DE_output could not be matched to sinfo$Samples. ",
         "Please ensure sinfo$Samples matches the log2CPM column names.")
  
  sample_group_map <- setNames(group_labels, bare_col_names)
  cat("Groups detected:", paste(sort(unique(group_labels)), collapse = ", "), "\n")
  
  unique_groups <- unique(unname(sample_group_map))
  
  # for each comparison, store which groups/samples belong to A and B
  cmp_to_sample_groups <- lapply(comparisons, function(cmp) {
    groups <- comparison_info[[cmp]]
    if (is_timepoint_cmp(cmp)) {
      valid <- groups[groups %in% unique_groups]
      return(list(A = valid, B = character(), label_A = cmp, label_B = "", timepoint = TRUE))
    } else {
      return(list(A = groups[1], B = groups[2], label_A = groups[1], label_B = groups[2], timepoint = FALSE))
    }
  })
  names(cmp_to_sample_groups) <- comparisons
  
  # build MA and Volcano data for each comparison
  # stat columns + all annotation columns added dynamically
  MA_list      <- list()
  Volcano_list <- list()
  
  for (cmp in comparisons) {
    logFC_col  <- paste0(cmp, "_logFC")
    FDR_col    <- paste0(cmp, "_FDR")
    logCPM_col <- paste0(cmp, "_logCPM")
    
    ma_df <- data.frame(
      avg_logCPM = df[[logCPM_col]],   # x-axis of MA plot — edgeR average logCPM
      logFC      = df[[logFC_col]],
      FDR        = df[[FDR_col]],
      pvalue_raw = df[[paste0(cmp, "_PValue")]],
      stringsAsFactors = FALSE)
    for (col in all_annot_cols) ma_df[[col]] <- df[[col]]
    MA_list[[cmp]] <- ma_df
    
    vol_df <- data.frame(
      logFC      = df[[logFC_col]],
      FDR        = df[[FDR_col]],
      pvalue_raw = df[[paste0(cmp, "_PValue")]],
      negLogFDR  = -log10(pmax(df[[FDR_col]], 1e-300)),  # y-axis, capped to avoid Inf
      stringsAsFactors = FALSE)
    for (col in all_annot_cols) vol_df[[col]] <- df[[col]]
    Volcano_list[[cmp]] <- vol_df
  }
  
  # expression table — GeneID + GeneSymbol + all sample log2CPM values
  # used for expression plot when user clicks a gene
  if (!is.null(gene_id_col) && !is.null(gene_symbol_col) && gene_id_col != gene_symbol_col) {
    expression_log2CPM <- df %>%
      dplyr::select(dplyr::all_of(c(gene_id_col, gene_symbol_col)), dplyr::all_of(log2CPM_cols))
    colnames(expression_log2CPM)[1:2] <- c("GeneID", "GeneSymbol")
  } else {
    id_col <- if (!is.null(gene_id_col)) gene_id_col else gene_symbol_col
    expression_log2CPM <- df %>% dplyr::select(dplyr::all_of(id_col), dplyr::all_of(log2CPM_cols))
    expression_log2CPM <- cbind(GeneID     = expression_log2CPM[[1]],
                                GeneSymbol = expression_log2CPM[[1]],
                                expression_log2CPM[, -1, drop = FALSE])
  }
  
  # MDS and PCA on all samples
  logCPMmat  <- as.matrix(expression_log2CPM[, -(1:2)])
  bare_names <- gsub("_log2CPM$", "", colnames(logCPMmat), ignore.case = TRUE)
  
  mds    <- limma::plotMDS(logCPMmat, plot = FALSE)
  mds_df <- data.frame(Sample = bare_names, Dim1 = mds$x, Dim2 = mds$y,
                       Group  = unname(sample_group_map[bare_names]), stringsAsFactors = FALSE)
  pca    <- prcomp(t(logCPMmat), center = TRUE, scale. = FALSE)
  pca_df <- data.frame(Sample = bare_names, PC1 = pca$x[,1], PC2 = pca$x[,2],
                       Group  = unname(sample_group_map[bare_names]), stringsAsFactors = FALSE)
  
  # for timepoint comparisons recompute MDS/PCA on the relevant subset only
  tp_comparisons <- comparisons[grepl("_Timepoint$", comparisons)]
  mds_subsets    <- list()
  pca_subsets    <- list()
  
  for (tp_cmp in tp_comparisons) {
    tp_groups  <- get_timepoint_groups(tp_cmp, sample_group_map, comparison_info)
    keep_idx   <- which(unname(sample_group_map[bare_names]) %in% tp_groups)
    if (length(keep_idx) < 3) next
    sub_mat    <- logCPMmat[, keep_idx, drop = FALSE]
    sub_bare   <- bare_names[keep_idx]
    sub_groups <- unname(sample_group_map[sub_bare])
    tryCatch({
      sub_mds <- limma::plotMDS(sub_mat, plot = FALSE)
      mds_subsets[[tp_cmp]] <- data.frame(
        Sample = sub_bare, Dim1 = sub_mds$x, Dim2 = sub_mds$y,
        Group  = sub_groups, stringsAsFactors = FALSE)
    }, error = function(e) cat("Subset MDS failed for", tp_cmp, ":", e$message, "\n"))
    tryCatch({
      sub_pca <- prcomp(t(sub_mat), center = TRUE, scale. = FALSE)
      pca_subsets[[tp_cmp]] <- data.frame(
        Sample = sub_bare, PC1 = sub_pca$x[,1], PC2 = sub_pca$x[,2],
        Group  = sub_groups, stringsAsFactors = FALSE)
    }, error = function(e) cat("Subset PCA failed for", tp_cmp, ":", e$message, "\n"))
  }
  
  list(
    metadata = list(
      comparisons          = comparisons,
      gene_cols            = gene_cols,           # full list: id, symbol, gene_name, entrez, all_annot
      sample_group_map     = sample_group_map,
      cmp_to_sample_groups = cmp_to_sample_groups,
      DE_df                = df,
      species              = species,
      gsea_seed            = gsea_seed,
      comparison_info      = comparison_info
    ),
    plot_data  = list(MA = MA_list, Volcano = Volcano_list,
                      MDS = mds_df, PCA = pca_df,
                      MDS_subsets = mds_subsets, PCA_subsets = pca_subsets),
    expression = list(log2CPM = expression_log2CPM)
  )
}


# UI

ui <- fluidPage(
  tags$head(
    tags$link(rel = "preconnect", href = "https://fonts.googleapis.com"),
    tags$link(rel = "stylesheet",
              href = "https://fonts.googleapis.com/css2?family=IBM+Plex+Sans:wght@300;400;500;600&family=IBM+Plex+Mono:wght@400;500&display=swap"),
    tags$style(HTML("
      body { font-family: 'IBM Plex Sans', sans-serif; font-size: 13px; background: #f7f8fa; color: #1a1a2e; }
      .app-header { background: linear-gradient(135deg, #1a1a2e 0%, #16213e 60%, #0f3460 100%); color: white; padding: 14px 28px 12px; margin-bottom: 0; display: flex; align-items: baseline; gap: 14px; border-bottom: 2px solid #e94560; }
      .app-header h2 { font-size: 20px; font-weight: 600; margin: 0; letter-spacing: -0.3px; color: white; }
      .app-header span { font-size: 11px; color: #a0aec0; font-weight: 300; letter-spacing: 0.8px; text-transform: uppercase; }
      .well { background: #ffffff; border: 1px solid #e2e8f0; border-radius: 8px; padding: 14px 16px; box-shadow: 0 1px 4px rgba(0,0,0,0.06); }
      .shiny-input-container > label { font-size: 11px; font-weight: 600; text-transform: uppercase; letter-spacing: 0.6px; color: #718096; margin-bottom: 4px; }
      .btn-file { font-size: 12px !important; padding: 4px 10px !important; background: #1a1a2e !important; color: white !important; border: none !important; border-radius: 4px !important; }
      .form-control[type=text] { font-size: 11px !important; height: 28px !important; padding: 3px 8px !important; border-radius: 4px !important; border: 1px solid #e2e8f0 !important; background: #f7f8fa !important; font-family: 'IBM Plex Mono', monospace !important; }
      .progress-bar { background-color: #e94560 !important; }
      .project-badge { background: linear-gradient(135deg, #0f3460, #16213e); color: white; border-radius: 6px; padding: 8px 12px; margin: 8px 0 6px; }
      .project-badge .proj-cutoffs { font-size: 10px; color: #a0c4ff; font-family: 'IBM Plex Mono', monospace; margin: 0; }
      .selectize-input { font-size: 12px !important; border-radius: 5px !important; border: 1px solid #cbd5e0 !important; min-height: 32px !important; white-space: nowrap !important; overflow: hidden !important; text-overflow: ellipsis !important; }
      .selectize-dropdown { width: auto !important; min-width: 100% !important; max-width: 420px !important; font-size: 12px !important; }
      .selectize-dropdown-content .option { white-space: nowrap !important; overflow: visible !important; }
      pre.shiny-text-output { font-family: 'IBM Plex Mono', monospace !important; font-size: 10.5px !important; line-height: 1.6 !important; background: #f0f4f8 !important; border: 1px solid #e2e8f0 !important; border-radius: 5px !important; padding: 8px 10px !important; color: #2d3748 !important; }
      .sidebar-section-label { font-size: 9.5px; font-weight: 700; text-transform: uppercase; letter-spacing: 0.9px; color: #a0aec0; margin: 10px 0 4px; }
      .main-panel-inner { background: #ffffff; border-radius: 8px; border: 1px solid #e2e8f0; padding: 16px 20px; box-shadow: 0 1px 4px rgba(0,0,0,0.05); }
      .nav-tabs { border-bottom: 2px solid #e2e8f0; margin-bottom: 16px; }
      .nav-tabs > li > a { font-size: 12px !important; font-weight: 500 !important; color: #718096 !important; padding: 7px 14px !important; border: none !important; border-radius: 0 !important; background: transparent !important; text-transform: uppercase; letter-spacing: 0.4px; }
      .nav-tabs > li.active > a, .nav-tabs > li.active > a:focus, .nav-tabs > li.active > a:hover { color: #1a1a2e !important; border-bottom: 2px solid #e94560 !important; background: transparent !important; font-weight: 600 !important; }
      .nav-tabs > li > a:hover { color: #1a1a2e !important; }
      .tab-section-title { font-size: 11px; font-weight: 700; text-transform: uppercase; letter-spacing: 0.8px; color: #718096; margin: 16px 0 8px; padding-bottom: 4px; border-bottom: 1px solid #edf2f7; }
      .btn-dl { font-size: 11px !important; padding: 4px 10px !important; border-radius: 4px !important; font-weight: 500 !important; margin-right: 5px !important; margin-bottom: 8px !important; }
      .btn-dl-pdf { background: #1a1a2e !important; color: white !important; border: none !important; }
      .btn-dl-png { background: white !important; color: #1a1a2e !important; border: 1px solid #cbd5e0 !important; }
      .btn-dl-csv { background: #276749 !important; color: white !important; border: none !important; }
      .stats-bar { display: flex; align-items: center; gap: 16px; padding: 5px 10px; background: #f7f8fa; border: 1px solid #e2e8f0; border-radius: 5px; margin-bottom: 4px; font-size: 12px; }
      .stats-bar .stat-up   { color: #e94560; font-weight: 600; }
      .stats-bar .stat-down { color: #0f3460; font-weight: 600; }
      .stats-bar .stat-label { color: #a0aec0; font-size: 10px; text-transform: uppercase; letter-spacing: 0.5px; }
      .legend-bar { display: flex; align-items: center; gap: 18px; padding: 5px 10px; background: #f7f8fa; border: 1px solid #e2e8f0; border-radius: 5px; margin-top: 4px; font-size: 12px; font-weight: 500; }
      .legend-dot { display: inline-block; width: 11px; height: 11px; border-radius: 50%; margin-right: 4px; vertical-align: middle; }
      .no-sig-msg { color: #718096; font-style: italic; font-size: 13px; padding: 16px 0; }
      hr { border-top: 1px solid #edf2f7; margin: 14px 0; }
      .dataTables_wrapper { font-size: 12px !important; }
      table.dataTable thead th { font-size: 11px !important; font-weight: 600 !important; text-transform: uppercase !important; letter-spacing: 0.4px !important; background: #f7f8fa !important; color: #4a5568 !important; border-bottom: 2px solid #e2e8f0 !important; }
      table.dataTable tbody tr:hover { background: #ebf8ff !important; cursor: pointer; }
      .table-click-hint { font-size: 11px; color: #a0aec0; font-style: italic; margin: 2px 0 6px; }
      .slider-section { background: #fff; border: 1px solid #e2e8f0; border-radius: 6px; padding: 6px 10px 2px; margin: 4px 0 8px; }
      .irs--shiny .irs-bar { background: #0f3460; border-top: 1px solid #0f3460; border-bottom: 1px solid #0f3460; }
      .irs--shiny .irs-handle { background: #e94560 !important; border: 2px solid #e94560 !important; }
      .irs--shiny .irs-from, .irs--shiny .irs-to, .irs--shiny .irs-single { background: #1a1a2e; }
      .species-prompt { font-size: 10px; color: #718096; font-style: italic; margin: 2px 0 6px; }
      .gsea-running-msg { background: #fffbeb; border: 1px solid #f6e05e; border-radius: 6px; padding: 8px 12px; font-size: 11px; color: #744210; margin-bottom: 8px; }
      select:disabled, .selectize-input.disabled { opacity: 0.5 !important; pointer-events: none !important; }
    "))
  ),
  
  div(class = "app-header",
      h2("RNA-seq DEG Explorer"),
      span("with GSEA \u00b7 Gang Lab")),
  
  sidebarLayout(
    sidebarPanel(width = 2,
                 div(class = "sidebar-section-label", "Species (for GSEA)"),
                 uiOutput("species_ui"),
                 div(class = "sidebar-section-label", "Data"),
                 p("1. Please refresh the page before uploading next RData file.",
                   style = "font-size:10px; color:black; font-style:italic; margin:2px 0 6px;"),
                 p("2. Please allow a few seconds for the plots to load.",
                   style = "font-size:10px; color:black; font-style:italic; margin:2px 0 6px;"),
                 fileInput("rdata", "", buttonLabel = "Browse RData\u2026",
                           placeholder = "No file selected", accept = ".RData"),
                 uiOutput("project_ui"),
                 div(class = "sidebar-section-label", "Comparison"),
                 uiOutput("cmp_ui"),
                 div(class = "sidebar-section-label", "Significance Metric"),
                 div(class = "slider-section",
                     radioButtons("sig_metric", label = NULL,
                                  choices  = c("FDR (adjusted p-value)" = "fdr",
                                               "P-value (raw)"          = "pvalue"),
                                  selected = "fdr")),
                 div(class = "sidebar-section-label", "FDR / P-value Cutoff"),
                 div(class = "slider-section",
                     sliderInput("fdr_cutoff", label = NULL,
                                 min = 0, max = 0.1, value = 0.05, step = 0.01, width = "100%")),
                 div(class = "sidebar-section-label", "logFC Cutoff"),
                 div(class = "slider-section",
                     sliderInput("logfc_cutoff", label = NULL,
                                 min = 0, max = 2, value = 1, step = 0.5, width = "100%")),
                 div(class = "sidebar-section-label", "Dataset Info"),
                 verbatimTextOutput("dataset_info")
    ),
    
    mainPanel(width = 10,
              div(class = "main-panel-inner",
                  tabsetPanel(id = "main_tabs",
                              
                              tabPanel("MDS / PCA",
                                       plotlyOutput("mdsPlot", height = "580px"), hr(),
                                       plotlyOutput("pcaPlot", height = "580px")),
                              
                              tabPanel("MA Plot",
                                       uiOutput("maStatsBar"),
                                       fluidRow(
                                         column(7, plotlyOutput("maPlot", height = "420px"), uiOutput("maLegendBar")),
                                         column(5, plotlyOutput("exprPlot", height = "420px"))
                                       ),
                                       hr(),
                                       div(class = "tab-section-title", "DEG Table \u2014 ordered by |logFC|"),
                                       div(class = "table-click-hint", "Click any row to highlight that gene on the plots above"),
                                       downloadButton("downloadMATableCSV", "\u2193 CSV (with log2CPM)", class = "btn-dl btn-dl-csv"),
                                       DTOutput("degTableMA")),
                              
                              tabPanel("Volcano Plot",
                                       uiOutput("volcanoStatsBar"),
                                       fluidRow(
                                         column(7, plotlyOutput("volcanoPlot", height = "420px"), uiOutput("volcanoLegendBar")),
                                         column(5, plotlyOutput("exprPlot2", height = "420px"))
                                       ),
                                       hr(),
                                       div(class = "tab-section-title", "DEG Table \u2014 ordered by |logFC|"),
                                       div(class = "table-click-hint", "Click any row to highlight that gene on the plots above"),
                                       downloadButton("downloadVolcanoTableCSV", "\u2193 CSV", class = "btn-dl btn-dl-csv"),
                                       DTOutput("degTableVolcano")),
                              
                              tabPanel("Heatmap",
                                       div(class = "tab-section-title", "Top 50 DEGs \u00b7 Euclidean / Complete \u00b7 Z-score"),
                                       fluidRow(column(12,
                                                       downloadButton("downloadHeatmapPDF", "\u2193 PDF", class = "btn-dl btn-dl-pdf"),
                                                       downloadButton("downloadHeatmapPNG", "\u2193 PNG", class = "btn-dl btn-dl-png")
                                       )),
                                       div(style = "width:100%; overflow-x:auto;",
                                           plotOutput("heatmapPlot", width = "100%", height = "1000px")),
                                       hr()),
                              
                              tabPanel("GSEA - GO BP",  uiOutput("gsea_go_content")),
                              tabPanel("GSEA - KEGG",   uiOutput("gsea_kegg_content"))
                  )
              )
    )
  )
)


# SERVER

server <- function(input, output, session) {
  
  rv <- reactiveValues(
    data       = NULL,
    gene       = NULL,
    gsea_cache = list(),
    gsea_running = FALSE   # tracks if GSEA is currently running
  )
  
  logfc_cut  <- reactive({ as.numeric(input$logfc_cutoff %||% 1) })
  fdr_cut    <- reactive({ as.numeric(input$fdr_cutoff   %||% 0.05) })
  sig_metric <- reactive({ input$sig_metric %||% "fdr" })
  
  classify_deg <- function(logFC, FDR, pvalue_raw = NULL, fdr_cutoff, lfc_cutoff) {
    threshold <- if (sig_metric() == "fdr") FDR else if (!is.null(pvalue_raw)) pvalue_raw else FDR
    factor(dplyr::case_when(
      threshold < fdr_cutoff & logFC >  lfc_cutoff ~ "upregulated",
      threshold < fdr_cutoff & logFC < -lfc_cutoff ~ "downregulated",
      TRUE ~ "nonDE"
    ), levels = c("nonDE", "downregulated", "upregulated"))
  }
  
  make_hover <- function(gene_id, gene_symbol, logFC, avg_logCPM = NULL, FDR) {
    base <- paste0("<b>Gene ID:</b> ", gene_id, "<br><b>Symbol:</b> ", gene_symbol,
                   "<br><b>logFC:</b> ", round(logFC, 3), "<br>")
    if (!is.null(avg_logCPM))
      base <- paste0(base, "<b>Avg logCPM:</b> ", round(avg_logCPM, 3), "<br>")
    paste0(base, "<b>FDR:</b> ", signif(FDR, 3))
  }
  
  # ── Load RData ──────────────────────────────────────────────────────────────
  observeEvent(input$rdata, {
    req(input$rdata)
    withProgress(message = "Loading RData file...", value = 0, {
      incProgress(0.3, detail = "Reading file")
      tryCatch({
        env <- new.env()
        load(input$rdata$datapath, envir = env)
        incProgress(0.5, detail = "Processing data")
        if (!exists("sinfo",          envir = env)) stop("RData must contain 'sinfo'")
        if (!exists("comparison_info",envir = env)) stop("RData must contain 'comparison_info'")
        if (!exists("DE_output",      envir = env)) stop("RData must contain 'DE_output'")
        
        species   <- if (exists("species",  envir = env)) env$species  else NULL
        gsea_seed <- if (exists("set_seed", envir = env)) env$set_seed else 123456
        cat("Species from RData:", ifelse(is.null(species), "not found", species), "\n")
        cat("GSEA seed:", gsea_seed, "\n")
        
        rv$data        <- process_rdata(env$DE_output, env$sinfo, env$comparison_info,
                                        species = species, gsea_seed = gsea_seed)
        rv$gsea_cache  <- list()
        rv$gsea_running <- FALSE
        incProgress(1, detail = "Complete!")
        showNotification("RData loaded successfully!", type = "message", duration = 3)
      }, error = function(e) {
        showNotification(paste("Error loading file:", e$message), type = "error", duration = 10)
      })
    })
  })
  
  # ── Sidebar UI ──────────────────────────────────────────────────────────────
  output$species_ui <- renderUI({
    if (!is.null(rv$data) && !is.null(rv$data$metadata$species))
      div(p(style = "font-size:11px;color:#276749;font-weight:600;margin:2px 0 6px;",
            paste0("\u2713 ", rv$data$metadata$species)))
    else
      p(class = "species-prompt", "Species will be loaded from RData automatically.")
  })
  
  output$project_ui <- renderUI({
    req(rv$data); m <- rv$data$metadata
    metric_label <- if (sig_metric() == "fdr") "FDR" else "P-value"
    div(class = "project-badge",
        p(class = "proj-cutoffs", paste0(metric_label, " < ", fdr_cut(), "  \u00b7  logFC > ", logfc_cut())),
        if (!is.null(m$species))
          p(class = "proj-cutoffs", paste0("Species: ", m$species)))
  })
  
  # comparison dropdown — disabled when GSEA is running
  output$cmp_ui <- renderUI({
    req(rv$data)
    is_running <- isTRUE(rv$gsea_running)
    tagList(
      if (is_running)
        div(class = "gsea-running-msg",
            "\u23f3 GSEA running \u2014 comparison locked. You can still switch tabs."),
      selectInput("cmp", "Select comparison", rv$data$metadata$comparisons,
                  selected = isolate(input$cmp)),
      # disable the select when GSEA is running via JS
      if (is_running) tags$script("$('#cmp')[0].selectize.disable();")
      else tags$script("if($('#cmp')[0].selectize) $('#cmp')[0].selectize.enable();")
    )
  })
  
  output$dataset_info <- renderText({
    req(rv$data); m <- rv$data$metadata; sg <- m$sample_group_map
    paste0("FDR cutoff  : ", fdr_cut(),                                          "\n",
           "logFC cutoff: ", logfc_cut(),                                         "\n",
           "Comparisons : ", length(m$comparisons),                               "\n",
           "Samples     : ", length(sg),                                          "\n",
           "Species     : ", ifelse(is.null(m$species), "not found", m$species),  "\n",
           "Groups      : ", paste(sort(unique(sg)), collapse = ", "))
  })
  
  get_group_reactive <- function(sample_names) {
    req(rv$data); sg <- rv$data$metadata$sample_group_map
    vapply(sample_names, function(s) if (s %in% names(sg)) sg[[s]] else s, character(1))
  }
  
  # ── Gene selection ──────────────────────────────────────────────────────────
  observeEvent(event_data("plotly_click", source = "ma"),
               { rv$gene <- event_data("plotly_click", source = "ma")$key })
  observeEvent(event_data("plotly_click", source = "volcano"),
               { rv$gene <- event_data("plotly_click", source = "volcano")$key })
  observeEvent(input$degTableMA_clicked_gene,
               { req(input$degTableMA_clicked_gene); rv$gene <- input$degTableMA_clicked_gene })
  observeEvent(input$degTableVolcano_clicked_gene,
               { req(input$degTableVolcano_clicked_gene); rv$gene <- input$degTableVolcano_clicked_gene })
  
  # ── Table data (display — no log2CPM, no genomic coord cols) ───────────────
  # cols we never want in the display table — genomic coordinates not useful
  EXCLUDE_TABLE_COLS <- c("Chr", "Start", "End", "Strand", "Length",
                          "Chr_x", "Chr_y", "Start_x", "Start_y")
  
  ma_table_data <- reactive({
    req(rv$data, input$cmp)
    df  <- rv$data$plot_data$MA[[input$cmp]]
    lfc <- logfc_cut(); fdr <- fdr_cut()
    gc  <- rv$data$metadata$gene_cols
    
    tbl <- df %>%
      dplyr::mutate(DE_Status  = classify_deg(logFC, FDR, pvalue_raw, fdr, lfc),
                    avg_logCPM = round(avg_logCPM, 3),
                    logFC      = round(logFC, 3),
                    AdjPValue  = signif(FDR, 3)) %>%
      dplyr::arrange(desc(abs(logFC)))
    
    # annotation cols — exclude genomic coordinate cols
    annot_display <- gc$all_annot[!gc$all_annot %in% EXCLUDE_TABLE_COLS]
    keep_cols     <- c(annot_display, "DE_Status", "avg_logCPM", "logFC", "AdjPValue")
    keep_cols     <- keep_cols[keep_cols %in% colnames(tbl)]
    tbl[, keep_cols, drop = FALSE]
  })
  
  volcano_table_data <- reactive({
    req(rv$data, input$cmp)
    df  <- rv$data$plot_data$Volcano[[input$cmp]]
    lfc <- logfc_cut(); fdr <- fdr_cut()
    gc  <- rv$data$metadata$gene_cols
    
    tbl <- df %>%
      dplyr::mutate(DE_Status     = classify_deg(logFC, FDR, pvalue_raw, fdr, lfc),
                    logFC         = round(logFC, 3),
                    AdjPValue     = signif(FDR, 3),
                    neg_log10_FDR = round(negLogFDR, 2)) %>%
      dplyr::arrange(desc(abs(logFC)))
    
    annot_display <- gc$all_annot[!gc$all_annot %in% EXCLUDE_TABLE_COLS]
    keep_cols     <- c(annot_display, "DE_Status", "logFC", "AdjPValue", "neg_log10_FDR")
    keep_cols     <- keep_cols[keep_cols %in% colnames(tbl)]
    tbl[, keep_cols, drop = FALSE]
  })
  
  # ── Download data — MA gets log2CPM, Volcano does not (no avg_logCPM either) ─
  append_log2cpm <- function(tbl, gc) {
    id_col      <- if (!is.null(gc$id)) gc$id else gc$symbol
    expr        <- rv$data$expression$log2CPM
    id_match    <- if (id_col %in% colnames(tbl)) tbl[[id_col]] else tbl[[gc$all_annot[1]]]
    expr_subset <- expr[match(id_match, expr$GeneID), -(1:2), drop = FALSE]
    colnames(expr_subset) <- gsub("_log2CPM$", "", colnames(expr_subset))
    cbind(tbl, expr_subset)
  }
  
  ma_download_data <- reactive({
    req(rv$data, input$cmp)
    append_log2cpm(ma_table_data(), rv$data$metadata$gene_cols)
  })
  
  # volcano download — no log2CPM, just the stat table
  volcano_download_data <- reactive({
    req(rv$data, input$cmp)
    volcano_table_data()  # no log2CPM appended for volcano
  })
  
  # ── MDS / PCA ───────────────────────────────────────────────────────────────
  filtered_mds_pca <- reactive({
    req(rv$data, input$cmp)
    if (is_timepoint_cmp(input$cmp)) {
      sub_mds <- rv$data$plot_data$MDS_subsets[[input$cmp]]
      sub_pca <- rv$data$plot_data$PCA_subsets[[input$cmp]]
      if (!is.null(sub_mds) && nrow(sub_mds) > 0)
        return(list(mds = sub_mds, pca = sub_pca, is_subset = TRUE))
      cmp_map      <- rv$data$metadata$cmp_to_sample_groups[[input$cmp]]
      valid_groups <- cmp_map$A
      mds <- rv$data$plot_data$MDS[rv$data$plot_data$MDS$Group %in% valid_groups, ]
      pca <- rv$data$plot_data$PCA[rv$data$plot_data$PCA$Group %in% valid_groups, ]
      list(mds = mds, pca = pca, is_subset = FALSE)
    } else {
      list(mds = rv$data$plot_data$MDS, pca = rv$data$plot_data$PCA, is_subset = FALSE)
    }
  })
  
  output$mdsPlot <- renderPlotly({
    req(rv$data, input$cmp)
    fp <- filtered_mds_pca(); df <- fp$mds
    validate(need(nrow(df) > 0, "Not enough samples for MDS"))
    title_txt <- if (is_timepoint_cmp(input$cmp)) {
      suffix <- if (fp$is_subset) " (subset \u2014 recomputed)" else " (filtered)"
      paste0("MDS Plot \u2014 ", input$cmp, suffix)
    } else "MDS Plot \u2014 all samples"
    df <- df %>% mutate(HoverText = paste0(
      "<b>Sample</b>: ", Sample, "<br><b>Group</b>: ", Group, "<br>",
      "<b>MDS1</b>: ", round(Dim1,4), "<br><b>MDS2</b>: ", round(Dim2,4)))
    plot_ly(df, x=~Dim1, y=~Dim2, color=~Group, symbol=~Group,
            text=~HoverText, hoverinfo="text",
            type="scatter", mode="markers", marker=list(size=9)) %>%
      layout(title=title_txt, xaxis=list(title="MDS1"), yaxis=list(title="MDS2"))
  })
  
  output$pcaPlot <- renderPlotly({
    req(rv$data, input$cmp)
    fp <- filtered_mds_pca(); df <- fp$pca
    validate(need(nrow(df) > 0, "Not enough samples for PCA"))
    title_txt <- if (is_timepoint_cmp(input$cmp)) {
      suffix <- if (fp$is_subset) " (subset \u2014 recomputed)" else " (filtered)"
      paste0("PCA Plot \u2014 ", input$cmp, suffix)
    } else "PCA Plot \u2014 all samples"
    df <- df %>% mutate(HoverText = paste0(
      "<b>Sample</b>: ", Sample, "<br><b>Group</b>: ", Group, "<br>",
      "<b>PC1</b>: ", round(PC1,4), "<br><b>PC2</b>: ", round(PC2,4)))
    plot_ly(df, x=~PC1, y=~PC2, color=~Group, symbol=~Group,
            text=~HoverText, hoverinfo="text",
            type="scatter", mode="markers", marker=list(size=9)) %>%
      layout(title=title_txt, xaxis=list(title="PC1"), yaxis=list(title="PC2"))
  })
  
  # ── Stats bars ──────────────────────────────────────────────────────────────
  make_stats_bar <- function(df_col, fdr, lfc) {
    deg  <- classify_deg(df_col$logFC, df_col$FDR, df_col$pvalue_raw, fdr, lfc)
    n_up <- sum(deg == "upregulated",   na.rm = TRUE)
    n_dn <- sum(deg == "downregulated", na.rm = TRUE)
    metric_label <- if (sig_metric() == "fdr") "FDR" else "P-value"
    div(class = "stats-bar",
        span(class="stat-up",   HTML(paste0("<b>\u2191 Upregulated:</b> ", metric_label, " < ", fdr, " & logFC > ", lfc,
                                            " \u2192 <b style='color:#e94560;font-size:14px;'>", n_up, "</b> genes"))),
        span(class="stat-label","\u00b7"),
        span(class="stat-down", HTML(paste0("<b>\u2193 Downregulated:</b> ", metric_label, " < ", fdr, " & logFC < \u2212", lfc,
                                            " \u2192 <b style='color:#0f3460;font-size:14px;'>", n_dn, "</b> genes"))))
  }
  
  output$maStatsBar      <- renderUI({ req(rv$data,input$cmp); make_stats_bar(rv$data$plot_data$MA[[input$cmp]],      fdr_cut(), logfc_cut()) })
  output$volcanoStatsBar <- renderUI({ req(rv$data,input$cmp); make_stats_bar(rv$data$plot_data$Volcano[[input$cmp]], fdr_cut(), logfc_cut()) })
  
  legend_bar_html <- function() {
    div(class = "legend-bar",
        span(tags$span(class="legend-dot", style="background:#e94560;"), "upregulated"),
        span(tags$span(class="legend-dot", style="background:#0f3460;"), "downregulated"),
        span(tags$span(class="legend-dot", style="background:#CBD5E0;"), "nonDE"),
        span(tags$span(class="legend-dot", style="background:#f5a623; border:2px solid #333;"), "selected"))
  }
  output$maLegendBar      <- renderUI({ req(rv$data); legend_bar_html() })
  output$volcanoLegendBar <- renderUI({ req(rv$data); legend_bar_html() })
  
  # ── MA Plot ─────────────────────────────────────────────────────────────────
  output$maPlot <- renderPlotly({
    req(rv$data, input$cmp)
    gc  <- rv$data$metadata$gene_cols
    raw <- rv$data$plot_data$MA[[input$cmp]]
    fdr <- fdr_cut(); lfc <- logfc_cut()
    id_col  <- if (!is.null(gc$id)) gc$id else gc$symbol
    sym_col <- if (!is.null(gc$symbol)) gc$symbol else id_col
    df <- raw %>% mutate(
      DEG   = classify_deg(logFC, FDR, pvalue_raw, fdr, lfc),
      hover = make_hover(.data[[id_col]], .data[[sym_col]], logFC, avg_logCPM, FDR))
    p <- plot_ly(df, x=~avg_logCPM, y=~logFC, key=~.data[[id_col]], source="ma",
                 color=~DEG,
                 colors=c("upregulated"="#e94560","downregulated"="#0f3460","nonDE"="#CBD5E0"),
                 text=~hover, hoverinfo="text", type="scatter", mode="markers",
                 marker=list(size=7, opacity=0.75)) %>%
      layout(title=list(text=paste0("<b>MA Plot</b>  <span style='font-size:12px;color:#718096'>",
                                    input$cmp,"</span>"), font=list(size=14), x=0),
             xaxis=list(title="Average log CPM",  titlefont=list(size=12)),
             yaxis=list(title="log2 Fold Change", titlefont=list(size=12)),
             showlegend=FALSE, margin=list(t=50))
    if (!is.null(rv$gene) && rv$gene != "") {
      sel <- df[df[[id_col]] == rv$gene, ]
      if (nrow(sel) > 0)
        p <- p %>% add_trace(data=sel, x=~avg_logCPM, y=~logFC,
                             type="scatter", mode="markers+text",
                             marker=list(size=15, color="#f5a623", line=list(color="#333333",width=2)),
                             text=sel[[sym_col]], textposition="top center",
                             textfont=list(size=11, color="#333333"),
                             name="Selected", inherit=FALSE, showlegend=FALSE)
    }
    event_register(p, "plotly_click"); p
  })
  
  # ── Volcano Plot ─────────────────────────────────────────────────────────────
  output$volcanoPlot <- renderPlotly({
    req(rv$data, input$cmp)
    gc  <- rv$data$metadata$gene_cols
    raw <- rv$data$plot_data$Volcano[[input$cmp]]
    fdr <- fdr_cut(); lfc <- logfc_cut()
    use_fdr <- sig_metric() == "fdr"
    id_col  <- if (!is.null(gc$id)) gc$id else gc$symbol
    sym_col <- if (!is.null(gc$symbol)) gc$symbol else id_col
    
    df <- raw %>% mutate(
      thresh    = if (use_fdr) FDR else pvalue_raw,
      negLogVal = -log10(pmax(thresh, 1e-300)),
      DEG   = classify_deg(logFC, FDR, pvalue_raw, fdr, lfc),
      hover = paste0("<b>Gene ID:</b> ", .data[[id_col]], "<br>",
                     "<b>Symbol:</b> ",  .data[[sym_col]], "<br>",
                     "<b>logFC:</b> ",   round(logFC, 3), "<br>",
                     if (use_fdr) paste0("<b>FDR:</b> ", signif(FDR, 3))
                     else paste0("<b>P-value:</b> ", signif(pvalue_raw, 3))))
    
    ylab <- if (use_fdr) "-log10(FDR)" else "-log10(P-value)"
    
    # fix y-axis range so top-right labels don't overlap
    y_max <- max(df$negLogVal, na.rm = TRUE)
    y_pad <- y_max * 0.12
    
    p <- plot_ly(df, x=~logFC, y=~negLogVal, key=~.data[[id_col]], source="volcano",
                 color=~DEG,
                 colors=c("upregulated"="#e94560","downregulated"="#0f3460","nonDE"="#CBD5E0"),
                 text=~hover, hoverinfo="text", type="scatter", mode="markers",
                 marker=list(size=5, opacity=0.75)) %>%
      layout(title=list(text=paste0("<b>Volcano Plot</b>  <span style='font-size:12px;color:#718096'>",
                                    input$cmp,"</span>"), font=list(size=14), x=0),
             xaxis=list(title="log2 Fold Change", titlefont=list(size=12)),
             yaxis=list(title=ylab, titlefont=list(size=12),
                        range=c(0, y_max + y_pad)),  # extra padding at top to prevent overlap
             showlegend=FALSE,
             shapes=list(
               list(type="line", x0=lfc,  x1=lfc,  y0=0, y1=1, yref="paper",
                    line=list(color="#0f3460", width=1, dash="dot")),
               list(type="line", x0=-lfc, x1=-lfc, y0=0, y1=1, yref="paper",
                    line=list(color="#0f3460", width=1, dash="dot")),
               list(type="line", x0=0, x1=1, xref="paper",
                    y0=-log10(fdr), y1=-log10(fdr),
                    line=list(color="#e94560", width=1, dash="dot"))),
             margin=list(t=50))
    
    if (!is.null(rv$gene) && rv$gene != "") {
      sel <- df[df[[id_col]] == rv$gene, ]
      if (nrow(sel) > 0)
        p <- p %>% add_trace(data=sel, x=~logFC, y=~negLogVal,
                             type="scatter", mode="markers+text",
                             marker=list(size=15, color="#f5a623", line=list(color="#333333",width=2)),
                             text=sel[[sym_col]], textposition="top center",
                             textfont=list(size=11, color="#333333"),
                             name="Selected", inherit=FALSE, showlegend=FALSE)
    }
    event_register(p, "plotly_click"); p
  })
  
  # ── Expression plot ──────────────────────────────────────────────────────────
  render_expr_plot <- reactive({
    req(rv$data, rv$gene, input$cmp)
    gc      <- rv$data$metadata$gene_cols
    id_col  <- if (!is.null(gc$id)) gc$id else gc$symbol
    cmp_map <- rv$data$metadata$cmp_to_sample_groups[[input$cmp]]
    if (is.null(cmp_map)) return(plotly_empty())
    set.seed(42)
    
    make_expr_hover <- function(sample, group, expr) {
      paste0("<b>Sample</b>: ", sample, "<br><b>Group</b>: ", group,
             "<br><b>Expression</b>: ", round(expr, 3))
    }
    
    if (isTRUE(cmp_map$timepoint)) {
      valid_groups <- cmp_map$A
      expr <- rv$data$expression$log2CPM %>%
        dplyr::filter(GeneID == rv$gene) %>%
        tidyr::pivot_longer(-c(GeneID, GeneSymbol), names_to="SampleCol", values_to="Expr") %>%
        dplyr::mutate(BareSample = gsub("_log2CPM$","",SampleCol),
                      Group      = get_group_reactive(BareSample)) %>%
        dplyr::filter(Group %in% valid_groups) %>%
        dplyr::mutate(Group     = factor(Group, levels=valid_groups),
                      # wider jitter spread so points don't cluster on top of each other
                      JitterX   = as.numeric(Group) + runif(dplyr::n(), -0.3, 0.3),
                      HoverText = make_expr_hover(BareSample, Group, Expr))
      gene_label <- if (!is.na(expr$GeneSymbol[1]) && expr$GeneSymbol[1] != "")
        paste0(expr$GeneSymbol[1], " (", rv$gene, ")") else rv$gene
      plot_ly(expr, x=~JitterX, y=~Expr, color=~Group,
              text=~HoverText, hoverinfo="text",
              hoveron="points",
              type="scatter", mode="markers",
              marker=list(size=8, opacity=0.8, line=list(width=0.5, color="white"))) %>%
        layout(title=list(text=paste("Expression:", gene_label), font=list(size=13)),
               xaxis=list(title="Timepoint", tickmode="array",
                          tickvals=seq_along(levels(expr$Group)),
                          ticktext=levels(expr$Group), tickangle=90,
                          zeroline=FALSE, range=c(0.3, length(levels(expr$Group))+0.7)),
               yaxis=list(title="log2 CPM"), showlegend=TRUE,
               hovermode="closest",
               hoverlabel=list(namelength=0, align="left"))
    } else {
      valid_groups <- c(cmp_map$A, cmp_map$B)
      expr <- rv$data$expression$log2CPM %>%
        dplyr::filter(GeneID == rv$gene) %>%
        tidyr::pivot_longer(-c(GeneID, GeneSymbol), names_to="SampleCol", values_to="Expr") %>%
        dplyr::mutate(BareSample = gsub("_log2CPM$","",SampleCol),
                      Group      = get_group_reactive(BareSample)) %>%
        dplyr::filter(Group %in% valid_groups) %>%
        dplyr::mutate(CompGroup = factor(Group, levels=valid_groups),
                      JitterNum = as.numeric(CompGroup) + runif(dplyr::n(), -0.3, 0.3),
                      HoverText = make_expr_hover(BareSample, Group, Expr))
      if (nrow(expr) == 0) {
        expr <- rv$data$expression$log2CPM %>%
          dplyr::filter(GeneID == rv$gene) %>%
          tidyr::pivot_longer(-c(GeneID, GeneSymbol), names_to="SampleCol", values_to="Expr") %>%
          dplyr::mutate(BareSample = gsub("_log2CPM$","",SampleCol),
                        Group      = get_group_reactive(BareSample),
                        CompGroup  = factor(Group),
                        JitterNum  = as.numeric(factor(Group)) + runif(dplyr::n(), -0.3, 0.3),
                        HoverText  = make_expr_hover(BareSample, Group, Expr))
      }
      gene_label <- if (!is.na(expr$GeneSymbol[1]) && expr$GeneSymbol[1] != "")
        paste0(expr$GeneSymbol[1], " (", rv$gene, ")") else rv$gene
      plot_ly(expr, x=~JitterNum, y=~Expr, color=~CompGroup,
              text=~HoverText, hoverinfo="text",
              hoveron="points",           # hover only on actual dot
              type="scatter", mode="markers",
              marker=list(size=8, opacity=0.8, line=list(width=0.5, color="white"))) %>%
        layout(title=list(text=paste("Expression:", gene_label), font=list(size=13)),
               xaxis=list(title="Group", tickmode="array",
                          tickvals=seq_along(levels(expr$CompGroup)),
                          ticktext=levels(expr$CompGroup), tickangle=90,
                          zeroline=FALSE, range=c(0.3, length(levels(expr$CompGroup))+0.7)),
               yaxis=list(title="log2 CPM"), showlegend=TRUE,
               hovermode="closest",
               hoverlabel=list(namelength=0, align="left"))
    }
  })
  
  output$exprPlot  <- renderPlotly({ render_expr_plot() })
  output$exprPlot2 <- renderPlotly({ render_expr_plot() })
  
  # ── DEG Tables ──────────────────────────────────────────────────────────────
  # helper — finds which column index in the display table holds the gene ID
  # so the JS click always sends the right ID regardless of column order
  get_id_col_index <- function(tbl, gc) {
    id_col <- if (!is.null(gc$id)) gc$id else gc$symbol
    idx    <- which(colnames(tbl) == id_col) - 1  # JS is 0-indexed
    if (length(idx) == 0) idx <- 0  # fallback to first column
    idx[1]
  }
  
  output$degTableMA <- renderDT({
    tbl <- ma_table_data()
    gc  <- rv$data$metadata$gene_cols
    idx <- get_id_col_index(tbl, gc)
    datatable(tbl, rownames=FALSE, selection="none",
              options=list(scrollX=TRUE, scrollY="300px", paging=TRUE, pageLength=10,
                           initComplete=JS(sprintf(
                             "function(settings){var api=this.api();api.on('click','tbody tr',function(){var d=api.row(this).data();if(d)Shiny.setInputValue('degTableMA_clicked_gene',d[%d],{priority:'event'});});}", idx))),
              filter="top")
  })
  
  output$degTableVolcano <- renderDT({
    tbl <- volcano_table_data()
    gc  <- rv$data$metadata$gene_cols
    idx <- get_id_col_index(tbl, gc)
    datatable(tbl, rownames=FALSE, selection="none",
              options=list(scrollX=TRUE, scrollY="300px", paging=TRUE, pageLength=10,
                           initComplete=JS(sprintf(
                             "function(settings){var api=this.api();api.on('click','tbody tr',function(){var d=api.row(this).data();if(d)Shiny.setInputValue('degTableVolcano_clicked_gene',d[%d],{priority:'event'});});}", idx))),
              filter="top")
  })
  
  # downloads — MA includes log2CPM, Volcano is stat table only
  output$downloadMATableCSV <- downloadHandler(
    filename = function() paste0(input$cmp, "_MA_DEG_table_with_log2CPM.csv"),
    content  = function(file) write.csv(ma_download_data(), file, row.names=FALSE))
  output$downloadVolcanoTableCSV <- downloadHandler(
    filename = function() paste0(input$cmp, "_Volcano_DEG_table.csv"),
    content  = function(file) write.csv(volcano_download_data(), file, row.names=FALSE))
  
  # ── Heatmap ─────────────────────────────────────────────────────────────────
  heatmap_data <- reactive({
    req(rv$data, input$cmp)
    fdr <- fdr_cut(); lfc <- logfc_cut()
    gc  <- rv$data$metadata$gene_cols
    id_col <- if (!is.null(gc$id)) gc$id else gc$symbol
    
    df_raw <- rv$data$plot_data$MA[[input$cmp]]
    degs   <- df_raw %>%
      dplyr::mutate(DEG = classify_deg(logFC, FDR, pvalue_raw, fdr, lfc)) %>%
      dplyr::filter(DEG != "nonDE" & abs(logFC) > lfc) %>%
      dplyr::arrange(FDR) %>% head(50)
    validate(need(nrow(degs) > 0, paste("No significant DEGs found for:", input$cmp)))
    
    cmp_map <- rv$data$metadata$cmp_to_sample_groups[[input$cmp]]
    sg      <- rv$data$metadata$sample_group_map
    
    valid_groups <- if (isTRUE(cmp_map$timepoint)) cmp_map$A else c(cmp_map$A, cmp_map$B)
    keep_samples <- names(sg)[sg %in% valid_groups]
    
    sym_col <- if (!is.null(gc$symbol)) gc$symbol else id_col
    
    expr <- rv$data$expression$log2CPM %>%
      dplyr::filter(GeneID %in% degs[[id_col]]) %>%
      dplyr::mutate(GeneLabel = ifelse(is.na(GeneSymbol) | GeneSymbol=="" | GeneSymbol==GeneID,
                                       GeneID, paste0(GeneSymbol," (",GeneID,")"))) %>%
      tibble::column_to_rownames("GeneLabel") %>%
      dplyr::select(-GeneID, -GeneSymbol)
    
    colnames(expr) <- gsub("_log2CPM$","", colnames(expr))
    keep_cols <- colnames(expr)[colnames(expr) %in% keep_samples]
    expr      <- expr[, keep_cols, drop=FALSE]
    
    mat     <- as.matrix(expr)
    row_var <- apply(mat, 1, var, na.rm = TRUE)
    mat     <- mat[!is.na(row_var) & row_var > 0, , drop = FALSE]
    validate(need(nrow(mat) > 0, "No genes with sufficient variance for heatmap."))
    mat_scaled <- t(scale(t(mat)))
    mat_scaled[!is.finite(mat_scaled)] <- 0
    mat_scaled[mat_scaled >  3] <-  3
    mat_scaled[mat_scaled < -3] <- -3
    
    col_annotation <- data.frame(Group = unname(sg[colnames(mat_scaled)]),
                                 row.names = colnames(mat_scaled))
    list(mat=mat_scaled, annotation=col_annotation,
         title=paste0("Top ",nrow(degs)," DEGs (FDR<",fdr," & |logFC|>",lfc,") \u2014 z-score: ",input$cmp))
  })
  
  draw_pheatmap <- function(hd) {
    row_hc <- hclust(dist(hd$mat,    method="euclidean"), method="complete")
    col_hc <- hclust(dist(t(hd$mat), method="euclidean"), method="complete")
    pheatmap(hd$mat, cluster_rows=row_hc, cluster_cols=col_hc,
             clustering_distance_rows="euclidean", clustering_distance_cols="euclidean",
             clustering_method="complete", scale="none",
             color=colorRampPalette(c("#0f3460","white","#e94560"))(100),
             annotation_col=hd$annotation, fontsize_row=8, fontsize_col=10, fontsize=10,
             angle_col=45, main=hd$title, border_color=NA,
             treeheight_row=120, treeheight_col=80, annotation_legend=TRUE)
  }
  
  output$heatmapPlot <- renderPlot({ draw_pheatmap(heatmap_data()) }, height=1000)
  output$downloadHeatmapPDF <- downloadHandler(
    filename = function() paste0(input$cmp,"_heatmap.pdf"),
    content  = function(file) { hd <- isolate(heatmap_data()); pdf(file,width=12,height=14); draw_pheatmap(hd); dev.off() })
  output$downloadHeatmapPNG <- downloadHandler(
    filename = function() paste0(input$cmp,"_heatmap.png"),
    content  = function(file) { hd <- isolate(heatmap_data()); png(file,width=12,height=14,units="in",res=300); draw_pheatmap(hd); dev.off() })
  
  # ── GSEA — on-demand, cached per comparison ──────────────────────────────────
  gsea_triggered <- reactiveVal(NULL)
  
  observeEvent(input$main_tabs, {
    if (input$main_tabs %in% c("GSEA - GO BP", "GSEA - KEGG"))
      gsea_triggered(input$cmp)
    else
      gsea_triggered(NULL)
  })
  
  # only reset trigger when NOT on GSEA tab
  observeEvent(input$cmp, {
    if (!is.null(input$main_tabs) && !input$main_tabs %in% c("GSEA - GO BP", "GSEA - KEGG"))
      gsea_triggered(NULL)
  })
  
  gsea_for_cmp <- reactive({
    cmp_trigger <- gsea_triggered()
    req(rv$data, cmp_trigger)
    cmp      <- isolate(input$cmp)
    if (!is.null(rv$gsea_cache[[cmp]])) return(rv$gsea_cache[[cmp]])
    
    species   <- rv$data$metadata$species
    gsea_seed <- rv$data$metadata$gsea_seed %||% 123456
    if (is.null(species) || !species %in% names(SPECIES_MAP)) {
      showNotification("Species not found or not supported. Check the RData file.",
                       type="error", duration=8)
      return(list(GO_BP = NULL, KEGG = NULL))
    }
    species_info <- SPECIES_MAP[[species]]
    
    # lock comparison selector during GSEA
    rv$gsea_running <- TRUE
    on.exit(rv$gsea_running <- FALSE)
    
    withProgress(message=paste("Running GSEA for", cmp, "..."), value=0.1, {
      incProgress(0.4, detail="GO BP + KEGG analysis (~1-4 min)")
      result <- run_gsea_for_cmp(
        cmp             = cmp,
        df              = rv$data$metadata$DE_df,
        gene_symbol_col = rv$data$metadata$gene_cols$symbol,
        OrgDb           = species_info$OrgDb,
        KEGGOrg         = species_info$KEGGOrg,
        gsea_seed       = gsea_seed)
      incProgress(0.5, detail="Done")
    })
    
    rv$gsea_cache[[cmp]] <- result
    result
  })
  
  # ── GSEA helpers ─────────────────────────────────────────────────────────────
  make_gsea_table <- function(result_df, direction="up") {
    if (is.null(result_df)) return(data.frame(Message="No data available"))
    filtered <- if (direction=="up") result_df %>% dplyr::filter(NES > 0)
    else result_df %>% dplyr::filter(NES < 0)
    filtered <- filtered %>% dplyr::arrange(p.adjust)
    sig_col  <- if (sig_metric() == "fdr") filtered$p.adjust else filtered$pvalue
    sig      <- filtered[sig_col < fdr_cut(), ]
    non_sig  <- filtered[sig_col >= fdr_cut(), ]
    display  <- dplyr::bind_rows(sig, head(non_sig, max(0, 20 - nrow(sig))))
    display %>% transmute(Description,
                          `P value`    = signif(pvalue,   4),
                          `P adjusted` = signif(p.adjust, 4),
                          Significant  = ifelse(if (sig_metric()=="fdr") p.adjust < fdr_cut() else pvalue < fdr_cut(), "\u2713",""),
                          ES  = round(enrichmentScore, 3),
                          NES = round(NES, 3))
  }
  
  filter_to_genesets <- function(obj) {
    if (is.null(obj) || is.null(obj$gsea_object)) return(NULL)
    obj$result[obj$result$ID %in% names(obj$gsea_object@geneSets), ]
  }
  
  get_sig_gsea_obj <- function(gsea_obj) {
    if (is.null(gsea_obj)) return(NULL)
    sig_col    <- if (sig_metric() == "fdr") gsea_obj@result$p.adjust else gsea_obj@result$pvalue
    sig_result <- gsea_obj@result[sig_col < fdr_cut(), ]
    if (nrow(sig_result) == 0) return(NULL)
    valid_ids  <- intersect(sig_result$ID, names(gsea_obj@geneSets))
    sig_result <- sig_result[sig_result$ID %in% valid_ids, ]
    if (nrow(sig_result) == 0) return(NULL)
    slim <- gsea_obj
    slim@result   <- sig_result
    slim@geneSets <- gsea_obj@geneSets[valid_ids]
    slim
  }
  
  no_sig_plot <- function(msg="No significant pathways (FDR < 0.05) \u2014 no figure shown") {
    plot.new(); text(0.5, 0.5, msg, cex=1.2, col="#718096", font=3)
  }
  
  make_dotplot <- function(data_obj) {
    if (is.null(data_obj) || is.null(data_obj$gsea_object)) return(NULL)
    sig_obj <- get_sig_gsea_obj(data_obj$gsea_object)
    if (is.null(sig_obj)) return(NULL)
    dotplot(sig_obj, showCategory=10, split=".sign", font.size=11) + facet_grid(.~.sign)
  }
  
  render_enrichment_ui <- function(data_obj, prefix) {
    if (is.null(data_obj) || is.null(data_obj$gsea_object))
      return(p(class="no-sig-msg", "No GSEA data available for this comparison."))
    sig_obj <- get_sig_gsea_obj(data_obj$gsea_object)
    if (is.null(sig_obj))
      return(p(class="no-sig-msg", "No significant pathways (FDR < 0.05) \u2014 enrichment score plots not shown."))
    result_df     <- sig_obj@result
    upregulated   <- result_df %>% dplyr::filter(NES > 0) %>% dplyr::arrange(p.adjust)
    downregulated <- result_df %>% dplyr::filter(NES < 0) %>% dplyr::arrange(p.adjust)
    n_up   <- min(5, nrow(upregulated))
    n_down <- min(5, nrow(downregulated))
    ui_elements <- list()
    if (n_up > 0)
      ui_elements <- c(ui_elements, list(
        h4(paste0("Top ", n_up, " Upregulated (FDR < 0.05)")),
        lapply(seq_len(n_up), function(i) {
          pname <- paste0(prefix,"Up",i)
          tagList(fluidRow(column(12,
                                  downloadButton(paste0("dl_",pname,"_pdf"),"PDF",style="margin-bottom:4px;margin-right:4px;"),
                                  downloadButton(paste0("dl_",pname,"_png"),"PNG",style="margin-bottom:4px;")
          )), plotOutput(pname, height="400px"))
        })))
    if (n_down > 0)
      ui_elements <- c(ui_elements, list(
        hr(), h4(paste0("Top ", n_down, " Downregulated (FDR < 0.05)")),
        lapply(seq_len(n_down), function(i) {
          pname <- paste0(prefix,"Down",i)
          tagList(fluidRow(column(12,
                                  downloadButton(paste0("dl_",pname,"_pdf"),"PDF",style="margin-bottom:4px;margin-right:4px;"),
                                  downloadButton(paste0("dl_",pname,"_png"),"PNG",style="margin-bottom:4px;")
          )), plotOutput(pname, height="400px"))
        })))
    if (length(ui_elements) == 0) return(p(class="no-sig-msg","No terms to plot."))
    do.call(tagList, ui_elements)
  }
  
  render_enrichment_plots <- function(data_obj, prefix, output, cmp_label) {
    if (is.null(data_obj) || is.null(data_obj$gsea_object)) return()
    sig_obj <- get_sig_gsea_obj(data_obj$gsea_object)
    if (is.null(sig_obj)) return()
    result_df <- sig_obj@result
    up_df <- result_df %>% dplyr::filter(NES > 0) %>% dplyr::arrange(p.adjust)
    dn_df <- result_df %>% dplyr::filter(NES < 0) %>% dplyr::arrange(p.adjust)
    register_enrichment <- function(pname, term_id, title_text) {
      output[[pname]] <- renderPlot({ gseaplot(sig_obj,by="all",title=title_text,geneSetID=term_id,font.size=14) })
      output[[paste0("dl_",pname,"_pdf")]] <- downloadHandler(
        filename=function() paste0(cmp_label,"_",pname,".pdf"),
        content =function(file) { pdf(file,width=10,height=6); print(gseaplot(sig_obj,by="all",title=title_text,geneSetID=term_id,font.size=14)); dev.off() })
      output[[paste0("dl_",pname,"_png")]] <- downloadHandler(
        filename=function() paste0(cmp_label,"_",pname,".png"),
        content =function(file) { png(file,width=10,height=6,units="in",res=300); print(gseaplot(sig_obj,by="all",title=title_text,geneSetID=term_id,font.size=14)); dev.off() })
    }
    for (i in seq_len(min(5,nrow(up_df)))) {
      local({ my_i <- i; term_id <- which(result_df$ID==up_df$ID[my_i]); pname <- paste0(prefix,"Up",my_i)
      register_enrichment(pname, term_id, paste("UP:", result_df$Description[term_id])) })
    }
    for (i in seq_len(min(5,nrow(dn_df)))) {
      local({ my_i <- i; term_id <- which(result_df$ID==dn_df$ID[my_i]); pname <- paste0(prefix,"Down",my_i)
      register_enrichment(pname, term_id, paste("DOWN:", result_df$Description[term_id])) })
    }
  }
  
  # ── GO BP tab ────────────────────────────────────────────────────────────────
  output$gsea_go_content <- renderUI({
    req(rv$data, input$cmp)
    if (is.null(rv$data$metadata$species) || !rv$data$metadata$species %in% names(SPECIES_MAP))
      return(div(style="padding:20px;",
                 p(style="color:#e94560;font-weight:600;font-size:14px;",
                   "\u26a0 Species not found in RData. Please ensure the RData file contains a valid species variable.")))
    tagList(
      h4(textOutput("goTitle"), style="margin-top:8px;"), hr(),
      uiOutput("goUpTitle"),
      downloadButton("downloadGoUpCSV","\u2193 CSV",class="btn-dl btn-dl-csv"),
      DTOutput("goUpregulatedTable"), hr(),
      uiOutput("goDownTitle"),
      downloadButton("downloadGoDownCSV","\u2193 CSV",class="btn-dl btn-dl-csv"),
      DTOutput("goDownregulatedTable"), hr(),
      div(class="tab-section-title","Dot Plot \u2014 significant terms only"),
      fluidRow(column(12,
                      downloadButton("downloadGoDotPDF","\u2193 PDF",class="btn-dl btn-dl-pdf"),
                      downloadButton("downloadGoDotPNG","\u2193 PNG",class="btn-dl btn-dl-png")
      )),
      plotOutput("goDotPlot", height="720px"), hr(),
      div(class="tab-section-title","Enrichment Score Plots \u2014 significant terms only"),
      uiOutput("goEnrichmentPlots")
    )
  })
  
  # ── KEGG tab ─────────────────────────────────────────────────────────────────
  output$gsea_kegg_content <- renderUI({
    req(rv$data, input$cmp)
    if (is.null(rv$data$metadata$species) || !rv$data$metadata$species %in% names(SPECIES_MAP))
      return(div(style="padding:20px;",
                 p(style="color:#e94560;font-weight:600;font-size:14px;",
                   "\u26a0 Species not found in RData. Please ensure the RData file contains a valid species variable.")))
    tagList(
      h4(textOutput("keggTitle"), style="margin-top:8px;"), hr(),
      uiOutput("keggUpTitle"),
      downloadButton("downloadKeggUpCSV","\u2193 CSV",class="btn-dl btn-dl-csv"),
      DTOutput("keggUpregulatedTable"), hr(),
      uiOutput("keggDownTitle"),
      downloadButton("downloadKeggDownCSV","\u2193 CSV",class="btn-dl btn-dl-csv"),
      DTOutput("keggDownregulatedTable"), hr(),
      div(class="tab-section-title","Dot Plot \u2014 significant terms only"),
      fluidRow(column(12,
                      downloadButton("downloadKeggDotPDF","\u2193 PDF",class="btn-dl btn-dl-pdf"),
                      downloadButton("downloadKeggDotPNG","\u2193 PNG",class="btn-dl btn-dl-png")
      )),
      plotOutput("keggDotPlot", height="720px"), hr(),
      div(class="tab-section-title","Enrichment Score Plots \u2014 significant terms only"),
      uiOutput("keggEnrichmentPlots")
    )
  })
  
  gsea_title_ui <- function(direction, type) {
    renderUI({
      metric <- if (sig_metric() == "fdr") "FDR" else "P-value"
      div(class = "tab-section-title",
          paste0(direction, " ", type, " \u2014 ", metric, " < ", fdr_cut(),
                 " (all significant shown; remaining filled to 20)"))
    })
  }
  output$goUpTitle     <- gsea_title_ui("Upregulated",   "GO BP Terms")
  output$goDownTitle   <- gsea_title_ui("Downregulated", "GO BP Terms")
  output$keggUpTitle   <- gsea_title_ui("Upregulated",   "KEGG Pathways")
  output$keggDownTitle <- gsea_title_ui("Downregulated", "KEGG Pathways")
  
  # ── GO BP outputs ─────────────────────────────────────────────────────────────
  output$goTitle <- renderText({
    req(rv$data, input$cmp)
    groups <- extract_groups_from_comparison(input$cmp, rv$data$metadata$comparison_info)
    if (!is.null(groups) && length(groups) >= 2) paste("GO BP GSEA:", groups[1], "vs", groups[2])
    else paste("GO BP GSEA:", input$cmp)
  })
  
  go_up_data   <- reactive({ req(rv$data,input$cmp); make_gsea_table(filter_to_genesets(gsea_for_cmp()$GO_BP),  "up") })
  go_down_data <- reactive({ req(rv$data,input$cmp); make_gsea_table(filter_to_genesets(gsea_for_cmp()$GO_BP),  "down") })
  
  output$goUpregulatedTable <- renderDT({
    tbl <- go_up_data()
    dt  <- datatable(tbl, rownames=FALSE, options=list(scrollX=TRUE,pageLength=20,order=list(list(2,"asc"))))
    if ("P value" %in% colnames(tbl)) dt <- dt %>% formatSignif(columns=c("P value","P adjusted"),digits=3)
    dt
  })
  output$goDownregulatedTable <- renderDT({
    tbl <- go_down_data()
    dt  <- datatable(tbl, rownames=FALSE, options=list(scrollX=TRUE,pageLength=20,order=list(list(2,"asc"))))
    if ("P value" %in% colnames(tbl)) dt <- dt %>% formatSignif(columns=c("P value","P adjusted"),digits=3)
    dt
  })
  output$downloadGoUpCSV   <- downloadHandler(filename=function() paste0(input$cmp,"_GO_up.csv"),   content=function(file) write.csv(go_up_data(),  file,row.names=FALSE))
  output$downloadGoDownCSV <- downloadHandler(filename=function() paste0(input$cmp,"_GO_down.csv"), content=function(file) write.csv(go_down_data(),file,row.names=FALSE))
  output$goDotPlot <- renderPlot({
    req(rv$data,input$cmp); p <- make_dotplot(gsea_for_cmp()$GO_BP)
    if (is.null(p)) no_sig_plot() else print(p)
  })
  output$downloadGoDotPDF <- downloadHandler(filename=function() paste0(input$cmp,"_GO_dotplot.pdf"), content=function(file) { p <- make_dotplot(gsea_for_cmp()$GO_BP); req(!is.null(p)); ggplot2::ggsave(file,plot=p,device="pdf",width=10,height=7) })
  output$downloadGoDotPNG <- downloadHandler(filename=function() paste0(input$cmp,"_GO_dotplot.png"), content=function(file) { p <- make_dotplot(gsea_for_cmp()$GO_BP); req(!is.null(p)); ggplot2::ggsave(file,plot=p,device="png",width=10,height=7,dpi=300) })
  output$goEnrichmentPlots <- renderUI({ req(rv$data,input$cmp); render_enrichment_ui(gsea_for_cmp()$GO_BP,"goEnrichPlot") })
  observe({ req(rv$data,input$cmp); render_enrichment_plots(gsea_for_cmp()$GO_BP,"goEnrichPlot",output,paste0(input$cmp,"_GO")) })
  
  # ── KEGG outputs ──────────────────────────────────────────────────────────────
  output$keggTitle <- renderText({
    req(rv$data, input$cmp)
    groups <- extract_groups_from_comparison(input$cmp, rv$data$metadata$comparison_info)
    if (!is.null(groups) && length(groups) >= 2) paste("KEGG GSEA:", groups[1], "vs", groups[2])
    else paste("KEGG GSEA:", input$cmp)
  })
  
  kegg_up_data   <- reactive({ req(rv$data,input$cmp); make_gsea_table(filter_to_genesets(gsea_for_cmp()$KEGG), "up") })
  kegg_down_data <- reactive({ req(rv$data,input$cmp); make_gsea_table(filter_to_genesets(gsea_for_cmp()$KEGG), "down") })
  
  output$keggUpregulatedTable <- renderDT({
    tbl <- kegg_up_data()
    dt  <- datatable(tbl, rownames=FALSE, options=list(scrollX=TRUE,pageLength=20,order=list(list(2,"asc"))))
    if ("P value" %in% colnames(tbl)) dt <- dt %>% formatSignif(columns=c("P value","P adjusted"),digits=3)
    dt
  })
  output$keggDownregulatedTable <- renderDT({
    tbl <- kegg_down_data()
    dt  <- datatable(tbl, rownames=FALSE, options=list(scrollX=TRUE,pageLength=20,order=list(list(2,"asc"))))
    if ("P value" %in% colnames(tbl)) dt <- dt %>% formatSignif(columns=c("P value","P adjusted"),digits=3)
    dt
  })
  output$downloadKeggUpCSV   <- downloadHandler(filename=function() paste0(input$cmp,"_KEGG_up.csv"),   content=function(file) write.csv(kegg_up_data(),  file,row.names=FALSE))
  output$downloadKeggDownCSV <- downloadHandler(filename=function() paste0(input$cmp,"_KEGG_down.csv"), content=function(file) write.csv(kegg_down_data(),file,row.names=FALSE))
  output$keggDotPlot <- renderPlot({
    req(rv$data,input$cmp); p <- make_dotplot(gsea_for_cmp()$KEGG)
    if (is.null(p)) no_sig_plot() else print(p)
  })
  output$downloadKeggDotPDF <- downloadHandler(filename=function() paste0(input$cmp,"_KEGG_dotplot.pdf"), content=function(file) { p <- make_dotplot(gsea_for_cmp()$KEGG); req(!is.null(p)); ggplot2::ggsave(file,plot=p,device="pdf",width=10,height=7) })
  output$downloadKeggDotPNG <- downloadHandler(filename=function() paste0(input$cmp,"_KEGG_dotplot.png"), content=function(file) { p <- make_dotplot(gsea_for_cmp()$KEGG); req(!is.null(p)); ggplot2::ggsave(file,plot=p,device="png",width=10,height=7,dpi=300) })
  output$keggEnrichmentPlots <- renderUI({ req(rv$data,input$cmp); render_enrichment_ui(gsea_for_cmp()$KEGG,"keggEnrichPlot") })
  observe({ req(rv$data,input$cmp); render_enrichment_plots(gsea_for_cmp()$KEGG,"keggEnrichPlot",output,paste0(input$cmp,"_KEGG")) })
}

shinyApp(ui, server)
