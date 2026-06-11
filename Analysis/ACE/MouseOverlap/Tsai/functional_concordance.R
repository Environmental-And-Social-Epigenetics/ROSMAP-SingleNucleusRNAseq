#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(WebGestaltR)
  library(ggplot2)
})

split_csv <- function(x) {
  x <- trimws(unlist(strsplit(x, ",", fixed = TRUE)))
  x[nzchar(x)]
}

parse_args <- function(argv) {
  out <- list(
    mouse_deg = NULL,
    human_deg = NULL,
    ortholog_table = NULL,
    output_dir = NULL,
    figures_dir = NULL,
    databases = "geneontology_Biological_Process_noRedundant,pathway_Reactome",
    human_targets = "Exc,Inh,In-PV_Basket,In-PV_Chandelier",
    phenotypes = "tot_adverse_exp,ace_aggregate,early_hh_ses",
    integrations = "derived_batch,projid",
    sexes = "Male,Female",
    custom_gmt = "",
    permutations = 1000L,
    min_size = 10L,
    max_size = 500L,
    fdr_threshold = 1,
    alpha = 0.05,
    suggestive_alpha = 0.2,
    threads = 1L,
    seed = 1L,
    skip_existing = FALSE
  )

  key_map <- c(
    "--mouse-deg" = "mouse_deg",
    "--human-deg" = "human_deg",
    "--ortholog-table" = "ortholog_table",
    "--output-dir" = "output_dir",
    "--figures-dir" = "figures_dir",
    "--databases" = "databases",
    "--human-targets" = "human_targets",
    "--phenotypes" = "phenotypes",
    "--integrations" = "integrations",
    "--sexes" = "sexes",
    "--custom-gmt" = "custom_gmt",
    "--permutations" = "permutations",
    "--min-size" = "min_size",
    "--max-size" = "max_size",
    "--fdr-threshold" = "fdr_threshold",
    "--alpha" = "alpha",
    "--suggestive-alpha" = "suggestive_alpha",
    "--threads" = "threads",
    "--seed" = "seed"
  )

  i <- 1
  while (i <= length(argv)) {
    arg <- argv[[i]]
    if (arg == "--skip-existing") {
      out$skip_existing <- TRUE
      i <- i + 1
    } else if (arg %in% names(key_map)) {
      if (i == length(argv)) stop("Missing value for ", arg)
      out[[key_map[[arg]]]] <- argv[[i + 1]]
      i <- i + 2
    } else {
      stop("Unknown argument: ", arg)
    }
  }

  integer_fields <- c("permutations", "min_size", "max_size", "threads", "seed")
  for (field in integer_fields) out[[field]] <- as.integer(out[[field]])
  numeric_fields <- c("fdr_threshold", "alpha", "suggestive_alpha")
  for (field in numeric_fields) out[[field]] <- as.numeric(out[[field]])

  required <- c("mouse_deg", "human_deg", "ortholog_table", "output_dir", "figures_dir")
  missing <- required[vapply(required, function(x) is.null(out[[x]]) || !nzchar(out[[x]]), logical(1))]
  if (length(missing)) stop("Missing required argument(s): ", paste(missing, collapse = ", "))
  out
}

db_short_name <- function(db) {
  db <- sub("^geneontology_", "", db)
  db <- sub("^pathway_", "", db)
  db <- sub("^network_", "", db)
  db
}

sanitize <- function(x) {
  x <- gsub("[^A-Za-z0-9._-]+", "_", x)
  gsub("_+", "_", x)
}

read_table <- function(path) {
  read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}

num <- function(x) suppressWarnings(as.numeric(x))

safe_signed_rank <- function(stat, lfc, pvalue) {
  stat <- num(stat)
  lfc <- num(lfc)
  pvalue <- num(pvalue)
  pvalue[!is.finite(pvalue) | pvalue <= 0] <- .Machine$double.xmin
  fallback <- sign(lfc) * -log10(pvalue)
  rank <- ifelse(is.finite(stat), stat, fallback)
  rank[!is.finite(rank)] <- NA_real_
  rank
}

deduplicate_best <- function(df, group_cols, padj_col, rank_col) {
  if (!nrow(df)) return(df)
  key <- do.call(paste, c(df[group_cols], sep = "\r"))
  padj <- num(df[[padj_col]])
  padj[!is.finite(padj)] <- Inf
  abs_rank <- abs(num(df[[rank_col]]))
  abs_rank[!is.finite(abs_rank)] <- -Inf
  df$.__key <- key
  df$.__padj_sort <- padj
  df$.__abs_rank_sort <- abs_rank
  df <- df[order(df$.__key, df$.__padj_sort, -df$.__abs_rank_sort), , drop = FALSE]
  df <- df[!duplicated(df$.__key), , drop = FALSE]
  df[, setdiff(names(df), c(".__key", ".__padj_sort", ".__abs_rank_sort")), drop = FALSE]
}

empty_gsea <- function() {
  data.frame(
    species = character(), analysis_group_id = character(),
    mouse_dataset = character(), mouse_population = character(),
    mouse_condition = character(), human_integration = character(),
    human_phenotype = character(), human_cell_type = character(),
    human_sex = character(), database = character(), database_source = character(),
    geneSet = character(), description = character(), enrichmentScore = numeric(),
    NES = numeric(), pValue = numeric(), FDR = numeric(), size = numeric(),
    leadingEdgeNum = numeric(), leadingEdgeId = character(), userId = character(),
    n_ranked_genes = integer(), stringsAsFactors = FALSE
  )
}

empty_summary <- function() {
  data.frame(
    comparison_id = character(), comparison_role = character(),
    mouse_dataset = character(), mouse_population = character(),
    mouse_condition = character(), human_integration = character(),
    human_phenotype = character(), human_cell_type = character(),
    human_sex = character(), shared_pathway_n = integer(),
    shared_go_bp_pathway_n = integer(), shared_reactome_pathway_n = integer(),
    spearman_r = numeric(), spearman_pvalue = numeric(), pearson_r = numeric(),
    pearson_pvalue = numeric(), sign_concordance = numeric(),
    concordant_sig_n = integer(), discordant_sig_n = integer(),
    concordant_suggestive_n = integer(), discordant_suggestive_n = integer(),
    lnb_up_ace_up_pathway_n = integer(), lnb_down_ace_down_pathway_n = integer(),
    lnb_up_ace_down_pathway_n = integer(), lnb_down_ace_up_pathway_n = integer(),
    lnb_up_ace_up_fisher_pvalue = numeric(),
    lnb_down_ace_down_fisher_pvalue = numeric(),
    lnb_up_ace_down_fisher_pvalue = numeric(),
    lnb_down_ace_up_fisher_pvalue = numeric(),
    best_concordant_fisher_pvalue = numeric(),
    best_discordant_fisher_pvalue = numeric(),
    stringsAsFactors = FALSE
  )
}

empty_pathways <- function() {
  data.frame(
    comparison_id = character(), comparison_role = character(),
    mouse_dataset = character(), mouse_population = character(),
    mouse_condition = character(), human_integration = character(),
    human_phenotype = character(), human_cell_type = character(),
    human_sex = character(), database = character(), geneSet = character(),
    description = character(), mouse_NES = numeric(), human_NES = numeric(),
    mouse_FDR = numeric(), human_FDR = numeric(),
    mouse_pValue = numeric(), human_pValue = numeric(),
    mouse_size = numeric(), human_size = numeric(),
    mouse_leadingEdgeNum = numeric(), human_leadingEdgeNum = numeric(),
    mouse_userId = character(), human_userId = character(),
    mouse_pathway_direction = character(), human_pathway_direction = character(),
    is_concordant = logical(), is_discordant = logical(),
    mouse_sig = logical(), human_sig = logical(), both_sig = logical(),
    both_suggestive = logical(), pathway_status = character(),
    mean_abs_NES = numeric(), max_pair_FDR = numeric(), stringsAsFactors = FALSE
  )
}

prepare_ranked_genes <- function(args) {
  mouse <- read_table(args$mouse_deg)
  human <- read_table(args$human_deg)
  orth <- read_table(args$ortholog_table)

  required_mouse <- c("mouse_dataset", "mouse_population", "condition", "mouse_symbol",
                      "mouse_stat", "mouse_lnb_log2fc", "mouse_pvalue", "mouse_padj")
  required_human <- c("integration", "phenotype", "cell_type", "sex", "gene_symbol",
                      "stat", "log2FoldChange", "pvalue", "padj")
  required_orth <- c("mouse_symbol", "human_symbol", "is_one_to_one")
  for (col in required_mouse) if (!col %in% names(mouse)) stop("Mouse table missing column: ", col)
  for (col in required_human) if (!col %in% names(human)) stop("Human table missing column: ", col)
  for (col in required_orth) if (!col %in% names(orth)) stop("Ortholog table missing column: ", col)

  orth <- orth[tolower(as.character(orth$is_one_to_one)) %in% c("true", "1"), , drop = FALSE]
  orth <- orth[!duplicated(orth[, c("mouse_symbol", "human_symbol")]), , drop = FALSE]

  mouse <- merge(mouse, orth[, c("mouse_symbol", "human_symbol")], by = "mouse_symbol", all = FALSE)
  mouse$rank <- safe_signed_rank(-num(mouse$mouse_stat), mouse$mouse_lnb_log2fc, mouse$mouse_pvalue)
  mouse$source_stat <- num(mouse$mouse_stat)
  mouse$signed_stat <- -num(mouse$mouse_stat)
  mouse <- mouse[is.finite(mouse$rank), , drop = FALSE]
  mouse <- deduplicate_best(mouse, c("mouse_dataset", "human_symbol"), "mouse_padj", "rank")

  mouse_ranked <- data.frame(
    species = "mouse",
    analysis_group_id = paste("mouse", mouse$mouse_dataset, sep = "|"),
    mouse_dataset = mouse$mouse_dataset,
    mouse_population = mouse$mouse_population,
    mouse_condition = mouse$condition,
    human_integration = NA_character_,
    human_phenotype = NA_character_,
    human_cell_type = NA_character_,
    human_sex = NA_character_,
    human_symbol = mouse$human_symbol,
    source_symbol = mouse$mouse_symbol,
    rank = mouse$rank,
    source_stat = mouse$source_stat,
    signed_stat = mouse$signed_stat,
    log2fc = num(mouse$mouse_lnb_log2fc),
    pvalue = num(mouse$mouse_pvalue),
    padj = num(mouse$mouse_padj),
    stringsAsFactors = FALSE
  )

  human$human_sex <- ifelse(human$sex %in% c("Fem", "Female"), "Female", human$sex)
  human <- human[human$cell_type %in% split_csv(args$human_targets), , drop = FALSE]
  human <- human[human$phenotype %in% split_csv(args$phenotypes), , drop = FALSE]
  human <- human[human$integration %in% split_csv(args$integrations), , drop = FALSE]
  human <- human[human$human_sex %in% split_csv(args$sexes), , drop = FALSE]
  human$human_symbol <- human$gene_symbol
  human$rank <- safe_signed_rank(human$stat, human$log2FoldChange, human$pvalue)
  human$source_stat <- num(human$stat)
  human$signed_stat <- num(human$stat)
  human <- human[is.finite(human$rank), , drop = FALSE]
  human <- deduplicate_best(
    human,
    c("integration", "phenotype", "cell_type", "human_sex", "human_symbol"),
    "padj",
    "rank"
  )

  human_ranked <- data.frame(
    species = "human",
    analysis_group_id = paste("human", human$integration, human$phenotype,
                              human$cell_type, human$human_sex, sep = "|"),
    mouse_dataset = NA_character_,
    mouse_population = NA_character_,
    mouse_condition = NA_character_,
    human_integration = human$integration,
    human_phenotype = human$phenotype,
    human_cell_type = human$cell_type,
    human_sex = human$human_sex,
    human_symbol = human$human_symbol,
    source_symbol = human$gene_symbol,
    rank = human$rank,
    source_stat = human$source_stat,
    signed_stat = human$signed_stat,
    log2fc = num(human$log2FoldChange),
    pvalue = num(human$pvalue),
    padj = num(human$padj),
    stringsAsFactors = FALSE
  )

  ranked <- rbind(mouse_ranked, human_ranked)
  ranked <- ranked[order(ranked$species, ranked$analysis_group_id, -ranked$rank), , drop = FALSE]
  ranked
}

database_table <- function(args) {
  if (nzchar(args$custom_gmt)) {
    return(data.frame(database = "custom", database_source = args$custom_gmt, stringsAsFactors = FALSE))
  }
  db <- split_csv(args$databases)
  data.frame(database = db_short_name(db), database_source = db, stringsAsFactors = FALSE)
}

normalize_webgestalt_result <- function(result) {
  if (is.null(result) || !nrow(result)) return(data.frame())
  get_col <- function(name, default = NA) {
    if (name %in% names(result)) result[[name]] else rep(default, nrow(result))
  }
  data.frame(
    geneSet = as.character(get_col("geneSet")),
    description = as.character(get_col("description")),
    enrichmentScore = num(get_col("enrichmentScore")),
    NES = num(get_col("normalizedEnrichmentScore")),
    pValue = num(get_col("pValue")),
    FDR = num(get_col("FDR")),
    size = num(get_col("size")),
    leadingEdgeNum = num(get_col("leadingEdgeNum")),
    leadingEdgeId = as.character(get_col("leadingEdgeId", "")),
    userId = as.character(get_col("userId", "")),
    stringsAsFactors = FALSE
  )
}

run_gsea_for_group <- function(args, rank_group, db_row, cache_dir) {
  group_id <- unique(rank_group$analysis_group_id)
  if (length(group_id) != 1) stop("Rank group has multiple IDs")
  cache_path <- file.path(cache_dir, paste0(sanitize(group_id), "__", sanitize(db_row$database), ".rds"))

  if (args$skip_existing && file.exists(cache_path) && file.info(cache_path)$size > 0) {
    result <- readRDS(cache_path)
  } else {
    interest <- data.frame(gene = rank_group$human_symbol, rank = num(rank_group$rank), stringsAsFactors = FALSE)
    interest <- interest[is.finite(interest$rank) & nzchar(interest$gene), , drop = FALSE]
    interest <- interest[order(-interest$rank), , drop = FALSE]

    if (nrow(interest) < args$min_size) {
      result <- data.frame()
    } else {
      result <- tryCatch({
        if (nzchar(args$custom_gmt)) {
          WebGestaltR(
            enrichMethod = "GSEA",
            organism = "others",
            enrichDatabaseFile = db_row$database_source,
            interestGene = interest,
            minNum = args$min_size,
            maxNum = args$max_size,
            sigMethod = "fdr",
            fdrThr = args$fdr_threshold,
            perNum = args$permutations,
            isOutput = FALSE,
            nThreads = args$threads
          )
        } else {
          WebGestaltR(
            enrichMethod = "GSEA",
            organism = "hsapiens",
            enrichDatabase = db_row$database_source,
            enrichDatabaseType = "genesymbol",
            interestGene = interest,
            interestGeneType = "genesymbol",
            minNum = args$min_size,
            maxNum = args$max_size,
            sigMethod = "fdr",
            fdrThr = args$fdr_threshold,
            perNum = args$permutations,
            isOutput = FALSE,
            nThreads = args$threads,
            cache = file.path(args$output_dir, "webgestalt_cache")
          )
        }
      }, error = function(e) {
        warning("WebGestaltR failed for ", group_id, " / ", db_row$database, ": ", conditionMessage(e))
        data.frame()
      })
    }
    saveRDS(result, cache_path)
  }

  result <- normalize_webgestalt_result(result)
  if (!nrow(result)) return(empty_gsea())

  meta <- rank_group[1, c("species", "analysis_group_id", "mouse_dataset", "mouse_population",
                         "mouse_condition", "human_integration", "human_phenotype",
                         "human_cell_type", "human_sex"), drop = FALSE]
  meta <- meta[rep(1, nrow(result)), , drop = FALSE]
  result <- cbind(meta, database = db_row$database, database_source = db_row$database_source,
                  result, n_ranked_genes = length(unique(rank_group$human_symbol)))
  result
}

run_all_gsea <- function(args, ranked) {
  cache_dir <- file.path(args$output_dir, "gsea_cache")
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(args$output_dir, "webgestalt_cache"), recursive = TRUE, showWarnings = FALSE)

  dbs <- database_table(args)
  groups <- split(ranked, ranked$analysis_group_id)
  out <- list()
  idx <- 1L
  for (group_id in names(groups)) {
    cat("GSEA group:", group_id, "(", nrow(groups[[group_id]]), "ranked genes )\n")
    for (i in seq_len(nrow(dbs))) {
      cat("  database:", dbs$database[[i]], "\n")
      out[[idx]] <- run_gsea_for_group(args, groups[[group_id]], dbs[i, ], cache_dir)
      idx <- idx + 1L
    }
  }
  out <- out[vapply(out, nrow, integer(1)) > 0]
  if (!length(out)) return(empty_gsea())
  do.call(rbind, out)
}

target_match <- function(mouse_population, human_cell_type) {
  if (mouse_population == "NeuN") return(human_cell_type == "Exc")
  if (mouse_population == "PV") return(human_cell_type %in% c("Inh", "In-PV_Basket", "In-PV_Chandelier"))
  FALSE
}

comparison_role <- function(mouse_population, human_cell_type, human_integration, human_phenotype, human_sex) {
  first_line <- human_integration == "derived_batch" && human_phenotype == "tot_adverse_exp" && human_sex == "Male"
  if (first_line && mouse_population == "NeuN" && human_cell_type == "Exc") return("primary")
  if (first_line && mouse_population == "PV" && human_cell_type == "Inh") return("primary")
  if (first_line && mouse_population == "PV" && human_cell_type %in% c("In-PV_Basket", "In-PV_Chandelier")) {
    return("pv_sensitivity")
  }
  if (human_sex == "Female") return("female_sensitivity")
  "secondary_sensitivity"
}

direction_label <- function(nes, up_label, down_label) {
  ifelse(nes > 0, up_label, ifelse(nes < 0, down_label, "zero"))
}

fisher_direction <- function(df, mouse_direction, human_direction, alpha) {
  universe <- unique(df$geneSet)
  mouse_set <- unique(df$geneSet[df$mouse_sig & df$mouse_pathway_direction == mouse_direction])
  human_set <- unique(df$geneSet[df$human_sig & df$human_pathway_direction == human_direction])
  a <- length(intersect(mouse_set, human_set))
  b <- length(setdiff(mouse_set, human_set))
  c <- length(setdiff(human_set, mouse_set))
  d <- max(length(universe) - a - b - c, 0)
  p <- if (length(universe) > 0) stats::fisher.test(matrix(c(a, b, c, d), nrow = 2), alternative = "greater")$p.value else NA_real_
  c(n = a, p = p)
}

cor_or_na <- function(x, y, method) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 3) return(c(r = NA_real_, p = NA_real_))
  test <- suppressWarnings(stats::cor.test(x[ok], y[ok], method = method, exact = FALSE))
  c(r = unname(test$estimate), p = test$p.value)
}

build_concordance <- function(args, gsea) {
  if (!nrow(gsea)) return(list(summary = empty_summary(), pathways = empty_pathways()))

  mouse_meta <- unique(gsea[gsea$species == "mouse", c("analysis_group_id", "mouse_dataset",
                                                       "mouse_population", "mouse_condition"), drop = FALSE])
  human_meta <- unique(gsea[gsea$species == "human", c("analysis_group_id", "human_integration",
                                                       "human_phenotype", "human_cell_type", "human_sex"), drop = FALSE])
  summary_rows <- list()
  pathway_rows <- list()
  si <- 1L
  pi <- 1L

  for (m_i in seq_len(nrow(mouse_meta))) {
    m <- mouse_meta[m_i, , drop = FALSE]
    for (h_i in seq_len(nrow(human_meta))) {
      h <- human_meta[h_i, , drop = FALSE]
      if (!target_match(m$mouse_population, h$human_cell_type)) next

      mres <- gsea[gsea$analysis_group_id == m$analysis_group_id, , drop = FALSE]
      hres <- gsea[gsea$analysis_group_id == h$analysis_group_id, , drop = FALSE]
      joined <- merge(mres, hres, by = c("database", "geneSet"), suffixes = c("_mouse", "_human"))
      joined <- joined[is.finite(joined$NES_mouse) & is.finite(joined$NES_human), , drop = FALSE]
      if (!nrow(joined)) next

      role <- comparison_role(m$mouse_population, h$human_cell_type, h$human_integration, h$human_phenotype, h$human_sex)
      comparison_id <- paste(m$mouse_dataset, h$human_integration, h$human_phenotype,
                             h$human_cell_type, h$human_sex, sep = "|")

      pathway_df <- data.frame(
        comparison_id = comparison_id,
        comparison_role = role,
        mouse_dataset = m$mouse_dataset,
        mouse_population = m$mouse_population,
        mouse_condition = m$mouse_condition,
        human_integration = h$human_integration,
        human_phenotype = h$human_phenotype,
        human_cell_type = h$human_cell_type,
        human_sex = h$human_sex,
        database = joined$database,
        geneSet = joined$geneSet,
        description = ifelse(nzchar(joined$description_mouse), joined$description_mouse, joined$description_human),
        mouse_NES = joined$NES_mouse,
        human_NES = joined$NES_human,
        mouse_FDR = joined$FDR_mouse,
        human_FDR = joined$FDR_human,
        mouse_pValue = joined$pValue_mouse,
        human_pValue = joined$pValue_human,
        mouse_size = joined$size_mouse,
        human_size = joined$size_human,
        mouse_leadingEdgeNum = joined$leadingEdgeNum_mouse,
        human_leadingEdgeNum = joined$leadingEdgeNum_human,
        mouse_userId = joined$userId_mouse,
        human_userId = joined$userId_human,
        stringsAsFactors = FALSE
      )
      pathway_df$mouse_pathway_direction <- direction_label(pathway_df$mouse_NES, "LNB_up", "LNB_down")
      pathway_df$human_pathway_direction <- direction_label(pathway_df$human_NES, "ACE_up", "ACE_down")
      pathway_df$is_concordant <- pathway_df$mouse_pathway_direction %in% c("LNB_up", "LNB_down") &
        pathway_df$human_pathway_direction %in% c("ACE_up", "ACE_down") &
        ((pathway_df$mouse_pathway_direction == "LNB_up" & pathway_df$human_pathway_direction == "ACE_up") |
           (pathway_df$mouse_pathway_direction == "LNB_down" & pathway_df$human_pathway_direction == "ACE_down"))
      pathway_df$is_discordant <- pathway_df$mouse_pathway_direction %in% c("LNB_up", "LNB_down") &
        pathway_df$human_pathway_direction %in% c("ACE_up", "ACE_down") & !pathway_df$is_concordant
      pathway_df$mouse_sig <- is.finite(pathway_df$mouse_FDR) & pathway_df$mouse_FDR < args$alpha
      pathway_df$human_sig <- is.finite(pathway_df$human_FDR) & pathway_df$human_FDR < args$alpha
      pathway_df$both_sig <- pathway_df$mouse_sig & pathway_df$human_sig
      pathway_df$both_suggestive <- is.finite(pathway_df$mouse_FDR) & is.finite(pathway_df$human_FDR) &
        pathway_df$mouse_FDR < args$suggestive_alpha & pathway_df$human_FDR < args$suggestive_alpha
      pathway_df$mean_abs_NES <- (abs(pathway_df$mouse_NES) + abs(pathway_df$human_NES)) / 2
      pathway_df$max_pair_FDR <- pmax(pathway_df$mouse_FDR, pathway_df$human_FDR)
      pathway_df$pathway_status <- ifelse(pathway_df$both_sig & pathway_df$is_concordant, "concordant_sig",
        ifelse(pathway_df$both_sig & pathway_df$is_discordant, "discordant_sig",
          ifelse(pathway_df$both_suggestive & pathway_df$is_concordant, "concordant_suggestive",
            ifelse(pathway_df$both_suggestive & pathway_df$is_discordant, "discordant_suggestive", "other"))))

      spear <- cor_or_na(pathway_df$mouse_NES, pathway_df$human_NES, "spearman")
      pear <- cor_or_na(pathway_df$mouse_NES, pathway_df$human_NES, "pearson")
      sign_conc <- mean(pathway_df$is_concordant, na.rm = TRUE)
      up_up <- fisher_direction(pathway_df, "LNB_up", "ACE_up", args$alpha)
      down_down <- fisher_direction(pathway_df, "LNB_down", "ACE_down", args$alpha)
      up_down <- fisher_direction(pathway_df, "LNB_up", "ACE_down", args$alpha)
      down_up <- fisher_direction(pathway_df, "LNB_down", "ACE_up", args$alpha)

      summary_rows[[si]] <- data.frame(
        comparison_id = comparison_id,
        comparison_role = role,
        mouse_dataset = m$mouse_dataset,
        mouse_population = m$mouse_population,
        mouse_condition = m$mouse_condition,
        human_integration = h$human_integration,
        human_phenotype = h$human_phenotype,
        human_cell_type = h$human_cell_type,
        human_sex = h$human_sex,
        shared_pathway_n = nrow(pathway_df),
        shared_go_bp_pathway_n = sum(pathway_df$database == "Biological_Process_noRedundant"),
        shared_reactome_pathway_n = sum(pathway_df$database == "Reactome"),
        spearman_r = spear[["r"]],
        spearman_pvalue = spear[["p"]],
        pearson_r = pear[["r"]],
        pearson_pvalue = pear[["p"]],
        sign_concordance = sign_conc,
        concordant_sig_n = sum(pathway_df$both_sig & pathway_df$is_concordant),
        discordant_sig_n = sum(pathway_df$both_sig & pathway_df$is_discordant),
        concordant_suggestive_n = sum(pathway_df$both_suggestive & pathway_df$is_concordant),
        discordant_suggestive_n = sum(pathway_df$both_suggestive & pathway_df$is_discordant),
        lnb_up_ace_up_pathway_n = up_up[["n"]],
        lnb_down_ace_down_pathway_n = down_down[["n"]],
        lnb_up_ace_down_pathway_n = up_down[["n"]],
        lnb_down_ace_up_pathway_n = down_up[["n"]],
        lnb_up_ace_up_fisher_pvalue = up_up[["p"]],
        lnb_down_ace_down_fisher_pvalue = down_down[["p"]],
        lnb_up_ace_down_fisher_pvalue = up_down[["p"]],
        lnb_down_ace_up_fisher_pvalue = down_up[["p"]],
        best_concordant_fisher_pvalue = min(up_up[["p"]], down_down[["p"]], na.rm = TRUE),
        best_discordant_fisher_pvalue = min(up_down[["p"]], down_up[["p"]], na.rm = TRUE),
        stringsAsFactors = FALSE
      )
      pathway_rows[[pi]] <- pathway_df
      si <- si + 1L
      pi <- pi + 1L
    }
  }

  summary <- if (length(summary_rows)) do.call(rbind, summary_rows) else empty_summary()
  pathways <- if (length(pathway_rows)) do.call(rbind, pathway_rows) else empty_pathways()

  if (nrow(summary)) {
    pcols <- c("spearman_pvalue", "pearson_pvalue",
               "lnb_up_ace_up_fisher_pvalue", "lnb_down_ace_down_fisher_pvalue",
               "lnb_up_ace_down_fisher_pvalue", "lnb_down_ace_up_fisher_pvalue",
               "best_concordant_fisher_pvalue", "best_discordant_fisher_pvalue")
    for (pcol in pcols) {
      summary[[sub("_pvalue$", "_padj", pcol)]] <- p.adjust(summary[[pcol]], method = "BH")
    }
  }
  list(summary = summary, pathways = pathways)
}

parse_gene_list <- function(x) {
  if (is.na(x) || !nzchar(x)) return(character())
  genes <- trimws(unlist(strsplit(x, "[;,[:space:]]+")))
  unique(genes[nzchar(genes)])
}

build_leading_edge_overlap <- function(pathways) {
  if (!nrow(pathways)) {
    return(data.frame(
      comparison_id = character(), comparison_role = character(), database = character(),
      geneSet = character(), description = character(), mouse_NES = numeric(),
      human_NES = numeric(), mouse_FDR = numeric(), human_FDR = numeric(),
      mouse_leading_edge_n = integer(), human_leading_edge_n = integer(),
      shared_leading_edge_n = integer(), shared_leading_edge_genes = character(),
      stringsAsFactors = FALSE
    ))
  }
  keep <- pathways$both_suggestive & pathways$is_concordant
  pathways <- pathways[keep, , drop = FALSE]
  if (!nrow(pathways)) {
    return(data.frame(
      comparison_id = character(), comparison_role = character(), database = character(),
      geneSet = character(), description = character(), mouse_NES = numeric(),
      human_NES = numeric(), mouse_FDR = numeric(), human_FDR = numeric(),
      mouse_leading_edge_n = integer(), human_leading_edge_n = integer(),
      shared_leading_edge_n = integer(), shared_leading_edge_genes = character(),
      stringsAsFactors = FALSE
    ))
  }
  rows <- vector("list", nrow(pathways))
  for (i in seq_len(nrow(pathways))) {
    mg <- parse_gene_list(pathways$mouse_userId[[i]])
    hg <- parse_gene_list(pathways$human_userId[[i]])
    shared <- intersect(mg, hg)
    rows[[i]] <- data.frame(
      comparison_id = pathways$comparison_id[[i]],
      comparison_role = pathways$comparison_role[[i]],
      database = pathways$database[[i]],
      geneSet = pathways$geneSet[[i]],
      description = pathways$description[[i]],
      mouse_NES = pathways$mouse_NES[[i]],
      human_NES = pathways$human_NES[[i]],
      mouse_FDR = pathways$mouse_FDR[[i]],
      human_FDR = pathways$human_FDR[[i]],
      mouse_leading_edge_n = length(mg),
      human_leading_edge_n = length(hg),
      shared_leading_edge_n = length(shared),
      shared_leading_edge_genes = paste(shared, collapse = ";"),
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows)
}

blank_plot <- function(message) {
  ggplot() +
    annotate("text", x = 0, y = 0, label = message, size = 4) +
    theme_void()
}

save_plot_pair <- function(plot, path_base, width = 10, height = 7) {
  ggsave(paste0(path_base, ".png"), plot, width = width, height = height, dpi = 300)
  ggsave(paste0(path_base, ".pdf"), plot, width = width, height = height)
}

plot_outputs <- function(summary, pathways, figures_dir) {
  dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

  if (nrow(summary)) {
    heat_df <- summary[summary$human_integration == "derived_batch" &
                         summary$human_phenotype == "tot_adverse_exp" &
                         summary$comparison_role %in% c("primary", "pv_sensitivity", "female_sensitivity"), , drop = FALSE]
    if (!nrow(heat_df)) heat_df <- summary
    heat_df$comparison_label <- paste(heat_df$mouse_dataset, "vs", heat_df$human_cell_type,
                                      heat_df$human_sex, sep = " | ")
    heat <- ggplot(heat_df, aes(x = comparison_role, y = reorder(comparison_label, spearman_r), fill = spearman_r)) +
      geom_tile(color = "white", linewidth = 0.2) +
      scale_fill_gradient2(low = "#3B6FB6", mid = "white", high = "#B63B3B", midpoint = 0, na.value = "grey85") +
      labs(x = NULL, y = NULL, fill = "Spearman r", title = "Functional NES Concordance") +
      theme_minimal(base_size = 11) +
      theme(axis.text.x = element_text(angle = 30, hjust = 1))
  } else {
    heat <- blank_plot("No functional concordance summary rows")
  }
  save_plot_pair(heat, file.path(figures_dir, "functional_nes_correlation_heatmap"), 9, 7)

  primary <- pathways[pathways$comparison_role == "primary" &
                        pathways$human_integration == "derived_batch" &
                        pathways$human_phenotype == "tot_adverse_exp" &
                        pathways$human_sex == "Male", , drop = FALSE]
  if (nrow(primary)) {
    primary$facet_label <- paste(primary$mouse_dataset, "vs", primary$human_cell_type)
    scatter <- ggplot(primary, aes(mouse_NES, human_NES, color = pathway_status)) +
      geom_hline(yintercept = 0, color = "grey70", linewidth = 0.3) +
      geom_vline(xintercept = 0, color = "grey70", linewidth = 0.3) +
      geom_point(alpha = 0.55, size = 1.4) +
      facet_grid(facet_label ~ database, scales = "free") +
      scale_color_manual(values = c(
        concordant_sig = "#157F3B",
        discordant_sig = "#B53A30",
        concordant_suggestive = "#63A85B",
        discordant_suggestive = "#D77A61",
        other = "grey60"
      ), drop = FALSE) +
      labs(x = "Mouse LNB NES", y = "Human ACE NES", color = NULL,
           title = "Primary Male Mouse-Human Pathway Concordance") +
      theme_minimal(base_size = 10)
  } else {
    scatter <- blank_plot("No primary male pathway rows")
  }
  save_plot_pair(scatter, file.path(figures_dir, "functional_nes_scatter_primary_male"), 12, 8)

  first_line <- pathways[pathways$human_integration == "derived_batch" &
                           pathways$human_phenotype == "tot_adverse_exp" &
                           pathways$human_sex == "Male", , drop = FALSE]
  concordant <- first_line[first_line$is_concordant & first_line$both_suggestive, , drop = FALSE]
  if (!nrow(concordant)) concordant <- first_line[first_line$is_concordant, , drop = FALSE]
  if (nrow(concordant)) {
    concordant <- concordant[order(concordant$max_pair_FDR, -concordant$mean_abs_NES), , drop = FALSE]
    concordant <- concordant[seq_len(min(30, nrow(concordant))), , drop = FALSE]
    concordant$comparison_label <- paste(concordant$mouse_dataset, "vs", concordant$human_cell_type)
    concordant$pathway_label <- paste(concordant$database, concordant$description, sep = ": ")
    dot <- ggplot(concordant, aes(x = comparison_label, y = reorder(pathway_label, mean_abs_NES))) +
      geom_point(aes(size = pmin(-log10(pmax(max_pair_FDR, 1e-300)), 20),
                     color = (mouse_NES + human_NES) / 2)) +
      scale_color_gradient2(low = "#3B6FB6", mid = "white", high = "#B63B3B", midpoint = 0) +
      labs(x = NULL, y = NULL, size = "-log10 max FDR", color = "Mean NES",
           title = "Top Concordant Functional Pathways") +
      theme_minimal(base_size = 10) +
      theme(axis.text.x = element_text(angle = 30, hjust = 1))
  } else {
    dot <- blank_plot("No concordant pathway rows")
  }
  save_plot_pair(dot, file.path(figures_dir, "functional_concordant_pathways_dotplot"), 12, 8)

  discordant <- first_line[first_line$is_discordant & first_line$both_suggestive, , drop = FALSE]
  if (!nrow(discordant)) discordant <- first_line[first_line$is_discordant, , drop = FALSE]
  if (nrow(discordant)) {
    discordant <- discordant[order(discordant$max_pair_FDR, -discordant$mean_abs_NES), , drop = FALSE]
    discordant <- discordant[seq_len(min(25, nrow(discordant))), , drop = FALSE]
    discordant$comparison_label <- paste(discordant$mouse_dataset, "vs", discordant$human_cell_type)
    discordant$pathway_label <- paste(discordant$database, discordant$description, sep = ": ")
    dis <- ggplot(discordant, aes(x = comparison_label, y = reorder(pathway_label, mean_abs_NES))) +
      geom_point(aes(size = pmin(-log10(pmax(max_pair_FDR, 1e-300)), 20),
                     color = mouse_NES - human_NES)) +
      scale_color_gradient2(low = "#3B6FB6", mid = "white", high = "#B63B3B", midpoint = 0) +
      labs(x = NULL, y = NULL, size = "-log10 max FDR", color = "Mouse - human NES",
           title = "Top Discordant Functional Pathways") +
      theme_minimal(base_size = 10) +
      theme(axis.text.x = element_text(angle = 30, hjust = 1))
  } else {
    dis <- blank_plot("No discordant pathway rows")
  }
  save_plot_pair(dis, file.path(figures_dir, "functional_discordant_pathways_dotplot"), 12, 8)
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
dir.create(args$output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(args$figures_dir, recursive = TRUE, showWarnings = FALSE)
set.seed(args$seed)

cat("Preparing ranked genes...\n")
ranked <- prepare_ranked_genes(args)
write.csv(ranked, file.path(args$output_dir, "functional_ranked_genes.csv"), row.names = FALSE)
cat("Ranked rows:", nrow(ranked), "\n")

cat("Running ranked GSEA...\n")
gsea <- run_all_gsea(args, ranked)
write.csv(gsea, file.path(args$output_dir, "functional_gsea_results.csv"), row.names = FALSE)
cat("GSEA rows:", nrow(gsea), "\n")

cat("Computing pathway concordance...\n")
concordance <- build_concordance(args, gsea)
write.csv(concordance$summary, file.path(args$output_dir, "functional_concordance_summary.csv"), row.names = FALSE)
write.csv(concordance$pathways, file.path(args$output_dir, "functional_concordance_pathways.csv"), row.names = FALSE)
leading <- build_leading_edge_overlap(concordance$pathways)
write.csv(leading, file.path(args$output_dir, "functional_leading_edge_overlap.csv"), row.names = FALSE)

cat("Writing figures...\n")
plot_outputs(concordance$summary, concordance$pathways, args$figures_dir)

cat("Functional concordance complete.\n")
