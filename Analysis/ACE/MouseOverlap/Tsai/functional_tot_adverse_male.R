#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)

suppressPackageStartupMessages({
  library(WebGestaltR)
  library(ggplot2)
})

GSEA_DATABASES <- c(
  "geneontology_Biological_Process_noRedundant",
  "geneontology_Cellular_Component_noRedundant",
  "geneontology_Molecular_Function_noRedundant",
  "pathway_KEGG",
  "pathway_Panther",
  "pathway_Reactome",
  "pathway_Wikipathway",
  "network_Transcription_Factor_target"
)

GSEA_LIBRARY_NAMES <- c(
  "Biological_Process_noRedundant",
  "Cellular_Component_noRedundant",
  "Molecular_Function_noRedundant",
  "KEGG",
  "Panther",
  "Reactome",
  "Wikipathway",
  "Transcription_Factor_target"
)

EX_SUBTYPES <- c("Ex-L2_3", "Ex-L4", "Ex-L4_5", "Ex-L5", "Ex-L5_6", "Ex-L5_6-CC", "Ex-NRGN")
PV_SENSITIVITY <- c("In-PV_Basket", "In-PV_Chandelier")

parse_args <- function(argv) {
  out <- list(
    mouse_deg = NULL,
    ortholog_table = NULL,
    human_gsea_dir = NULL,
    output_dir = NULL,
    figures_dir = NULL,
    custom_gmt = "",
    permutations = 1000L,
    min_size = 10L,
    max_size = 500L,
    fdr_threshold = 0.2,
    alpha = 0.05,
    suggestive_alpha = 0.2,
    threads = 1L,
    seed = 1L,
    skip_existing = FALSE
  )
  key_map <- c(
    "--mouse-deg" = "mouse_deg",
    "--ortholog-table" = "ortholog_table",
    "--human-gsea-dir" = "human_gsea_dir",
    "--output-dir" = "output_dir",
    "--figures-dir" = "figures_dir",
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
  for (field in c("permutations", "min_size", "max_size", "threads", "seed")) {
    out[[field]] <- as.integer(out[[field]])
  }
  for (field in c("fdr_threshold", "alpha", "suggestive_alpha")) {
    out[[field]] <- as.numeric(out[[field]])
  }
  required <- c("mouse_deg", "ortholog_table", "human_gsea_dir", "output_dir", "figures_dir")
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

num <- function(x) suppressWarnings(as.numeric(x))

read_csv_safe <- function(path) {
  read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}

write_csv <- function(df, path) {
  write.csv(df, path, row.names = FALSE, na = "")
}

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

parse_library_from_file <- function(path, libraries = c(GSEA_LIBRARY_NAMES, "custom")) {
  stem <- sub("\\.rds$", "", basename(path))
  for (library in libraries) {
    suffix <- paste0("_", library)
    if (endsWith(stem, suffix)) {
      return(list(cell_type = sub(paste0("_", library, "$"), "", stem), library = library))
    }
  }
  NULL
}

empty_gsea <- function() {
  data.frame(
    species = character(), mouse_dataset = character(), mouse_population = character(),
    mouse_condition = character(), human_cell_type = character(), human_reference_type = character(),
    library = character(), geneSet = character(), term = character(), NES = numeric(),
    pvalue = numeric(), FDR = numeric(), size = numeric(), leadingEdgeNum = numeric(),
    userId = character(), source_file = character(), stringsAsFactors = FALSE
  )
}

normalize_gsea_result <- function(x) {
  if (is.null(x)) return(data.frame())
  df <- tryCatch(as.data.frame(x), error = function(e) data.frame())
  if (!nrow(df)) return(data.frame())
  get_col <- function(names, default = NA) {
    hit <- intersect(names, colnames(df))[1]
    if (is.na(hit)) rep(default, nrow(df)) else df[[hit]]
  }
  data.frame(
    geneSet = as.character(get_col(c("geneSet", "ID", "term"))),
    term = as.character(get_col(c("description", "Description", "term", "geneSet"))),
    enrichmentScore = num(get_col(c("enrichmentScore"))),
    NES = num(get_col(c("normalizedEnrichmentScore", "NES", "enrichmentScore"))),
    pvalue = num(get_col(c("pValue", "pvalue", "p.value"))),
    FDR = num(get_col(c("FDR", "padj", "p.adjust", "qvalue"))),
    size = num(get_col(c("size", "setSize"))),
    leadingEdgeNum = num(get_col(c("leadingEdgeNum"))),
    leadingEdgeId = as.character(get_col(c("leadingEdgeId"), "")),
    userId = as.character(get_col(c("userId", "geneID"), "")),
    stringsAsFactors = FALSE
  )
}

flatten_human_gsea <- function(human_gsea_dir) {
  files <- list.files(human_gsea_dir, pattern = "\\.rds$", full.names = TRUE)
  rows <- list()
  for (path in files) {
    parsed <- parse_library_from_file(path)
    if (is.null(parsed)) next
    x <- tryCatch(readRDS(path), error = function(e) NULL)
    df <- normalize_gsea_result(x)
    if (!nrow(df)) next
    rows[[length(rows) + 1]] <- data.frame(
      species = "human",
      mouse_dataset = NA_character_,
      mouse_population = NA_character_,
      mouse_condition = NA_character_,
      human_cell_type = parsed$cell_type,
      human_reference_type = "single_cell_type",
      library = parsed$library,
      df[, c("geneSet", "term", "NES", "pvalue", "FDR", "size", "leadingEdgeNum", "userId")],
      source_file = path,
      stringsAsFactors = FALSE
    )
  }
  if (!length(rows)) return(empty_gsea())
  out <- do.call(rbind, rows)
  out <- out[order(out$human_cell_type, out$library, out$FDR, out$pvalue), , drop = FALSE]
  rownames(out) <- NULL
  out
}

parse_gene_list <- function(x) {
  if (is.na(x) || !nzchar(x)) return(character())
  genes <- trimws(unlist(strsplit(x, "[;,/[:space:]]+")))
  unique(genes[nzchar(genes)])
}

collapse_gene_lists <- function(x) {
  genes <- unique(unlist(lapply(x, parse_gene_list), use.names = FALSE))
  paste(sort(genes), collapse = ";")
}

aggregate_ex_reference <- function(human_gsea) {
  ex <- human_gsea[human_gsea$human_cell_type %in% EX_SUBTYPES, , drop = FALSE]
  if (!nrow(ex)) return(empty_gsea())
  keys <- unique(ex[, c("library", "geneSet"), drop = FALSE])
  rows <- vector("list", nrow(keys))
  for (i in seq_len(nrow(keys))) {
    sub <- ex[ex$library == keys$library[[i]] & ex$geneSet == keys$geneSet[[i]], , drop = FALSE]
    signs <- sign(num(sub$NES))
    signs <- signs[signs != 0 & is.finite(signs)]
    positive_n <- sum(num(sub$NES) > 0, na.rm = TRUE)
    negative_n <- sum(num(sub$NES) < 0, na.rm = TRUE)
    support_n <- length(unique(sub$human_cell_type))
    sign_consistency <- if (length(signs)) max(positive_n, negative_n) / length(signs) else NA_real_
    rows[[i]] <- data.frame(
      species = "human",
      mouse_dataset = NA_character_,
      mouse_population = NA_character_,
      mouse_condition = NA_character_,
      human_cell_type = "Ex_aggregate",
      human_reference_type = "ex_subtype_aggregate",
      library = keys$library[[i]],
      geneSet = keys$geneSet[[i]],
      term = sub$term[[which.min(num(sub$FDR))]],
      NES = stats::median(num(sub$NES), na.rm = TRUE),
      pvalue = min(num(sub$pvalue), na.rm = TRUE),
      FDR = min(num(sub$FDR), na.rm = TRUE),
      size = stats::median(num(sub$size), na.rm = TRUE),
      leadingEdgeNum = stats::median(num(sub$leadingEdgeNum), na.rm = TRUE),
      userId = collapse_gene_lists(sub$userId),
      source_file = paste(unique(sub$source_file), collapse = ";"),
      subtype_support_n = support_n,
      subtype_support_fraction = support_n / length(EX_SUBTYPES),
      subtype_positive_n = positive_n,
      subtype_negative_n = negative_n,
      subtype_sign_consistency = sign_consistency,
      subtype_list = paste(sort(unique(sub$human_cell_type)), collapse = ";"),
      stringsAsFactors = FALSE
    )
  }
  out <- do.call(rbind, rows)
  out <- out[order(out$library, out$FDR, out$pvalue), , drop = FALSE]
  rownames(out) <- NULL
  out
}

prepare_mouse_ranked <- function(mouse_deg_path, ortholog_path) {
  mouse <- read_csv_safe(mouse_deg_path)
  orth <- read_csv_safe(ortholog_path)
  required_mouse <- c("mouse_dataset", "mouse_population", "condition", "mouse_symbol",
                      "mouse_stat", "mouse_lnb_log2fc", "mouse_pvalue", "mouse_padj")
  required_orth <- c("mouse_symbol", "human_symbol", "is_one_to_one")
  for (col in required_mouse) if (!col %in% names(mouse)) stop("Mouse table missing column: ", col)
  for (col in required_orth) if (!col %in% names(orth)) stop("Ortholog table missing column: ", col)

  orth <- orth[tolower(as.character(orth$is_one_to_one)) %in% c("true", "1"), , drop = FALSE]
  orth <- orth[!duplicated(orth[, c("mouse_symbol", "human_symbol")]), , drop = FALSE]
  mouse <- merge(mouse, orth[, c("mouse_symbol", "human_symbol")], by = "mouse_symbol", all = FALSE)

  mouse$mouse_lnb_stat <- safe_signed_rank(-num(mouse$mouse_stat), mouse$mouse_lnb_log2fc, mouse$mouse_pvalue)
  mouse <- mouse[is.finite(mouse$mouse_lnb_stat) & nzchar(mouse$human_symbol), , drop = FALSE]
  mouse <- deduplicate_best(mouse, c("mouse_dataset", "human_symbol"), "mouse_padj", "mouse_lnb_stat")
  out <- data.frame(
    mouse_dataset = mouse$mouse_dataset,
    mouse_population = mouse$mouse_population,
    mouse_condition = mouse$condition,
    human_symbol = mouse$human_symbol,
    mouse_symbol = mouse$mouse_symbol,
    rank = mouse$mouse_lnb_stat,
    mouse_stat_original = num(mouse$mouse_stat),
    mouse_lnb_stat = mouse$mouse_lnb_stat,
    mouse_lnb_log2fc = num(mouse$mouse_lnb_log2fc),
    mouse_pvalue = num(mouse$mouse_pvalue),
    mouse_padj = num(mouse$mouse_padj),
    stringsAsFactors = FALSE
  )
  out[order(out$mouse_dataset, -out$rank), , drop = FALSE]
}

database_table <- function(custom_gmt = "") {
  if (nzchar(custom_gmt)) {
    return(data.frame(library = "custom", database_source = custom_gmt, stringsAsFactors = FALSE))
  }
  data.frame(library = db_short_name(GSEA_DATABASES), database_source = GSEA_DATABASES, stringsAsFactors = FALSE)
}

run_mouse_gsea_one <- function(args, mouse_group, db_row, cache_dir) {
  dataset <- unique(mouse_group$mouse_dataset)
  if (length(dataset) != 1) stop("Mouse group has multiple datasets")
  cache_path <- file.path(cache_dir, paste0(sanitize(dataset), "__", sanitize(db_row$library), ".rds"))

  if (args$skip_existing && file.exists(cache_path) && file.info(cache_path)$size > 0) {
    result <- readRDS(cache_path)
  } else {
    interest <- data.frame(gene = mouse_group$human_symbol, rank = num(mouse_group$rank), stringsAsFactors = FALSE)
    interest <- interest[is.finite(interest$rank) & nzchar(interest$gene), , drop = FALSE]
    interest <- interest[order(-interest$rank), , drop = FALSE]

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
      warning("WebGestaltR failed for ", dataset, " / ", db_row$library, ": ", conditionMessage(e))
      data.frame()
    })
    saveRDS(result, cache_path)
  }

  result <- normalize_gsea_result(result)
  if (!nrow(result)) return(empty_gsea())
  meta <- mouse_group[1, c("mouse_dataset", "mouse_population", "mouse_condition"), drop = FALSE]
  data.frame(
    species = "mouse",
    mouse_dataset = meta$mouse_dataset,
    mouse_population = meta$mouse_population,
    mouse_condition = meta$mouse_condition,
    human_cell_type = NA_character_,
    human_reference_type = NA_character_,
    library = db_row$library,
    result[, c("geneSet", "term", "NES", "pvalue", "FDR", "size", "leadingEdgeNum", "userId")],
    source_file = cache_path,
    stringsAsFactors = FALSE
  )
}

run_mouse_gsea <- function(args, ranked) {
  cache_dir <- file.path(args$output_dir, "mouse_gsea_cache")
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(args$output_dir, "webgestalt_cache"), recursive = TRUE, showWarnings = FALSE)
  dbs <- database_table(args$custom_gmt)
  groups <- split(ranked, ranked$mouse_dataset)
  rows <- list()
  idx <- 1L
  for (dataset in names(groups)) {
    cat("Mouse GSEA:", dataset, "(", nrow(groups[[dataset]]), "ranked genes )\n")
    for (i in seq_len(nrow(dbs))) {
      cat("  library:", dbs$library[[i]], "\n")
      rows[[idx]] <- run_mouse_gsea_one(args, groups[[dataset]], dbs[i, ], cache_dir)
      idx <- idx + 1L
    }
  }
  rows <- rows[vapply(rows, nrow, integer(1)) > 0]
  if (!length(rows)) return(empty_gsea())
  out <- do.call(rbind, rows)
  out[order(out$mouse_dataset, out$library, out$FDR, out$pvalue), , drop = FALSE]
}

build_comparison_table <- function(mouse_gsea, human_gsea, ex_aggregate) {
  mouse_meta <- unique(mouse_gsea[, c("mouse_dataset", "mouse_population", "mouse_condition"), drop = FALSE])
  rows <- list()
  i <- 1L
  for (m in seq_len(nrow(mouse_meta))) {
    mm <- mouse_meta[m, , drop = FALSE]
    if (mm$mouse_population == "NeuN") {
      rows[[i]] <- data.frame(mm, human_cell_type = "Ex_aggregate", human_reference_type = "ex_subtype_aggregate", comparison_role = "primary")
      i <- i + 1L
      for (ct in intersect(EX_SUBTYPES, unique(human_gsea$human_cell_type))) {
        rows[[i]] <- data.frame(mm, human_cell_type = ct, human_reference_type = "single_cell_type", comparison_role = "ex_subtype_sensitivity")
        i <- i + 1L
      }
    } else if (mm$mouse_population == "PV") {
      if ("Inh" %in% unique(human_gsea$human_cell_type)) {
        rows[[i]] <- data.frame(mm, human_cell_type = "Inh", human_reference_type = "single_cell_type", comparison_role = "primary")
        i <- i + 1L
      }
      for (ct in intersect(PV_SENSITIVITY, unique(human_gsea$human_cell_type))) {
        rows[[i]] <- data.frame(mm, human_cell_type = ct, human_reference_type = "single_cell_type", comparison_role = "pv_sensitivity")
        i <- i + 1L
      }
    }
  }
  if (!length(rows)) return(data.frame())
  out <- do.call(rbind, rows)
  out$comparison_id <- paste(out$mouse_dataset, out$human_cell_type, out$comparison_role, sep = "|")
  out
}

direction <- function(nes, up, down) {
  ifelse(nes > 0, up, ifelse(nes < 0, down, "zero"))
}

cor_or_na <- function(x, y, method) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 3) return(c(r = NA_real_, p = NA_real_))
  test <- suppressWarnings(stats::cor.test(x[ok], y[ok], method = method, exact = FALSE))
  c(r = unname(test$estimate), p = test$p.value)
}

fisher_direction <- function(df, mouse_dir, human_dir, mouse_sig_col, human_sig_col) {
  universe <- unique(df$pathway_key)
  mouse_set <- unique(df$pathway_key[df[[mouse_sig_col]] & df$mouse_direction == mouse_dir])
  human_set <- unique(df$pathway_key[df[[human_sig_col]] & df$human_direction == human_dir])
  a <- length(intersect(mouse_set, human_set))
  b <- length(setdiff(mouse_set, human_set))
  c <- length(setdiff(human_set, mouse_set))
  d <- max(length(universe) - a - b - c, 0)
  p <- if (length(universe)) stats::fisher.test(matrix(c(a, b, c, d), nrow = 2), alternative = "greater")$p.value else NA_real_
  c(n = a, p = p)
}

compare_gsea <- function(args, mouse_gsea, human_gsea, ex_aggregate) {
  human_refs <- rbind(
    human_gsea,
    ex_aggregate[, intersect(names(human_gsea), names(ex_aggregate)), drop = FALSE]
  )
  comparisons <- build_comparison_table(mouse_gsea, human_gsea, ex_aggregate)
  pathway_rows <- list()
  summary_rows <- list()
  pi <- 1L
  si <- 1L

  for (i in seq_len(nrow(comparisons))) {
    comp <- comparisons[i, , drop = FALSE]
    m <- mouse_gsea[mouse_gsea$mouse_dataset == comp$mouse_dataset, , drop = FALSE]
    h <- human_refs[human_refs$human_cell_type == comp$human_cell_type &
                      human_refs$human_reference_type == comp$human_reference_type, , drop = FALSE]
    if (!nrow(m) || !nrow(h)) next
    joined <- merge(m, h, by = c("library", "geneSet"), suffixes = c("_mouse", "_human"))
    joined <- joined[is.finite(joined$NES_mouse) & is.finite(joined$NES_human), , drop = FALSE]
    if (!nrow(joined)) next
    pathway_key <- paste(joined$library, joined$geneSet, sep = "|")
    mouse_sig <- is.finite(joined$FDR_mouse) & joined$FDR_mouse < args$alpha
    human_sig <- is.finite(joined$FDR_human) & joined$FDR_human < args$alpha
    mouse_suggestive <- is.finite(joined$FDR_mouse) & joined$FDR_mouse < args$suggestive_alpha
    human_suggestive <- is.finite(joined$FDR_human) & joined$FDR_human < args$suggestive_alpha
    mouse_dir <- direction(joined$NES_mouse, "LNB_up", "LNB_down")
    human_dir <- direction(joined$NES_human, "ACE_up", "ACE_down")
    concordant <- (mouse_dir == "LNB_up" & human_dir == "ACE_up") |
      (mouse_dir == "LNB_down" & human_dir == "ACE_down")
    discordant <- (mouse_dir %in% c("LNB_up", "LNB_down")) &
      (human_dir %in% c("ACE_up", "ACE_down")) & !concordant
    both_sig <- mouse_sig & human_sig
    both_suggestive <- mouse_suggestive & human_suggestive
    mean_abs_nes <- (abs(joined$NES_mouse) + abs(joined$NES_human)) / 2
    max_fdr <- pmax(joined$FDR_mouse, joined$FDR_human)
    pathway_status <- ifelse(both_sig & concordant, "concordant_sig",
      ifelse(both_sig & discordant, "discordant_sig",
        ifelse(both_suggestive & concordant, "concordant_suggestive",
          ifelse(both_suggestive & discordant, "discordant_suggestive", "other"))))

    subtype_cols <- c("subtype_support_n", "subtype_support_fraction", "subtype_positive_n",
                      "subtype_negative_n", "subtype_sign_consistency", "subtype_list")
    subtype_data <- data.frame(
      subtype_support_n = NA_real_, subtype_support_fraction = NA_real_,
      subtype_positive_n = NA_real_, subtype_negative_n = NA_real_,
      subtype_sign_consistency = NA_real_, subtype_list = NA_character_
    )
    if (comp$human_reference_type == "ex_subtype_aggregate") {
      agg_meta <- ex_aggregate[, c("library", "geneSet", subtype_cols), drop = FALSE]
      subtype_data <- merge(joined[, c("library", "geneSet"), drop = FALSE], agg_meta,
                            by = c("library", "geneSet"), all.x = TRUE, sort = FALSE)[, subtype_cols, drop = FALSE]
    } else {
      subtype_data <- subtype_data[rep(1, nrow(joined)), , drop = FALSE]
    }

    paths <- data.frame(
      comparison_id = comp$comparison_id,
      comparison_role = comp$comparison_role,
      mouse_dataset = comp$mouse_dataset,
      mouse_population = comp$mouse_population,
      mouse_condition = comp$mouse_condition,
      human_cell_type = comp$human_cell_type,
      human_reference_type = comp$human_reference_type,
      library = joined$library,
      geneSet = joined$geneSet,
      pathway_key = pathway_key,
      term = ifelse(nzchar(joined$term_mouse), joined$term_mouse, joined$term_human),
      mouse_NES = joined$NES_mouse,
      human_NES = joined$NES_human,
      mouse_pvalue = joined$pvalue_mouse,
      human_pvalue = joined$pvalue_human,
      mouse_FDR = joined$FDR_mouse,
      human_FDR = joined$FDR_human,
      mouse_size = joined$size_mouse,
      human_size = joined$size_human,
      mouse_leadingEdgeNum = joined$leadingEdgeNum_mouse,
      human_leadingEdgeNum = joined$leadingEdgeNum_human,
      mouse_userId = joined$userId_mouse,
      human_userId = joined$userId_human,
      mouse_direction = mouse_dir,
      human_direction = human_dir,
      mouse_sig = mouse_sig,
      human_sig = human_sig,
      both_sig = both_sig,
      mouse_suggestive = mouse_suggestive,
      human_suggestive = human_suggestive,
      both_suggestive = both_suggestive,
      is_concordant = concordant,
      is_discordant = discordant,
      pathway_status = pathway_status,
      mean_abs_NES = mean_abs_nes,
      max_pair_FDR = max_fdr,
      stringsAsFactors = FALSE
    )
    paths <- cbind(paths, subtype_data)
    pathway_rows[[pi]] <- paths
    pi <- pi + 1L

    spear <- cor_or_na(paths$mouse_NES, paths$human_NES, "spearman")
    pear <- cor_or_na(paths$mouse_NES, paths$human_NES, "pearson")
    up_up <- fisher_direction(paths, "LNB_up", "ACE_up", "mouse_sig", "human_sig")
    down_down <- fisher_direction(paths, "LNB_down", "ACE_down", "mouse_sig", "human_sig")
    up_down <- fisher_direction(paths, "LNB_up", "ACE_down", "mouse_sig", "human_sig")
    down_up <- fisher_direction(paths, "LNB_down", "ACE_up", "mouse_sig", "human_sig")

    summary_rows[[si]] <- data.frame(
      comparison_id = comp$comparison_id,
      comparison_role = comp$comparison_role,
      mouse_dataset = comp$mouse_dataset,
      mouse_population = comp$mouse_population,
      mouse_condition = comp$mouse_condition,
      human_cell_type = comp$human_cell_type,
      human_reference_type = comp$human_reference_type,
      shared_reported_pathway_n = nrow(paths),
      shared_library_n = length(unique(paths$library)),
      spearman_r = spear[["r"]],
      spearman_pvalue = spear[["p"]],
      pearson_r = pear[["r"]],
      pearson_pvalue = pear[["p"]],
      sign_concordance = mean(paths$is_concordant, na.rm = TRUE),
      concordant_sig_n = sum(paths$both_sig & paths$is_concordant),
      discordant_sig_n = sum(paths$both_sig & paths$is_discordant),
      concordant_suggestive_n = sum(paths$both_suggestive & paths$is_concordant),
      discordant_suggestive_n = sum(paths$both_suggestive & paths$is_discordant),
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
    si <- si + 1L
  }

  pathways <- if (length(pathway_rows)) do.call(rbind, pathway_rows) else data.frame()
  summary <- if (length(summary_rows)) do.call(rbind, summary_rows) else data.frame()
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

build_leading_edge_overlap <- function(pathways) {
  if (!nrow(pathways)) return(data.frame())
  keep <- pathways$both_suggestive & pathways$is_concordant
  pathways <- pathways[keep, , drop = FALSE]
  if (!nrow(pathways)) return(data.frame())
  rows <- vector("list", nrow(pathways))
  for (i in seq_len(nrow(pathways))) {
    mg <- parse_gene_list(pathways$mouse_userId[[i]])
    hg <- parse_gene_list(pathways$human_userId[[i]])
    shared <- intersect(mg, hg)
    rows[[i]] <- data.frame(
      comparison_id = pathways$comparison_id[[i]],
      comparison_role = pathways$comparison_role[[i]],
      mouse_dataset = pathways$mouse_dataset[[i]],
      human_cell_type = pathways$human_cell_type[[i]],
      library = pathways$library[[i]],
      geneSet = pathways$geneSet[[i]],
      term = pathways$term[[i]],
      mouse_NES = pathways$mouse_NES[[i]],
      human_NES = pathways$human_NES[[i]],
      mouse_FDR = pathways$mouse_FDR[[i]],
      human_FDR = pathways$human_FDR[[i]],
      mouse_leading_edge_n = length(mg),
      human_leading_edge_n = length(hg),
      shared_leading_edge_n = length(shared),
      shared_leading_edge_genes = paste(sort(shared), collapse = ";"),
      stringsAsFactors = FALSE
    )
  }
  out <- do.call(rbind, rows)
  out[order(-out$shared_leading_edge_n, out$mouse_FDR + out$human_FDR), , drop = FALSE]
}

blank_plot <- function(message) {
  ggplot() + annotate("text", x = 0, y = 0, label = message, size = 4) + theme_void()
}

save_plot_pair <- function(plot, path_base, width = 10, height = 7) {
  ggsave(paste0(path_base, ".png"), plot, width = width, height = height, dpi = 300)
  ggsave(paste0(path_base, ".pdf"), plot, width = width, height = height)
}

plot_outputs <- function(summary, pathways, ex_aggregate, figures_dir) {
  dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
  if (nrow(summary)) {
    heat_df <- summary[summary$comparison_role %in% c("primary", "pv_sensitivity"), , drop = FALSE]
    if (!nrow(heat_df)) heat_df <- summary
    heat_df$comparison_label <- paste(heat_df$mouse_dataset, "vs", heat_df$human_cell_type)
    heat <- ggplot(heat_df, aes(x = comparison_role, y = reorder(comparison_label, spearman_r), fill = spearman_r)) +
      geom_tile(color = "white", linewidth = 0.2) +
      scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0, na.value = "grey85") +
      labs(x = NULL, y = NULL, fill = "Spearman r", title = "Mouse LNB vs Male Human ACE Functional Concordance") +
      theme_minimal(base_size = 11) +
      theme(axis.text.x = element_text(angle = 30, hjust = 1))
  } else {
    heat <- blank_plot("No concordance summary rows")
  }
  save_plot_pair(heat, file.path(figures_dir, "nes_concordance_heatmap"), 9, 7)

  primary <- pathways[pathways$comparison_role == "primary", , drop = FALSE]
  if (nrow(primary)) {
    primary$facet_label <- paste(primary$mouse_dataset, "vs", primary$human_cell_type)
    scatter <- ggplot(primary, aes(mouse_NES, human_NES, color = pathway_status)) +
      geom_hline(yintercept = 0, color = "grey70", linewidth = 0.3) +
      geom_vline(xintercept = 0, color = "grey70", linewidth = 0.3) +
      geom_point(alpha = 0.55, size = 1.35) +
      facet_wrap(~facet_label, scales = "free") +
      scale_color_manual(values = c(
        concordant_sig = "#1B9E77",
        discordant_sig = "#D95F02",
        concordant_suggestive = "#66A61E",
        discordant_suggestive = "#E6AB02",
        other = "grey60"
      ), drop = FALSE) +
      labs(x = "Mouse LNB NES", y = "Human ACE NES", color = NULL,
           title = "Primary Pathway NES Concordance") +
      theme_minimal(base_size = 10)
  } else {
    scatter <- blank_plot("No primary pathway rows")
  }
  save_plot_pair(scatter, file.path(figures_dir, "primary_nes_scatterplots"), 11, 8)

  concordant <- pathways[pathways$is_concordant & pathways$both_suggestive, , drop = FALSE]
  if (!nrow(concordant)) concordant <- pathways[pathways$is_concordant, , drop = FALSE]
  if (nrow(concordant)) {
    concordant <- concordant[order(concordant$max_pair_FDR, -concordant$mean_abs_NES), , drop = FALSE]
    concordant <- concordant[seq_len(min(35, nrow(concordant))), , drop = FALSE]
    concordant$comparison_label <- paste(concordant$mouse_dataset, "vs", concordant$human_cell_type)
    concordant$term_short <- ifelse(nchar(concordant$term) > 75, paste0(substr(concordant$term, 1, 72), "..."), concordant$term)
    dot <- ggplot(concordant, aes(x = comparison_label, y = reorder(term_short, mean_abs_NES))) +
      geom_point(aes(size = pmin(-log10(pmax(max_pair_FDR, 1e-300)), 20),
                     color = (mouse_NES + human_NES) / 2)) +
      scale_color_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0) +
      labs(x = NULL, y = NULL, size = "-log10 max FDR", color = "Mean NES",
           title = "Top Concordant Functional Terms") +
      theme_minimal(base_size = 10) +
      theme(axis.text.x = element_text(angle = 30, hjust = 1))
  } else {
    dot <- blank_plot("No concordant terms")
  }
  save_plot_pair(dot, file.path(figures_dir, "top_concordant_terms_dotplot"), 12, 8)

  discordant <- pathways[pathways$is_discordant & pathways$both_suggestive, , drop = FALSE]
  if (!nrow(discordant)) discordant <- pathways[pathways$is_discordant, , drop = FALSE]
  if (nrow(discordant)) {
    discordant <- discordant[order(discordant$max_pair_FDR, -discordant$mean_abs_NES), , drop = FALSE]
    discordant <- discordant[seq_len(min(35, nrow(discordant))), , drop = FALSE]
    discordant$comparison_label <- paste(discordant$mouse_dataset, "vs", discordant$human_cell_type)
    discordant$term_short <- ifelse(nchar(discordant$term) > 75, paste0(substr(discordant$term, 1, 72), "..."), discordant$term)
    dis <- ggplot(discordant, aes(x = comparison_label, y = reorder(term_short, mean_abs_NES))) +
      geom_point(aes(size = pmin(-log10(pmax(max_pair_FDR, 1e-300)), 20),
                     color = mouse_NES - human_NES)) +
      scale_color_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0) +
      labs(x = NULL, y = NULL, size = "-log10 max FDR", color = "Mouse - human NES",
           title = "Top Discordant Functional Terms") +
      theme_minimal(base_size = 10) +
      theme(axis.text.x = element_text(angle = 30, hjust = 1))
  } else {
    dis <- blank_plot("No discordant terms")
  }
  save_plot_pair(dis, file.path(figures_dir, "top_discordant_terms_dotplot"), 12, 8)

  ex_paths <- pathways[pathways$human_reference_type == "ex_subtype_aggregate" &
                         is.finite(pathways$subtype_support_n), , drop = FALSE]
  if (nrow(ex_paths)) {
    ex_paths <- ex_paths[order(-ex_paths$subtype_support_n, -ex_paths$subtype_sign_consistency, ex_paths$max_pair_FDR), , drop = FALSE]
    ex_paths <- ex_paths[seq_len(min(35, nrow(ex_paths))), , drop = FALSE]
    ex_paths$comparison_label <- paste(ex_paths$mouse_dataset, "vs Ex aggregate")
    ex_paths$term_short <- ifelse(nchar(ex_paths$term) > 75, paste0(substr(ex_paths$term, 1, 72), "..."), ex_paths$term)
    ex_plot <- ggplot(ex_paths, aes(x = comparison_label, y = reorder(term_short, subtype_support_n))) +
      geom_point(aes(size = subtype_support_n, color = subtype_sign_consistency)) +
      scale_color_gradient(low = "#7570B3", high = "#1B9E77", na.value = "grey70") +
      labs(x = NULL, y = NULL, size = "Ex subtype support", color = "Sign consistency",
           title = "Human Excitatory Subtype Support for NeuN Comparisons") +
      theme_minimal(base_size = 10) +
      theme(axis.text.x = element_text(angle = 30, hjust = 1))
  } else {
    ex_plot <- blank_plot("No Ex aggregate pathway rows")
  }
  save_plot_pair(ex_plot, file.path(figures_dir, "ex_subtype_support_plot"), 12, 8)
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
dir.create(args$output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(args$figures_dir, recursive = TRUE, showWarnings = FALSE)
set.seed(args$seed)

cat("Flattening human male GSEA reference...\n")
human_gsea <- flatten_human_gsea(args$human_gsea_dir)
write_csv(human_gsea, file.path(args$output_dir, "human_gsea_reference.csv"))
cat("Human reference rows:", nrow(human_gsea), "\n")

cat("Aggregating human excitatory subtype reference...\n")
ex_aggregate <- aggregate_ex_reference(human_gsea)
write_csv(ex_aggregate, file.path(args$output_dir, "human_ex_aggregate_reference.csv"))
cat("Human Ex aggregate rows:", nrow(ex_aggregate), "\n")

cat("Preparing mouse ranked gene lists...\n")
mouse_ranked <- prepare_mouse_ranked(args$mouse_deg, args$ortholog_table)
write_csv(mouse_ranked, file.path(args$output_dir, "mouse_ranked_genes.csv"))
cat("Mouse ranked rows:", nrow(mouse_ranked), "\n")

cat("Running mouse WebGestaltR GSEA...\n")
mouse_gsea <- run_mouse_gsea(args, mouse_ranked)
write_csv(mouse_gsea, file.path(args$output_dir, "mouse_gsea_results.csv"))
cat("Mouse GSEA rows:", nrow(mouse_gsea), "\n")

cat("Computing functional concordance...\n")
concordance <- compare_gsea(args, mouse_gsea, human_gsea, ex_aggregate)
write_csv(concordance$summary, file.path(args$output_dir, "functional_concordance_summary.csv"))
write_csv(concordance$pathways, file.path(args$output_dir, "functional_concordance_pathways.csv"))

concordant_terms <- concordance$pathways[concordance$pathways$is_concordant, , drop = FALSE]
if (nrow(concordant_terms)) {
  concordant_terms <- concordant_terms[order(!concordant_terms$both_sig, !concordant_terms$both_suggestive,
                                             concordant_terms$max_pair_FDR, -concordant_terms$mean_abs_NES), , drop = FALSE]
}
discordant_terms <- concordance$pathways[concordance$pathways$is_discordant, , drop = FALSE]
if (nrow(discordant_terms)) {
  discordant_terms <- discordant_terms[order(!discordant_terms$both_sig, !discordant_terms$both_suggestive,
                                             discordant_terms$max_pair_FDR, -discordant_terms$mean_abs_NES), , drop = FALSE]
}
write_csv(concordant_terms, file.path(args$output_dir, "functional_concordant_terms.csv"))
write_csv(discordant_terms, file.path(args$output_dir, "functional_discordant_terms.csv"))

leading <- build_leading_edge_overlap(concordance$pathways)
write_csv(leading, file.path(args$output_dir, "functional_leading_edge_overlap.csv"))

cat("Writing figures...\n")
plot_outputs(concordance$summary, concordance$pathways, ex_aggregate, args$figures_dir)

cat("Functional concordance complete.\n")
