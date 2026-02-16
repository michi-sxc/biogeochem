# ==============================================================================
# FORMATTER v4 — Function Edition
# Michael Schneider, 2026
# https://github.com/michi-sxc
# Universal tool to convert sink() output from R to formatted markdown and PDF.
#
# USAGE:
#   source("format.R")
#   format_report("data/globocarb_r.txt")                     # minimal
#   format_report("data/globocarb_r.txt",                      # full options
#                 output_pdf  = "my_report.pdf",
#                 title       = "GLOBOCARB Analysis Report",
#                 subtitle    = "Biogeochemistry Seminar 2026",
#                 toc_only_first_page = TRUE,
#                 page_break_sections = TRUE,
#                 min_orphan_lines    = 4)
# ==============================================================================


# ══════════════════════════════════════════════════════════════════════════════
# 0. DEPENDENCY BOOTSTRAP (runs once when the file is sourced)
# ══════════════════════════════════════════════════════════════════════════════

local({
  packages <- c("rmarkdown", "here", "tinytex", "knitr", "kableExtra", "stringr")
  for (pkg in packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
  }
  if (!tinytex::is_tinytex()) tinytex::install_tinytex()
})


# ══════════════════════════════════════════════════════════════════════════════
# 1. INTERNAL HELPERS (not exported — used by format_report)
# ══════════════════════════════════════════════════════════════════════════════

# ── Line classifiers ─────────────────────────────────────────────────────────

.is_separator <- function(l) {
  grepl("^={4,}", l) || (grepl("^--- ", l) && grepl(" ---$", l))
}

.is_col_header <- function(l) {
  stripped <- trimws(l)
  if (nchar(stripped) < 3) return(FALSE)
  if (!grepl("^[A-Za-z._<]", stripped)) return(FALSE)
  tokens <- strsplit(stripped, "\\s+")[[1]]
  if (length(tokens) < 2) return(FALSE)
  non_num <- sum(is.na(suppressWarnings(as.numeric(tokens))))
  return(non_num / length(tokens) >= 0.4)
}

.is_data_row <- function(l) {
  grepl("^\\s*\\d+\\s+\\S", l)
}

.is_stats_line <- function(l) {
  tl <- trimws(l)
  grepl(paste0(
    "^(Call:|Residuals:|Coefficients:|Welch Two|Welch's|Paired|",
    "One-way analysis|Linear mixed model|Formula:|Data:|REML criterion|",
    "Scaled residuals:|Random effects:|Fixed effects:|Correlation of Fixed|",
    "Number of obs|Groups:|Signif\\. codes|Residual standard|",
    "Multiple R-squared|Adjusted R-squared|F-statistic|AIC|BIC|logLik|",
    "alternative hypothesis|95 percent confidence|sample estimates:|",
    "mean of|conf\\.level|data:|Spearman|Shapiro-Wilk|Levene)"
  ), tl)
}

.is_stats_continuation <- function(l) {
  tl <- trimws(l)
  grepl(paste0(
    "^(Min|1Q|Median|3Q|Max|Estimate|Std\\.|---$|\\(Intercept\\)|",
    "Decade|Region|Biome|Signif|Residual|Multiple|Adjusted|F-stat|",
    "Number|Groups:|Correlation|AIC|mean of|rho|data:|S =|",
    "alternative|true|95 per|sample|conf)"
  ), tl) ||
    grepl("\\*\\*\\*$|\\*\\*$|\\*$|\\.$", tl) ||
    (grepl("^-?[0-9]", tl) && nchar(tl) < 80)
}


# ── Console table parser (positional column detection) ────────────────────────

.parse_console_table <- function(lines_block) {
  lines_block <- lines_block[nchar(trimws(lines_block)) > 0]
  if (length(lines_block) < 2) return(NULL)
  
  # Find header row
  header_idx <- NA
  for (j in seq_along(lines_block)) {
    if (.is_col_header(lines_block[j])) { header_idx <- j; break }
  }
  if (is.na(header_idx)) return(NULL)
  
  header_line <- lines_block[header_idx]
  data_lines  <- lines_block[(header_idx + 1):length(lines_block)]
  data_lines  <- data_lines[nchar(trimws(data_lines)) > 0]
  if (length(data_lines) < 1) return(NULL)
  
  # Column positions from header tokens
  hdr_matches <- gregexpr("\\S+", header_line)[[1]]
  col_starts  <- as.integer(hdr_matches)
  col_widths  <- attr(hdr_matches, "match.length")
  col_names   <- mapply(function(s, w) substring(header_line, s, s + w - 1),
                        col_starts, col_widths, USE.NAMES = FALSE)
  
  n_cols <- length(col_names)
  if (n_cols < 2) return(NULL)
  
  col_bounds <- data.frame(
    start = col_starts,
    end   = c(col_starts[-1] - 1, 999)
  )
  
  # Positional row parser
  parse_row <- function(line) {
    idx_match <- regexpr("^\\s*\\d+\\s+", line)
    if (idx_match > 0) {
      idx_len <- attr(idx_match, "match.length")
      line_noindex <- paste0(strrep(" ", idx_len), substring(line, idx_len + 1))
    } else {
      line_noindex <- line
    }
    line_pad <- sprintf("%-999s", line_noindex)
    vals <- character(n_cols)
    for (k in seq_len(n_cols)) {
      vals[k] <- trimws(substring(line_pad, col_bounds$start[k], col_bounds$end[k]))
    }
    return(vals)
  }
  
  valid_lines <- data_lines[grepl("^\\s*\\d+\\s", data_lines)]
  if (length(valid_lines) < 1) return(NULL)
  
  rows <- lapply(valid_lines, parse_row)
  df   <- as.data.frame(do.call(rbind, rows), stringsAsFactors = FALSE)
  if (ncol(df) != n_cols) return(NULL)
  colnames(df) <- col_names
  
  # Coerce numeric columns
  for (col in seq_len(ncol(df))) {
    vals <- df[[col]]
    vals_clean <- vals[nchar(vals) > 0]
    if (length(vals_clean) == 0) next
    num_vals <- suppressWarnings(as.numeric(vals_clean))
    if (sum(!is.na(num_vals)) / length(vals_clean) > 0.7) {
      df[[col]] <- suppressWarnings(as.numeric(df[[col]]))
    }
  }
  
  # Sanity: reject if >40% empty
  total_cells <- nrow(df) * ncol(df)
  empty_cells <- sum(is.na(df) | df == "" | df == "NA")
  if (empty_cells / total_cells > 0.4) return(NULL)
  
  return(df)
}


# ── Render data.frame → kableExtra LaTeX ─────────────────────────────────────

.df_to_kable <- function(df, font_sz = 8, header_color = "#3c5488") {
  for (col in seq_len(ncol(df))) {
    if (is.numeric(df[[col]])) {
      max_val <- max(abs(df[[col]]), na.rm = TRUE)
      if (is.na(max_val) || is.infinite(max_val)) next
      dp <- ifelse(max_val < 1, 4, ifelse(max_val < 100, 3, 2))
      df[[col]] <- round(df[[col]], dp)
    }
  }
  
  k <- kableExtra::kable_styling(
    knitr::kable(df,
                 format    = "latex",
                 booktabs  = TRUE,
                 row.names = FALSE,
                 linesep   = ""),
    latex_options = c("striped", "scale_down", "HOLD_position"),
    font_size     = font_sz,
    stripe_color  = "gray!8"
  )
  k <- kableExtra::row_spec(k, 0, bold = TRUE, color = "white", background = header_color)
  return(as.character(k))
}


# ── Styled boxes ─────────────────────────────────────────────────────────────

.styled_mono_box <- function(lines_block, color = "gray") {
  paste0(
    "\n\n\\begin{tcolorbox}[",
    "colback=", color, "!3, colframe=", color, "!40, ",
    "boxrule=0.3pt, arc=2pt, left=5pt, right=5pt, top=3pt, bottom=3pt, ",
    "breakable]\n",
    "\\begin{verbatim}\n",
    paste(lines_block, collapse = "\n"), "\n",
    "\\end{verbatim}\n",
    "\\end{tcolorbox}\n\n"
  )
}

.styled_stats_box <- function(lines_block) {
  .styled_mono_box(lines_block, color = "blue")
}


# ── Key findings formatter ───────────────────────────────────────────────────

.format_key_findings <- function(lines_block) {
  joined <- paste(lines_block, collapse = " ")
  joined <- gsub("\\s+", " ", trimws(joined))
  
  parts <- strsplit(joined, "(?<=\\.)\\s+(?=\\d+\\.\\s)", perl = TRUE)[[1]]
  all_items <- c()
  for (part in parts) {
    sub_parts <- strsplit(part, "(?<=\\.)\\s*(?=\\d+\\.\\s)", perl = TRUE)[[1]]
    all_items <- c(all_items, sub_parts)
  }
  parts <- trimws(all_items)
  parts <- parts[nchar(parts) > 0]
  
  rmd <- c()
  for (part in parts) {
    if (grepl("^\\d+\\.", part)) {
      num  <- sub("^(\\d+)\\..+", "\\1", part)
      rest <- sub("^\\d+\\.\\s*", "", part)
      if (grepl(":", rest)) {
        label <- sub("^([^:]+):.*", "\\1", rest)
        body  <- sub("^[^:]+:\\s*", "", rest)
        rmd <- c(rmd, paste0(num, ". **", label, ":** ", body), "")
      } else {
        rmd <- c(rmd, paste0(num, ". ", rest), "")
      }
    } else {
      rmd <- c(rmd, paste0("*", part, "*"), "")
    }
  }
  return(rmd)
}


# ══════════════════════════════════════════════════════════════════════════════
# 2. MAIN FUNCTION
# ══════════════════════════════════════════════════════════════════════════════

#' Convert sink() console output to a formatted PDF report.
#'
#' @param input_file Path to the .txt file produced by sink().
#' @param output_pdf Path (or just filename) for the rendered PDF.
#'   If NULL, uses the same name/location as input_file but with .pdf extension.
#' @param title Report title (appears on cover page).
#' @param subtitle Report subtitle.
#' @param toc_only_first_page If TRUE (default), insert \\newpage after TOC.
#' @param page_break_sections If TRUE (default), major sections (# level)
#'   start on new pages.
#' @param min_orphan_lines Minimum lines that must stay with a heading at the
#'   bottom of a page. Enforced via LaTeX \\widowpenalty and \\clubpenalty.
#'   Default 4. Range 2-10.
#' @param header_color Hex color for table header bars. Default "#3c5488".
#' @param keep_rmd If TRUE, don't delete the intermediate .Rmd file.
#' @return Invisible path to the generated PDF.

format_report <- function(
    input_file,
    output_pdf          = NULL,
    title               = "Analysis Report",
    subtitle            = NULL,
    toc_only_first_page = TRUE,
    page_break_sections = TRUE,
    min_orphan_lines    = 4,
    header_color        = "#3c5488",
    keep_rmd            = FALSE
) {
  
  # ── Resolve paths ──────────────────────────────────────────────────────────
  if (!file.exists(input_file)) {
    # Try relative to here()
    alt <- here::here(input_file)
    if (file.exists(alt)) input_file <- alt
    else stop("Input file not found: ", input_file)
  }
  
  if (is.null(output_pdf)) {
    output_pdf <- sub("\\.[^.]+$", ".pdf", basename(input_file))
  }
  
  output_rmd <- sub("\\.pdf$", ".Rmd", output_pdf)
  # If output_pdf is just a filename (no dir), put it alongside the input
  if (dirname(output_pdf) == ".") {
    output_dir <- dirname(input_file)
    output_pdf <- file.path(output_dir, output_pdf)
    output_rmd <- file.path(output_dir, output_rmd)
  }
  
  # ── Build subtitle line ────────────────────────────────────────────────────
  subtitle_yaml <- if (!is.null(subtitle)) {
    paste0('subtitle: "', subtitle, '"')
  } else {
    ""
  }
  
  # ── YAML header ────────────────────────────────────────────────────────────
  rmd_lines <- c(
    "---",
    paste0('title: "', title, '"'),
    subtitle_yaml,
    "date: \"`r format(Sys.Date(), '%B %d, %Y')`\"",
    "output:",
    "  pdf_document:",
    "    latex_engine: xelatex",
    "    toc: true",
    "    toc_depth: 3",
    "    number_sections: false",
    "    highlight: tango",
    "geometry: margin=0.85in",
    "header-includes:",
    # Core table packages
    "  - \\usepackage{booktabs}",
    "  - \\usepackage{longtable}",
    "  - \\usepackage{array}",
    "  - \\usepackage{multirow}",
    "  - \\usepackage{float}",
    "  - \\usepackage{colortbl}",
    "  - \\usepackage{tabu}",
    "  - \\usepackage[table]{xcolor}",
    # Styled stat boxes
    "  - \\usepackage[breakable]{tcolorbox}",
    # Typography & spacing
    "  - \\usepackage{setspace}",
    "  - \\onehalfspacing",
    # Orphan / widow control
    paste0("  - \\widowpenalty=", min_orphan_lines * 2500),
    paste0("  - \\clubpenalty=",  min_orphan_lines * 2500),
    "  - \\brokenpenalty=10000",
    "  - \\predisplaypenalty=10000",
    # Float tuning
    "  - \\renewcommand{\\topfraction}{0.9}",
    "  - \\renewcommand{\\bottomfraction}{0.9}",
    "  - \\renewcommand{\\textfraction}{0.05}",
    "  - \\renewcommand{\\floatpagefraction}{0.8}",
    # Section styling
    "  - \\usepackage{titlesec}",
    paste0("  - \\titleformat{\\section}{\\Large\\bfseries\\color[HTML]{",
           gsub("#", "", header_color),
           "}}{\\thesection}{1em}{}[\\vspace{2pt}\\titlerule]"),
    "  - \\titleformat{\\subsection}{\\large\\bfseries\\color[HTML]{4dbbd5}}{\\thesubsection}{1em}{}",
    "  - \\titleformat{\\subsubsection}{\\normalsize\\bfseries\\color[HTML]{00a087}}{\\thesubsubsection}{1em}{}",
    "  - \\titlespacing*{\\section}{0pt}{20pt plus 4pt minus 2pt}{10pt plus 2pt minus 2pt}",
    "  - \\titlespacing*{\\subsection}{0pt}{16pt plus 3pt minus 2pt}{6pt plus 2pt minus 1pt}",
    "  - \\raggedbottom",
    # Header / footer
    "  - \\usepackage{fancyhdr}",
    "  - \\pagestyle{fancy}",
    "  - \\fancyhf{}",
    paste0("  - \\fancyhead[L]{\\small\\textcolor{gray}{", title, "}}"),
    if (!is.null(subtitle)) paste0("  - \\fancyhead[R]{\\small\\textcolor{gray}{", subtitle, "}}") else NULL,
    "  - \\fancyfoot[C]{\\thepage}",
    "  - \\renewcommand{\\headrulewidth}{0.4pt}",
    "  - \\fancypagestyle{plain}{\\fancyhf{}\\fancyfoot[C]{\\thepage}\\renewcommand{\\headrulewidth}{0pt}}",
    "---",
    "",
    "```{r setup, include=FALSE}",
    "knitr::opts_chunk$set(echo = FALSE, warning = FALSE, message = FALSE)",
    "library(kableExtra)",
    "```",
    ""
  )
  
  # Remove NULLs from conditional header lines
  rmd_lines <- rmd_lines[!sapply(rmd_lines, is.null)]
  
  if (toc_only_first_page) {
    rmd_lines <- c(rmd_lines, "\\newpage", "")
  }
  
  
  # ── Read input ─────────────────────────────────────────────────────────────
  raw_lines <- readLines(input_file)
  n_lines   <- length(raw_lines)
  
  
  # ── State machine ──────────────────────────────────────────────────────────
  # Using an environment so flush helpers can modify buffers without <<-
  env <- new.env(parent = emptyenv())
  env$text_buf  <- c()
  env$table_buf <- c()
  env$stats_buf <- c()
  env$kf_buf    <- c()
  env$mode      <- "text"
  
  section_count <- 0
  
  # Flush helpers (read/write via env$ — no <<- needed)
  flush_text <- function() {
    if (length(env$text_buf) == 0) return(character(0))
    out <- c(env$text_buf, "")
    env$text_buf <- c()
    return(out)
  }
  
  flush_table <- function() {
    if (length(env$table_buf) == 0) return(character(0))
    buf <- env$table_buf
    env$table_buf <- c()
    df <- tryCatch(.parse_console_table(buf), error = function(e) NULL)
    if (!is.null(df) && ncol(df) >= 2 && nrow(df) >= 1) {
      return(c("", .df_to_kable(df, header_color = header_color), ""))
    }
    return(.styled_mono_box(buf))
  }
  
  flush_stats <- function() {
    if (length(env$stats_buf) == 0) return(character(0))
    buf <- env$stats_buf
    env$stats_buf <- c()
    return(.styled_stats_box(buf))
  }
  
  flush_keyfindings <- function() {
    if (length(env$kf_buf) == 0) return(character(0))
    buf <- env$kf_buf
    env$kf_buf <- c()
    return(.format_key_findings(buf))
  }
  
  flush_all <- function() {
    c(flush_text(), flush_table(), flush_stats(), flush_keyfindings())
  }
  
  # Look-ahead helper
  peek_table <- function(start_idx) {
    count <- 0
    j <- start_idx
    while (j <= n_lines && j < start_idx + 5) {
      tl <- trimws(raw_lines[j])
      if (nchar(tl) == 0) { j <- j + 1; next }
      if (.is_data_row(raw_lines[j]) || .is_col_header(raw_lines[j])) count <- count + 1
      j <- j + 1
    }
    return(count >= 2)
  }
  
  
  # ════════════════════════════════════════════════════════════════════════════
  # MAIN LOOP
  # ════════════════════════════════════════════════════════════════════════════
  i <- 1
  while (i <= n_lines) {
    line      <- raw_lines[i]
    trim_line <- trimws(line)
    
    # ── A) Section separators ──────────────────────────────────────────────
    if (.is_separator(line)) {
      rmd_lines <- c(rmd_lines, flush_all())
      env$mode  <- "text"
      clean_header <- trimws(gsub("[=-]", "", line))
      
      if (grepl("^(A tibble|Signif)", clean_header, ignore.case = TRUE)) {
        i <- i + 1; next
      }
      if (nchar(clean_header) < 2) { i <- i + 1; next }
      
      section_count <- section_count + 1
      level <- if (grepl("^=", line)) "#" else "##"
      
      if (page_break_sections && level == "#" && section_count > 1) {
        rmd_lines <- c(rmd_lines, "", "\\newpage", "")
      }
      
      rmd_lines <- c(rmd_lines, paste(level, clean_header), "")
      
      if (grepl("KEY FINDINGS", clean_header, ignore.case = TRUE)) {
        env$mode <- "keyfindings"
      }
      i <- i + 1; next
    }
    
    # ── B) "Effect Size" / "A tibble" lines ────────────────────────────────
    if (env$mode == "text" && grepl("^(Effect Size|A tibble)", trim_line)) {
      rmd_lines <- c(rmd_lines, flush_all())
      rmd_lines <- c(rmd_lines, paste0("**", trim_line, "**"), "")
      i <- i + 1; next
    }
    
    # ── C) Empty lines ────────────────────────────────────────────────────
    if (trim_line == "") {
      if (env$mode == "table") {
        if (i + 1 <= n_lines && .is_data_row(raw_lines[i + 1])) {
          env$table_buf <- c(env$table_buf, line)
          i <- i + 1; next
        }
        rmd_lines <- c(rmd_lines, flush_table())
        env$mode <- "text"
      } else if (env$mode == "stats") {
        if (i + 1 <= n_lines) {
          next_trim <- trimws(raw_lines[i + 1])
          if (nchar(next_trim) > 0 &&
              (.is_stats_line(raw_lines[i + 1]) ||
               .is_stats_continuation(raw_lines[i + 1]))) {
            env$stats_buf <- c(env$stats_buf, "")
            i <- i + 1; next
          }
        }
        rmd_lines <- c(rmd_lines, flush_stats())
        env$mode <- "text"
      } else if (env$mode == "keyfindings") {
        env$kf_buf <- c(env$kf_buf, " ")
      } else {
        if (length(env$text_buf) > 0) rmd_lines <- c(rmd_lines, flush_text())
      }
      i <- i + 1; next
    }
    
    # ── D) Key findings accumulator ───────────────────────────────────────
    if (env$mode == "keyfindings") {
      env$kf_buf <- c(env$kf_buf, line)
      i <- i + 1; next
    }
    
    # ── E) Stats detection & continuation ─────────────────────────────────
    if (env$mode == "text" && .is_stats_line(line)) {
      rmd_lines <- c(rmd_lines, flush_text())
      env$stats_buf <- c(line)
      env$mode <- "stats"
      i <- i + 1; next
    }
    if (env$mode == "stats") {
      env$stats_buf <- c(env$stats_buf, line)
      i <- i + 1; next
    }
    
    # ── F) Table detection ────────────────────────────────────────────────
    if (env$mode == "text") {
      start_table <- FALSE
      
      if (.is_col_header(line) && i + 1 <= n_lines && .is_data_row(raw_lines[i + 1])) {
        start_table <- TRUE
      }
      if (!start_table && .is_data_row(line) && peek_table(i)) {
        for (back in (i - 1):max(1, i - 3)) {
          bl <- trimws(raw_lines[back])
          if (nchar(bl) > 0) {
            if (.is_col_header(raw_lines[back])) {
              rmd_lines <- c(rmd_lines, flush_text())
              env$table_buf <- c(raw_lines[back])
            }
            break
          }
        }
        start_table <- TRUE
      }
      
      if (start_table) {
        rmd_lines <- c(rmd_lines, flush_text())
        env$table_buf <- c(env$table_buf, line)
        env$mode <- "table"
        i <- i + 1; next
      }
    }
    
    if (env$mode == "table") {
      if (.is_data_row(line) || .is_col_header(line)) {
        env$table_buf <- c(env$table_buf, line)
      } else {
        rmd_lines <- c(rmd_lines, flush_table())
        env$mode <- "text"
        next  # re-process line
      }
      i <- i + 1; next
    }
    
    # ── G) summary() output ───────────────────────────────────────────────
    if (env$mode == "text" && grepl("Min\\.", line) && grepl("(Median|1st Qu)", line)) {
      rmd_lines <- c(rmd_lines, flush_text())
      summary_buf <- c(line)
      j <- i + 1
      while (j <= n_lines && nchar(trimws(raw_lines[j])) > 0) {
        summary_buf <- c(summary_buf, raw_lines[j])
        j <- j + 1
      }
      rmd_lines <- c(rmd_lines, .styled_mono_box(summary_buf))
      i <- j; next
    }
    
    # ── H) str() output ──────────────────────────────────────────────────
    if (env$mode == "text" && grepl("^'data\\.frame'", trim_line)) {
      rmd_lines <- c(rmd_lines, flush_text())
      str_buf <- c(line)
      j <- i + 1
      while (j <= n_lines && (grepl("^\\s*\\$", raw_lines[j]) || nchar(trimws(raw_lines[j])) == 0)) {
        if (nchar(trimws(raw_lines[j])) == 0) break
        str_buf <- c(str_buf, raw_lines[j])
        j <- j + 1
      }
      rmd_lines <- c(rmd_lines, .styled_mono_box(str_buf, "teal"))
      i <- j; next
    }
    
    # ── I) Plain text (default) ──────────────────────────────────────────
    formatted <- gsub("^\\[1\\]\\s*", "", line)
    env$text_buf <- c(env$text_buf, formatted, "")
    i <- i + 1
  }
  
  # Final flush
  rmd_lines <- c(rmd_lines, flush_all())
  
  
  # ── Post-processing: collapse excessive blank lines ────────────────────────
  clean_rmd   <- c()
  blank_count <- 0
  for (line in rmd_lines) {
    if (is.null(line) || length(line) == 0) next
    if (trimws(line) == "") {
      blank_count <- blank_count + 1
      if (blank_count <= 2) clean_rmd <- c(clean_rmd, line)
    } else {
      blank_count <- 0
      clean_rmd <- c(clean_rmd, line)
    }
  }
  rmd_lines <- clean_rmd
  
  
  # ── Write & render ─────────────────────────────────────────────────────────
  writeLines(rmd_lines, output_rmd)
  message("\u2713 RMarkdown generated: ", output_rmd)
  message("  Rendering PDF ...")
  
  result <- tryCatch({
    rmarkdown::render(output_rmd, output_file = basename(output_pdf),
                      output_dir = dirname(output_pdf), quiet = TRUE)
    message("\n\u2713 SUCCESS \u2014 Report saved to: ", output_pdf)
    output_pdf
  }, error = function(e) {
    message("\n\u2717 Rendering Error: ", e$message)
    message("  Tip: tinytex::tlmgr_install('tcolorbox')")
    NULL
  })
  
  # Cleanup
  if (!keep_rmd && file.exists(output_rmd)) {
    file.remove(output_rmd)
  }
  
  return(invisible(result))
}