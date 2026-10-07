#!/usr/bin/env Rscript
## Re-apply the §19 peak rules to a finished run, without re-running it.
##   Rscript tools/apply_peak_rules.R <run_dir> <cache_prefix> [driver]
## e.g.  Rscript tools/apply_peak_rules.R ~/Desktop/AFxQTL/hoffmanDrop/notAF/v8b v8b_
## The rule code and settings are read FROM THE DRIVER (default
## pipeline/xQTL_analysis8b.R), so there is one copy of them, not two.
## Writes <run_dir>/peakrules/.
a = commandArgs(TRUE)
run_dir = normalizePath(a[1]); prefix = a[2]
here = dirname(normalizePath(sub('--file=', '', grep('--file=', commandArgs(FALSE), value = TRUE))))
driver = if (length(a) >= 3) a[3] else file.path(here, '..', 'pipeline', 'xQTL_analysis8b.R')
ex = parse(driver)
lhs = function(e) if (is.call(e) && (identical(e[[1]], as.name('=')) || identical(e[[1]], as.name('<-'))) &&
                      is.name(e[[2]])) as.character(e[[2]]) else NA_character_
for (e in ex) { n = lhs(e); if (!is.na(n) && (n %in% c('qtl_peak_rules_run', 'LOD_threshold') || grepl('^rule_', n))) eval(e, globalenv()) }
qtl_peak_rules_run(RUN_DIR = run_dir, CACHE_PREFIX = prefix, OUT = file.path(run_dir, 'peakrules'),
                   LOD_THR = LOD_threshold, MERGE_VALLEY = rule_merge_valley,
                   DEAD_MIN_READS = rule_dead_min_reads, DEAD_WINDOW = rule_dead_window, DEAD_FRAC = rule_dead_frac,
                   FAR_AF = rule_far_af, FAR_EXCESS = rule_far_excess, SIG_FRAC = rule_sig_frac)
