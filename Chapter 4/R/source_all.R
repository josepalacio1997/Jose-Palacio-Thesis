# source_all.R ------------------------------------------------------------
# Source every R file that defines the public API. After this, you have:
#   buffer_capture_analysis()    -- the orchestrator
#   save_buffer_analysis()       -- writes figures + RDS
#   default_colors() / default_labels()
#   generate_test_data() / make_nesting_outer()
#   plus the lower-level pieces in core.R
# -------------------------------------------------------------------------

for (f in c("R/core.R", "R/plot.R", "R/analysis.R", "R/generate_test_data.R")) {
  print(f)
  source(f)
}
