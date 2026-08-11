# =============================================================================
# house_style.R — Rice + HHD branded, COLORBLIND-FRIENDLY palette
# -----------------------------------------------------------------------------
# Diseno: la combinacion blue + orange es la mas robusta bajo TODOS los tipos
# de daltonismo (deuteranopia, protanopia, tritanopia) porque las longitudes
# de onda son maximamente separadas. Para palettes de >2 categorias se usa
# luminance contrast y/o linetype para distinguir niveles dentro de cada
# familia de color en vez de hues similares.
#
# Para mas de 6 categorias se anade el Okabe-Ito palette (CB-safe estandar).
#
# Uso:
#   source(file.path(HERE, "R", "house_style.R"))
#   ggplot(...) + scale_color_house() + theme_house()
#
# Para distinguir 3 niveles dentro de un paradigma sin perder CB-safety, usa
#   geom_line(aes(color = paradigm, linetype = model))
# en vez de 6 hues distintos.
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(scales)
})

# ---- Anclas institucionales (CB-safe blue/orange pair) ---------------------
RICE_BLUE       <- "#00205B"   # Rice navy primary (deep blue)
HHD_ORANGE      <- "#F58025"   # HHD orange primary

# Variantes con LUMINANCE alta (NO solo hue) para diferencias detectables CB:
# Mas blanco -> menos saturado -> mas claro. Esto produce diferencia
# perceivable bajo cualquier tipo de daltonismo.
RICE_BLUE_LIGHT <- "#7398C2"   # luminance ~ +30, hue conservado
RICE_BLUE_DARK  <- "#001640"   # muy oscuro, alto contraste con LIGHT

HHD_ORANGE_LIGHT<- "#FBB979"   # luminance ~ +30
HHD_ORANGE_DARK <- "#8C4D12"   # muy oscuro

# Neutros para overlays, NA fill, grids
NEUTRAL_GRAY    <- "#737373"
NEUTRAL_LIGHT   <- "#D9D9D9"
NEUTRAL_DARK    <- "#252525"

# ---- Okabe-Ito CB-safe palette (8 colores, todos distinguibles bajo CB) ----
# Reference: Okabe, M. & Ito, K. (2008). "Color Universal Design".
# Estos 8 colores son los recomendados internacionalmente para
# CB-accessible publishing. RICE_BLUE y HHD_ORANGE estan alineados
# con los anchors blue (#0072B2) y vermillion (#D55E00) del set.
OKABE_ITO <- c(
  "#000000",  # black
  "#E69F00",  # orange (similar to HHD)
  "#56B4E9",  # sky blue
  "#009E73",  # bluish green
  "#F0E442",  # yellow
  "#0072B2",  # blue (similar to Rice)
  "#D55E00",  # vermillion (orange-red)
  "#CC79A7"   # reddish purple
)

# ---- Palettes graduados por tamano -----------------------------------------
# 2 categorias: Rice blue + HHD orange (maxima CB-safety)
PALETTE_2 <- c(RICE_BLUE, HHD_ORANGE)

# 3 categorias: + neutro oscuro (gris). Evita anadir un 3er hue confundible.
PALETTE_3 <- c(RICE_BLUE, HHD_ORANGE, NEUTRAL_DARK)

# 4 categorias: + verde Okabe (bluish-green) que es distinguible bajo CB
PALETTE_4 <- c(RICE_BLUE, HHD_ORANGE, OKABE_ITO[4], NEUTRAL_DARK)

# 4 categorias para A/C x MARSS/FSV (Model B descartado 2026-06-07):
# blue (RICE) para A, orange (HHD) para C; luminance distinguishable
# para MARSS vs FSV dentro de cada modelo.
PALETTE_4_FITS <- c(
  RICE_BLUE,        # A_marss: navy
  RICE_BLUE_LIGHT,  # A_fsv:   light navy (luminance distinguishable)
  HHD_ORANGE_DARK,  # C_marss: dark orange
  HHD_ORANGE        # C_fsv:   orange
)

# Named para los 4 fits del comparativo de modelos (A vs C, ambos paradigmas)
PALETTE_FITS <- c(
  A_marss = RICE_BLUE,
  A_fsv   = RICE_BLUE_LIGHT,
  C_marss = HHD_ORANGE_DARK,
  C_fsv   = HHD_ORANGE
)

# Linetypes recomendados para distinguir A/C cuando se mantiene 1 color
# por paradigma (mejor opcion CB de todas)
LINETYPES_MODELS <- c(A = "solid", C = "dashed")

# ---- Scales -----------------------------------------------------------------
scale_color_house <- function(n = NULL, ...) {
  cols <- if (is.null(n))     PALETTE_6
          else if (n <= 2)    PALETTE_2
          else if (n <= 3)    PALETTE_3
          else if (n <= 4)    PALETTE_4
          else                PALETTE_6
  ggplot2::scale_color_manual(values = cols, ...)
}
scale_fill_house <- function(n = NULL, ...) {
  cols <- if (is.null(n))     PALETTE_6
          else if (n <= 2)    PALETTE_2
          else if (n <= 3)    PALETTE_3
          else if (n <= 4)    PALETTE_4
          else                PALETTE_6
  ggplot2::scale_fill_manual(values = cols, ...)
}
scale_color_fits <- function(...) {
  ggplot2::scale_color_manual(values = PALETTE_FITS, ...)
}
scale_fill_fits <- function(...) {
  ggplot2::scale_fill_manual(values = PALETTE_FITS, ...)
}
scale_color_okabe <- function(...) {
  ggplot2::scale_color_manual(values = OKABE_ITO, ...)
}
scale_fill_okabe <- function(...) {
  ggplot2::scale_fill_manual(values = OKABE_ITO, ...)
}

# Linetype para A/B/C (preferred CB-strategy)
scale_linetype_models <- function(...) {
  ggplot2::scale_linetype_manual(values = LINETYPES_MODELS, ...)
}

# Continuous: divergente CB-safe (blue-white-orange)
scale_fill_house_gradient <- function(midpoint = 0, ...) {
  ggplot2::scale_fill_gradient2(low = RICE_BLUE, mid = "white",
                                 high = HHD_ORANGE, midpoint = midpoint, ...)
}
# Sequential: blue or orange (one direction only)
scale_fill_house_seq_blue <- function(...) {
  ggplot2::scale_fill_gradient(low = "white", high = RICE_BLUE, ...)
}
scale_fill_house_seq_orange <- function(...) {
  ggplot2::scale_fill_gradient(low = "white", high = HHD_ORANGE, ...)
}

# Viridis tambien es CB-safe y bueno para heatmaps:
scale_fill_viridis_cb <- function(...) {
  ggplot2::scale_fill_viridis_c(option = "viridis", ...)
}

# ---- Theme ------------------------------------------------------------------
theme_house <- function(base_size = 11) {
  ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      plot.title       = element_text(face = "bold", color = RICE_BLUE_DARK),
      plot.subtitle    = element_text(color = NEUTRAL_GRAY),
      strip.background = element_rect(fill = "#F0F2F5", color = NA),
      strip.text       = element_text(face = "bold", color = RICE_BLUE_DARK),
      panel.grid.minor = element_blank(),
      legend.position  = "bottom",
      legend.title     = element_text(face = "bold")
    )
}
