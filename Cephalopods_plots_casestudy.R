###############################################################################
# North Sea, 2016-2025: Loligo vulgaris and Illex coindetii in DATRAS
#   1. Presence/absence maps (rows = quarter, columns = 5-year block)
#   2. Length-frequency plot (to check/adjust the length classes)
#   3. Length-class maps: CPUE (n/hour) per length class (rows = quarter,
#      columns = length class), all years pooled
#
# Self-contained: it downloads HH (haul info) and HL (length data) straight
# from DATRAS and caches them, so the haul universe for presence/absence and
# for the length maps is exactly the same.
###############################################################################

# install.packages(c("icesDatras", "worrms", "rnaturalearth", "rnaturalearthdata"))
library(tidyverse)
library(sf)
library(rnaturalearth)

## Parameters (REVIEW) --------------------------------------------------------
years        <- 2016:2025
surveys      <- "BTS"          # add "BTS", "SNS" if wanted (different gears -> CPUE not comparable)
quarters     <- 3                # quarters without data are skipped automatically
species      <- c("Loligo vulgaris", "Illex coindetii", "Loligo forbesii") 
period_width <- 5                  # years per column in the presence maps
ns_xlim      <- c(-4, 9)           # North Sea window
ns_ylim      <- c(51, 62)

# Length classes = dorsal mantle length in cm, lower bound inclusive.
# Defaults: <10 cm (small/recruits), 10-19 cm (medium), >=20 cm (large).
# Check the length-frequency plot (step 2) and adjust per species.
len_breaks <- list(
  "Loligo vulgaris" = c(0, 11, 20, Inf),
  "Illex coindetii" = c(0, 8, 17, Inf), 
  "Loligo forbesii" = c(0, 14, 20, Inf)
)

if (!exists("out_dir"))   out_dir   <- "./ADHOC/WGNSNETSEA/output/casestudy"
cache_dir <- "./ADHOC/WGNSNETSEA/cache_casestudy"
plot_dir  <- file.path(out_dir, "plots", "northsea_squid")
for (d in c(cache_dir, plot_dir)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

land <- ne_countries(scale = "medium", returnclass = "sf")

## Species -> AphiaID -----------------------------------------------------------
aphia <- map_dbl(species, ~ as.numeric(worrms::wm_name2id(.x)))
names(aphia) <- species
print(aphia)   # sanity check

# aphia <- c("Loligo vulgaris" = 140281, "Illex coindetii" = 140301, "Loligo forbesii" = 140270)  # fill in from WoRMS

## Download + cache -------------------------------------------------------------
num_cols <- c("Quarter", "Year", "HaulNo", "ShootLat", "ShootLong", "HaulDur", "SubFactor",
              "LngtClass", "HLNoAtLngt", "TotalNo", "Valid_Aphia")
chr_cols <- c("Survey", "Country", "Ship", "Gear", "StNo", "HaulVal", "DataType", "LngtCode")

get_cached <- function(record, survey, year, quarter) {
  f <- file.path(cache_dir, sprintf("%s_%s_%d_Q%d%s.rds", record, survey, year, quarter,
                                    if (record == "HL") "_ceph" else ""))
  if (file.exists(f)) return(readRDS(f))
  x <- tryCatch(icesDatras::getDATRAS(record, survey, year, quarter),
                error = function(e) NULL)
  if (is.null(x)) { warning("download failed: ", record, " ", survey, " ", year, " Q", quarter); return(tibble()) }
  x <- if (is.data.frame(x)) as_tibble(x) else tibble()      # no data for this survey/year/quarter
  if (nrow(x) > 0) {
    x <- x |>
      mutate(across(any_of(num_cols), as.numeric), across(any_of(chr_cols), as.character),
             StNo = as.character(StNo)) |>
      select(any_of(c(chr_cols, num_cols, "StNo")))
    if (record == "HL") x <- filter(x, Valid_Aphia %in% aphia)  # keep only our species
  }
  saveRDS(x, f); x
}

grid <- expand_grid(survey = surveys, year = years, quarter = quarters)
message("Downloading/reading ", nrow(grid), " survey-year-quarter combinations (cached after first run)")
hh_raw <- pmap(grid, \(survey, year, quarter) get_cached("HH", survey, year, quarter)) |> bind_rows()
hl_raw <- pmap(grid, \(survey, year, quarter) get_cached("HL", survey, year, quarter)) |> bind_rows()

mk_id <- function(d) with(d, paste(Survey, Year, Quarter, Country, Ship, Gear, StNo, HaulNo, sep = "_"))

## Hauls (the denominator for presence/absence) --------------------------------
hh <- hh_raw |>
  mutate(haul_id = mk_id(pick(everything()))) |>
  filter(HaulVal == "V", HaulDur > 0,
         between(ShootLat, ns_ylim[1], ns_ylim[2]),
         between(ShootLong, ns_xlim[1], ns_xlim[2]),
         !(ShootLong > 7 & ShootLat > 57.2)) |>      # drop Skagerrak
  distinct(haul_id, .keep_all = TRUE) |>
  transmute(haul_id, Survey, Year, lon = ShootLong, lat = ShootLat, HaulDur, DataType,
            quarter = paste0("Q", Quarter),
            p_start = years[1] + floor((Year - years[1]) / period_width) * period_width,
            period  = paste0(p_start, "-", pmin(p_start + period_width - 1, max(years))))
message(nrow(hh), " valid North Sea hauls")

hl <- hl_raw |>
  mutate(haul_id = mk_id(pick(everything()))) |>
  filter(haul_id %in% hh$haul_id) |>
  mutate(taxon = names(aphia)[match(Valid_Aphia, aphia)])

## 1. Presence/absence maps ---------------------------------------------------------
pres <- hl |> distinct(haul_id, taxon) |> mutate(present = TRUE)
pa <- expand_grid(hh, taxon = species) |>
  left_join(pres, by = c("haul_id", "taxon")) |>
  mutate(present = replace_na(present, FALSE))

map_layers <- list(
  geom_sf(data = land, fill = "grey85", colour = "grey50", linewidth = 0.1),
  coord_sf(xlim = ns_xlim, ylim = ns_ylim, expand = FALSE),
  theme_bw(base_size = 9),
  theme(axis.title = element_blank(), axis.text = element_text(size = 6))
)
safe_name <- function(x) gsub("[^A-Za-z0-9]+", "_", x)

for (tx in species) {
  d <- filter(pa, taxon == tx)
  p <- ggplot() +
    geom_point(data = filter(d, !present), aes(lon, lat), colour = "grey65", size = 0.3, alpha = 0.6) +
    geom_point(data = filter(d, present),  aes(lon, lat), colour = "red3", size = 1.8) +
    map_layers +
    facet_grid(quarter ~ period) +
    labs(title = paste(tx, "- presence/absence"),
         subtitle = sprintf("Present in %d of %d valid hauls (%s)", sum(d$present), nrow(d),
                            paste(surveys, collapse = ", ")),
         caption = "Grey = haul without the species; red = haul with the species.")
  ggsave(file.path(plot_dir, paste0(safe_name(tx), "_presence.png")), p,
         width = 3.2 * n_distinct(d$period) + 1, height = 3.6 * n_distinct(d$quarter) + 1, dpi = 150)
}

## 2. Length data -------------------------------------------------------------------
# Mantle length in cm. DATRAS: LngtCode "." and "0" are stored in mm, "1"/"2" in cm.
# CPUE (n/hour): DataType "C" is already per hour; otherwise raise by SubFactor and scale to 60 min.
len <- hl |>
  filter(HLNoAtLngt >= 0, LngtClass > 0) |>
  left_join(select(hh, haul_id, HaulDur, DataType, lon, lat, quarter, Year), by = "haul_id") |>
  mutate(len_cm = if_else(LngtCode %in% c(".", "0"), LngtClass / 10, LngtClass),
         n_h    = HLNoAtLngt * if_else(DataType == "C", 1, SubFactor * 60 / HaulDur))

n_nolen <- pres |> anti_join(distinct(len, haul_id, taxon), by = c("haul_id", "taxon")) |> count(taxon)
message("Hauls with the species present but no valid length data (excluded from length maps):")
print(n_nolen)

# length-frequency plot (use it to check the class limits)
lf <- len |> mutate(len_bin = floor(len_cm)) |> count(taxon, quarter, len_bin, wt = n_h, name = "n")
brk <- imap_dfr(len_breaks, \(b, nm) tibble(taxon = nm, x = b[is.finite(b) & b > 0]))
p_lf <- ggplot(lf, aes(len_bin, n)) +
  geom_col(width = 1, just = 0, fill = "grey40")+
  geom_vline(data = brk, aes(xintercept = x), colour = "red3", linetype = "dashed") +
  facet_grid(quarter ~ taxon, scales = "free_y") +
  labs(x = "Mantle length (cm)", y = "Sum of CPUE (n/hour) over hauls",
       title = "Length frequency, North Sea", caption = "Dashed lines = length-class limits.") +
  theme_bw(base_size = 9)
ggsave(file.path(plot_dir, "length_frequency.png"), p_lf, width = 9, height = 6, dpi = 150)

## 3. Length-class maps ---------------------------------------------------------------
class_labels <- function(b) {
  n <- length(b) - 1
  map_chr(seq_len(n), \(i) if (i == 1) paste0("<", b[2], " cm")
          else if (i == n) paste0("\u2265", b[i], " cm")
          else paste0(b[i], "-", b[i + 1] - 1, " cm"))
}

for (tx in species) {
  b   <- len_breaks[[tx]]
  lab <- class_labels(b)
  ld  <- filter(len, taxon == tx) |>
    mutate(len_class = cut(len_cm, breaks = b, labels = lab, right = FALSE))
  if (nrow(ld) == 0) { message("No length data for ", tx); next }
  
  pos <- ld |> group_by(haul_id, len_class) |> summarise(n_h = sum(n_h), .groups = "drop")
  all <- expand_grid(hh, len_class = factor(lab, levels = lab)) |>
    left_join(pos, by = c("haul_id", "len_class"))
  
  n_pos <- pos |> count(len_class) |> mutate(txt = paste0(len_class, ": ", n, " hauls"))
  message(tx, " - hauls per class: ", paste(n_pos$txt, collapse = "; "))
  
  p <- ggplot() +
    geom_point(data = filter(all, is.na(n_h)), aes(lon, lat), colour = "grey70", size = 0.3, alpha = 0.6) +
    geom_point(data = arrange(filter(all, !is.na(n_h)), n_h), aes(lon, lat, colour = n_h), size = 1.8) +
    map_layers +
    scale_colour_viridis_c(name = "CPUE\n(n/hour)", trans = "log10", option = "magma", end = 0.9) +
    facet_grid(quarter ~ len_class) +
    labs(title = paste(tx, "- CPUE by mantle length class"),
         subtitle = paste0(min(years), "-", max(years), ", ", paste(surveys, collapse = ", ")),
         caption = "Grey = haul without that length class. Hauls with no usable length data are not shown as present.")
  ggsave(file.path(plot_dir, paste0(safe_name(tx), "_length_classes.png")), p,
         width = 3.4 * length(lab) + 1.5, height = 3.6 * n_distinct(all$quarter) + 1, dpi = 150)
}

message("Plots saved in ", plot_dir)
