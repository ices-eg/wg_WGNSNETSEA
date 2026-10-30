#################################################
# Check age-length structure for plaice: BTS Q3 and NS-IBTS Q1 and Q3
# Change surv / qtr below to rerun for another survey x quarter
#################################################

library(icesDatras); library(dplyr); library(tidyr); library(purrr); library(ggplot2)

# 1. Data (2016 to 2025) -------------------------------------------------
years <- 2016:2025
surv  <- "BTS"      # "BTS" or "NS-IBTS"
qtr   <- 3          # BTS: 3 ; NS-IBTS: 1 or 3

ca_raw <- getDATRAS("CA", survey = surv, years = years, quarters = qtr)
head(ca_raw)
table(ca_raw$CANoAtLngt)     # check: mostly 1 = one aged fish per row
unique(ca_raw$LngtCode)      # check length units

rect_to_lonlat <- function(code) {
  letter <- substr(code, 3, 3)
  lon0 <- ifelse(letter == "A", -44,
                 -40 + 10 * (match(letter, c("B","C","D","E","F","G","H","J","K","L","M")) - 1))
  tibble(lon = lon0 + as.numeric(substr(code, 4, 4)) + 0.5,
         lat = 35.75 + 0.5 * as.numeric(substr(code, 1, 2)))
}

ca <- ca_raw |>
  filter(Valid_Aphia == 127143, AreaType == "0", Age >= 0, CANoAtLngt > 0, !is.na(LngtClass)) |>
  mutate(length_cm = case_when(LngtCode %in% c(".", "0") ~ LngtClass / 10,   # mm -> cm
                               TRUE                       ~ LngtClass),  # codes 1, 2, 5 are already in cm
         age     = Age,
         NoAtALK = CANoAtLngt,                    # number of aged fish in this row
         rect    = AreaCode,
         sex     = ifelse(Sex %in% c("F", "M"), Sex, "U"),
         haul_id = paste(Survey, Quarter, Year, Country, Ship, StNo, HaulNo, sep = "_"))

ca <- ca |>
  bind_cols(rect_to_lonlat(ca$rect)) |>
  # crude study-area filter: drop Kattegat (verify against an ICES area map)
  filter(!(lat < 57.5 & lon > 10), lat > 51, lat < 62, lon > -4, lon < 13)

ca_sub <- ca |> filter(Survey == surv, Quarter == qtr)
table(ca_sub$Survey, ca_sub$Quarter)   # sanity check
nrow(ca_sub)

# 2. Group neighbouring rectangles until each has enough fish -------------

# 2a. The function ---------------------------------------------------------
# Merges neighbouring ICES rectangles until every region has at least
# `min_per_age` aged fish in EACH age from `min_age` to `plus_age`
# (all fish older than plus_age are pooled into a plus group).
# `max_rect` optionally stops regions growing beyond that many rectangles
# (only respected when a merge partner under the cap exists).

group_rectangles <- function(ca_sub, min_age = 1, plus_age = 10,
                             min_per_age = 30, max_rect = Inf) {
  
  ages_req <- min_age:plus_age
  
  # aged fish per rectangle x age (plus group applied here)
  cnt <- ca_sub |>
    mutate(age = factor(pmin(age, plus_age), levels = ages_req)) |>
    filter(!is.na(age)) |>
    count(rect, lon, lat, age, wt = NoAtALK) |>
    pivot_wider(names_from = age, values_from = n,
                values_fill = 0, names_expand = TRUE) |>
    arrange(rect)
  M <- as.matrix(cnt[, as.character(ages_req)])
  
  # queen adjacency; lat x2 so rectangles (1 deg lon x 0.5 deg lat) are ~square
  xy  <- cbind(cnt$lon, cnt$lat * 2)
  d   <- as.matrix(dist(xy))
  adj <- d <= 1.5 & d > 0
  
  grp <- seq_len(nrow(cnt))            # start: every rectangle is its own group
  
  repeat {
    gt <- rowsum(M, grp)                               # group x age counts
    ok <- apply(gt >= min_per_age, 1, all)             # group passes if ALL ages >= min
    if (all(ok)) break
    if (nrow(gt) == 1) { warning("Not enough data even when pooled"); break }
    
    bad <- as.integer(rownames(gt))[!ok]
    g   <- bad[which.min(rowSums(gt)[as.character(bad)])]   # smallest failing group
    mem <- which(grp == g)                                  # its rectangles
    oth <- which(grp != g)                                  # all other rectangles
    
    # rectangles that touch the group (or the nearest one if it is isolated)
    nb <- oth[colSums(adj[mem, oth, drop = FALSE]) > 0]
    if (length(nb) == 0)
      nb <- oth[which.min(apply(d[mem, oth, drop = FALSE], 2, min))]
    
    cand <- unique(grp[nb])                                 # neighbouring groups
    
    # optional size cap (falls back to all neighbours if none would stay under it)
    sz    <- as.integer(table(grp)[as.character(cand)])
    ok_sz <- cand[sz + length(mem) <= max_rect]
    if (length(ok_sz) > 0) cand <- ok_sz
    
    tgt <- cand[which.min(rowSums(gt)[as.character(cand)])] # merge with smallest neighbour
    grp[grp == g] <- tgt
  }
  
  cnt |> mutate(region = match(grp, unique(grp))) |> select(rect, lon, lat, region)
}

# 2b. How many aged fish per age? (helps choose plus_age) -------------------
ca_sub |> count(age, wt = NoAtALK) |> as_tibble() |> print(n = 30)

# 2c. Sensitivity: how many regions for different settings? ----------------
sens <- expand.grid(plus_age = c(5, 8, 10), min_n = c(20, 30, 50)) |>
  pmap_dfr(function(plus_age, min_n) {
    rm <- group_rectangles(ca_sub, min_age = 1, plus_age = plus_age, min_per_age = min_n)
    tibble(plus_age = plus_age, min_n = min_n, n_regions = n_distinct(rm$region))
  })
sens

# 2d. Final regions (chosen settings) ----------------------------------------
plus_age    <- 10
min_per_age <- 30

reg_map <- group_rectangles(ca_sub, min_age = 1, plus_age = plus_age,
                            min_per_age = min_per_age)
# Optional: limit region size, e.g. max_rect = 10
# reg_map <- group_rectangles(ca_sub, 1, plus_age, min_per_age, max_rect = 10)

n_distinct(reg_map$region)

# 2e. Checks ------------------------------------------------------------------
# Fish per region x age: every cell for ages 1..plus_age should be >= min_per_age
cell_counts <- ca_sub |>
  mutate(age = pmin(age, plus_age)) |>
  inner_join(reg_map |> select(rect, region), by = "rect") |>
  count(region, age, wt = NoAtALK) |>
  pivot_wider(names_from = age, values_from = n, values_fill = 0)
print(cell_counts, n = Inf)

cell_counts |> select(all_of(as.character(1:plus_age))) |> min()   # should be >= min_per_age

# Region sizes (number of rectangles, centre of region)
reg_map |>
  group_by(region) |>
  summarise(n_rect = n(), lat = mean(lat), lon = mean(lon)) |>
  arrange(desc(n_rect)) |>
  print(n = Inf)

# 2f. Map (shuffled colours + region numbers; needs rnaturalearth + rnaturalearthdata + sf)
land <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")

set.seed(1)
reg_plot <- reg_map |>
  mutate(region_f = factor(region, levels = sample(unique(region))))

ggplot() + geom_sf(data = land) +
  geom_tile(data = reg_plot, aes(lon, lat, fill = region_f),
            width = 1, height = 0.5, colour = "white") +
  geom_text(data = reg_plot, aes(lon, lat, label = region), size = 2.5) +
  coord_sf(xlim = c(-4, 13), ylim = c(51, 62)) +
  guides(fill = "none")

ggsave("./ADHOC/WGNSNETSEA/output_PLEotoliths/clusters_geogrsamplesize.png", width = 11, height = 5, dpi = 300)






# 3. Test age-length variability among regions ----------------------------
library(lme4); library(lmerTest); library(emmeans)

fish <- ca_sub |>
  inner_join(reg_map |> select(rect, region), by = "rect") |>
  filter(age >= 1) |>                              # age 0 was not part of the grouping rule
  uncount(round(NoAtALK)) |>                       # one row per fish
  mutate(region = factor(region),
         age_f  = factor(pmin(age, plus_age), levels = 1:plus_age,
                         labels = c(1:(plus_age - 1), paste0(plus_age, "+"))),
         year_f = factor(Year)) |>
  droplevels()

table(fish$age_f, fish$region)                     # quick check of the cells

m_full  <- lmer(length_cm ~ age_f * region + sex + (1 | year_f) + (1 | haul_id),
                data = fish, REML = FALSE)
m_noint <- update(m_full, . ~ . - age_f:region)
m_noreg <- update(m_noint, . ~ . - region)

anova(m_noreg, m_noint, m_full)    # region main effect, then region x age interaction

emm_options(lmer.df = "asymptotic")                # avoids very slow df calculations on a large dataset
emm <- emmeans(m_full, ~ region | age_f)
plot(emm)
contrast(emm, method = "eff")                      # each region vs. the mean of all regions, per age


eff <- as.data.frame(contrast(emm, method = "eff")) |>
  mutate(region = as.integer(gsub("region| effect", "", contrast)))

map_eff <- reg_map |>
  left_join(eff |> filter(age_f %in% c("3", "6", "10+")), by = "region")

ggplot() + geom_sf(data = land) +
  geom_tile(data = map_eff, aes(lon, lat, fill = estimate), width = 1, height = 0.5) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", name = "cm vs mean") +
  facet_wrap(~ age_f) +
  coord_sf(xlim = c(-4, 13), ylim = c(51, 62))

# Red: fish in that area are longer than average at that age (up to about +4 cm in the 10+ panel).
# Blue: fish are shorter than average (down to about −4 cm).
# White or pale: about average, so no real difference.

# difference in mean length at that age, in cm, between a region and the average of all 11 regions
ggplot() + geom_sf(data = land) +
  geom_tile(data = map_eff, aes(lon, lat, fill = estimate), width = 1, height = 0.5) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                       limits = c(-5, 5), oob = scales::squish,
                       name = "Length-at-age\ndifference (cm)\nvs. average region") +
  facet_wrap(~ age_f, ncol = 3) +
  coord_sf(xlim = c(-4, 13), ylim = c(51, 62)) +
  labs(title = "Plaice length-at-age by region, relative to the average region",
       subtitle = "Red = longer than average, blue = shorter")
ggsave("./ADHOC/WGNSNETSEA/output_PLEotoliths/diffmLatA.png", width = 8,height = 6,dpi = 300)


head(map_eff)
table(map_eff$region)

#identify rectangles

# 4. Growth curves per region (von Bertalanffy, Kimura LRT) ---------------
fish_vb <- fish |>
  filter(age <= 20) |>                         # a handful of very old fish carry almost no data
  mutate(reg = as.integer(region)) |>
  droplevels()
nr <- nlevels(fish_vb$region)

vb_common <- nls(length_cm ~ Linf * (1 - exp(-K * (age - t0))),
                 data = fish_vb, start = list(Linf = 50, K = 0.2, t0 = 0))

vb_full <- nls(length_cm ~ Linf[reg] * (1 - exp(-K[reg] * (age - t0[reg]))),
               data = fish_vb,
               start = list(Linf = rep(coef(vb_common)["Linf"], nr),
                            K    = rep(coef(vb_common)["K"],    nr),
                            t0   = rep(coef(vb_common)["t0"],   nr)))

anova(vb_common, vb_full)

# Letting each region have its own curve reduces the residual sum of squares from 541,834 to 481,830, and the F-test is highly significant (F = 118.7, p < 2.2e-16). So the regions differ in their growth curves, which agrees with the mixed-model result in section 3.
# So most of the variation in length-at-age is individual variation within regions and ages, not a regional effect. Regional shifts of 2 to 4 cm, as seen in the maps, are real but small compared with the spread among individual fish of the same age (SD about 4 cm).

# Tidy table of parameters per region
co <- summary(vb_full)$coefficients
vb_tab <- tibble(term = rownames(co), est = co[, 1], se = co[, 2]) |>
  mutate(par    = sub("\\[.*", "", term),
         reg    = as.integer(sub(".*\\[(\\d+)\\]", "\\1", term)),
         region = levels(fish_vb$region)[reg]) |>
  select(region, par, est, se) |>
  pivot_wider(names_from = par, values_from = c(est, se))
print(vb_tab, n = Inf)


st <- as.list(coef(vb_common))

m_Linf <- nls(length_cm ~ Linf[reg] * (1 - exp(-K * (age - t0))),
              data = fish_vb, start = list(Linf = rep(st$Linf, nr), K = st$K, t0 = st$t0))
m_K    <- nls(length_cm ~ Linf * (1 - exp(-K[reg] * (age - t0))),
              data = fish_vb, start = list(Linf = st$Linf, K = rep(st$K, nr), t0 = st$t0))
m_t0   <- nls(length_cm ~ Linf * (1 - exp(-K * (age - t0[reg]))),
              data = fish_vb, start = list(Linf = st$Linf, K = st$K, t0 = rep(st$t0, nr)))
m_LK   <- nls(length_cm ~ Linf[reg] * (1 - exp(-K[reg] * (age - t0))),
              data = fish_vb, start = list(Linf = rep(st$Linf, nr), K = rep(st$K, nr), t0 = st$t0))

anova(vb_common, m_Linf, vb_full)
anova(vb_common, m_K,    vb_full)
anova(vb_common, m_t0,   vb_full)
anova(vb_common, m_LK,   vb_full)    # does letting both Linf and K vary capture most of it?

AIC(vb_common, m_Linf, m_K, m_t0, m_LK, vb_full)   # parsimonious comparison


newd <- expand.grid(reg = 1:nr, age = seq(1, 20, 0.25)) |>
  mutate(region = levels(fish_vb$region)[reg])
newd$pred <- predict(vb_full, newdata = newd)

ggplot(newd, aes(age, pred, colour = region)) + geom_line() +
  coord_cartesian(xlim = c(1, 15)) + labs(x = "Age", y = "Length (cm)")

# 
# vb_by_region <- fish_vb |>
#   group_split(region) |>
#   map(~ try(nls(length_cm ~ Linf * (1 - exp(-K * (age - t0))),
#                 data = .x, start = as.list(coef(vb_common))), silent = TRUE))


#### identify clusters: 
# 1. Cluster the regions on their profile across ages
library(cluster)

# 'eff' = the region x age contrast table made earlier
# (as.data.frame(contrast(emm, method = "eff")) with 'region' parsed from the contrast name)
ages_use <- c("2","3","4","5","6","7","8","9","10+")     # age 1 left out: region 8 has only 30 fish

eff_wide <- eff |>
  filter(age_f %in% ages_use) |>
  select(region, age_f, estimate) |>
  pivot_wider(names_from = age_f, values_from = estimate) |>
  arrange(region)

X <- as.matrix(eff_wide[, ages_use])
rownames(X) <- eff_wide$region

d  <- dist(X)                          # Euclidean distance between profiles, in cm
hc <- hclust(d, method = "ward.D2")
plot(hc)                               # look for natural breaks

# silhouette width for k = 2 to 6 (noisy with only 11 objects, use as a guide)
sapply(2:6, function(k) mean(silhouette(cutree(hc, k), d)[, 3]))



# 2. Assign groups to regions, then to rectangles
k <- 3                                               # choose from dendrogram / silhouette / the map

grp_tab <- tibble(region = as.integer(rownames(X)), cl = cutree(hc, k))

# Name groups by their mean effect (low -> high) so labels are ordered, not arbitrary
ord <- grp_tab |>
  left_join(tibble(region = as.integer(rownames(X)), m = rowMeans(X)), by = "region") |>
  group_by(cl) |> summarise(m = mean(m)) |> arrange(m) |>
  mutate(group = factor(seq_len(n()), labels = c("Low", "Mid", "High")[seq_len(n())]))
# for k != 3 use e.g. labels = paste0("G", seq_len(k))

rect_groups <- reg_map |>
  left_join(grp_tab, by = "region") |>
  left_join(ord |> select(cl, group), by = "cl") |>
  select(rect, lon, lat, region, group)

table(rect_groups$group)                   # number of rectangles per group
write.csv(rect_groups, "./ADHOC/WGNSNETSEA/output_PLEotoliths/plaice_rectangle_groups.csv", row.names = FALSE)


# 3. Check the result
# Map of the groups
ggplot() + geom_sf(data = land) +
  geom_tile(data = rect_groups, aes(lon, lat, fill = group), width = 1, height = 0.5, colour = "white") +
  scale_fill_manual(values = c(Low = "#4575b4", Mid = "grey85", High = "#d73027")) +
  coord_sf(xlim = c(-4, 13), ylim = c(51, 62))

# Profiles: thin lines = regions, thick = group mean
eff |>
  filter(age_f %in% ages_use) |>
  mutate(region = as.integer(region)) |>
  left_join(rect_groups |> distinct(region, group), by = "region") |>
  ggplot(aes(age_f, estimate, group = region, colour = group)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_line(alpha = 0.6) +
  stat_summary(aes(group = group), fun = mean, geom = "line", linewidth = 1.5) +
  scale_colour_manual(values = c(Low = "#4575b4", Mid = "grey50", High = "#d73027")) +
  labs(x = "Age", y = "Difference vs. average region (cm)")


fish <- fish |> left_join(rect_groups |> select(rect, group), by = "rect")

m_grp <- lmer(length_cm ~ age_f * group + sex + (1 | year_f) + (1 | haul_id),
              data = fish, REML = FALSE)
emmeans(m_grp, ~ group | age_f)



cols <- c(Low = "#4575b4", Mid = "grey85", High = "#d73027")

p_groups <- ggplot() +
  geom_sf(data = land, fill = "grey95", colour = "grey60") +
  geom_tile(data = rect_groups, aes(lon, lat, fill = group),
            width = 1, height = 0.5, colour = "white", linewidth = 0.2) +
  scale_fill_manual(values = cols, name = "Length-at-age\ngroup") +
  coord_sf(xlim = c(-4, 13), ylim = c(51, 62)) +
  labs(x = "Longitude", y = "Latitude",
       title = "Plaice length-at-age groups (BTS Q3, 2016-2025)") +
  theme_bw()

p_groups

# # Optional: add the number of each original region
# lab <- rect_groups |>
#   group_by(region, group) |>
#   summarise(lon = mean(lon), lat = mean(lat), .groups = "drop")
# 
# p_groups + geom_text(data = lab, aes(lon, lat, label = region), size = 3)

ggsave("./ADHOC/WGNSNETSEA/output_PLEotoliths/plaice_groups_map.png", p_groups, width = 7, height = 6, dpi = 300)


emm_g <- as.data.frame(emmeans(m_grp, ~ group | age_f)) |>
  filter(age_f %in% c("3", "6", "10+")) |>
  mutate(age_f = factor(age_f, levels = c("3", "6", "10+")))

map_len <- rect_groups |> left_join(emm_g, by = "group")

p_len <- ggplot() +
  geom_sf(data = land, fill = "grey95", colour = "grey60") +
  geom_tile(data = map_len, aes(lon, lat, fill = emmean),
            width = 1, height = 0.5, colour = "white", linewidth = 0.2) +
  scale_fill_viridis_c(name = "Mean length\n(cm)", option = "C") +
  facet_wrap(~ age_f, labeller = labeller(age_f = function(x) paste("Age", x))) +
  coord_sf(xlim = c(-4, 13), ylim = c(51, 62)) +
  labs(x = "Longitude", y = "Latitude") +
  theme_bw()

p_len
ggsave("./ADHOC/WGNSNETSEA/output_PLEotoliths/plaice_length_at_age_map.png", p_len, width = 11, height = 5, dpi = 300)






# library(ClustGeo)
# # install.packages("ClustGeo")
# 
# cen <- reg_map |> group_by(region) |> summarise(lon = mean(lon), lat = mean(lat)) |> arrange(region)
# D0  <- dist(cbind(cen$lon, cen$lat * 2))
# cr  <- choicealpha(D0, d, range.alpha = seq(0, 1, 0.1), K = k)   # pick alpha from the plot
# hc2 <- hclustgeo(d, D0, alpha = 0.3)                              # use cutree(hc2, k) instead of hc




# # 5. Optional: age composition via age-length keys -------------------------
# key <- c("Survey","Quarter","Country","Ship","Gear","StNo","HaulNo","Year")
# 
# hh <- getDATRAS("HH", survey = surv, years = years, quarters = qtr)
# hl <- getDATRAS("HL", survey = surv, years = years, quarters = qtr)
# 
# table(hl$DataType)   # check: "R" = raw per haul, "C" = already per hour
# 
# lf <- hl |>
#   filter(Valid_Aphia == 127143, !is.na(LngtClass)) |>
#   inner_join(hh |> select(all_of(key), StatRec, HaulDur), by = key) |>
#   mutate(len  = floor(ifelse(LngtCode == ".", LngtClass / 10, LngtClass)),
#          cpue = case_when(DataType == "C" ~ HLNoAtLngt,
#                           DataType == "R" ~ HLNoAtLngt * SubFactor * 60 / HaulDur)) |>
#   inner_join(reg_map |> select(rect, region), by = c("StatRec" = "rect")) |>
#   group_by(region, len) |> summarise(N = sum(cpue, na.rm = TRUE), .groups = "drop")
# 
# alk <- ca_sub |>
#   inner_join(reg_map |> select(rect, region), by = "rect") |>
#   mutate(len = floor(length_cm)) |>
#   group_by(region, len, age) |> summarise(n = sum(NoAtALK), .groups = "drop") |>
#   group_by(region, len) |> mutate(p = n / sum(n)) |> ungroup()
# 
# age_comp <- lf |>
#   inner_join(alk, by = c("region", "len"), relationship = "many-to-many") |>
#   group_by(region, age) |> summarise(N = sum(N * p), .groups = "drop") |>
#   group_by(region) |> mutate(prop = N / sum(N)) |> ungroup()
# 
# # Share of the catch whose length classes have no key
# unkeyed <- lf |>
#   anti_join(alk |> distinct(region, len), by = c("region", "len")) |>
#   summarise(x = sum(N)) |> pull(x)
# unkeyed / sum(lf$N)