## This script generates supplementary figures not produced by previous scripts




##### Load Data #####

## Load species data (nutrient concentrations)
species <- read_csv("data/clean-data/02_SpeciesData.csv")

## Load catch data
fish <- read_csv("data/clean-data/03_FishData.csv")

## Load DRIs
dri <- read_csv("data/clean-data/02_DietaryReferenceIntakes.csv")

## Load EEZs
eez <- readRDS("data/clean-data/04_NationalCatchCompositions.rds")





##### Supplementary plots #####

# Nutrient concentrations by species size and feeding path

# Vector of species from eez data
eez <- do.call(rbind, eez)
captured.species <- unique(eez$valid_name)

# Add species from detailed catch data
captured.species <- append(captured.species, unique(fish$species))
captured.species <- unique(captured.species)

# Subset of species actually caught
species2 <- filter(species, species %in% captured.species)

# Capitalize feeding path
species2$feeding.path <- str_to_title(species2$feeding.path)

# Plot nutrient concentrations against size

# Calcium
ca.lmat <- glmmTMB(calcium_mg.100g ~ Lmat_cm * feeding.path,
  data = species2)
summary(ca.lmat)
pred.ca.lmat <- ggpredict(ca.lmat, terms = c("Lmat_cm", "feeding.path[Benthic, Pelagic]"))
a <- ggplot(pred.ca.lmat, aes(x = x, y = predicted, color = group)) +
  geom_point(data = species2,
    aes(x = Lmat_cm, y = calcium_mg.100g, color = feeding.path),
    alpha = 0.5) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, x = x),
    fill = "gray", linetype = 0, alpha = 0.5) +
  geom_line(linewidth = 1) +
  geom_hline(yintercept = dri$dri[dri$nutrients == "calcium_mg"],
    linetype = "dashed") +
  labs(title = "", color = "",
    x = expression(paste("Length at Maturity (", L[mat], ") (cm)", sep = "")),
    y = expression(paste("Calcium Concentration ", bgroup("(", frac('mg', '100g'), ")"), sep = ""))) +
  theme_pubr() +
  coord_cartesian(xlim = c(0, 100), ylim = c(0, 800)) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_color_npg() +
  theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
    text = element_text(size = 5))


# Iron
fe.lmat <- glmmTMB(iron_mg.100g ~ Lmat_cm * feeding.path,
  data = species2)
summary(fe.lmat)
pred.fe.lmat <- ggpredict(fe.lmat)[[2]]
b <- ggplot(species2, aes(x = feeding.path, y = iron_mg.100g, color = feeding.path)) +
  geom_violin(fill = "gray", alpha = 0.5) +
  geom_point(aes(x = x, y = predicted, color = x),
    data = pred.fe.lmat) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high, x = x, color = x),
    data = pred.fe.lmat,
    width = 0.2,
    inherit.aes = FALSE) +
  scale_color_npg() +
  geom_hline(yintercept = dri$dri[dri$nutrients == "iron_mg"],
    linetype = "dashed") +
  labs(title = "", color = "",
    x = "",
    y = expression(paste("Iron Concentration ", bgroup("(", frac('mg', '100g'), ")"), sep = ""))) +
  coord_cartesian(ylim = c(0, 8.5)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_pubr() +
  theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
    text = element_text(size = 5))


# Omega 3
o3.lmat <- glmmTMB(omega3_g.100g ~ Lmat_cm * feeding.path,
  data = species2)
summary(o3.lmat)
pred.o3.lmat <- ggpredict(o3.lmat, terms = c("Lmat_cm", "feeding.path[Benthic, Pelagic]"))
c <- ggplot(pred.o3.lmat, aes(x = x, y = predicted, color = group)) +
  geom_point(data = species2,
    aes(x = Lmat_cm, y = omega3_g.100g, color = feeding.path),
    alpha = 0.5) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, x = x),
    fill = "gray", linetype = 0, alpha = 0.5) +
  geom_line(linewidth = 1) +
  geom_hline(yintercept = dri$dri[dri$nutrients == "omega3_g"],
    linetype = "dashed") +
  labs(title = "", color = "",
    x = expression(paste("Length at Maturity (", L[mat], ") (cm)", sep = "")),
    y = expression(paste("Omega-3 Concentration ", bgroup("(", frac('g', '100g'), ")"), sep = ""))) +
  theme_pubr() +
  coord_cartesian(xlim = c(0, 100), ylim = c(0, 1)) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_color_npg() +
  theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
    text = element_text(size = 5))


# Selenium
se.lmat <- glmmTMB(selenium_ug.100g ~ Lmat_cm * feeding.path,
  data = species2)
summary(se.lmat)
pred.se.lmat <- ggpredict(se.lmat, terms = c("Lmat_cm", "feeding.path[Benthic, Pelagic]"))
d <- ggplot(pred.se.lmat, aes(x = x, y = predicted, color = group)) +
  geom_point(data = species2,
    aes(x = Lmat_cm, y = selenium_ug.100g, color = feeding.path),
    alpha = 0.5) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, x = x),
    fill = "gray", linetype = 0, alpha = 0.5) +
  geom_line(linewidth = 1) +
  geom_hline(yintercept = dri$dri[dri$nutrients == "selenium_ug"],
    linetype = "dashed") +
  labs(title = "",
    x = expression(paste("Length at Maturity (", L[mat], ") (cm)", sep = "")),
    y = expression(paste("Selenium Concentration ", bgroup("(", frac(paste("\u00b5", g, sep = ""), '100g'), ")"), sep = "")),
    color = NULL) +
  theme_pubr() +
  coord_cartesian(xlim = c(0, 100), ylim = c(0, 150)) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_color_npg() +
  theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
    text = element_text(size = 5))

# Vitamin A
va.lmat <- glmmTMB(vitamin.a_ug.100g ~ Lmat_cm * feeding.path,
  data = species2)
summary(va.lmat)
pred.va.lmat <- ggpredict(va.lmat)[[2]]
e <- ggplot(species2, aes(x = feeding.path, y = vitamin.a_ug.100g, color = feeding.path)) +
  geom_violin(fill = "gray", alpha = 0.5) +
  geom_point(aes(x = x, y = predicted, color = x),
    data = pred.va.lmat) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high, x = x, color = x),
    data = pred.va.lmat,
    width = 0.2,
    inherit.aes = FALSE) +
  scale_color_npg() +
  geom_hline(yintercept = dri$dri[dri$nutrients == "vitamina_ug"],
    linetype = "dashed") +
  labs(title = "", color = "",
    x = "",
    y = expression(paste("Vitamin A Concentration ", bgroup("(", frac(paste("\u00b5", g, sep = ""), '100g'), ")"), sep = ""))) +
  coord_cartesian(ylim = c(0, 500)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_pubr() +
  theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
    text = element_text(size = 5))


# Zinc
zn.lmat <- glmmTMB(zinc_mg.100g ~ Lmat_cm * feeding.path,
  data = species2)
summary(zn.lmat)
pred.zn.lmat <- ggpredict(zn.lmat, terms = c("Lmat_cm", "feeding.path[Benthic, Pelagic]"))
f <- ggplot(pred.zn.lmat, aes(x = x, y = predicted, color = group)) +
  geom_point(data = species2,
    aes(x = Lmat_cm, y = zinc_mg.100g, color = feeding.path),
    alpha = 0.5) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, x = x),
    fill = "gray", linetype = 0, alpha = 0.5) +
  geom_line(linewidth = 1) +
  geom_hline(yintercept = dri$dri[dri$nutrients == "zinc_mg"],
    linetype = "dashed") +
  labs(title = "", color = "",
    x = expression(paste("Length at Maturity (", L[mat], ") (cm)", sep = "")),
    y = expression(paste("Zinc Concentration ", bgroup("(", frac('mg', '100g'), ")"), sep = ""))) +
  theme_pubr() +
  coord_cartesian(xlim = c(0, 100), ylim = c(0, 4)) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_color_npg() +
  theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
    text = element_text(size = 5))


## Combine plots and save

# Combine
ggarrange(a, b, c, d, e, f,
  ncol = 2, nrow = 3,
  common.legend = TRUE, legend = "bottom")

# Save
ggsave("figures/07_SmallPelagicSpecies.pdf",
  width = 190, height = 270,
  units = "mm",
  dpi = 300)



