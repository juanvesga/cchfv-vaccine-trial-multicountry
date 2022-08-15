
# 12/08/2022
# JF Vesga - London School of Hygiene and Tropical Medicine
# This script combines, plots and saves manuscript plots 2 to 4



rm(list = ls())
graphics.off()
library(ggplot2)
library(gridExtra)
library(here)


# Create user-defined function, which extracts legends from ggplots
extract_legend <- function(my_ggp) {
  step1 <- ggplot_gtable(ggplot_build(my_ggp))
  step2 <- which(sapply(step1$grobs, function(x) x$name) == "guide-box")
  step3 <- step1$grobs[[step2]]
  return(step3)
}



country <- "AFG"
load(here("output", country, "endPFig.rdata"))
af.end <- qtls_plot + ylim(0, 170)
load(here("output", country, "attackRFig.rdata"))
af.ar <- attackR
load(here("output", country, "cumCasesFig.rdata"))
af.c <- cumcases +
  xlim(0, 35)
load(here("output", country, "heat_cases.rdata"))
af.heatc <- heat_cases
load(here("output", country, "heat_inf.rdata"))
af.heati <- heat_inf
load(here("output", country, "pmonth_Fig.rdata"))
af.pmonth <- qtls_plot


country <- "SA"
load(here("output", country, "endPFig.rdata"))
sa.end <- qtls_plot + ylim(0, 170)
load(here("output", country, "attackRFig.rdata"))
sa.ar <- attackR
load(here("output", country, "cumCasesFig.rdata"))
sa.c <- cumcases +
  xlim(0, 35)
load(here("output", country, "heat_cases.rdata"))
sa.heatc <- heat_cases
load(here("output", country, "heat_inf.rdata"))
sa.heati <- heat_inf
load(here("output", country, "pmonth_Fig.rdata"))
sa.pmonth <- qtls_plot

country <- "TKY"
load(here("output", country, "endPFig.rdata"))
tk.end <- qtls_plot + ylim(0, 170)
load(here("output", country, "attackRFig.rdata"))
tk.ar <- attackR
load(here("output", country, "cumCasesFig.rdata"))
tk.c <- cumcases +
  xlim(0, 35)
load(here("output", country, "heat_cases.rdata"))
tk.heatc <- heat_cases
load(here("output", country, "heat_inf.rdata"))
tk.heati <- heat_inf
load(here("output", country, "pmonth_Fig.rdata"))
tk.pmonth <- qtls_plot


###################
## Figure 2: Incidence trajectories and fits to sleected years
###################

country <- "AFG"
load(here("output", country, "CasesFit.rdata"))
af <- hinc_mo2 +
  labs(y = "Monthly\n incidence")

shared_cases <- extract_legend(af)

country <- "SA"
load(here("output", country, "CasesFit.rdata"))
sa <- hinc_mo2 +
  ylim(0, 25) +
  labs(y = " ")

country <- "TKY"
load(here("output", country, "CasesFit.rdata"))
tk <- hinc_mo2 +
  ylim(0, 200) +
  labs(y = " ")

country <- "AFG"
load(here("output", country, "stack_incidence.rdata"))
af.stack <- stack +
  labs(y = "Monthly\n incidence")

shared_stack <- extract_legend(af.stack)

country <- "SA"
load(here("output", country, "stack_incidence.rdata"))
sa.stack <- stack +
  labs(y = " ")

country <- "TKY"
load(here("output", country, "stack_incidence.rdata"))
tk.stack <- stack +
  labs(y = " ")

# Plot
windows()
gridExtra::grid.arrange(gridExtra::arrangeGrob(af.stack + theme(legend.position = "none"),
  tk.stack + theme(legend.position = "none"),
  sa.stack + theme(legend.position = "none"),
  shared_stack,
  af + theme(legend.position = "none"),
  tk + theme(legend.position = "none"),
  sa + theme(legend.position = "none"),
  shared_cases,
  ncol = 4
))


# Save figure 2
ggsave(path = here("output"), filename = "Fig2.tiff", device = "tiff")


###################
## Figure 3 : patient months
###################

# Apply user-defined function to extract legend
shared_legend <- extract_legend(qtls_plot)

windows()
gridExtra::grid.arrange(gridExtra::arrangeGrob(af.pmonth + theme(legend.position = "none") +
  labs(title = "Afghanistan", x = "", y = "Person-months"),
tk.pmonth + theme(legend.position = "none") +
  labs(title = "Turkey", x = "Cases", y = " "),
sa.pmonth + theme(legend.position = "none") +
  labs(title = "South Africa", x = "", y = " "),
ncol = 3, nrow = 1
), shared_legend, heights = c(10, 1))

ggsave(path = here("output"), filename = "Fig3.tiff", device = "tiff")


###############
## Figure 4 : Heat countour plots
##############

shared_c <- extract_legend(tk.heatc)
shared_i <- extract_legend(tk.heati)


windows()
gridExtra::grid.arrange(gridExtra::arrangeGrob(af.heatc + theme(legend.position = "none"),
  tk.heatc + theme(legend.position = "none"),
  sa.heatc + theme(legend.position = "none"),
  shared_c,
  af.heati + labs(title = "") +
    theme(legend.position = "none"),
  tk.heati + labs(title = "") +
    theme(legend.position = "none"),
  sa.heati + labs(title = "") +
    theme(legend.position = "none"),
  shared_i,
  ncol = 4
))


ggsave(path = here("output"), filename = "Fig4.tiff", device = "tiff")

########
## End of code
#########