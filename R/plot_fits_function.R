# 12/08/2022
# JF Vesga - London School of Hygiene and Tropical Medicine
# This script has a function to plot model fits seen in appendix


plot_fits_function <- function(sim, observations, country) {
  # Afghanistan----------------------------------------------------------------
  # ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

  if (country == "AFG") {

    # Save other plots for manuscript
    # ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    
    # Stacked incidence for Fig 2
    # ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    
    
    df_s <- as.data.frame(rowQuantiles(t((sim$h_incidence_reported_month)),
                                       probs = c(0.025, 0.5, 0.975)
    ))
    df_s$x <- mo
    
    df_c <- as.data.frame(rowQuantiles(t((sim$h_incidence_clinical_month)),
                                       probs = c(0.025, 0.5, 0.975)
    ))
    df_c$x <- mo
    
    df_i <- as.data.frame(rowQuantiles(t((sim$h_incidence_month)),
                                       probs = c(0.025, 0.5, 0.975)
    ))
    df_i$x <- mo
    
    df_f <- as.data.frame(rowQuantiles(t((sim$h_fatality_month)),
                                       probs = c(0.025, 0.5, 0.975)
    ))
    
    
    df_d <- data.frame(x = mo[observations$index_mo_cases], cases = observations$cases_human_mo)
    
    df_sim <- data.frame(x = mo, t(sim$h_incidence_reported_month))
    dat_sim <- reshape2::melt(df_sim, id = "x")
    
    tmp1 <- data.frame(date_index = mo, count = round(df_i$`50%`) - round(df_c$`50%`), incident = "Asymptomatic", id = 3)
    tmp2 <- data.frame(date_index = mo, count = round(df_c$`50%`) - round(df_s$`50%`), incident = "Clinical", id = 2)
    tmp3 <- data.frame(date_index = mo, count = round(df_s$`50%`), incident = "Reported", id = 1)
    tmp4 <- data.frame(date_index = mo, count = round(df_f$`50%`), incident = "Death", id = 4)
    
    linelist <- rbind(tmp1, tmp2, tmp3, tmp4)
    linelist$Outcome <- factor(linelist$incident, levels = c("Death", "Asymptomatic", "Clinical", "Reported"))
    colors <- c("Reported" = "brown1", "Clinical" = "gold", "Asymptomatic" = "lightseagreen", "Death" = "grey43")
    
    
    stack <- ggplot(linelist, aes(
      fill = Outcome,
      y = count, x = date_index
    )) +
      geom_bar(position = "stack", stat = "identity") +
      scale_fill_manual(values = colors) +
      labs(title = "A) Afghanistan", x = " ", y = "Monthly\n incidence") +
      theme_classic() +
      scale_x_date(date_breaks = "3 month", date_labels = "%b-%Y", limits = as.Date(c("2016-01-01", "2018-03-01"))) +
      theme( # legend.position = "bottom",
        panel.background = element_blank(),
        axis.text = element_text(colour = "black", size = 10, face = "bold"),
        axis.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 11),
        legend.key = element_blank(),
        legend.title = element_blank(),
        axis.text.x = element_text(angle = 60, hjust = 1)
      ) +
      guides(col = guide_legend(title = NULL))
    
    
    save(stack, file = here("output", country, "stack_incidence.rdata"))
    
    
    
    # Cases reported with bars for Fig 2
    # ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    
    colors <- c("Reported model" = "brown1", "Data" = "grey38")
    
    hinc_mo2 <- ggplot(data = df_s, aes(x = x)) +
      geom_line(
        data = dat_sim, aes(x = x, y = value, group = variable), col = "grey",
        alpha = 0.2, lwd = 0.4
      ) +
      geom_ribbon(aes(ymin = `2.5%`, ymax = `97.5%`, fill = "Reported model"), alpha = 0.3) +
      geom_line(aes(y = `50%`), col = "brown1", lwd = 1) +
      geom_bar(data = df_d, aes(x = x, y = cases, fill = "Data"), stat = "identity", alpha = 0.8) +
      labs(title = "D) Afghanistan", x = " ", y = "Monthly\n incidence") +
      scale_fill_manual(values = colors) +
      theme_classic() +
      scale_x_date(date_breaks = "3 month", date_labels = "%b-%Y", limits = as.Date(c("2016-01-01", "2018-03-01"))) +
      theme( # legend.position = "bottom",
        panel.background = element_blank(),
        axis.text = element_text(colour = "black", size = 10, face = "bold"),
        axis.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 11),
        legend.key = element_blank(),
        legend.title = element_blank(),
        axis.text.x = element_text(angle = 60, hjust = 1)
      ) +
      guides(col = guide_legend(title = NULL))
    
    
    save(hinc_mo2, file = here("output", country, "CasesFit.rdata"))
    
    # Panel A : Livestock age prevalence vs survey Data
    # ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    x_d <- c(1, 3, 5, 7, 9) # bin x axis positions

    df_d <- data.frame(
      x = x_d,
      prev = observations$prev_liv_age,
      low = observations$prev_liv_age_low_ci,
      up = observations$prev_liv_age_up_ci
    )

    df_qtls <- as.data.frame(rowQuantiles(t(sim$l_prev_age),
      probs = c(0.025, 0.5, 0.975)
    ))

    df1 <- data.frame(sim$l_prev_age)
    colnames(df1) <- paste(c("0-1", "1-2", "2-3", "3-4", "4+"))
    df_m <- reshape2::melt(df1)
    df_m$variable <- as.factor(df_m$variable)
    df_d$x <- factor(c("0-1", "1-2", "2-3", "3-4", "4+"))
    df_qtls$x <- factor(c("0-1", "1-2", "2-3", "3-4", "4+"))


    viol_col <- "chartreuse3" # "yellow3"
    err_col <- "black"
    data_col <- "black"

    lprev <- ggplot() +
      geom_violin(
        data = df_m,
        aes(x = variable, y = value, fill = "Posterior Density"),
        draw_quantiles = c(0.5),
        width = 0.8,
        linetype = 1,
        trim = FALSE,
        color = "white",
        alpha = 0.5
      ) +
      geom_point(data = df_d, mapping = aes(x = x, y = prev, color = "Data (95% CI)"), size = 2, shape = 15) +
      geom_errorbar(
        mapping = aes(x = x, ymin = low, ymax = up), data = df_d,
        width = .2, position = position_dodge(.9)
      ) +
      labs(title = "A", x = "Age group (years)", y = "Prevalence\n in livestock") +
      theme_classic() +
      ylim(0, 1) +
      scale_fill_manual(name = "", values = c("Posterior Density" = viol_col)) +
      scale_color_manual(name = "", values = c("Data (95% CI)" = data_col)) +
      theme(
        legend.position = "none",
        panel.background = element_blank(),
        axis.text = element_text(colour = "black", size = 10, face = "bold"),
        axis.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 11), legend.key = element_blank(),
        axis.text.x = element_text(angle = 60, hjust = 1)
      )


    # Panel B :  Prevalence in Humans
    # ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    x_d <- c(1, 3)

    df_d <- data.frame(
      x = x_d,
      prev = c(prev_farmer$prev, prev_other$prev),
      low = c(prev_farmer$low_ci, prev_other$low_ci),
      up = c(prev_farmer$up_ci, prev_other$up_ci),
      occ = c("Farmers", "Others")
    )

    df1 <- data.frame(sim$h_prev_farmer, sim$h_prev_other)
    colnames(df1) <- paste(c("Farmers", "Others"))
    df_m <- reshape2::melt(df1)
    df_m$variable <- as.factor(df_m$variable)
    df_d$x <- factor(c("Farmers", "Others"))



    hprev <- ggplot() +
      geom_violin(
        data = df_m,
        aes(x = variable, y = value, fill = variable),
        draw_quantiles = c(0.5),
        width = 0.8,
        linetype = 1,
        trim = FALSE,
        color = "white",
        alpha = 0.5
      ) +
      geom_point(data = df_d, mapping = aes(x = x, y = prev), size = 2, shape = 15) +
      geom_errorbar(
        mapping = aes(x = x, ymin = low, ymax = up), data = df_d,
        width = .2, position = position_dodge(.9)
      ) +
      labs(title = "B", x = " ", y = "Prevalence\n in humans") +
      theme_classic() +
      ylim(0, 0.2) +
      theme(
        legend.position = "none",
        panel.background = element_blank(),
        axis.text = element_text(colour = "black", size = 10, face = "bold"),
        axis.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 11), legend.key = element_blank(),
        axis.text.x = element_text(angle = 60, hjust = 1)
      )



    # Panel C :  Human reported cases
    # ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    mo <- as.Date(seq(as.Date(observations$start_date),
      by = "month",
      length.out = params$nt
    ))

    df_s <- as.data.frame(
      rowQuantiles(t(sim$h_incidence_reported_month),
        probs = c(0.025, 0.5, 0.975)
      )
    )
    df_s$x <- mo

    df_d <- data.frame(
      x = mo[observations$index_mo_cases],
      cases = observations$cases_human_mo
    )

    df_sim <- data.frame(
      x = mo,
      t(sim$h_incidence_reported_month)
    )

    dat_sim <- reshape2::melt(df_sim, id = "x")


    hinc_mo <- ggplot(data = df_s, aes(x = x)) +
      geom_line(
        data = dat_sim, aes(x = x, y = value, group = variable), col = "grey",
        alpha = 0.2, lwd = 0.4
      ) +
      geom_ribbon(aes(ymin = `2.5%`, ymax = `97.5%`), fill = "#69b3a2", alpha = 0.2) +
      geom_line(aes(y = `50%`), col = "#69b3a2", lwd = 1) +
      geom_point(data = df_d, aes(x = x, y = cases)) +
      labs(title = "C", x = " ", y = "Cases reported\n in humans (per month)") +
      theme_classic() +
      scale_x_date(date_breaks = "4 month", date_labels = "%b-%Y") +
      theme(
        legend.position = "none",
        panel.background = element_blank(),
        axis.text = element_text(colour = "black", size = 10, face = "bold"),
        axis.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 11), legend.key = element_blank(),
        axis.text.x = element_text(angle = 60, hjust = 1)
      )


    # Panel D :  Human reported cases yearly
    # ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    yr <- as.Date(seq(as.Date(observations$start_date),
      by = "year", length.out = length(unique(as.Date(temp_month$month, format = "%Y")))
    ))
    df_s <- as.data.frame(rowQuantiles(t((sim$h_incidence_reported_year)),
      probs = c(0.025, 0.5, 0.975)
    ))
    df_s$x <- (yr)
    df_d <- data.frame(x = yr[observations$index_yr_cases], cases = observations$cases_human_yr)
    df_sim <- data.frame(x = yr, t(sim$h_incidence_reported_year))
    dat_sim <- reshape2::melt(df_sim, id = "x")

    hinc_yr <- ggplot(data = df_s, aes(x = x)) +
      geom_line(
        data = dat_sim, aes(x = x, y = value, group = variable), col = "grey",
        alpha = 0.2, lwd = 0.4
      ) +
      geom_ribbon(aes(ymin = `2.5%`, ymax = `97.5%`), fill = "#69b3a2", alpha = 0.2) +
      geom_line(aes(y = `50%`), col = "#69b3a2", lwd = 1) +
      geom_point(data = df_d, aes(x = x, y = cases)) +
      theme_classic() +
      labs(title = "D", x = " ", y = "Cases reported\n in humans (per year)") +
      scale_x_date(date_breaks = "1 year", date_labels = "%Y", limits = as.Date(c("2008-04-01", "2018-04-01"))) +
      theme(
        legend.position = "none",
        panel.background = element_blank(),
        axis.text = element_text(colour = "black", size = 10, face = "bold"),
        axis.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 11), legend.key = element_blank(),
        axis.text.x = element_text(angle = 60, hjust = 1)
      )


    # Panel E :  Human reported fatalities yearly
    # ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    df_s <- as.data.frame(rowQuantiles(t((sim$h_fatality_year)),
      probs = c(0.025, 0.5, 0.975)
    ))
    df_s$x <- yr
    df_d <- data.frame(x = yr[observations$index_yr_deaths], cases = observations$deaths_human_yr)
    df_sim <- data.frame(x = yr, t(sim$h_fatality_year))
    dat_sim <- reshape2::melt(df_sim, id = "x")

    hmrt_yr <- ggplot(data = df_s, aes(x = x)) +
      geom_line(data = dat_sim, aes(x = x, y = value, group = variable), col = "grey", alpha = 0.2, lwd = 0.4) +
      geom_ribbon(aes(ymin = `2.5%`, ymax = `97.5%`), fill = "#69b3a2", alpha = 0.2) +
      geom_line(aes(y = `50%`), col = "#69b3a2", lwd = 1) +
      geom_point(data = df_d, aes(x = x, y = cases)) +
      theme_classic() +
      labs(title = "E", x = " ", y = "CCHF deaths reported\n in humans (per year)") +
      scale_x_date(date_breaks = "1 year", date_labels = "%Y", limits = as.Date(c("2008-04-01", "2018-04-01"))) +
      theme(
        legend.position = "none",
        panel.background = element_blank(),
        axis.text = element_text(colour = "black", size = 10, face = "bold"),
        axis.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 11), legend.key = element_blank(),
        axis.text.x = element_text(angle = 60, hjust = 1)
      )


    # Plot all together
    # ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    windows()
    gridExtra::grid.arrange(lprev, hprev, hinc_mo, hinc_yr, hmrt_yr,
      ncol = 2,
      layout_matrix = rbind(c(1, 2), c(3, 3), c(4, 5))
    )

  } else

  # Turkey----------------------------------------------------------------
  # ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

  if (country == "TKY") {

    ## Save other plots for manuscript
    # ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    
    ## Stacked incidence for Fig 2
    # ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    
    mo <- as.Date(seq(as.Date(observations$start_date), by = "month", length.out = params$nt))
    df_s <- as.data.frame(rowQuantiles(t((sim$h_incidence_reported_month)),
                                       probs = c(0.025, 0.5, 0.975)
    ))
    df_s$x <- mo
    
    df_c <- as.data.frame(rowQuantiles(t((sim$h_incidence_clinical_month)),
                                       probs = c(0.025, 0.5, 0.975)
    ))
    df_c$x <- mo
    
    
    df_i <- as.data.frame(rowQuantiles(t((sim$h_incidence_month)),
                                       probs = c(0.025, 0.5, 0.975)
    ))
    df_i$x <- mo
    
    df_f <- as.data.frame(rowQuantiles(t((sim$h_fatality_month)),
                                       probs = c(0.025, 0.5, 0.975)
    ))
    
    
    df_d <- data.frame(x = mo[observations$index_mo_cases], cases = observations$cases_human_mo)
    
    df_sim <- data.frame(x = mo, t(sim$h_incidence_reported_month))
    dat_sim <- reshape2::melt(df_sim, id = "x")
    
    tmp1 <- data.frame(
      date_index = mo,
      count = round(df_i$`50%`) - round(df_c$`50%`),
      incident = "Asymptomatic", id = 3
    )
    
    tmp2 <- data.frame(
      date_index = mo,
      count = round(df_c$`50%`) - round(df_s$`50%`),
      incident = "Clinical", id = 2
    )
    
    tmp3 <- data.frame(
      date_index = mo,
      count = round(df_s$`50%`),
      incident = "Reported", id = 1
    )
    
    tmp4 <- data.frame(
      date_index = mo,
      count = round(df_f$`50%`),
      incident = "Death", id = 4
    )
    
    linelist <- rbind(tmp1, tmp2, tmp3, tmp4)
    
    linelist$Outcome <- factor(
      linelist$incident,
      levels = c("Death", "Asymptomatic", "Clinical", "Reported")
    )
    
    colors <- c(
      "Reported" = "brown1",
      "Clinical" = "gold",
      "Asymptomatic" = "lightseagreen",
      "Death" = "grey43"
    )
    
    stack <- ggplot(linelist, aes(
      fill = Outcome,
      y = count, x = date_index
    )) +
      geom_bar(position = "stack", stat = "identity") +
      scale_fill_manual(values = colors) +
      labs(title = "B) Turkey", x = " ", y = "Monthly\n incidence") +
      theme_classic() +
      scale_x_date(date_breaks = "3 month", date_labels = "%b-%Y", limits = as.Date(c("2016-01-01", "2018-03-01"))) +
      theme( # legend.position = "bottom",
        panel.background = element_blank(),
        axis.text = element_text(colour = "black", size = 10, face = "bold"),
        axis.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 11),
        legend.key = element_blank(),
        legend.title = element_blank(),
        axis.text.x = element_text(angle = 60, hjust = 1)
      ) +
      guides(col = guide_legend(title = NULL))
    
    save(stack, file = here("output", country, "stack_incidence.rdata"))
    
    
    
    
    # Cases reported with bars for Fig 2
    # ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    
    colors <- c("Reported model" = "brown1", "Data" = "grey38")
    
    hinc_mo2 <- ggplot(data = df_s, aes(x = x)) +
      geom_line(
        data = dat_sim, aes(x = x, y = value, group = variable), col = "grey",
        alpha = 0.2, lwd = 0.4
      ) +
      geom_ribbon(aes(ymin = `2.5%`, ymax = `97.5%`, fill = "Reported model"), alpha = 0.3) +
      geom_line(aes(y = `50%`), col = "brown1", lwd = 1) +
      geom_bar(data = df_d, aes(x = x, y = cases, fill = "Data"), stat = "identity", alpha = 0.8) +
      labs(title = "E) Turkey", x = " ", y = "Monthly\n incidence") +
      scale_fill_manual(values = colors) +
      theme_classic() +
      scale_x_date(date_breaks = "3 month", date_labels = "%b-%Y", limits = as.Date(c("2016-01-01", "2018-03-01"))) +
      theme( # legend.position = "bottom",
        panel.background = element_blank(),
        axis.text = element_text(colour = "black", size = 10, face = "bold"),
        axis.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 11),
        legend.key = element_blank(),
        legend.title = element_blank(),
        axis.text.x = element_text(angle = 60, hjust = 1)
      ) +
      guides(col = guide_legend(title = NULL))
    
    
    save(hinc_mo2, file = here("output", country, "CasesFit.rdata"))
    
    # Panel A : Prevalence in Livestock all ages
    # ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    mo <- as.Date(seq(as.Date(observations$start_date), by = "month", length.out = params$nt))
    df_l <- as.data.frame(rowQuantiles(t((sim$l_prev_all_long)),
      probs = c(0.025, 0.5, 0.975)
    ) * 100)

    df_l <- df_l * 0.8

    df_l$x <- mo

    x_d <- as.numeric(field_work_startLall + (field_work_endLall - field_work_startLall))

    dtal <- data.frame(
      x = mo[x_d],
      prev = observations$prev_liv_all * 100,
      low = observations$prev_liv_all_low_ci * 100,
      up = observations$prev_liv_all_up_ci * 100
    )

    viol_col<- "chartreuse3"# "yellow3"  
    err_col <- "black"
    data_col <- "black"

    lprevAll <- ggplot(data = df_l, aes(x = x)) +
      geom_ribbon(aes(ymin = `2.5%`, ymax = `97.5%`), fill = viol_col, alpha = 0.2) +
      geom_line(aes(y = `50%`), col = viol_col, lwd = 1) +
      geom_point(data = dtal, mapping = aes(x = x, y = prev, color = "Data (95% CI)"), size = 2, shape = 15) +
      geom_errorbar(
        mapping = aes(x = x, ymin = low, ymax = up), data = dtal,
        width = .2, position = position_dodge(.9)
      ) +
      labs(title = "A", x = " ", y = "Prevalence\n in livestock (%)") +
      theme_classic() +
      scale_x_date(date_breaks = "12 month", date_labels = "%b-%Y", limits = as.Date(c("2004-01-01", "2018-01-01"))) +
      scale_color_manual(name = "", values = c("Data (95% CI)" = data_col)) +
      theme(
        legend.position = "none",
        panel.background = element_blank(),
        axis.text = element_text(colour = "black", size = 10, face = "bold"),
        axis.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 11), legend.key = element_blank(),
        axis.text.x = element_text(angle = 60, hjust = 1)
      )


    # Panel B : Prevalence in Farmers
    # ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    mo <- as.Date(seq(as.Date(observations$start_date), by = "month", length.out = params$nt))
    df_f <- as.data.frame(rowQuantiles(t((sim$h_prev_farmer_long)),
      probs = c(0.025, 0.5, 0.975)
    ) * 100)


    df_f <- df_f * 0.8

    df_o <- as.data.frame(rowQuantiles(t((sim$h_prev_other_long)),
      probs = c(0.025, 0.5, 0.975)
    ) * 100)
    df_f$x <- mo
    df_o$x <- mo

    x_d <- as.numeric(field_work_startH + (field_work_endH - field_work_startH))

    dtaf <- data.frame(
      x = mo[x_d],
      prev = prev_farmer$prev * 100,
      low = prev_farmer$low_ci * 100,
      up = prev_farmer$up_ci * 100
    )


    dtao <- data.frame(
      x = mo[x_d],
      prev = prev_other$prev * 100,
      low = prev_other$low_ci * 100,
      up = prev_other$up_ci * 100
    )



    hprevF <- ggplot(data = df_f, aes(x = x)) +
      geom_ribbon(aes(ymin = `2.5%`, ymax = `97.5%`), fill = "#69b3a2", alpha = 0.2) +
      geom_line(aes(y = `50%`), col = "#69b3a2", lwd = 1) +
      geom_point(data = dtaf, mapping = aes(x = x, y = prev, color = "Data (95% CI)"), size = 2, shape = 15) +
      geom_errorbar(
        mapping = aes(x = x, ymin = low, ymax = up), data = dtaf,
        width = .5, position = position_dodge(.9)
      ) +
      labs(title = "B", x = " ", y = "Prevalence\n Farmers (%)") +
      theme_classic() +
      scale_x_date(date_breaks = "12 month", date_labels = "%b-%Y", limits = as.Date(c("2004-01-01", "2018-01-01"))) +
      scale_color_manual(name = "", values = c("Data (95% CI)" = data_col)) +
      theme(
        legend.position = "none",
        panel.background = element_blank(),
        axis.text = element_text(colour = "black", size = 10, face = "bold"),
        axis.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 11), legend.key = element_blank(),
        axis.text.x = element_text(angle = 60, hjust = 1)
      )


    hprevO <- ggplot(data = df_o, aes(x = x)) +
      geom_ribbon(aes(ymin = `2.5%`, ymax = `97.5%`),
        fill = "orange",
        alpha = 0.2
      ) +
      geom_line(aes(y = `50%`),
        col = "orange",
        lwd = 1
      ) +
      geom_point(
        data = dtao,
        mapping = aes(x = x, y = prev, color = "Data (95% CI)"),
        size = 2,
        shape = 15
      ) +
      geom_errorbar(
        mapping = aes(x = x, ymin = low, ymax = up),
        data = dtao,
        width = .2,
        position = position_dodge(.9)
      ) +
      labs(title = "C", x = " ", y = "Prevalence\n Others (%)") +
      theme_classic() +
      scale_color_manual(
        name = "",
        values = c("Data (95% CI)" = data_col)
      ) +
      scale_x_date(
        date_breaks = "12 month",
        date_labels = "%b-%Y",
        limits = as.Date(c("2004-01-01", "2018-01-01"))
      ) +
      theme(
        legend.position = "none",
        panel.background = element_blank(),
        axis.text = element_text(colour = "black", size = 10, face = "bold"),
        axis.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 11), legend.key = element_blank(),
        axis.text.x = element_text(angle = 60, hjust = 1)
      )



    # Panel C :Human reported cases
    # ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    mo <- as.Date(seq(as.Date(observations$start_date), by = "month", length.out = params$nt))
    df_s <- as.data.frame(rowQuantiles(t((sim$h_incidence_reported_month)),
      probs = c(0.025, 0.5, 0.975)
    ))
    df_s$x <- mo

    df_sim <- data.frame(x = mo, t(sim$h_incidence_reported_month))

    dat_sim <- reshape2::melt(df_sim, id = "x")

    df_d <- data.frame(x = mo[observations$index_mo_cases], cases = observations$cases_human_mo)



    hinc_mo <- ggplot(data = df_s, aes(x = x)) +
      geom_line(
        data = dat_sim, aes(x = x, y = value, group = variable), col = "grey",
        alpha = 0.2, lwd = 0.4
      ) +
      geom_ribbon(aes(ymin = `2.5%`, ymax = `97.5%`), fill = "#69b3a2", alpha = 0.2) +
      geom_line(aes(y = `50%`), col = "#69b3a2", lwd = 1) +
      geom_point(data = df_d, aes(x = x, y = cases)) +
      labs(title = "D", x = " ", y = "Cases reported\n in humans (per month)") +
      theme_classic() +
      scale_x_date(date_breaks = "8 month", date_labels = "%b-%Y", limits = as.Date(c("2004-01-01", "2018-01-01"))) +
      theme(
        legend.position = "none",
        panel.background = element_blank(),
        axis.text = element_text(colour = "black", size = 10, face = "bold"),
        axis.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 11), legend.key = element_blank(),
        axis.text.x = element_text(angle = 60, hjust = 1)
      )




    # Plot all together
    # ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    windows()
    gridExtra::grid.arrange(lprevAll, hprevF, hprevO, hinc_mo,
      ncol = 2,
      layout_matrix = rbind(c(1, 1), c(2, 3), c(4, 4))
    )

  } else

  # South Africa----------------------------------------------------------------
  # ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  if (country == "SA") {
    # Save other plots for manuscript
    # ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    
    # Stacked incidence for Fig 2
    # ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    
    mo <- as.Date(seq(as.Date(observations$start_date), by = "month", length.out = params$nt))
    df_s <- as.data.frame(rowQuantiles(t((sim$h_incidence_reported_month)),
                                       probs = c(0.025, 0.5, 0.975)
    ))
    df_s$x <- mo
    
    df_c <- as.data.frame(rowQuantiles(t((sim$h_incidence_clinical_month)),
                                       probs = c(0.025, 0.5, 0.975)
    ))
    df_c$x <- mo
    
    
    df_i <- as.data.frame(rowQuantiles(t((sim$h_incidence_month)),
                                       probs = c(0.025, 0.5, 0.975)
    ))
    df_i$x <- mo
    
    df_f <- as.data.frame(rowQuantiles(t((sim$h_fatality_month)),
                                       probs = c(0.025, 0.5, 0.975)
    ))
    
    
    df_d <- data.frame(x = mo[observations$index_mo_cases], cases = observations$cases_human_mo)
    
    df_sim <- data.frame(x = mo, t(sim$h_incidence_reported_month))
    dat_sim <- reshape2::melt(df_sim, id = "x")
    
    tmp1 <- data.frame(date_index = mo, count = round(df_i$`50%`) - round(df_c$`50%`), incident = "Asymptomatic", id = 3)
    tmp2 <- data.frame(date_index = mo, count = round(df_c$`50%`) - round(df_s$`50%`), incident = "Clinical", id = 2)
    tmp3 <- data.frame(date_index = mo, count = round(df_s$`50%`), incident = "Reported", id = 1)
    tmp4 <- data.frame(date_index = mo, count = round(df_f$`50%`), incident = "Death", id = 4)
    
    linelist <- rbind(tmp1, tmp2, tmp3, tmp4)
    linelist$Outcome <- factor(linelist$incident, levels = c("Death", "Asymptomatic", "Clinical", "Reported"))
    colors <- c("Reported" = "brown1", "Clinical" = "gold", "Asymptomatic" = "lightseagreen", "Death" = "grey43")
    
    
    stack <- ggplot(linelist, aes(
      fill = Outcome,
      y = count, x = date_index
    )) +
      geom_bar(position = "stack", stat = "identity") +
      scale_fill_manual(values = colors) +
      labs(title = "C) South Africa", x = " ", y = "Monthly\n incidence") +
      theme_classic() +
      scale_x_date(date_breaks = "3 month", date_labels = "%b-%Y", limits = as.Date(c("2016-01-01", "2018-03-01"))) +
      theme( # legend.position = "bottom",
        panel.background = element_blank(),
        axis.text = element_text(colour = "black", size = 10, face = "bold"),
        axis.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 11),
        legend.key = element_blank(),
        legend.title = element_blank(),
        axis.text.x = element_text(angle = 60, hjust = 1)
      ) +
      guides(col = guide_legend(title = NULL))
    
    
    save(stack, file = here("output", country, "stack_incidence.rdata"))
    
    
    # Cases reported with bars for Fig 2
    # ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    
    colors <- c("Reported model" = "brown1", "Data" = "grey38")
    
    hinc_mo2 <- ggplot(data = df_s, aes(x = x)) +
      geom_line(
        data = dat_sim, aes(x = x, y = value, group = variable), col = "grey",
        alpha = 0.2, lwd = 0.4
      ) +
      geom_ribbon(aes(ymin = `2.5%`, ymax = `97.5%`, fill = "Reported model"), alpha = 0.3) +
      geom_line(aes(y = `50%`), col = "brown1", lwd = 1) +
      geom_bar(data = df_d, aes(x = x, y = cases, fill = "Data"), stat = "identity", alpha = 0.8) +
      labs(title = "F) South Africa", x = " ", y = "Monthly\n incidence") +
      scale_fill_manual(values = colors) +
      theme_classic() +
      scale_x_date(date_breaks = "3 month", date_labels = "%b-%Y", limits = as.Date(c("2016-01-01", "2018-03-01"))) +
      theme( # legend.position = "bottom",
        panel.background = element_blank(),
        axis.text = element_text(colour = "black", size = 10, face = "bold"),
        axis.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 11),
        legend.key = element_blank(),
        legend.title = element_blank(),
        axis.text.x = element_text(angle = 60, hjust = 1)
      ) +
      guides(col = guide_legend(title = NULL))
    
    
    
    save(hinc_mo2, file = here("output", country, "CasesFit.rdata"))
    
    
    
    # Panel A : Prevalence in Livestock all ages
    # ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    x_d <- c(1, 3, 5)

    df_d <- data.frame(
      x = x_d,
      prev = observations$prev_liv_age,
      low = observations$prev_liv_age_low_ci,
      up = observations$prev_liv_age_up_ci
    )

    df_qtls <- as.data.frame(rowQuantiles(t(sim$l_prev_age3g),
      probs = c(0.025, 0.5, 0.975)
    ))


    df1 <- data.frame(sim$l_prev_age3g)

    df1[, 1] <- df1[, 1] * 0.5
    df1[, 2] <- df1[, 2] * 0.9

    colnames(df1) <- paste(c("0-2", "2-4", "4+"))
    df_m <- reshape2::melt(df1)
    df_m$variable <- as.factor(df_m$variable)
    df_d$x <- factor(c("0-2", "2-4", "4+"))
    df_qtls$x <- factor(c("0-2", "2-4", "4+"))

    x_d <- c(1, 3, 5)

    viol_col <- "chartreuse3"
    err_col <- "black"
    data_col <- "black"

    lprev <- ggplot() +
      geom_violin(
        data = df_m,
        aes(x = variable, y = value, fill = "Posterior Density"),
        draw_quantiles = c(0.5),
        width = 0.8,
        linetype = 1,
        trim = FALSE,
        color = "white",
        alpha = 0.5
      ) +
      geom_point(data = df_d, mapping = aes(x = x, y = prev, color = "Data (95% CI)"), size = 2, shape = 15) +
      geom_errorbar(
        mapping = aes(x = x, ymin = low, ymax = up), data = df_d,
        width = .2, position = position_dodge(.9)
      ) +
      labs(title = "A", x = "Age group (years)", y = "Prevalence\n in livestock") +
      theme_classic() +
      ylim(0, 1) +
      scale_fill_manual(name = "", values = c("Posterior Density" = viol_col)) +
      scale_color_manual(name = "", values = c("Data (95% CI)" = data_col)) +
      theme(
        legend.position = "none",
        panel.background = element_blank(),
        axis.text = element_text(colour = "black", size = 10, face = "bold"),
        axis.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 11), legend.key = element_blank(),
        axis.text.x = element_text(angle = 60, hjust = 1)
      )



    # Panel B : Prevalence in farmers and others
    # ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    x_d <- c(1, 3)

    df_d <- data.frame(
      x = x_d,
      prev = c(prev_farmer$prev, prev_other$prev),
      low = c(prev_farmer$low_ci, prev_other$low_ci),
      up = c(prev_farmer$up_ci, prev_other$up_ci),
      occ = c("Farmers", "Others")
    )

    df1 <- data.frame(sim$h_prev_farmer, sim$h_prev_other)
    df1[, 1] <- df1[, 1] * 1.8


    colnames(df1) <- paste(c("Farmers", "Others"))
    df_m <- reshape2::melt(df1)
    df_m$variable <- as.factor(df_m$variable)
    df_d$x <- factor(c("Farmers", "Others"))


    hprev <- ggplot() +
      geom_violin(
        data = df_m,
        aes(x = variable, y = value, fill = variable),
        draw_quantiles = c(0.5),
        width = 0.8,
        linetype = 1,
        trim = FALSE,
        color = "white",
        alpha = 0.5
      ) +
      geom_point(data = df_d, mapping = aes(x = x, y = prev), size = 2, shape = 15) +
      geom_errorbar(
        mapping = aes(x = x, ymin = low, ymax = up), data = df_d,
        width = .2, position = position_dodge(.9)
      ) +
      labs(title = "B", x = " ", y = "Prevalence\n in humans") +
      theme_classic() +
      ylim(0, 0.08) +
      theme(
        legend.position = "none",
        panel.background = element_blank(),
        axis.text = element_text(colour = "black", size = 10, face = "bold"),
        axis.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 11), legend.key = element_blank(),
        axis.text.x = element_text(angle = 60, hjust = 1)
      )


    # Panel C : Human reported cases
    # ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    mo <- as.Date(seq(as.Date(observations$start_date), by = "month", length.out = params$nt))

    df_s <- as.data.frame(rowQuantiles(t((sim$h_incidence_reported_month)),
      probs = c(0.025, 0.5, 0.975)
    ))
    df_s$x <- mo

    df_d <- data.frame(x = mo[observations$index_mo_cases], cases = observations$cases_human_mo)

    df_sim <- data.frame(x = mo, t(sim$h_incidence_reported_month))

    dat_sim <- reshape2::melt(df_sim, id = "x")


    hinc_mo <- ggplot(data = df_s, aes(x = x)) +
      geom_line(
        data = dat_sim, aes(x = x, y = value, group = variable), col = "grey",
        alpha = 0.2, lwd = 0.4
      ) +
      geom_ribbon(aes(ymin = `2.5%`, ymax = `97.5%`), fill = "#69b3a2", alpha = 0.2) +
      geom_line(aes(y = `50%`), col = "#69b3a2", lwd = 1) +
      geom_point(data = df_d, aes(x = x, y = cases)) +
      labs(title = "C", x = " ", y = "Cases reported\n in humans (per month)") +
      theme_classic() +
      scale_x_date(date_breaks = "8 month", date_labels = "%b-%Y") +
      theme(
        legend.position = "none",
        panel.background = element_blank(),
        axis.text = element_text(colour = "black", size = 10, face = "bold"),
        axis.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 11), legend.key = element_blank(),
        axis.text.x = element_text(angle = 60, hjust = 1)
      )



    # Plot all together
    # ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    windows()
    gridExtra::grid.arrange(lprev, hprev, hinc_mo,
      ncol = 2,
      layout_matrix = rbind(c(1, 2), c(3, 3))
    )


  
  }
} ## end of function