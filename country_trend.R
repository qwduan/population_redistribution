
#### The script includes functions for statistical analysis to identify trends and seasonal patterns in urban user proportions using time-series decomposition and breakpoint analysis.
#### The 'lulc_night_update.csv' is Facebook daily data during nighttime, which has undergone data imputation and normalisation, and is stored in files named by each country.
#### The functions 'DTECTslope_ro' and 'DTECTslope_0_ro' for detecting robust linear regressions should be run first.

library(dplyr); library(ggplot2); library(readr);library(tidyr)
library(lubridate);library(stringr);library(ISOweek);library(sfsmisc)
library(strucchange);library(bfast);library(stlplus)
library(MASS);library(car);library(relaimpo)

# Set the directory path where the FB data for each country is stored
directory_path <- "*/imputation_urban"
# List all country directories in the root directory
country_dirs <- list.dirs(directory_path, full.names = TRUE, recursive = FALSE)

# Function to process each country's data
country_data <- function(country_dir) {
  country_name <- basename(country_dir)
  print(country_name)
  # Construct the path to the night.csv file within the country's directory
  night_data <- file.path(country_dir, "lulc_night_update.csv")
  
  # Check if the night.csv file exists
  if(!file.exists(night_data)) {
    warning(paste("night.csv not found in", country_dir))
    return(NULL)
  } else {
    pop <- read_csv(night_data)
    pop$est[pop$est < 0] <- 0
    
    # Define workdays depending on the country
    if (country_name %in% c('Algeria', 'Egypt', 'Jordan', 'Qatar', 'Saudi Arabia')) {
      workday <- c('Sun', 'Mon', 'Tue', 'Wed', 'Thu')
    } else if (country_name == 'Nepal') {
      workday <- c('Sun', 'Mon', 'Tue', 'Wed', 'Thu', 'Fri')
    } else if (country_name == 'United Arab Emirates') {
      workday <- c('Mon', 'Tue', 'Wed', 'Thu')
    } else {
      workday <- c('Mon', 'Tue', 'Wed', 'Thu', 'Fri')  # Default case for other countries
    }
    
    pop_week <- pop %>% filter(week %in% workday)
    
    # Re-categorise urban areas into new categories for analysis
    # 1 is urban centre, 2 is suburban, 3 is rural areas
    # 21 is urban areas, 23 is rural areas
    pop_week <- pop_week %>%
      mutate(urban_category2 = case_when(
        urban_category %in% c(1, 2) ~ 21,
        urban_category == 3 ~ 23
      ))
    
    proportions <- pop_week %>%
      dplyr::group_by(Date, urban_category2) %>%
      dplyr::summarise(category_sum = sum(est)) %>%
      left_join(pop_week %>% group_by(Date) %>%
                  dplyr::summarise(total_num = sum(est)), by = "Date") %>%
      mutate(proportion = category_sum / total_num * 100)
    
    # Adjust the year 2020 to 52 weeks for weekly matching.
    proportions <- proportions %>%
      mutate(
        Year = isoyear(Date),
        Week = isoweek(Date),
        adj_Week = if_else(Year == 2020, Week - 1, Week),
        weekdate = Date - wday(Date, week_start = 1) + 1   # Assuming the week starts on Monday
      )
    proportions$weekdate <- as.Date(proportions$weekdate)

    weekly_pro <- proportions %>%
      dplyr::group_by(urban_category2, weekdate, Year, Week, adj_Week) %>%
      dplyr::summarise(ave = mean(proportion , na.rm = TRUE),
                       std_error = sd(proportion , na.rm = TRUE) / sqrt(n()),
                       .groups = 'drop') 
    
    # Filter for urban areas
    weekly_pro_3 <- weekly_pro %>% filter(urban_category2==21)%>%
      arrange(weekdate)
    
    
    complete_adj_week <- function(data) {
      years <- unique(data$Year)
      list_data <- lapply(years, function(y) {
        year_data <- data[data$Year == y, ]
        start_week <- min(year_data$adj_Week)
        if (y == 2020 || y == 2021) {
          end_week <- 52  
        } else {
          end_week <- max(year_data$adj_Week) 
        } 
        full_weeks <- data.frame(adj_Week = start_week:end_week)
        full_data <- merge(full_weeks, year_data, by = "adj_Week", all.x = TRUE)
        full_data$Year <- y  
        return(full_data)
      })
      
      full_data <- do.call(rbind, list_data)
      return(full_data)
    }
    

    full_week_pro_3 <- complete_adj_week(weekly_pro_3)
    full_week_pro_3$Week = if_else(full_week_pro_3$Year == 2020, full_week_pro_3$adj_Week + 1, full_week_pro_3$adj_Week)
  
    final_df <- full_week_pro_3 %>%
      arrange(Year, adj_Week, Week) %>%
      group_by(Year) %>%
      mutate(
        weekdate = ifelse(is.na(weekdate), ISOweek2date(sprintf("%d-W%02d-%d", Year, Week, 1)), weekdate),
        urban_category2 = ifelse(is.na(urban_category2), 21, urban_category2),
        ave = ifelse(is.na(ave), NA, ave)
      ) %>%
      ungroup()
    
    # Create a time series object for further analysis
    final_df$weekdate = as.Date(final_df$weekdate)
    final_df$number <- as.numeric(as.factor((final_df$weekdate)))
    final_df$ts_data <- ts(final_df$ave, frequency = 52, start = c(2020, final_df$adj_Week[1]))
    
    # Perform STL decomposition to extract trend and seasonal components
    if (length(final_df$ts_data) <2 * frequency(final_df$ts_data)) {
      NULL
    } else {
      stlv <- stlplus(final_df$ts_data, s.window = "periodic")
      trend_vector <- as.vector(stlv[["data"]][["trend"]])
      season_vector <- as.vector(stlv[["data"]][["seasonal"]])
      time_vector <- as.vector(time(final_df$ts_data))
      
      plot(stlv)
      final_df$trend <- trend_vector
      final_df$season <- season_vector
      # Deseasonalised data
      final_df$deseason <- final_df$ts_data - final_df$season
      final_df$ts_time <- time_vector

      # Apply the bfast method to detect structural breaks in the time series
      bf_result <- bfast(final_df$ts_data , h = 0.4,  season = 'dummy', decomp='stlplus', max.iter = 1, level = 0.05)
      plot(bf_result)
      
      # Extract breakpoints and trend components from the bfast results
      bps <- bf_result[["output"]][[1]][["Vt.bp"]]
      Tt <- bf_result[["output"]][[1]][["Tt"]]
      final_df$Tt <- as.vector(Tt)
      
      # If breakpoints are detected, calcualte trend
      # should run DTECTslope_ro function firstly for regression analysis before and after the breakpoint
      if (bps > 0) {
        bp_date<- final_df$weekdate[final_df$number== bps]
        lower_ci <- bf_result[["output"]][[1]][["ci.Vt"]][["confint"]][1, "2.5 %"]         # Lower bound of the confidence interval
        lci_date <- final_df$weekdate[final_df$number == lower_ci]
        upper_ci <- bf_result[["output"]][[1]][["ci.Vt"]][["confint"]][1, "97.5 %"]  
        uci_date <- final_df$weekdate[final_df$number == upper_ci]

        final_df_filtered <- final_df[!is.na(final_df$ave) & !is.na(final_df$deseason), ]
        bptest <- DTECTslope_ro(final_df_filtered,bps)
        
        before_line_type <- ifelse(bptest$p1 < 0.05, "solid", "dashed")
        after_line_type <- ifelse(bptest$p2  < 0.05, "solid", "dashed")
        
        seg1_y <- bptest$s10 + bptest$s1 * final_df_filtered$ts_time[which.min(final_df_filtered$weekdate)]
        seg1_yend <- bptest$s10 + bptest$s1 * final_df_filtered$ts_time[final_df_filtered$weekdate == bp_date]
        
        seg2_y <- bptest$s20 + bptest$s2 * final_df_filtered$ts_time[final_df_filtered$weekdate == (bp_date + 7)]
        seg2_yend <- bptest$s20 + bptest$s2 * final_df_filtered$ts_time[which.max(final_df_filtered$weekdate)]
        
        plot1 <- ggplot(final_df_filtered , aes(x=weekdate))+
          geom_rect(aes(xmin = lci_date, xmax = uci_date, ymin = -Inf, ymax = Inf),
                    fill = "gray90", alpha = 0.5)+
          geom_ribbon(data = final_df_filtered[final_df_filtered$weekdate <= bp_date, ], aes(ymin = bptest$cil_1[, "lwr"], ymax = bptest$cil_1[, "upr"]), fill = "#FF00FF", alpha = 0.2) +
          geom_ribbon(data = final_df_filtered[final_df_filtered$weekdate > bp_date, ], aes(ymin = bptest$cil_2[, "lwr"], ymax = bptest$cil_2[, "upr"]), fill = "#FF00FF", alpha = 0.2) + 
          geom_line(aes(y = ave), color = "black", linewidth =0.5) +
          geom_line(aes(y = deseason), color = '#40ff00', linewidth=0.5) +
          geom_vline(xintercept = bp_date, linetype = "dashed", color = "gray45", linewidth = 0.5) +
          annotate("segment", x = min(final_df_filtered$weekdate), xend = bp_date, 
                   y = seg1_y, yend = seg1_yend,
                   color = "#FF00FF", linewidth = 0.5, linetype = before_line_type) +
          annotate("segment", x = bp_date + 7, xend = max(final_df_filtered$weekdate), 
                   y = seg2_y, yend = seg2_yend,
                   color = "#FF00FF", linewidth = 0.5, linetype = after_line_type) +
          theme_bw()+ 
          theme(panel.border = element_rect(linewidth = 0.25))+
          labs(x = "Year Month", y = "urban proportions")
        
        country_name <- basename(country_dir)
        #ggsave(paste0(country_dir, "/", country_name, "_bfast_plot_urban_rlm.pdf"), plot1, width = 4, height = 3)
        
        return(list(
          country_name = country_name,
          bp_num = bps,
          bp_date = bp_date,
          lower_ci = lci_date,
          upper_ci = uci_date,
          bp = bptest
        ))
      }
      else{
        # for those without breakpoint
        # should run the the DTECTslope_0_ro function firstly to calcualte the whole trend
        final_df_filtered <- final_df[!is.na(final_df$ave) & !is.na(final_df$deseason), ]
        bptest <- DTECTslope_0_ro(final_df_filtered)
        
        line_type <- ifelse(bptest$pall  < 0.05, "solid", "dashed")
        segall_y <- bptest$sa0 + bptest$sall * final_df_filtered$ts_time[which.min(final_df_filtered$weekdate)]
        segall_yend <- bptest$sa0 + bptest$sall * final_df_filtered$ts_time[which.max(final_df_filtered$weekdate)]
        
        plot1 <- ggplot(final_df_filtered, aes(x=weekdate))+
          geom_ribbon(aes(ymin = bptest$cil_all[, "lwr"], ymax = bptest$cil_all[, "upr"]), fill = "#FF00FF", alpha = 0.2) + 
          geom_line(aes(y = ave), color = "black", linewidth =0.5) +
          geom_line(aes(y = deseason), color = '#40ff00', linewidth=0.5) +
          annotate("segment", x = min(final_df_filtered$weekdate), xend = max(final_df_filtered$weekdate), 
                   y = segall_y, yend = segall_yend,
                   color = "#FF00FF", linewidth = 0.5, linetype = after_line_type) +
          theme_bw()+
          theme(panel.border = element_rect(linewidth = 0.25))+
          labs(x = "Year Month", y = "urban proportions")
        
        country_name <- basename(country_dir)
        #ggsave(paste0(country_dir, "/", country_name, "_bfast_plot_urban_rlm.pdf"), plot1, width = 4, height = 3)
        
        return(list(
          country_name = country_name,
          bp_num = NA,
          bp_date = NA,
          lower_ci = NA,
          upper_ci =NA,
          bp = bptest
        ))
      }
    }
  }
}

cc <- lapply(country_dirs, country_data)

cc_list <- vector("list", length(cc))

# Loop through each item in the results to extract and store specific details
for (i in seq_along(cc)) {
  item <- cc[[i]]
  if (!is.null(item$country_name)) { 
    cc_list[[i]] <- list(
      country = item$country_name,
      bp_num = item$bp_num,
      bp_date = item$bp_date,
      upper_ci = item$upper_ci,
      lower_ci = item$lower_ci,
      s1 = item$bp$s1,
      p1 = item$bp$p1,
      s10 = item$bp$s10,
      s2 = item$bp$s2,
      p2 = item$bp$p2,
      s20 = item$bp$s20,
      sall = item$bp$sall,
      sa0 = item$bp$sa0,
      pall = item$bp$pall,
      sd1 = item$bp$sd1,
      sd2 = item$bp$sd2,
      sdall = item$bp$sdall
    )
  }
}
cc_list <- cc_list[!sapply(cc_list, is.null)]
cc_df <- do.call(rbind, lapply(cc_list, data.frame, stringsAsFactors = FALSE))

write.csv(cc_df, '*/country_bp_urban_rlm.csv')  


######################################################

# name include country's social, economic and environmental characteristics
name <- read_csv('*/country_name.csv')
# Update bp dataframe handling no breakpoint countries, using the whole trend and siginficant to represent their trend in phase 1 and 2.
bp <- bp %>%
  mutate(s1 = ifelse(is.na(bp_num) & is.na(s1), sall, s1),
         s2 = ifelse(is.na(bp_num) & is.na(s2), sall, s2))
bp <- bp %>%
  mutate(p1 = ifelse(is.na(bp_num) & is.na(p1), pall, p1),
         p2 = ifelse(is.na(bp_num) & is.na(p2), pall, p2))

bp_link <- bp %>%
  left_join(name , by ='country')

lm1 <- lm(s1 ~ hdi + hdi_grow + gdp_grow + gdp_indu + tourism + 
            protect + pmdif1 + greendif, data = bp_link)
summary(lm1)
#check for multicollinearity
vif(lm1)
# Extract BIC and AIC
extractAIC(lm1, k=log(nrow(bp_link)))
extractAIC(lm1)

# Stepwise regression using BIC/AIC criteria
selected_model <- MASS::stepAIC(lm1, direction="both", k=log(nrow(bp_link)), trace=TRUE)
selected_mode2 <- MASS::stepAIC(lm1, direction="both", trace=TRUE) 

summary(selected_model)
summary(selected_mode2)

# Calculate the relative importance of variables in the selected model
rel_importance <- calc.relimp(selected_mode2, type = "lmg", rela = TRUE)
importance_values <- rel_importance@lmg # Convert to percentage
Variable = names(importance_values)
coe1 <- rel_importance@ave.coeffs[, "1X"]

# Plot the relative importance
importance_df <- data.frame(
  Variable = Variable,
  Importance = as.numeric(importance_values)*100,
  coe1 <- rel_importance @ave.coeffs[, "1X"],
  p <- as.numeric(summary(selected_mode2)$coefficients[, "Pr(>|t|)"][-1] )
)
importance_df$sig <- ifelse(importance_df$p < 0.001, "***",
                            ifelse(importance_df$p < 0.01, "**",
                                   ifelse(importance_df$p< 0.05, "*",
                                          ifelse(importance_df$p < 0.1, "-", " "))))
ggplot(importance_df, aes(x = reorder(Variable, -Importance), y = Importance, fill = coe1 < 0)) +
  geom_col(width=0.5) +
  scale_fill_manual(values = c("#f678a7", "#62aec5")) +
  geom_text(aes(label = sig, y = Importance + 0.05), vjust = 0) + 
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) + 
  labs(title = "Relative Importance of Variales",
       x = "Variable",
       y = "Importance (%)") +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), 
    axis.line = element_line(colour = "black")
  )

lm2 <- lm(s2 ~ hdi + hdi_grow + gdp_grow + gdp_indu + tourism + 
            protect + pmdif1 + greendif, data = bp_link)
summary(lm2)
vif(lm2)
extractAIC(lm2, k=log(nrow(bp_link)))
extractAIC(lm2)

selected_mode3 <- MASS::stepAIC(lm2, direction="both", k=log(nrow(bp_link)), trace=TRUE)
selected_mode4 <- MASS::stepAIC(lm2, direction="both", trace=TRUE) 

summary(selected_mode3)
summary(selected_mode4)

rel_importance2 <- calc.relimp(selected_mode4, type = "lmg", rela = TRUE)
importance_values2 <- rel_importance2@lmg # Convert to percentage
Variable2 = names(importance_values2)
coe2 <- rel_importance2@ave.coeffs[, "1X"]

importance_df2 <- data.frame(
  Variable = Variable2,
  Importance = as.numeric(importance_values2)*100,
  coe <- rel_importance2 @ave.coeffs[, "1X"],
  p <- as.numeric(summary(selected_mode4)$coefficients[, "Pr(>|t|)"][-1] )
)
importance_df2$sig <- ifelse(importance_df2$p < 0.001, "***",
                             ifelse(importance_df2$p < 0.01, "**",
                                    ifelse(importance_df2$p< 0.05, "*",
                                           ifelse(importance_df2$p < 0.1, "-", " "))))

ggplot(importance_df2, aes(x = reorder(Variable, -Importance), y = Importance, fill = coe1 < 0)) +
  geom_col(width=0.5) +
  scale_fill_manual(values = c("#f678a7", "#62aec5")) +
  geom_text(aes(label = sig, y = Importance + 0.05), vjust = 0) + 
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) + 
  labs(title = "Relative Importance of Variales",
       x = "Variable",
       y = "Importance (%)") +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), 
    axis.line = element_line(colour = "black")
  )
