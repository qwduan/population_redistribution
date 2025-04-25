

library(dplyr); library(ggplot2); library(readr); library(tidyr); library(lubridate);
library(strucchange); library(bfast); library(segmented)
library(sp); library(sf); library(countrycode); library(r2country);
library(MASS); library(ISOweek); library(stlplus); library(sfsmisc)
library(spdep)
library(metaforest)
library(caret)
library(parallel)
library(doParallel)
# The functions 'DTECTslope_ro' and 'DTECTslope_0_ro' for detecting robust linear regressions should be run first.


######################################################
####Detect trend for each tile
####The breakpoint for each tile was derived from the breakpoint detected for this country.

directory_path <- "*"
country_dirs <- list.dirs(directory_path, full.names = TRUE, recursive = FALSE)

bp_data <- read_csv('*/country_bp_urban_rlm.csv')
record <- read_csv('*/country_name.csv')

# Define function to process data for each country directory
country_data <- function(country_dir) {
  country_name <- basename(country_dir)
  level <- record %>%
    filter(country == country_name) %>%
    pull(level)
  night_data <- file.path(country_dir, "lulc_night_update.csv")

  if(!file.exists(night_data )) {
    warning(paste("night.csv not found in", country_dir))
    return(NULL)
    } else {  
    pop <- read_csv(night_data)
    pop$est[pop$est < 0] <- 0
    
    # Define workdays based on country  
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
    pop_week <-  pop_week %>%
      mutate(  
        Year = isoyear(Date), 
        Week = isoweek(Date),   
        adj_Week = if_else(Year == 2020, Week - 1, Week),  
        weekdate = Date - lubridate::wday(Date, week_start = 1) + 1   # Assuming the week starts on Monday
      )
    pop_week$weekdate <- as.Date(pop_week$weekdate)
    pop_week_grouped <- pop_week %>%
      dplyr::group_by(lat, lon, urban_category, weekdate, urban, tree, crop, water, adj_Week, Week, Year) %>%
      dplyr::summarise(ave = mean(est, na.rm = TRUE)) %>%
      ungroup()
    
    tile_area <- calculate_tile_area(pop_week_grouped$lat, level)
    pop_week_grouped$area <- tile_area
    pop_week_grouped$den <-  pop_week_grouped$ave/pop_week_grouped$area
    pop_week_list <- pop_week_grouped %>%
        group_by(lat, lon, area, urban_category, urban, tree, crop, water) %>%
        group_split()
    pop_week_ts_list <- lapply(pop_week_list, function(df) {
        if (any(df$ave == 0)| length(df$den) <= 104) {
          return(NULL)
          }
      df$location_id = paste(df$lat[1], df$lon[1], sep = "_")
      df$lat = df$lat[1]
      df$lon = df$lon[1]
      df$area = df$area[1]
      df$urban_category = df$urban_category[1]
      df$urban = df$urban[1]
      df$crop = df$crop[1]
      df$tree = df$tree[1]
      df$water = df$water[1]
      df$weekdate = df$weekdate
      df$Week = df$Week
      df$adj_Week = df$adj_Week
      return(df)
      })
    pop_week_ts_list <- Filter(Negate(is.null), pop_week_ts_list)
    pop_week_ts_list <- lapply(pop_week_ts_list, function(list_item) {
      if (is.null(list_item)) return(NULL) 
      complete_adj_week <- function(data) {
        years <- unique(data$Year)
        list_data <- lapply(years, function(y) {
          year_data <- data[data$Year == y, ]
          start_week <- min(year_data$adj_Week)
          if (y == 2020 || y == 2021) { 
            end_week <- 52  # Fix the end week as 52 for the years 2020 and 2021
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
      completed_df = complete_adj_week(list_item) 
      completed_df$Week = if_else(completed_df$Year == 2020, completed_df$adj_Week + 1, completed_df$adj_Week)
      final_df <- completed_df  %>%
        arrange(Year, adj_Week, Week) %>%
        group_by(Year) %>%
        mutate(
          weekdate = ifelse(is.na(weekdate), ISOweek2date(sprintf("%d-W%02d-%d", Year, Week, 1)), weekdate),
          ave = ifelse(is.na(ave), NA, ave),
          den = ifelse(is.na(den), NA, den)
          ) %>%
        ungroup()
      final_df$weekdate = as.Date(final_df$weekdate)
      final_df$number <- as.numeric(as.factor((final_df$weekdate)))
      final_df$ts_data <- ts(final_df$den, frequency = 52, start = c(2020, final_df$adj_Week[1]))
      
      return(final_df)
        })
    pop_week_ts_list <- Filter(Negate(is.null), pop_week_ts_list)    
    stl_results <- list() 
    #Use breakpoint detected for this country before
    bps <- bp_data %>%
      filter(country == country_name) %>%
      pull(bp_num)
    check_missing_subseries <- function(ts_data) {
      subseries <- split(ts_data, cycle(ts_data))
      all_missing <- sapply(subseries, function(sub) all(is.na(sub)))
      return(any(all_missing))
      }
    for (i in seq_along(pop_week_ts_list)) {
      print(paste("Processing index:", i))  # Print the current index
      df <- pop_week_ts_list[[i]]
      if (!is.null(df) && length(df$ts_data) > 2 * frequency(df$ts_data)) {
        if (check_missing_subseries(df$ts_data)) {
          print(paste("Skipping index:", i, "due to all values missing in at least one subseries"))
          stl_results[[i]] <- NULL
          next
          } 
        df <- as.list(df)
        stlv <- stlplus(df$ts_data, s.window = "periodic")
        trend_vector <- as.vector(stlv[["data"]][["trend"]])
        season_vector <- as.vector(stlv[["data"]][["seasonal"]])
        time_vector <- as.vector(time(df$ts_data))
        df$trend <- trend_vector
        df$season <- season_vector
        df$deseason <- df$ts_data - df$season
        df$ts_time <- time_vector
        
        df_filtered <-  as.data.frame(df) %>% filter(!is.na(ts_data))
          if (!is.na(bps)) {
            bptest <- DTECTslope_ro(df_filtered, bps)
            } else {
              bptest <- DTECTslope_0_ro(df_filtered)
              }  
        df$bp <- bptest
        stl_results[[i]] <- df 
      } else {  
        stl_results[[i]] <- NULL  # Or handle short series differently
      }
    }
nonull_count <- sum(!sapply(stl_results, is.null))
print(nonull_count)

s1 <- numeric()
s2 <- numeric()
      sall <- numeric()
      p1 <- numeric()
      p2 <- numeric()
      pall <- numeric()
      lat <- numeric()
      lon <- numeric()
      urban <- numeric()
      tree <- numeric()
      crop <- numeric()
      water <- numeric()
      
      for (result in stl_results) {
        if (!is.null(result)) {
          lat <- c(lat, result$lat[1])
          lon <- c(lon, result$lon[1])
          urban <- c(urban, result$urban[1])
          water <- c(water, result$water[1])
          tree <- c(tree, result$tree[1])
          crop <- c(crop, result$crop[1])
          bp <- result$bp
        
          s1 <- c(s1, ifelse(!is.null(bp$s1), bp$s1, NA))
          s2 <- c(s2, ifelse(!is.null(bp$s2), bp$s2, NA))
          sall <-  c(sall, ifelse(!is.null(bp$sall), bp$sall, NA))
          p1 <- c(p1, ifelse(!is.null(bp$p1), bp$p1, NA))
          p2 <- c(p2, ifelse(!is.null(bp$p2), bp$p2, NA))
          pall <- c(pall, ifelse(!is.null(bp$pall), bp$pall, NA))
        }
        }
      result_df <- data.frame(
      lat = lat,
      lon = lon,
      urban= urban,
      tree =tree,
      crop = crop,
      water = water,
      s1 = s1,
      s2 = s2,
      sall = sall,
      p1 = p1,
      p2 = p2,
      pall = pall
      )
    
    write.csv(result_df, paste0(country_dir, "/", country_name, "_slopes_den.csv"))
    }
}


###################################################
####This function provides correlations between trend and built-up density across urban areas for each country in the list.

directory_path <- "*"

country_dirs <- list.dirs(directory_path, full.names = TRUE, recursive = FALSE)

country_corr <- function(country_dir_list) {
  result_list <- list()
  for (country_dir in country_dir_list) {
    # Construct the path to the country_name_slopes.csv file within the country's directory
    country_name <- basename(country_dir)
    grid_file <- file.path(country_dir, paste0(country_name, "_slopes_den.csv"))
    
    # Check if the grid_file exists
    if (!file.exists(grid_file)) {
      warning(paste(country_name, "_slopes.csv not found in", country_dir))
      return(NULL)
    } else {
      # Read the data and filter
      s <- read_csv(grid_file)
      
      s <- s %>% filter(urban>0.25)
      if (country_name %in% c('Poland', 'Turkey', 'Egypt', 'Italy', 'Nepal')) {
        s$s1 <- s$sall
        s$s2 <- s$sall
        s$p1 <- s$pall
        s$p2 <- s$pall
      }
     
      s$urban_rank <- rank(s$urban)
      s$s1_rank <- rank(s$s1)
      s$s2_rank <- rank(s$s2)
      
      # Compute Modified t-Test for S1 vs. Urban
      result_s1 <- modified.ttest(
        x = s$urban_rank, 
        y = s$s1_rank, 
        coords = s[, c("lon", "lat")]
      )
      
      # Compute Modified t-Test for S2 vs. Urban
      result_s2 <- modified.ttest(
        x = s$urban_rank, 
        y = s$s2_rank, 
        coords = s[, c("lon", "lat")]
      )
      
      # Extract modified correlation and effective sample size
      corr1 <- result_s1$corr  # Adjusted correlation for S1
      p1 <- result_s1$p.value # Adjusted p-value for S1
      corr2 <- result_s2$corr  # Adjusted correlation for S2
      p2 <- result_s2$p.value  # Adjusted p-value for S2
      
      n <- nrow(s)
      
      df <- data.frame(country = country_name, n = n, corr1 = corr1, corr2 = corr2, 
                       p1 = p1, p2 = p2)
      
      result_list[[country_name]] <- df
      return(df)
    }
  }
  result_df <- do.call(rbind, result_list)
  
  # Return the combined dataframe
  return(result_df)
}

all_cor <- lapply(country_dirs, country_corr)
all_cor2 <- do.call(rbind, all_cor)
write.csv(all_cor2, '*/urban_grid_corre_spearman.csv')
