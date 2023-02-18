library(tidyverse)
library(tidydatras)  #devtools::install_github("fishvice/tidydatras")
library(icesDatras)

# if(fs::dir_exists("data-raw") & fs::file_exists("data-raw/cgfs.rds")) {
#   stop("You have what it takes")
# hh <- read_rds("data-raw/cgfs.rds")$hh
# hl <- read_rds("data-raw/cgfs.rds")$hl

# } else {
  fs::dir_create("data-raw")
  # datadownload -----------------------------------------------------------------

  years <- 1990:2022
  qs <- 4
  surveys <- "FR-CGFS"
  # res <- list()
  # hh <- hl <- tibble()
  # y = 1
  hh   <- 
    tidydatras::dr_getdata("HH", surveys, years, qs) %>% 
    tidydatras::dr_tidy() %>%                             # Make tidy with column specifications
    tidydatras::dr_idunite(., remove = FALSE)             # Add identifier
  hl   <- 
    tidydatras::dr_getdata("HL", surveys, years, qs) %>% 
    tidydatras::dr_tidy() %>%                             # Make tidy with column specifications
    tidydatras::dr_idunite(., remove = FALSE) %>%         # Add identifier
    tidydatras::dr_calccpue(hh) %>%                       # Calculated CPUE per hour
    left_join(tidydatras::aphia_latin) %>%                # Add scientific name
    left_join(tidydatras::asfis)                          # Add FAO species names and codes

  # results by year, length and station ------------------------------------------
  rbyls <-
    hl |>
    # filter(year >= 2000) |>
    rename(lon = shootlong,
           lat = shootlat)
  # get rid of species that have all zeros
  species <-
    rbyls |>
    group_by(species) |>
    summarise(n = sum(cpue_number_per_hour)) |>
    filter(n > 0) |>
    pull(species) |>
    sort()
  rbyls <-
    rbyls |>
    filter(species %in% species, cpue_number_per_hour > 0) |>
    # unite("id", survey, year, quarter, ship, gear, haulno, remove = FALSE) |>
    mutate(length = as.integer(floor(length ))) |>
    group_by(id, year, lon, lat, species, latin, english_name, length) |>
    summarise(n = sum(cpue_number_per_hour),
              .groups = "drop")

  # results by year and length ---------------------------------------------------
  # grid so that all lengths each year are filled. ensures correct estimates of
  #  the median by length and makes plotting easier
  g <- list()
  for(i in 1:length(species)) {
    length.range <-
      rbyls |>
      filter(species == species[i],
             latin   == latin[i],
             english_name == english_name[i]) |>
      summarise(min = min(length, na.rm=TRUE),
                max = max(length, na.rm=TRUE))
    g[[i]] <-
      expand_grid(year = as.character(years),
                  species = species[[i]],
                  length = length.range$min:length.range$max) |>
      left_join(dplyr::distinct(rbyls,
                                species, latin, english_name),
                by="species")
  }
  g <- bind_rows(g)
  
  rbyls <-
    rbyls |>
    right_join(g) |>
    mutate(n = replace_na(n, 0),
           b = n * 0.00001 * length^3) # kg)

  rbyl <-
    rbyls |>
    group_by(year, species, latin, english_name, length) |>
    # Mean catch over all the stations
    summarise(N = mean(n),
              B = mean(b),
              .groups = "drop") 

  rbys <-
    rbyls |>
    group_by(id, year, lon, lat, species, latin, english_name) |>
    summarise(N = sum(n),
              B = sum(b),
              .groups = "drop")

  my_boot = function(x, times = 100) {

    # Get column name from input object
    var = deparse(substitute(x))
    var = gsub("^\\.\\$","", var)

    # Bootstrap 95% CI
    cis = stats::quantile(replicate(times, mean(sample(x, replace=TRUE))), probs=c(0.025,0.975))

    # Return data frame of results
    data.frame(var, n=length(x), mean=mean(x), lower.ci=cis[1], upper.ci=cis[2])
  }

  print("Bootstrapping abundance:")

  boot.N <-
    rbys |>
    dplyr::group_by(species, latin, english_name, year) %>%
    dplyr::do(my_boot(.$N)) %>%
    dplyr::mutate(variable = "N",
                  var = as.character(variable))

  print("Bootstrapping biomass:")

  boot.B <-
    rbys %>%
    dplyr::group_by(species, latin, english_name, year) %>%
    dplyr::do(my_boot(.$B)) %>%
    dplyr::mutate(variable = "B",
                  var = as.character(var))

  boot <-
    bind_rows(boot.N,
              boot.B)


  # results by length, used in the length frequency plot ----------------------b--
  rbl <-
    rbyl |>
    group_by(species, latin, english_name, length) |>
    summarise(N = mean(N),
              B = mean(B),
              .groups = "drop")

  list(rbyl = rbyl, rbl = rbl, rbys = rbys, boot = boot, species = species, hh=hh, hl=hl) |>
    write_rds("data-raw/cgfs.rds")

# }
