library(tidyverse)
if(fs::dir_exists("data-raw") & fs::file_exists("data-raw/cgfs.rds")) {
  stop("You have what it takes")
} else {
  fs::dir_create("data-raw")
  # datadownload -----------------------------------------------------------------

  library(icesDatras)
  years <- 2010:2021
  qs <- 4
  surveys <- "FR-CGFS"
  # res <- list()
  for(y in 1:length(years)) {
    print(years[y])
    h   <- icesDatras::getHHdata(surveys, year = years[y], quarter = qs)
    h   <- icesDatras::getHLdata(surveys, year = years[y], quarter = qs)
  
  }
  raw <- bind_rows(res)

  # head(raw)
  # icesDatras::checkSurveyYearOK("FR-CGFS", 2021)
  # icesDatras::getHHdata("FR-CGFS", 2021, 4)
  # icesDatras::getCPUELength("FR-CGFS", 2021, 4)
  
  # results by year, length and station ------------------------------------------
  rbyls <-
    raw |>
    rename_all(tolower) |>
    filter(year >= 2000) |>
    rename(lon = shootlon,
           lat = shootlat)
  # get rid of species that have all zeros
  latin <-
    rbyls |>
    group_by(species) |>
    summarise(n = sum(cpue_number_per_hour)) |>
    filter(n > 0) |>
    pull(species) |>
    sort()
  rbyls <-
    rbyls |>
    filter(species %in% latin) |>
    unite("id", survey, year, quarter, ship, gear, haulno, remove = FALSE) |>
    mutate(length = as.integer(floor(lngtclas / 10))) |>
    group_by(id, year, lon, lat, species, length) |>
    summarise(n = sum(cpue_number_per_hour),
              .groups = "drop")

  # results by year and length ---------------------------------------------------
  # grid so that all lengths each year are filled. ensures correct estimates of
  #  the median by length and makes plotting easier
  g <- list()
  for(i in 1:length(latin)) {
    length.range <-
      rbyls |>
      filter(species == latin[i]) |>
      summarise(min = min(length),
                max = max(length))
    g[[i]] <-
      expand_grid(year = 2000:2022,
                  species = latin[[i]],
                  length = length.range$min:length.range$max)
  }
  g <- bind_rows(g)
  rbyls <-
    rbyls |>
    right_join(g) |>
    mutate(n = replace_na(n, 0),
           b = n * 0.00001 * length^3) # kg)

  rbyl <-
    rbyls |>
    group_by(year, species, length) |>
    # Mean catch over all the stations
    summarise(N = mean(n),
              B = mean(b),
              .groups = "drop") |>
    select(species, year, length, N, B)

  rbys <-
    rbyls |>
    group_by(id, year, lon, lat, species) |>
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
    dplyr::group_by(species, year) %>%
    dplyr::do(my_boot(.$N)) %>%
    dplyr::mutate(variable = "N",
                  var = as.character(variable))

  print("Bootstrapping biomass:")

  boot.B <-
    rbys %>%
    dplyr::group_by(species, year) %>%
    dplyr::do(my_boot(.$B)) %>%
    dplyr::mutate(variable = "B",
                  var = as.character(var))

  boot <-
    bind_rows(boot.N,
              boot.B)


  # results by length, used in the length frequency plot ----------------------b--
  rbl <-
    rbyl |>
    group_by(species, length) |>
    summarise(N = mean(N),
              B = mean(B),
              .groups = "drop")

  list(rbyl = rbyl, rbl = rbl, rbys = rbys, boot = boot, species = latin) |>
    write_rds("data-raw/cgfs.rds")

}
