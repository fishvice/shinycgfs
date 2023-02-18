scale_longitude_ices <- function(min = -44, max = 68.5, step = 1, ...) {
  breaks <- seq(min + 0.5, max - 0.5, step)
  labels <- geo::d2ir(60, breaks) %>% stringr::str_sub(3)
  return(ggplot2::scale_x_continuous(name = NULL, breaks = breaks, labels = labels, ...))
}
scale_latitude_ices <- function(min = 36, max = 84.5, step = 0.5, ...) {
  breaks <- seq(min + 0.25, max - 0.25, step)
  labels <- geo::d2ir(breaks, 0) %>% stringr::str_sub(1, 2)
  return(ggplot2::scale_y_continuous(name = NULL, breaks = breaks, labels = labels, ...))
}
gg_length <- function(median, by.year, SID, var, lab = "Fjöldi í hverju lengdarbili") {

  ggplot2::ggplot() +
    ggplot2::theme_bw(base_size = 14) +
    ggplot2::geom_ribbon(data = median %>%
                           dplyr::filter(sid == SID),
                         ggplot2::aes(length, ymax = {{ var }}, ymin = 0), fill = "grey") +
    ggplot2::geom_line(data = by.year  %>%
                         dplyr::filter(sid == SID),
                       ggplot2::aes(length, {{ var }})) +
    ggplot2::facet_wrap(~ year, dir = "v", ncol = 3, strip.position = "right") +
    ggplot2::labs(x = NULL, y = lab) +
    ggplot2::scale_x_continuous(breaks = seq(10, 200, by = 20)) +
    guides(x = guide_axis(n.dodge = 2))

}

gg_boot <- function(data, SID, ylab = "Fjöldi í togi") {
  data %>%
    dplyr::filter(sid == SID) %>%
    ggplot2::ggplot(ggplot2::aes(year, mean)) +
    ggplot2::theme_bw(base_size = 16) +
    ggplot2::geom_pointrange(ggplot2::aes(year, mean, ymin = lower.ci, ymax = upper.ci)) +
    ggplot2::scale_x_continuous(breaks = seq(1985, 2025, by = 5)) +
    ggplot2::expand_limits(y = 0) +
    ggplot2::labs(x = NULL, y = ylab)
}


#bb <- st_bbox(c(xmin = -29.6, xmax = -10, ymin = 62, ymax = 68.25), crs = st_crs(4326))
#fil <- "~/stasi/gis/GEBCO/data-raw/GEBCO_22_Dec_2022_ee31ae9610db/gebco_2022_n80.0_s30.0_w-70.0_e37.0.tif"


steps <- c(-Inf, -1200, -800, -400, -150, Inf)
p <-
  ggplot() +
  theme_void(base_size = 18) +
  #theme(legend.position = "none") +
  scale_x_continuous(expand = expansion(0), breaks = seq(-40, 40, by = 2)) +
  scale_y_continuous(expand = expansion(0), breaks = seq(50, 80, by = 1))

gg_bubble <- function(data, SID, var, lab = "Fjöldi") {
  data <-
    data %>%
    dplyr::filter(year %in% c(2000,
                              2005, 2010, 2014, 2015, 2016, 2017,
                              2018, 2019, 2020, 2021, 2022),
                  sid == SID)
  dummy <-
    data |>
    select(year) |>
    distinct()
    #ggplot2::ggplot() +
    #ggplot2::theme_void(base_size = 16) +
    #ggplot2::geom_sf(data = z, alpha = 0) +
    #ggplot2::geom_polygon(data = geo::island, ggplot2::aes(lon, lat), fill = "grey") +
  p +
    ggplot2::geom_point(data = data,
                        ggplot2::aes(lon, lat, size = {{ var }}),
                        alpha = 0.2, colour = "red") +
    ggplot2::geom_point(data = data,
                        ggplot2::aes(lon, lat),
                        size = 0.1, colour = "blue") +
    #geom_polygon(data = geo::island, aes(lon, lat), fill = "grey") +
    ggplot2::geom_text(data = dummy,
                        aes(label = year),
                        x = 10, y = 52,
                       size = 5,
                       colour = "black") +
    ggplot2::scale_size_area(max_size = 30) +
    ggplot2::labs(x = NULL, y = NULL, size = lab) +
    ggplot2::facet_wrap(~ year, nrow = 3, dir = "v") +
    theme(strip.text.x = element_blank()) +
    coord_quickmap()
}

gg_glyph <- function(data, SID, now.year, z) {

  # add a year before and after


  n.glyph <-
    data |>
    dplyr::filter(sid == SID) %>%
    GGally::glyphs(x_major = "lon",
                   y_major = "lat",
                   x_minor = "year",
                   y_minor = "Y",
                   width = 1,
                   height = 0.5)

  n.glyph %>%
    dplyr::mutate(years = ifelse(between(year, 2000, 2021), "history", "current"),
                  pos = ifelse(Y != 0, TRUE, FALSE),
                  base = lat - 0.25,
                  gy = ifelse(Y == 0 & between(year, 2000, 2022), gy + 0.005, gy)) %>%
    ggplot2::ggplot() +
    ggplot2::theme_bw() +
    geom_sf(data = z, alpha = 0, linewidth = 0.1) +
    geom_vline(xintercept = seq(-30, -8, by = 1), colour = "grey", linewidth = 0.1) +
    geom_hline(yintercept = seq(62, 68, by = 0.5), colour = "grey", linewidth = 0.1) +
    ggplot2::geom_linerange(ggplot2::aes(x = gx, ymin = base, ymax = gy,
                                         colour = years)) +
    #geom_point(data = d$rby |> filter(year == 2022, sid == 1), aes(lon, lat), colour = "black", size = 0.4) +
    ggplot2::geom_polygon(data = geo::island, ggplot2::aes(lon, lat), fill = "grey", alpha = 0.7) +
    ggplot2::scale_colour_manual(values = c("history" = "#377EB8", "current" = "#E41A1C")) +
    ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
                   panel.grid.minor = ggplot2::element_line(size = 1),
                   axis.ticks = ggplot2::element_blank(),
                   legend.position = "none") +
    labs(x = NULL, y = NULL)
}


