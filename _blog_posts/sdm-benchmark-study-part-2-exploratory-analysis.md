---
layout: post
title: "SDM Benchmark Study Part 2: Exploratory Analysis"
permalink: /blog_posts/sdm-benchmark-study-part-2-exploratory-analysis
---

## Introduction

The Species Distribution Model (SDM) Benchmark Study Part 2 presents data from eight bird species across Colorado, Oregon, North Carolina, and Vermont, using environmental rasters from Part 1. The analysis segments data distribution by state and species and evaluates explanatory rasters using methods like Principal Component Analysis and Factor Analysis for Mixed Data. Additionally, the study addresses pseudo-absence generation and conducts feature engineering/selection to prepare for subsequent modeling.

## Visualize Observations

The initial phase of the exploratory data analysis visualizes bird species observations across four U.S. states: Oregon, Vermont, Colorado, and North Carolina, covering the period from 2016 to 2019. These visual representations display the distribution of eight distinct bird species within each state while also highlighting major cities for geographical context. Such visualizations offer insights into spatial distribution patterns, potentially revealing habitat preferences, migration routes, or the influence of urban areas on these distributions. As the analysis progresses, understanding these distributions becomes pivotal, especially when addressing spatial autocorrelations and other features that influence the accuracy and validity of the Species Distribution Model (SDM) benchmarks.

```
data(us.cities)

# Get major cities for each sample region (state)
.states <- c("OR", "VT", "CO", "NC")
top.cities <- purrr::map_df(.states, function(s) {
  out <- us.cities %>% 
  filter(country.etc==s) %>%
  mutate(city = gsub(paste0(" ", s), "", name)) %>%
  arrange(-pop)
  if (s == "OR") {
    out <- out %>% 
      head() %>%
      filter(!(city %in% c("Gresham", "Hillsboro", "Corvallis",
                           "Beaverton", "Springfield")))
  } else if (s == "CO") {
    out <- out %>%
      head() %>%
      filter(!(city %in% c("Thornton", "Lakewood", "Aurora")))
  } else if (s == "NC") {
    out <- out %>%
      head() %>%
      filter(!(city %in% c("Greensboro", "Durham", "Fayetteville")))
  } else {
    out <- out %>% head()
  }
  out
})

# Load the map data
states <- map_data("state") %>% 
  filter(region %in% c("oregon", "north carolina", "colorado", "vermont"))

# Load your data
data.files <- list.files("data/final", full.names = T)

df <- purrr::map_df(data.files, readRDS) 

caps.after.ws <- function(string) {
  gsub("(?<=\\s)([a-z])", "\\U\\1", string, perl = T)
}

# Define a function to create a plot for each species
plot.for.species <- function(spec, st.abbr) {
  st <- case_when(st.abbr == "CO" ~ "colorado",
                  st.abbr == "NC" ~ "north carolina",
                  st.abbr == "VT" ~ "vermont",
                  st.abbr == "OR" ~ "oregon",
                  T ~ "")
  
  title <- caps.after.ws(paste(st.abbr, gsub("_", " ", spec), 
                             "Observations, 2016-2019"))
  
  p <- ggplot(data = states %>% filter(region == st)) +
    geom_polygon(aes(x = long, y = lat, group = group),
                 fill = "#989875", color = "black") +
    geom_point(data = df %>% filter(state == st.abbr & common.name == spec), 
               aes(x = lon, y = lat), 
               size=1, alpha=.5, fill = "red", shape=21) +
    geom_point(data = top.cities %>% filter(country.etc == st.abbr), 
               aes(x=long, y=lat),
               fill="gold", color="black", size=3.5, shape = 21) + 
    geom_text(data = top.cities %>% filter(country.etc == st.abbr), 
              aes(x=long, y=lat, label=city),
              color="white", hjust=case_when(st.abbr=="NC"~.2,
                                               st.abbr=="VT"~.65,
                                               T~.5),
              vjust=ifelse(st.abbr=="NC", -.65, 1.5),
              size=4) + 
    coord_map() +
    ggtitle(title) +
    theme_minimal() +
    theme(panel.background = element_blank(),
          axis.text = element_blank(),
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          panel.grid = element_blank())

  data.table(
    state=st.abbr,
    species=spec,
    plot=list(p)
  )
}

spec.state <- expand.grid(unique(df$common.name), unique(df$state)) %>%
  rename(spec=Var1, st.abbr=Var2) 

# Create a list of plots
plots <- purrr::map2_df(spec.state$spec, 
                        spec.state$st.abbr, 
                        ~plot.for.species(.x, .y))
```

```
# Plot Ruddy Duck plots
do.call(ggpubr::ggarrange, 
        c(plots[species == "Ruddy Duck"]$plot, 
          list(nrow=2, ncol=2)))
```

<img class="sdm-eda-img" src="{{ site.baseurl }}/assets/plots/sdm-eda-1.png">

```
# Plot Belted Kingfisher plots
do.call(ggpubr::ggarrange, 
        c(plots[species == "Belted Kingfisher"]$plot, 
          list(nrow=2, ncol=2)))
```

<img class="sdm-eda-img" src="{{ site.baseurl }}/assets/plots/sdm-eda-2.png">

```
# Plot Wild Turkey plots
do.call(ggpubr::ggarrange, 
        c(plots[species == "Wild Turkey"]$plot, 
          list(nrow=2, ncol=2)))
```

<img class="sdm-eda-img" src="{{ site.baseurl }}/assets/plots/sdm-eda-3.png">

```
# Plot Sharp-Shinned Hawk plots
do.call(ggpubr::ggarrange, 
        c(plots[species == "Sharp-shinned Hawk"]$plot, 
          list(nrow=2, ncol=2)))
```

<img class="sdm-eda-img" src="{{ site.baseurl }}/assets/plots/sdm-eda-4.png">

```
# Plot Downy Woodpecker Plots
do.call(ggpubr::ggarrange, 
        c(plots[species == "Downy Woodpecker"]$plot, 
          list(nrow=2, ncol=2)))
```

<img class="sdm-eda-img" src="{{ site.baseurl }}/assets/plots/sdm-eda-5.png">

```
# Plot Cedar Waxwing Plots
do.call(ggpubr::ggarrange, 
        c(plots[species == "Cedar Waxwing"]$plot, 
          list(nrow=2, ncol=2)))
```

<img class="sdm-eda-img" src="{{ site.baseurl }}/assets/plots/sdm-eda-6.png">

```
# Plot Sandhill Crane Plots
do.call(ggpubr::ggarrange, 
        c(plots[species == "Sandhill Crane"]$plot, 
          list(nrow=2, ncol=2)))
```

<img class="sdm-eda-img" src="{{ site.baseurl }}/assets/plots/sdm-eda-7.png">

```
# Plot Sanderling Plots
do.call(ggpubr::ggarrange, 
        c(plots[species == "Sanderling"]$plot, 
          list(nrow=2, ncol=2)))
```

<img class="sdm-eda-img" src="{{ site.baseurl }}/assets/plots/sdm-eda-8.png">

## Explore Explanatory Rasters

The environmental rasters from Colorado, North Carolina, Oregon, and Vermont used in this study are meant to provide insights into which factors might potentially influence bird species distributions. The rasters were prepared in Part 1 of this study, and saved into 4 state-wide multi-layer rasters, each with 1000x1000 meter resoultions.

```
# Load prepared explanatory rasters by state
states <- c("CO", "NC", "OR", "VT")
r.files <- paste0("data/final_rasters/", states, ".tif")
r.list <- purrr::map(r.files, rast)
names(r.list) <- states
```

### Land Cover Type Frequency

The table below shows the frequency distribution of various land cover types across the selected states. By evaluating this, the most common and rare habitats in each region can be ascertained. Moreover, understanding the prevalence of each land cover type can provide insights into potential hot-spots or habitats of interest for the bird species in question.

```
purrr::map_df(states, function(st) {
  freq.df <- terra::freq(r.list[[st]]$NLCD_Land) %>%
    mutate(state = st) %>%
    dplyr::select(state, value, count)
}) 
```

<div data-pagedtable="false" pagedtable-page="0" class="pagedtable-wrapper">
<script data-pagedtable-source="" type="application/json">
{"columns":[{"label":["state"],"name":[1],"type":["chr"],"align":["left"]},{"label":["value"],"name":[2],"type":["chr"],"align":["left"]},{"label":["count"],"name":[3],"type":["dbl"],"align":["right"]}],"data":[{"1":"CO","2":"Open Water","3":"1055"},{"1":"CO","2":"Perennial Snow/Ice","3":"229"},{"1":"CO","2":"Developed, Open Space","3":"3854"},{"1":"CO","2":"Developed, Low Intensity","3":"2080"},{"1":"CO","2":"Developed, Medium Intensity","3":"1403"},{"1":"CO","2":"Developed, High Intensity","3":"484"},{"1":"CO","2":"Barren Land","3":"1998"},{"1":"CO","2":"Deciduous Forest","3":"17525"},{"1":"CO","2":"Evergreen Forest","3":"54313"},{"1":"CO","2":"Mixed Forest","3":"2158"},{"1":"CO","2":"Shrub/Scrub","3":"64428"},{"1":"CO","2":"Herbaceous","3":"78109"},{"1":"CO","2":"Hay/Pasture","3":"3192"},{"1":"CO","2":"Cultivated Crops","3":"33453"},{"1":"CO","2":"Woody Wetlands","3":"2389"},{"1":"CO","2":"Emergent Herbaceous Wetlands","3":"2946"},{"1":"NC","2":"Open Water","3":"2140"},{"1":"NC","2":"Developed, Open Space","3":"9159"},{"1":"NC","2":"Developed, Low Intensity","3":"4655"},{"1":"NC","2":"Developed, Medium Intensity","3":"1963"},{"1":"NC","2":"Developed, High Intensity","3":"693"},{"1":"NC","2":"Barren Land","3":"283"},{"1":"NC","2":"Deciduous Forest","3":"23858"},{"1":"NC","2":"Evergreen Forest","3":"17709"},{"1":"NC","2":"Mixed Forest","3":"13270"},{"1":"NC","2":"Shrub/Scrub","3":"3577"},{"1":"NC","2":"Herbaceous","3":"2991"},{"1":"NC","2":"Hay/Pasture","3":"9863"},{"1":"NC","2":"Cultivated Crops","3":"17649"},{"1":"NC","2":"Woody Wetlands","3":"18251"},{"1":"NC","2":"Emergent Herbaceous Wetlands","3":"1889"},{"1":"OR","2":"Open Water","3":"2319"},{"1":"OR","2":"Perennial Snow/Ice","3":"36"},{"1":"OR","2":"Developed, Open Space","3":"3864"},{"1":"OR","2":"Developed, Low Intensity","3":"1795"},{"1":"OR","2":"Developed, Medium Intensity","3":"1195"},{"1":"OR","2":"Developed, High Intensity","3":"403"},{"1":"OR","2":"Barren Land","3":"1362"},{"1":"OR","2":"Deciduous Forest","3":"619"},{"1":"OR","2":"Evergreen Forest","3":"85874"},{"1":"OR","2":"Mixed Forest","3":"6006"},{"1":"OR","2":"Shrub/Scrub","3":"82498"},{"1":"OR","2":"Herbaceous","3":"40217"},{"1":"OR","2":"Hay/Pasture","3":"8049"},{"1":"OR","2":"Cultivated Crops","3":"11965"},{"1":"OR","2":"Woody Wetlands","3":"1590"},{"1":"OR","2":"Emergent Herbaceous Wetlands","3":"3327"},{"1":"VT","2":"Unclassified","3":"3"},{"1":"VT","2":"Open Water","3":"1003"},{"1":"VT","2":"Developed, Open Space","3":"759"},{"1":"VT","2":"Developed, Low Intensity","3":"519"},{"1":"VT","2":"Developed, Medium Intensity","3":"267"},{"1":"VT","2":"Developed, High Intensity","3":"68"},{"1":"VT","2":"Barren Land","3":"40"},{"1":"VT","2":"Deciduous Forest","3":"9084"},{"1":"VT","2":"Evergreen Forest","3":"3109"},{"1":"VT","2":"Mixed Forest","3":"5201"},{"1":"VT","2":"Shrub/Scrub","3":"279"},{"1":"VT","2":"Herbaceous","3":"118"},{"1":"VT","2":"Hay/Pasture","3":"2725"},{"1":"VT","2":"Cultivated Crops","3":"436"},{"1":"VT","2":"Woody Wetlands","3":"1136"},{"1":"VT","2":"Emergent Herbaceous Wetlands","3":"189"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
</script>
</div>

### Continuous Variable Distribution

Continuous rasters capture a gradient of values, often related to environmental gradients such as temperature, precipitation, or elevation. In this study, the distribution of such continuous variables for each state is illustrated using density plots. These plots provide a visual representation of the distribution patterns of each environmental variable, allowing for an assessment of the variations and patterns within each state. For instance, understanding the elevation range of Oregon compared to Colorado can offer insights into the altitudinal preferences of species. Additionally, spotting any anomalies or peaks in these distributions might suggest specific ecological or environmental significance, further aiding in the understanding and interpretation of model outputs.

```
# Plot densities of numeric rasters, by state
purrr::map(states, function(st) {
  r <- r.list[[st]]
  plots <- purrr::map(names(r)[names(r) != "NLCD_Land"], function(l) {
    ggplot(data=r[l] %>% as.data.frame()) + 
      geom_density(aes(x=!!sym(l))) + 
      theme_bw() +
      ggtitle(paste0(st, " Distribution for ", gsub(paste0("_", st), "", l)))
  })
  
  ggpubr::ggarrange(plotlist=plots, ncol=3, nrow=4) 
}) %>%
  ggpubr::ggarrange(plotlist=., ncol=1, nrow=4) + 
  ggtitle(paste0("Density Plots of Numeric Variables"))
```

<img class="sdm-eda-img" src="{{ site.baseurl }}/assets/plots/sdm-eda-9.png">

### Principal Component Analysis

In this section, the study employs Principal Component Analysis (PCA) to transform the high-dimensional data into a lower-dimensional form while retaining as much of the original variance as possible. The initial step involves merging data from different states into a single dataframe, followed by cleaning factor levels of the NLCD_Land column. Subsequently, the code converts factor columns into dummy variables and removes columns with only a single unique value. The PCA() function is used to perform the PCA, and the results are visualized using bar charts, displaying the importance of each variable for the first two dimensions. This visualization is important in understanding the significance of different environmental rasters and how they contribute to the variations captured by the principal components.

```
r.df <- map_df(states, function(s) {
  df <- r.list[[s]] %>% as.data.frame()
  names(df) <- names(df) %>% gsub(paste0("_", s), "", .)
  df %>% setDT()
  df[, state := factor(s, levels=states)]
  df[apply(df, 1, function(.x) !any(is.na(.x)))]
}) 

# Custom function to process factor levels
clean.levels <- function(x) {
  # Remove non-alphanumeric characters and replace with underscores
  x <- gsub("[^a-zA-Z0-9]", "_", x)
  # Convert to lowercase
  x <- tolower(x)
  # Remove any leading or trailing underscores
  x <- gsub("^_|_$", "", x)
  x <- gsub("__", "_", x)
  x <- gsub("NLCD_Land_", "", x)
  return(x)
}

r.df[, NLCD_Land := factor(clean.levels(levels(NLCD_Land))[NLCD_Land])]

# Convert factor columns to dummy variables
df.dummies <- data.table(model.matrix(~ . - 1, 
                                      data = r.df[, .(NLCD_Land, state)])) %>%
  cbind(r.df[, -which(names(r.df) %in% c("NLCD_Land", "state")), with=F]) 

names(df.dummies) <- gsub("NLCD_Land", "", names(df.dummies))

# Ensure that there is more than one value per column (remove otherwise)
uniq.1 <- t( df.dummies[, lapply(.SD, uniqueN)]) %>%
  as.data.frame() %>%
  filter(V1 == 1) %>%
  row.names()

if (length(uniq.1) >= 1) {
  df.dummies <- df.dummies[, -which(names(df.dummies) %in% uniq.1), with=F]
}

pca.fit <- PCA(df.dummies, graph=F)
plot.PCA(pca.fit, choix="var")
```

<img class="sdm-eda-img" src="{{ site.baseurl }}/assets/plots/sdm-eda-10.png">

```
res <- pca.fit$var$coord %>%
  as.data.frame() %>%
  mutate(var=as.factor(rownames(.))) %>%
  dplyr::select(var, everything()) %>%
  as_tibble()
rownames(res) <- NULL
  
p.d1 <- ggplot(res, aes(x = reorder(var, Dim.1), y = Dim.1)) +
  geom_bar(stat = "identity", fill="darkblue") +
  coord_flip() +  # Makes it a horizontal bar chart
  labs(title = "Variable importance for Dim.1", y = "Importance", x = "Variable") +
  theme_minimal()

p.d2 <- ggplot(res, aes(x = reorder(var, Dim.2), y = Dim.2)) +
  geom_bar(stat = "identity", fill="darkred") +
  coord_flip() +  # Makes it a horizontal bar chart
  labs(title = "Variable importance for Dim.2", y = "Importance", x = "Variable") +
  theme_minimal()

ggpubr::ggarrange(plotlist=list(p.d1, p.d2), nrow=2, ncol=1)
```

<img class="sdm-eda-img" src="{{ site.baseurl }}/assets/plots/sdm-eda-11.png">

### Factor Analysis for Mixed Data

Diving deeper into the structure of the data, Factor Analysis for Mixed Data (FAMD) is implemented next. FAMD is a unique technique that extends traditional Factor Analysis by handling both continuous and categorical data, which fits the nature of the dataset used in this study. The code initiates this process by applying the FAMD() function to the dataframe. The results are then visualized in three distinct plots that highlight the variance explained by different variables for both continuous and categorical data. These plots provide insights into the underlying structure of the data and the significance of each variable. By observing the variable importance for the first two dimensions, displayed through bar charts, it can be discern which variables are related. Later, multi-collinearity and autocorrelation will be addressed in greater depth.

```
famd.fit <- FAMD(r.df, graph=F)

ggpubr::ggarrange(plotlist=purrr::map(
  c("var", "quanti", "quali"), 
  ~plot.FAMD(famd.fit, choix=.x)),
  ncol=1, nrow=3)
```

<img class="sdm-eda-img" src="{{ site.baseurl }}/assets/plots/sdm-eda-12.png">

```
res <- famd.fit$var$coord %>%
  as.data.frame() %>%
  mutate(var=as.factor(rownames(.))) %>%
  dplyr::select(var, everything()) %>%
  as_tibble()
rownames(res) <- NULL
  
p.d1 <- ggplot(res, aes(x = reorder(var, Dim.1), y = Dim.1)) +
  geom_bar(stat = "identity", fill="darkblue") +
  coord_flip() +  # Makes it a horizontal bar chart
  labs(title = "Variable importance for Dim.1", y = "Importance", x = "Variable") +
  theme_minimal()

p.d2 <- ggplot(res, aes(x = reorder(var, Dim.2), y = Dim.2)) +
  geom_bar(stat = "identity", fill="darkred") +
  coord_flip() +  # Makes it a horizontal bar chart
  labs(title = "Variable importance for Dim.2", y = "Importance", x = "Variable") +
  theme_minimal()

ggpubr::ggarrange(plotlist=list(p.d1, p.d2), nrow=2, ncol=1)
```

<img class="sdm-eda-img" src="{{ site.baseurl }}/assets/plots/sdm-eda-13.png">

## Pseudo-Absence Generation

In species distribution modeling, the presence of a species in certain locations is often well-recorded, but the absence is typically under-reported or not reported at all. This creates a challenge when trying to understand the complete distribution of a species. To address this, the concept of pseudo-absence data is introduced. Pseudo-absence data are artificially generated points that represent locations where the species is assumed not to be present. In this section, pseudo-absence data points are generated for the eight bird species across the four states ensuring that these points do not overlap with observed presence data. The generated pseudo-absence data helps in providing a more comprehensive view of the species distribution and serves as an essential component for the upcoming modeling process.

### Buffering Raster Data

For the method used for generating pseudo-absence data in this analysis, initially a buffer is created around each observation point, each with a 5000 meter radius. This buffer essentially “blocks” the regions of the sample-area from where pseudo-absence points can be sampled. I.e., only the non-buffered zones can be sampled from.

```
# Set up output directory
output.dir <- "artifacts/masks_5k"
if (!dir.exists("artifacts")) dir.create("artifacts")
if (!dir.exists(output.dir)) dir.create(output.dir)

# Load Data

# Get raster data
states <- c("CO", "NC", "OR", "VT")
r.list <- purrr::map(paste0("data/final_rasters/", states, ".tif"), rast)
names(r.list) <- states

# Get observation data
obs.df <- list.files("data/final", full.names = T) %>%
  purrr::map_df(readRDS) %>%
  dplyr::select(state, common.name, observation.point=geometry)

# Buffering Raster Data

dist <- 5e3

mask.update <- function(i, mask.raster, obs.df, obs.field="observation.point",
                        dist=10000, u="m") {
  obs.pt <- st_transform(obs.df[i, "observation.point"], st_crs(mask.raster))
  # Create a buffer around the point, ensuring correct units
  buf <- st_buffer(obs.pt, dist=units::as_units(paste(dist, u)))
  return(terra::rasterize(buf, mask.raster, update=T, field=1))
}

# For each observation point, you can now create a distance 
# raster and then mask cells within the buffer distance
get.buffered.zones <- function(r, obs.df, obs.field="observation.point",
                               dist=10000, u="m") {
  # Create an empty raster with the same extent and resolution as r
  mask.raster <- terra::rast(r)
  # Recursively update mask.raster with additional buffered regions
  for(i in 1:nrow(obs.df)) {
    mask.raster <- mask.update(i, mask.raster, obs.df, 
                               obs.field=obs.field, 
                               dist=dist, u=u)
    gc()
  }
  return(mask.raster)
}

# Get masks by state, species
masks <- purrr::map(states, function(state) {
    specs <- sort(unique(obs.df$common.name))
    spec.masks <- purrr::map(specs, function(spec, st=state) {
      fname <- file.path(output.dir, paste0(st, "_", spec, ".tif"))
      if (file.exists(fname)) {
        cat("Reading", spec, st, "mask from", fname, "\n")
        r.mask <- rast(fname)
      } else {
        cat("Computing", spec, st, "mask, and saving to", fname, "\n")
        r.mask <- get.buffered.zones(r=r.list[[st]], 
                                     obs.df=filter(obs.df, state == st & 
                                                     common.name == spec),
                                     dist=dist)
        terra::writeRaster(r.mask, fname, overwrite=T)
      }
      gc()
      r.mask
    }, .progress=T)
    names(spec.masks) <- specs
    spec.masks
  })
names(masks) <- states
```

### Sampling Pseudo-Absences and Comparing Totals

This section identifies areas outside the buffers as potential zones for pseudo-absences. To ensure an accurate representation, a random sampling mechanism is implemented. For each bird species in a state, an equivalent number of pseudo-absence points (or a pre-defined minimum threshold), based on observed data are generated. The result is a set of coordinates representing regions where the bird species have not been observed.

```
# Set seed
set.seed(19)

# Function to sample n points from the non-masked parts
sample.inverse.mask <- function(r.original, r.mask, n, 
                                species, state,
                                sample.crs=4326, min.n=500,
                                output.dir="artifacts/pseudo_absence_regions") {
  if (!dir.exists(output.dir)) dir.create(output.dir)
  output.path <- file.path(output.dir,
                           gsub(" |\\-", "_", 
                                paste0(
                                  paste(state, tolower(species), sep="_"), 
                                  ".tif")
                           ))
  if (!file.exists(output.path)) {
    # Get inverse mask;
    # Set NA cells to 0, keep 0 cells as 0, change other cells to 1
    r.inverse <- terra::ifel(is.na(r.mask), 0, r.mask)
    # Set 0 to 1 and everything else to NA
    r.inverse <- terra::lapp(r.inverse, fun = function(x) ifelse(x == 0, 1, NA))
    # Crop so that anything outside of the state is NA
    r.cropped <- terra::crop(r.inverse, terra::ext(r.original))
    
    # Create a binary raster from r.original where valid values are 
    # set to 1 and NA values remain NA
    r.binary <- terra::lapp(r.original[[1]], 
                            fun = function(x) ifelse(!is.na(x), 1, NA))
    
    # Multiply the cropped raster by the binary raster to ensure 
    # outside values are set to NA
    r.final <- r.cropped * r.binary
    terra::writeRaster(r.final, output.path, overwrite=T)
  } else {
    r.final <- terra::rast(output.path)
  }
  
  # Convert the raster to SpatialPoints
  gdf <- terra::as.points(r.final) %>%
    st_as_sf() %>%
    st_transform(crs=sample.crs)
  if (nrow(gdf) > 0) {
    gdf <- gdf %>%
      filter(!is.na(layer)) %>%
      select(geometry)
  } else {
    return(gdf)
  }
  
  # Set to min.n size if n < min.n
  if (n < min.n) n <- min.n
  # Make sure there is sufficient available sample points
  if (n > nrow(gdf)) n <- nrow(gdf)
  
  # Randomly sample n points from the available (non-masked) space
  sample.idx <- sample(1:nrow(gdf), n)
  samples <- gdf[sample.idx,] %>%
    mutate(common.name = species, 
           state = state, 
           lon = NA, 
           lat = NA,
           observations=0)
  
  # Populate lon and lat values:
  coords <- st_coordinates(samples)
  samples$lon <- coords[, "X"]
  samples$lat <- coords[, "Y"]
  
  return(samples)
}

# Get totals by species and state
totals <- obs.df %>%
  as_tibble() %>%
  select(state, common.name) %>%
  group_by(state, common.name) %>%
  summarize(N=n(), .groups="keep")


if (!dir.exists(file.path("data", "final_pseudo_absence"))) {
  dir.create(file.path("data", "final_pseudo_absence"))
}

if (!all(
  file.exists(
    paste0(file.path("data", "final_pseudo_absence", 
                     paste0("pa_", states, ".rds")))
  )
)) {
  # Create a list of pseudo absence points, by species and state,
  # where the sample number `n` >= 500 | `n` == the total observed
  # for each respective state and species
  pseudo.absence.pts <- list()
  for (st in states) {
    r.original <- r.list[[st]]
    r.masks <- masks[[st]]
    pseudo.absence.pts[[st]] <- list()
    for (spec in names(r.masks)) {
      r.mask <- r.masks[[spec]]
      n <- totals %>% filter(state == st & common.name == spec) %>% pull(N)
      cat("Generating", n, "pseudo-absence points for the", spec, "in", st, "\n")
      pseudo.absence.pts[[st]][[spec]] <- sample.inverse.mask(
        r.original, r.mask, spec, st, n=n, sample.crs=4326)
      cat("\tGenerated", nrow(pseudo.absence.pts[[st]][[spec]]), "/", n, "points.\n")
    }
  }
  
  # Extract raster values for each point
  for (state in states) {
    out.file.all <- file.path("data", "final_pseudo_absence", paste0("pa_", state, ".rds"))
    if (!file.exists(out.file.all)) {
      r <- r.list[[state]]
      
      cat(sprintf("Extracting points to values for %s...\n", state))
      # Load observations shapefile
      geo.df <- pseudo.absence.pts[[state]] %>% do.call("rbind", .)
      rownames(geo.df) <- NULL
      
      geo.df.crs <- st_crs(geo.df)
      
      # Define target CRS and update
      target.crs <- st_crs(r)
      cat(sprintf("Updating CRS for %s...\n", state))
      geo.df <- st_transform(geo.df, target.crs)
      
      # Extract raster values
      for (r.name in names(r)) {
        cat("\tExtracting", r.name, "values for", state, "\n")
        x <- terra::extract(r[[r.name]], geo.df)[[r.name]]
        geo.df[[gsub(paste0("_", state), "", r.name)]] <- x
      }
      
      # Update crs back
      geo.df <- st_transform(geo.df, geo.df.crs)
      
      # Fix names; Filter NA values
      geo.df <- geo.df %>%
        filter(dplyr::if_all(names(.), ~!is.na(.))) %>%
        suppressWarnings() 
      
      saveRDS(geo.df, out.file.all)
      cat("--------------\n")
    }
  }
}
```

The table below compares the generated pseudo-absence data with the observed dataset. This comparison ensures that the number of pseudo-absence points matches the observed ones, balancing the dataset and making it ready for the next analysis steps.

```
# Get all pseudo-absence data
abs.df <- list.files(file.path("data", "final_pseudo_absence"), full.names = T) %>%
  purrr::map_df(readRDS) %>%
  select(state, common.name, observation.point=geometry)

# There might be some slight differences since there are occasionally NA values
abs.df %>%
  as_tibble() %>%
  select(state, common.name) %>%
  group_by(state, common.name) %>%
  summarize(psuedo.absence.N=n(), .groups="keep") %>%
  left_join(totals, by=c("state", "common.name")) %>%
  rename(observed.N = N) %>%
  knitr::kable()
```

<table style="min-width:320px; overflow-x:auto;">
    <thead>
        <tr>
            <th align="left">state</th>
            <th align="left">common.name</th>
            <th align="right">psuedo.absence.N</th>
            <th align="right">observed.N</th>
        </tr>
    </thead>
    <tbody>
        <tr>
        <td align="left">CO</td>
        <td align="left">Belted Kingfisher</td>
        <td align="right">4523</td>
        <td align="right">4551</td>
        </tr>
        <tr>
        <td align="left">CO</td>
        <td align="left">Cedar Waxwing</td>
        <td align="right">3431</td>
        <td align="right">3446</td>
        </tr>
        <tr>
        <td align="left">CO</td>
        <td align="left">Downy Woodpecker</td>
        <td align="right">7440</td>
        <td align="right">7511</td>
        </tr>
        <tr>
        <td align="left">CO</td>
        <td align="left">Ruddy Duck</td>
        <td align="right">1715</td>
        <td align="right">1726</td>
        </tr>
        <tr>
        <td align="left">CO</td>
        <td align="left">Sanderling</td>
        <td align="right">497</td>
        <td align="right">131</td>
        </tr>
        <tr>
        <td align="left">CO</td>
        <td align="left">Sandhill Crane</td>
        <td align="right">1524</td>
        <td align="right">1532</td>
        </tr>
        <tr>
        <td align="left">CO</td>
        <td align="left">Sharp-shinned Hawk</td>
        <td align="right">2241</td>
        <td align="right">2257</td>
        </tr>
        <tr>
        <td align="left">CO</td>
        <td align="left">Wild Turkey</td>
        <td align="right">2597</td>
        <td align="right">2611</td>
        </tr>
        <tr>
        <td align="left">NC</td>
        <td align="left">Belted Kingfisher</td>
        <td align="right">4042</td>
        <td align="right">4183</td>
        </tr>
        <tr>
        <td align="left">NC</td>
        <td align="left">Cedar Waxwing</td>
        <td align="right">4059</td>
        <td align="right">4195</td>
        </tr>
        <tr>
        <td align="left">NC</td>
        <td align="left">Downy Woodpecker</td>
        <td align="right">10415</td>
        <td align="right">10914</td>
        </tr>
        <tr>
        <td align="left">NC</td>
        <td align="left">Ruddy Duck</td>
        <td align="right">1076</td>
        <td align="right">1106</td>
        </tr>
        <tr>
        <td align="left">NC</td>
        <td align="left">Sanderling</td>
        <td align="right">488</td>
        <td align="right">311</td>
        </tr>
        <tr>
        <td align="left">NC</td>
        <td align="left">Sandhill Crane</td>
        <td align="right">489</td>
        <td align="right">118</td>
        </tr>
        <tr>
        <td align="left">NC</td>
        <td align="left">Sharp-shinned Hawk</td>
        <td align="right">1211</td>
        <td align="right">1254</td>
        </tr>
        <tr>
        <td align="left">NC</td>
        <td align="left">Wild Turkey</td>
        <td align="right">2261</td>
        <td align="right">2372</td>
        </tr>
        <tr>
        <td align="left">OR</td>
        <td align="left">Belted Kingfisher</td>
        <td align="right">5803</td>
        <td align="right">5837</td>
        </tr>
        <tr>
        <td align="left">OR</td>
        <td align="left">Cedar Waxwing</td>
        <td align="right">8405</td>
        <td align="right">8446</td>
        </tr>
        <tr>
        <td align="left">OR</td>
        <td align="left">Downy Woodpecker</td>
        <td align="right">8529</td>
        <td align="right">8576</td>
        </tr>
        <tr>
        <td align="left">OR</td>
        <td align="left">Ruddy Duck</td>
        <td align="right">1996</td>
        <td align="right">2010</td>
        </tr>
        <tr>
        <td align="left">OR</td>
        <td align="left">Sanderling</td>
        <td align="right">496</td>
        <td align="right">258</td>
        </tr>
        <tr>
        <td align="left">OR</td>
        <td align="left">Sandhill Crane</td>
        <td align="right">2443</td>
        <td align="right">2458</td>
        </tr>
        <tr>
        <td align="left">OR</td>
        <td align="left">Sharp-shinned Hawk</td>
        <td align="right">2696</td>
        <td align="right">2714</td>
        </tr>
        <tr>
        <td align="left">OR</td>
        <td align="left">Wild Turkey</td>
        <td align="right">2429</td>
        <td align="right">2440</td>
        </tr>
        <tr>
        <td align="left">VT</td>
        <td align="left">Belted Kingfisher</td>
        <td align="right">1956</td>
        <td align="right">2033</td>
        </tr>
        <tr>
        <td align="left">VT</td>
        <td align="left">Cedar Waxwing</td>
        <td align="right">1098</td>
        <td align="right">3938</td>
        </tr>
        <tr>
        <td align="left">VT</td>
        <td align="left">Downy Woodpecker</td>
        <td align="right">1598</td>
        <td align="right">4635</td>
        </tr>
        <tr>
        <td align="left">VT</td>
        <td align="left">Ruddy Duck</td>
        <td align="right">490</td>
        <td align="right">51</td>
        </tr>
        <tr>
        <td align="left">VT</td>
        <td align="left">Sanderling</td>
        <td align="right">493</td>
        <td align="right">39</td>
        </tr>
        <tr>
        <td align="left">VT</td>
        <td align="left">Sandhill Crane</td>
        <td align="right">492</td>
        <td align="right">76</td>
        </tr>
        <tr>
            <td align="left">VT</td>
            <td align="left">Sharp-shinned Hawk</td>
            <td align="right">730</td>
            <td align="right">748</td>
        </tr>
        <tr>
            <td align="left">VT</td>
            <td align="left">Wild Turkey</td>
            <td align="right">2090</td>
            <td align="right">2181</td>
        </tr>
    </tbody>
</table>

## Identifying Spatial Autocorrelation

### Moran’s I

This is a common measure of global spatial autocorrelation. A positive Moran’s I suggests clustering, a negative value suggests dispersion, and a value near zero suggests randomness. In other words, this is randomization approach to test for spatial autocorrelation. It is essentially checking if the data has a spatial pattern that is significantly different from what would be expected if the values were randomly distributed in space.

### Understanding the Results

Moran I statistic: This is the calculated Moran’s I value for the data, which quantifies the degree of spatial autocorrelation. A positive value indicates positive autocorrelation (similar values are closer together), while a negative value indicates negative autocorrelation (dissimilar values are closer together). A value close to zero indicates a random spatial pattern.

Moran I statistic standard deviate: This value represents the standardized Moran’s I value. The larger the absolute value of this statistic, the stronger the evidence against the null hypothesis (of no spatial autocorrelation). A positive value indicates positive spatial autocorrelation, suggesting clustering of similar values.

P-value: This is the probability of observing a Moran’s I value as extreme as, or more extreme than, the one computed from the data, assuming the null hypothesis of no spatial autocorrelation is true. A very small p-value provides strong evidence to reject the null hypothesis, indicating significant spatial autocorrelation in your data.

Expectation: This is the expected Moran’s I value under the null hypothesis of no spatial autocorrelation. For a large dataset, it’s typically close to zero.

Variance: This is the variance of the Moran’s I statistic under the null hypothesis.
Spatial Autocorrelation of Observations

```
if (!dir.exists("artifacts/obs_autocor_morans")) {
  dir.create("artifacts/obs_autocor_morans")
}

# Get observation data
df <- list.files(file.path("data", "final"), full.names = T) %>%
  purrr::map_df(readRDS) %>%
  select(state, common.name, observations, geometry) # We only interested in these data here

states <- sort(unique(df$state))
species <- sort(unique(df$common.name))

# Function to extract the desired results from output
extract.data <- function(st, spec, results) {
  data_frame(
    state = st,
    species = spec,
    statistic = results$statistic,
    p.value = results$p.value,
    moran.i.statistic = results$estimate['Moran I statistic'],
    expectation = results$estimate['Expectation'],
    variance = results$estimate['Variance']
  )
}

# Define the k for knn computation
k <- 5

if (!file.exists("artifacts/obs_autocor_morans/mi.results.rds")) {
  # Perform test by state, by species
  mi.results <- list()
  
  for (st in states) {
    mi.results[[st]] <- list()
    for (spec in species) {
      cat("Doing Moran's test for", spec, "in", st, "\n")
      # Filter data
      d <- df %>% filter(common.name == spec & state == st)
      mi.results[[st]][[spec]]$data <- d
      # Fit knn model
      knn <- d %>%
        knearneigh(k = k)
      mi.results[[st]][[spec]]$knn <- knn
      # Create a neighbor's list
      nb <- knn %>%
        knn2nb()
      mi.results[[st]][[spec]]$nb <- nb
      # Create a spatial weights matrix
      listw <- nb2listw(nb)
      mi.results[[st]][[spec]]$listw <- listw
      # Compute Moran's I
      results <- moran.test(d$observations, listw)
      mi.results[[st]][[spec]]$moran.test.results <- extract.data(st, spec, results)
    }
  }
  saveRDS(mi.results, "artifacts/obs_autocor_morans/mi.results.rds")
} else {
  mi.results <- readRDS("artifacts/obs_autocor_morans/mi.results.rds")
}

spec.stat <- expand.grid(species=species, state=states, stringsAsFactors=F)

mi.results.df <- purrr::map(1:nrow(spec.stat), function(i) {
    spec <- spec.stat[i, ]$species
    st <- spec.stat[i, ]$state
    mi.results[[st]][[spec]]$moran.test.results
  }) %>% 
  do.call("rbind", .) %>%
  as_tibble() %>%
  # Compute adjusted p-value, accounting for multiple comparisons
  mutate(adj.p.value = p.adjust(p.value, method="holm"))
```

```
# Filter to where the adjusted p-value <= 0.05; sort by moran i stat
mi.results.df %>% 
  filter(adj.p.value <= 0.05) %>%
  arrange(-moran.i.statistic)
```

<div data-pagedtable="false" pagedtable-page="0" class="pagedtable-wrapper">
<script data-pagedtable-source="" type="application/json">
{"columns":[{"label":["state"],"name":[1],"type":["chr"],"align":["left"]},{"label":["species"],"name":[2],"type":["chr"],"align":["left"]},{"label":["statistic"],"name":[3],"type":["dbl"],"align":["right"]},{"label":["p.value"],"name":[4],"type":["dbl"],"align":["right"]},{"label":["moran.i.statistic"],"name":[5],"type":["dbl"],"align":["right"]},{"label":["expectation"],"name":[6],"type":["dbl"],"align":["right"]},{"label":["variance"],"name":[7],"type":["dbl"],"align":["right"]},{"label":["adj.p.value"],"name":[8],"type":["dbl"],"align":["right"]}],"data":[{"1":"VT","2":"Sharp-shinned Hawk","3":"12.837294","4":"5.067812e-38","5":"0.12932887","6":"-0.0013386881","7":"1.036069e-04","8":"1.621700e-36"},{"1":"OR","2":"Belted Kingfisher","3":"5.686901","4":"6.468259e-09","5":"0.04018657","6":"-0.0001713502","7":"5.036238e-05","8":"2.005160e-07"},{"1":"OR","2":"Sharp-shinned Hawk","3":"3.493483","4":"2.383818e-04","5":"0.02846832","6":"-0.0003685957","7":"6.813658e-05","8":"6.913073e-03"},{"1":"OR","2":"Cedar Waxwing","3":"4.213068","4":"1.259624e-05","5":"0.02497281","6":"-0.0001184133","7":"3.546882e-05","8":"3.778871e-04"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
</script>
</div>

The Moran’s I statistic for the Belted Kingfisher and Cedar Waxwing are positive in Oregon, as well as the Sharp-shinned Hawk in both Vermont and Oregon. This suggests that locations with high observations of these species are close to other locations with high observations, and the same for low observations. In simple terms, the observation patterns of these species are clustered.

Below is an example of the actual observations of the Sharp-shinned Hawk in Vermont against the spatial lag (i.e., the weighted sum of the neighboring values at each point). Consider the following when interpreting the plot:

- First Quadrant (top-right): High values surrounded by high values (indicating clustering of high values).
- Second Quadrant (top-left): Low values surrounded by high values.
- Third Quadrant (bottom-left): Low values surrounded by low values (indicating clustering of low values).
- Fourth Quadrant (bottom-right): High values surrounded by low values.

```
# Get spatial weights list
vt.ssh <- mi.results[["VT"]][["Sharp-shinned Hawk"]]$data

# Calculate the spatial lag of the observations
vt.ssh$spatial.lag <- lag.listw(
  mi.results[["VT"]][["Sharp-shinned Hawk"]]$listw,
  vt.ssh$observations)

# Scatter plot for the Sharp-Shinned Hawk in Vermont
ggplot(vt.ssh, aes(x=observations, y=spatial.lag)) +
  geom_point() +
  geom_smooth(method="lm", col="red") + # Adds a linear regression line
  ggtitle("Moran Scatter Plot for Sharp-shinned Hawk in VT") +
  xlab("Observations (log10 Scale)") +
  ylab("Spatial Lag of Observations (log10 Scale)") +
  scale_y_log10() + scale_x_log10()
```

<img class="sdm-eda-img" src="{{ site.baseurl }}/assets/plots/sdm-eda-14.png">

### Spatial Autocorrelation of Environmental Factors

```
env.df.list <- list()
for (state in states) {
  r <- r.list[[state]]
  
  gdf <- terra::as.points(r) %>% 
    st_as_sf() %>%
    st_transform(crs=4326)
  # Fix names
  names(gdf) <- gsub(paste0("_", state, "$"), "", names(gdf))
  
  # Convert land cover to binary variables
  binary.lc <- caret::dummyVars(~NLCD_Land, data=gdf, sep=".") %>% 
    predict(., gdf) %>%
    as.data.frame()
  names(binary.lc) <- gsub("NLCD_Land", "", names(binary.lc)) %>% 
    gsub("\\/| |,", "_", .) %>% 
    gsub("__", "_", .) %>% 
    tolower()
  
  gdf <- gdf %>% select(-NLCD_Land) %>% cbind(binary.lc)
  
  env.df.list[[state]] <- gdf
}

# Perform test by state, by environmental variable
env.mi.results <- list()

output.dir <- file.path("artifacts", "env_autocor_morans")
if (!dir.exists(output.dir)) dir.create(output.dir)

for (st in states) {
  env.mi.results[[st]] <- list()
  gdf <- env.df.list[[st]]
  env.vars <- names(gdf)[names(gdf) != "geometry"]
  for (e in env.vars) {
    output.path <- file.path(output.dir, paste0(st, "_", e, ".rds"))
    if (!file.exists(output.path)) {
      cat("Doing Moran's test for", e, "in", st, "\n")
      # Filter data
      d <- gdf %>% select(!!sym(e))
      # env.mi.results[[st]][[e]]$data <- d
      n.distinct <- d %>% pull(!!sym(e)) %>% unique() %>% length()
      if (n.distinct > 1) {
        # Fit knn model
        knn <- d %>%
          knearneigh(k = k)
        # env.mi.results[[st]][[e]]$knn <- knn
        # Create a neighbor's list
        nb <- knn %>%
          knn2nb()
        # env.mi.results[[st]][[e]]$nb <- nb
        # Create a spatial weights matrix
        listw <- nb2listw(nb)
        # env.mi.results[[st]][[e]]$listw <- listw
        # Compute Moran's I
        results <- moran.test(d[[e]], listw)
        env.mi.results[[st]][[e]]$moran.test.results <- extract.data(st, e, results)
      }
      cat("\tSaving results...\n")
      saveRDS(env.mi.results[[st]][[e]]$moran.test.results, output.path)
    } else {
      cat("Getting Moran's test results for", e, "in", st, "\n")
      env.mi.results[[st]][[e]]$moran.test.results <- readRDS(output.path)
    }
  }
}


env.stat <- expand.grid(env=names(env.df.list$CO )[names(env.df.list$CO) != "geometry"], 
                        state=states, stringsAsFactors=F)

env.mi.results.df <- purrr::map(1:nrow(env.stat), function(i) {
    e <- env.stat[i, ]$env
    st <- env.stat[i, ]$state
    env.mi.results[[st]][[e]]$moran.test.results
  }) %>% 
  do.call("rbind", .) %>%
  as_tibble() %>%
  # Compute adjusted p-value, accounting for multiple comparisons
  mutate(adj.p.value = p.adjust(p.value, method="holm"))
```

```
# Filter to where the adjusted p-value <= 0.05; sort by moran i stat
env.mi.results.df %>% 
  filter(adj.p.value <= 0.05) %>%
  arrange(-moran.i.statistic)
```

<div data-pagedtable="false" pagedtable-page="0" class="pagedtable-wrapper">
<script data-pagedtable-source="" type="application/json">
{"columns":[{"label":["state"],"name":[1],"type":["chr"],"align":["left"]},{"label":["species"],"name":[2],"type":["chr"],"align":["left"]},{"label":["statistic"],"name":[3],"type":["dbl"],"align":["right"]},{"label":["p.value"],"name":[4],"type":["dbl"],"align":["right"]},{"label":["moran.i.statistic"],"name":[5],"type":["dbl"],"align":["right"]},{"label":["expectation"],"name":[6],"type":["dbl"],"align":["right"]},{"label":["variance"],"name":[7],"type":["dbl"],"align":["right"]},{"label":["adj.p.value"],"name":[8],"type":["dbl"],"align":["right"]}],"data":[{"1":"OR","2":"avg_prcp","3":"833.671668","4":"0.000000e+00","5":"0.99940437","6":"-3.989738e-06","7":"1.437129e-06","8":"0.000000e+00"},{"1":"OR","2":"tmax","3":"833.619806","4":"0.000000e+00","5":"0.99934156","6":"-3.989738e-06","7":"1.437128e-06","8":"0.000000e+00"},{"1":"OR","2":"tmin","3":"833.617716","4":"0.000000e+00","5":"0.99934108","6":"-3.989738e-06","7":"1.437133e-06","8":"0.000000e+00"},{"1":"NC","2":"tmin","3":"591.448645","4":"0.000000e+00","5":"0.99794424","6":"-7.905076e-06","7":"2.846986e-06","8":"0.000000e+00"},{"1":"NC","2":"avg_prcp","3":"591.224399","4":"0.000000e+00","5":"0.99756283","6":"-7.905076e-06","7":"2.846969e-06","8":"0.000000e+00"},{"1":"NC","2":"tmax","3":"590.924561","4":"0.000000e+00","5":"0.99706426","6":"-7.905076e-06","7":"2.847011e-06","8":"0.000000e+00"},{"1":"VT","2":"tmax","3":"261.784457","4":"0.000000e+00","5":"0.99702234","6":"-4.023659e-05","7":"1.450632e-05","8":"0.000000e+00"},{"1":"CO","2":"tmax","3":"861.805250","4":"0.000000e+00","5":"0.99644583","6":"-3.713345e-06","7":"1.336880e-06","8":"0.000000e+00"},{"1":"CO","2":"avg_prcp","3":"861.487288","4":"0.000000e+00","5":"0.99607690","6":"-3.713345e-06","7":"1.336876e-06","8":"0.000000e+00"},{"1":"VT","2":"avg_prcp","3":"261.492767","4":"0.000000e+00","5":"0.99591050","6":"-4.023659e-05","7":"1.450629e-05","8":"0.000000e+00"},{"1":"CO","2":"tmin","3":"861.240311","4":"0.000000e+00","5":"0.99578942","6":"-3.713345e-06","7":"1.336871e-06","8":"0.000000e+00"},{"1":"NC","2":"dem","3":"590.129109","4":"0.000000e+00","5":"0.99571749","6":"-7.905076e-06","7":"2.846984e-06","8":"0.000000e+00"},{"1":"VT","2":"tmin","3":"261.398001","4":"0.000000e+00","5":"0.99554869","6":"-4.023659e-05","7":"1.450626e-05","8":"0.000000e+00"},{"1":"CO","2":"dem","3":"860.345727","4":"0.000000e+00","5":"0.99476737","6":"-3.713345e-06","7":"1.336904e-06","8":"0.000000e+00"},{"1":"OR","2":"dem","3":"827.075401","4":"0.000000e+00","5":"0.99150303","6":"-3.989738e-06","7":"1.437148e-06","8":"0.000000e+00"},{"1":"OR","2":"Fall_NDVI","3":"822.520998","4":"0.000000e+00","5":"0.98604461","6":"-3.989738e-06","7":"1.437152e-06","8":"0.000000e+00"},{"1":"OR","2":"Spring_NDVI","3":"819.933014","4":"0.000000e+00","5":"0.98294123","6":"-3.989738e-06","7":"1.437149e-06","8":"0.000000e+00"},{"1":"OR","2":"Summer_NDVI","3":"819.118214","4":"0.000000e+00","5":"0.98196452","6":"-3.989738e-06","7":"1.437150e-06","8":"0.000000e+00"},{"1":"OR","2":"coastline","3":"818.942867","4":"0.000000e+00","5":"0.98157279","6":"-3.989738e-06","7":"1.436618e-06","8":"0.000000e+00"},{"1":"OR","2":"Winter_NDVI","3":"817.388495","4":"0.000000e+00","5":"0.97989088","6":"-3.989738e-06","7":"1.437149e-06","8":"0.000000e+00"},{"1":"NC","2":"waterbody","3":"580.494513","4":"0.000000e+00","5":"0.97947304","6":"-7.905076e-06","7":"2.847054e-06","8":"0.000000e+00"},{"1":"OR","2":"waterbody","3":"813.114548","4":"0.000000e+00","5":"0.97474214","6":"-3.989738e-06","7":"1.437075e-06","8":"0.000000e+00"},{"1":"VT","2":"waterbody","3":"254.230259","4":"0.000000e+00","5":"0.96821930","6":"-4.023659e-05","7":"1.450538e-05","8":"0.000000e+00"},{"1":"CO","2":"Fall_NDVI","3":"836.229712","4":"0.000000e+00","5":"0.96688386","6":"-3.713345e-06","7":"1.336906e-06","8":"0.000000e+00"},{"1":"CO","2":"waterbody","3":"834.947812","4":"0.000000e+00","5":"0.96534888","6":"-3.713345e-06","7":"1.336759e-06","8":"0.000000e+00"},{"1":"NC","2":"coastline","3":"568.891841","4":"0.000000e+00","5":"0.95983472","6":"-7.905076e-06","7":"2.846693e-06","8":"0.000000e+00"},{"1":"CO","2":"Winter_NDVI","3":"824.960759","4":"0.000000e+00","5":"0.95385255","6":"-3.713345e-06","7":"1.336901e-06","8":"0.000000e+00"},{"1":"VT","2":"dem","3":"250.389778","4":"0.000000e+00","5":"0.95362834","6":"-4.023659e-05","7":"1.450647e-05","8":"0.000000e+00"},{"1":"CO","2":"Spring_NDVI","3":"812.443689","4":"0.000000e+00","5":"0.93937746","6":"-3.713345e-06","7":"1.336895e-06","8":"0.000000e+00"},{"1":"CO","2":"Summer_NDVI","3":"807.996253","4":"0.000000e+00","5":"0.93423607","6":"-3.713345e-06","7":"1.336897e-06","8":"0.000000e+00"},{"1":"VT","2":"Fall_NDVI","3":"242.225424","4":"0.000000e+00","5":"0.92202577","6":"-4.023659e-05","7":"1.449054e-05","8":"0.000000e+00"},{"1":"CO","2":"urban_imperviousness","3":"792.248129","4":"0.000000e+00","5":"0.91587657","6":"-3.713345e-06","7":"1.336457e-06","8":"0.000000e+00"},{"1":"OR","2":"urban_imperviousness","3":"744.158947","4":"0.000000e+00","5":"0.89192159","6":"-3.989738e-06","7":"1.436567e-06","8":"0.000000e+00"},{"1":"VT","2":"Spring_NDVI","3":"232.522434","4":"0.000000e+00","5":"0.88556714","6":"-4.023659e-05","7":"1.450617e-05","8":"0.000000e+00"},{"1":"VT","2":"Summer_NDVI","3":"232.398958","4":"0.000000e+00","5":"0.88486681","6":"-4.023659e-05","7":"1.449863e-05","8":"0.000000e+00"},{"1":"NC","2":"urban_imperviousness","3":"520.976475","4":"0.000000e+00","5":"0.87897910","6":"-7.905076e-06","7":"2.846615e-06","8":"0.000000e+00"},{"1":"VT","2":"Winter_NDVI","3":"229.378529","4":"0.000000e+00","5":"0.87342710","6":"-4.023659e-05","7":"1.450066e-05","8":"0.000000e+00"},{"1":"NC","2":"Spring_NDVI","3":"505.494656","4":"0.000000e+00","5":"0.85289143","6":"-7.905076e-06","7":"2.846836e-06","8":"0.000000e+00"},{"1":"VT","2":"urban_imperviousness","3":"220.776087","4":"0.000000e+00","5":"0.84015815","6":"-4.023659e-05","7":"1.448304e-05","8":"0.000000e+00"},{"1":"NC","2":"Winter_NDVI","3":"489.419346","4":"0.000000e+00","5":"0.82576938","6":"-7.905076e-06","7":"2.846844e-06","8":"0.000000e+00"},{"1":"NC","2":"Fall_NDVI","3":"480.281737","4":"0.000000e+00","5":"0.81025929","6":"-7.905076e-06","7":"2.846193e-06","8":"0.000000e+00"},{"1":"NC","2":"Summer_NDVI","3":"479.963111","4":"0.000000e+00","5":"0.80982225","6":"-7.905076e-06","7":"2.846900e-06","8":"0.000000e+00"},{"1":"CO","2":"cultivated_crops","3":"535.590851","4":"0.000000e+00","5":"0.61926639","6":"-3.713345e-06","7":"1.336885e-06","8":"0.000000e+00"},{"1":"OR","2":"evergreen_forest","3":"503.986786","4":"0.000000e+00","5":"0.60418179","6":"-3.989738e-06","7":"1.437152e-06","8":"0.000000e+00"},{"1":"CO","2":"herbaceous","3":"521.029377","4":"0.000000e+00","5":"0.60243474","6":"-3.713345e-06","7":"1.336906e-06","8":"0.000000e+00"},{"1":"OR","2":"cultivated_crops","3":"499.636842","4":"0.000000e+00","5":"0.59894604","6":"-3.989738e-06","7":"1.437051e-06","8":"0.000000e+00"},{"1":"VT","2":"open_water","3":"153.955240","4":"0.000000e+00","5":"0.58608936","6":"-4.023659e-05","7":"1.449434e-05","8":"0.000000e+00"},{"1":"OR","2":"shrub_scrub","3":"451.974178","4":"0.000000e+00","5":"0.54182832","6":"-3.989738e-06","7":"1.437151e-06","8":"0.000000e+00"},{"1":"CO","2":"evergreen_forest","3":"423.401814","4":"0.000000e+00","5":"0.48955197","6":"-3.713345e-06","7":"1.336900e-06","8":"0.000000e+00"},{"1":"OR","2":"herbaceous","3":"402.811289","4":"0.000000e+00","5":"0.48288844","6":"-3.989738e-06","7":"1.437135e-06","8":"0.000000e+00"},{"1":"OR","2":"open_water","3":"390.048408","4":"0.000000e+00","5":"0.46749207","6":"-3.989738e-06","7":"1.436544e-06","8":"0.000000e+00"},{"1":"OR","2":"barren_land","3":"385.416929","4":"0.000000e+00","5":"0.46186882","6":"-3.989738e-06","7":"1.436095e-06","8":"0.000000e+00"},{"1":"CO","2":"shrub_scrub","3":"369.081308","4":"0.000000e+00","5":"0.42674479","6":"-3.713345e-06","7":"1.336903e-06","8":"0.000000e+00"},{"1":"OR","2":"hay_pasture","3":"355.094476","4":"0.000000e+00","5":"0.42566423","6":"-3.989738e-06","7":"1.436993e-06","8":"0.000000e+00"},{"1":"CO","2":"deciduous_forest","3":"361.799305","4":"0.000000e+00","5":"0.41831651","6":"-3.713345e-06","7":"1.336849e-06","8":"0.000000e+00"},{"1":"NC","2":"woody_wetlands","3":"222.964539","4":"0.000000e+00","5":"0.37620168","6":"-7.905076e-06","7":"2.847004e-06","8":"0.000000e+00"},{"1":"OR","2":"perennial_snow_ice","3":"305.420479","4":"0.000000e+00","5":"0.36101933","6":"-3.989738e-06","7":"1.397250e-06","8":"0.000000e+00"},{"1":"OR","2":"emergent_herbaceous_wetlands","3":"284.413660","4":"0.000000e+00","5":"0.34090477","6":"-3.989738e-06","7":"1.436731e-06","8":"0.000000e+00"},{"1":"NC","2":"deciduous_forest","3":"197.482089","4":"0.000000e+00","5":"0.33320718","6":"-7.905076e-06","7":"2.847042e-06","8":"0.000000e+00"},{"1":"NC","2":"cultivated_crops","3":"179.126922","4":"0.000000e+00","5":"0.30223401","6":"-7.905076e-06","7":"2.847001e-06","8":"0.000000e+00"},{"1":"NC","2":"open_water","3":"168.866418","4":"0.000000e+00","5":"0.28485705","6":"-7.905076e-06","7":"2.845715e-06","8":"0.000000e+00"},{"1":"CO","2":"barren_land","3":"215.824820","4":"0.000000e+00","5":"0.24948233","6":"-3.713345e-06","7":"1.336256e-06","8":"0.000000e+00"},{"1":"VT","2":"deciduous_forest","3":"62.008779","4":"0.000000e+00","5":"0.23614177","6":"-4.023659e-05","7":"1.450732e-05","8":"0.000000e+00"},{"1":"CO","2":"open_water","3":"198.327635","4":"0.000000e+00","5":"0.22920488","6":"-3.713345e-06","7":"1.335658e-06","8":"0.000000e+00"},{"1":"OR","2":"mixed_forest","3":"183.746168","4":"0.000000e+00","5":"0.22025636","6":"-3.989738e-06","7":"1.436932e-06","8":"0.000000e+00"},{"1":"OR","2":"developed_medium_intensity","3":"181.918214","4":"0.000000e+00","5":"0.21799103","6":"-3.989738e-06","7":"1.435955e-06","8":"0.000000e+00"},{"1":"CO","2":"hay_pasture","3":"173.966045","4":"0.000000e+00","5":"0.20111383","6":"-3.713345e-06","7":"1.336506e-06","8":"0.000000e+00"},{"1":"CO","2":"developed_medium_intensity","3":"169.039132","4":"0.000000e+00","5":"0.19537894","6":"-3.713345e-06","7":"1.335972e-06","8":"0.000000e+00"},{"1":"NC","2":"emergent_herbaceous_wetlands","3":"112.404069","4":"0.000000e+00","5":"0.18959654","6":"-7.905076e-06","7":"2.845332e-06","8":"0.000000e+00"},{"1":"VT","2":"hay_pasture","3":"46.528801","4":"0.000000e+00","5":"0.17715969","6":"-4.023659e-05","7":"1.450385e-05","8":"0.000000e+00"},{"1":"NC","2":"evergreen_forest","3":"104.004167","4":"0.000000e+00","5":"0.17547894","6":"-7.905076e-06","7":"2.847001e-06","8":"0.000000e+00"},{"1":"CO","2":"developed_low_intensity","3":"151.027766","4":"0.000000e+00","5":"0.17458088","6":"-3.713345e-06","7":"1.336282e-06","8":"0.000000e+00"},{"1":"NC","2":"hay_pasture","3":"92.256516","4":"0.000000e+00","5":"0.15565363","6":"-7.905076e-06","7":"2.846876e-06","8":"0.000000e+00"},{"1":"OR","2":"developed_high_intensity","3":"129.169402","4":"0.000000e+00","5":"0.15465331","6":"-3.989738e-06","7":"1.433579e-06","8":"0.000000e+00"},{"1":"CO","2":"emergent_herbaceous_wetlands","3":"126.563609","4":"0.000000e+00","5":"0.14631122","6":"-3.713345e-06","7":"1.336471e-06","8":"0.000000e+00"},{"1":"NC","2":"mixed_forest","3":"85.723667","4":"0.000000e+00","5":"0.14463281","6":"-7.905076e-06","7":"2.846949e-06","8":"0.000000e+00"},{"1":"VT","2":"cultivated_crops","3":"36.220428","4":"1.451912e-287","5":"0.13776798","6":"-4.023659e-05","7":"1.447582e-05","8":"2.177868e-286"},{"1":"CO","2":"developed_high_intensity","3":"118.506607","4":"0.000000e+00","5":"0.13687851","6":"-3.713345e-06","7":"1.334163e-06","8":"0.000000e+00"},{"1":"CO","2":"perennial_snow_ice","3":"109.122789","4":"0.000000e+00","5":"0.12589426","6":"-3.713345e-06","7":"1.331087e-06","8":"0.000000e+00"},{"1":"OR","2":"developed_low_intensity","3":"96.735223","4":"0.000000e+00","5":"0.11593157","6":"-3.989738e-06","7":"1.436362e-06","8":"0.000000e+00"},{"1":"NC","2":"developed_low_intensity","3":"68.596337","4":"0.000000e+00","5":"0.11572581","6":"-7.905076e-06","7":"2.846544e-06","8":"0.000000e+00"},{"1":"CO","2":"mixed_forest","3":"96.389391","4":"0.000000e+00","5":"0.11142114","6":"-3.713345e-06","7":"1.336305e-06","8":"0.000000e+00"},{"1":"OR","2":"woody_wetlands","3":"90.133377","4":"0.000000e+00","5":"0.10801537","6":"-3.989738e-06","7":"1.436256e-06","8":"0.000000e+00"},{"1":"NC","2":"barren_land","3":"58.193137","4":"0.000000e+00","5":"0.09796995","6":"-7.905076e-06","7":"2.834734e-06","8":"0.000000e+00"},{"1":"NC","2":"developed_medium_intensity","3":"57.250719","4":"0.000000e+00","5":"0.09656911","6":"-7.905076e-06","7":"2.845677e-06","8":"0.000000e+00"},{"1":"VT","2":"evergreen_forest","3":"25.241820","4":"6.963931e-141","5":"0.09609267","6":"-4.023659e-05","7":"1.450450e-05","8":"9.053110e-140"},{"1":"VT","2":"mixed_forest","3":"23.957639","4":"3.846267e-127","5":"0.09120757","6":"-4.023659e-05","7":"1.450630e-05","8":"4.615521e-126"},{"1":"VT","2":"woody_wetlands","3":"21.418530","4":"4.488921e-102","5":"0.08150907","6":"-4.023659e-05","7":"1.449643e-05","8":"4.937813e-101"},{"1":"NC","2":"developed_high_intensity","3":"46.820488","4":"0.000000e+00","5":"0.07893693","6":"-7.905076e-06","7":"2.842992e-06","8":"0.000000e+00"},{"1":"NC","2":"developed_open_space","3":"43.901404","4":"0.000000e+00","5":"0.07406525","6":"-7.905076e-06","7":"2.846852e-06","8":"0.000000e+00"},{"1":"VT","2":"developed_medium_intensity","3":"19.250367","4":"7.008266e-83","5":"0.07314891","6":"-4.023659e-05","7":"1.445491e-05","8":"6.307440e-82"},{"1":"CO","2":"developed_open_space","3":"47.545334","4":"0.000000e+00","5":"0.05496367","6":"-3.713345e-06","7":"1.336578e-06","8":"0.000000e+00"},{"1":"OR","2":"deciduous_forest","3":"43.441446","4":"0.000000e+00","5":"0.05203238","6":"-3.989738e-06","7":"1.434846e-06","8":"0.000000e+00"},{"1":"CO","2":"woody_wetlands","3":"44.409106","4":"0.000000e+00","5":"0.05133378","6":"-3.713345e-06","7":"1.336365e-06","8":"0.000000e+00"},{"1":"VT","2":"developed_high_intensity","3":"13.321124","4":"8.723086e-41","5":"0.05032680","6":"-4.023659e-05","7":"1.429588e-05","8":"6.106160e-40"},{"1":"VT","2":"emergent_herbaceous_wetlands","3":"12.888868","4":"2.600090e-38","5":"0.04892471","6":"-4.023659e-05","7":"1.443248e-05","8":"1.399931e-37"},{"1":"VT","2":"barren_land","3":"12.897218","4":"2.333219e-38","5":"0.04846861","6":"-4.023659e-05","7":"1.414653e-05","8":"1.399931e-37"},{"1":"OR","2":"developed_open_space","3":"35.819306","4":"2.764503e-281","5":"0.04293137","6":"-3.989738e-06","7":"1.436799e-06","8":"3.870304e-280"},{"1":"NC","2":"herbaceous","3":"19.329792","4":"1.508147e-83","5":"0.03260283","6":"-7.905076e-06","7":"2.846210e-06","8":"1.508147e-82"},{"1":"NC","2":"shrub_scrub","3":"18.847806","4":"1.531433e-79","5":"0.03179056","6":"-7.905076e-06","7":"2.846366e-06","8":"1.225147e-78"},{"1":"VT","2":"developed_low_intensity","3":"7.225024","4":"2.505060e-13","5":"0.02745403","6":"-4.023659e-05","7":"1.448123e-05","8":"1.002024e-12"},{"1":"VT","2":"developed_open_space","3":"5.398189","4":"3.365854e-08","5":"0.02050842","6":"-4.023659e-05","7":"1.449007e-05","8":"1.009756e-07"},{"1":"VT","2":"shrub_scrub","3":"5.043638","4":"2.283819e-07","5":"0.01913703","6":"-4.023659e-05","7":"1.445725e-05","8":"4.567639e-07"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
</script>
</div>

It is apparent that many of the environmental variables exhibit strong spatial clustering patterns. Notably, the weather layers have Moran’s I values nearing 1, signifying near-perfect spatial autocorrelation. Such patterns, while anticipated for certain types of variables, can introduce challenges in statistical analyses. Specifically, strong spatial autocorrelation often violates the assumption of observation independence fundamental to many traditional statistical models, potentially leading to biased parameter estimates and misleading significance levels. However, this spatial dependence can be accounted for using specialized spatial statistical methods or integrated as inherent structure in certain modeling techniques. It’s crucial to select an appropriate modeling method that accommodates or leverages this spatial structure, ensuring reliable and valid results.

### Mitigating Potential Problems Due to Spatial Autocorrelation

1. Multi-Collinearity
- Problem: When predictors are spatially autocorrelated, they may also be correlated with each other. This multi-collinearity can make it hard to determine the influence of individual predictors on the response variable.
- Solution: One method of detecting multi-collinearity is to use the condition index derived from the eigenvalues of a correlation matrix. Variables with high condition indices can be further examined through their contributing eigenvectors to identify key contributors. Furthermore, visualizing correlation matrices can help provide a clearer understanding. Based on these insights, dropping variables, combining variables, or applying dimension reduction techniques can help mitigate mutli-collinearity between features.
2. Spatial Autocorrelation in Residuals
- Problem: Spatially autocorrelated residuals indicate that the model hasn’t captured all spatial structures, leading to biased parameter estimates.
- Solution: Employ spatial regression models that can incorporate spatial structure, such as spatial lag or spatial error models. Alternatively, add spatially structured random effects to more standard models.
3. Overfitting
- Problem: Models that capture too much of the spatial structure in the predictors might perform exceptionally well on training data but poorly on validation or future data.
- Solution: Use spatial cross-validation to ensure models generalize well. Regularization techniques, such as Ridge or LASSO regression, can also help prevent overfitting.
4. Non-Stationarity
- Problem: The relationship between predictors and the response variable might change across space, which traditional models don’t account for.
- Solution: Geographically Weighted Regression (GWR) allows relationships to vary across space, accounting for non-stationarity. Another approach involves segmenting the study area into regions where relationships are stable and modeling each separately (this approach is used during this study, where study areas are separate states).

Multi-collinearity is addressed in the following section, while the other mitigation solutions will be implemented either during the modeling stage of the study (e.g., cross-validation), or else addressed during the feature engineering/selection portion.

#### Checking for Multi-Collinearity in Environmental Factors

```
all.columns <- unique(unlist(lapply(env.df.list, names)))
env.df <- lapply(env.df.list, function(df) {
  missing.cols <- setdiff(all.columns, names(df))
  for (col in missing.cols) {
    df[[col]] <- NA
  }
  return(df)
}) %>%
  do.call("rbind", .) %>%
  as.data.frame() %>%
  select(-geometry, -unclassified)

corr.matrix <- cor(env.df, use="na.or.complete")
eigenvalues <- eigen(corr.matrix)$values
ci.df <- tibble(
  variable=names(env.df),
  condition.index=sqrt(max(eigenvalues) / eigenvalues)
)

ci.df
```

<div data-pagedtable="false" pagedtable-page="0" class="pagedtable-wrapper">
<script data-pagedtable-source="" type="application/json">
{"columns":[{"label":["variable"],"name":[1],"type":["chr"],"align":["left"]},{"label":["condition.index"],"name":[2],"type":["dbl"],"align":["right"]}],"data":[{"1":"coastline","2":"1.000000e+00"},{"1":"dem","2":"1.814842e+00"},{"1":"Fall_NDVI","2":"1.960193e+00"},{"1":"Spring_NDVI","2":"2.004830e+00"},{"1":"Summer_NDVI","2":"2.052971e+00"},{"1":"Winter_NDVI","2":"2.109996e+00"},{"1":"avg_prcp","2":"2.253564e+00"},{"1":"tmax","2":"2.316360e+00"},{"1":"tmin","2":"2.355981e+00"},{"1":"urban_imperviousness","2":"2.380492e+00"},{"1":"waterbody","2":"2.386861e+00"},{"1":"barren_land","2":"2.395021e+00"},{"1":"cultivated_crops","2":"2.397248e+00"},{"1":"deciduous_forest","2":"2.401958e+00"},{"1":"developed_high_intensity","2":"2.403239e+00"},{"1":"developed_low_intensity","2":"2.527208e+00"},{"1":"developed_medium_intensity","2":"2.623054e+00"},{"1":"developed_open_space","2":"2.777975e+00"},{"1":"emergent_herbaceous_wetlands","2":"3.068005e+00"},{"1":"evergreen_forest","2":"3.524406e+00"},{"1":"hay_pasture","2":"3.724769e+00"},{"1":"herbaceous","2":"4.818900e+00"},{"1":"mixed_forest","2":"5.922854e+00"},{"1":"open_water","2":"7.705246e+00"},{"1":"perennial_snow_ice","2":"1.847924e+01"},{"1":"shrub_scrub","2":"6.530425e+01"},{"1":"woody_wetlands","2":"5.418423e+07"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
</script>
</div>

```
# Extract the eigenvectors
eigenvectors <- eigen(corr.matrix)$vectors

threshold <- 30
large.ci.indices <- which(sqrt(max(eigenvalues) / eigenvalues) > threshold)

# Examine the eigenvectors
for (index in large.ci.indices) {
  cat(paste("Eigenvalue:", eigenvalues[index]), "\n")
  cat(paste("Condition Index:", 
                  sqrt(max(eigenvalues) / eigenvalues[index])), "\n")
  
  # Sorting eigenvector components by absolute magnitude for clarity
  ev <- eigenvectors[, index]
  sorted.ev <- sort(abs(ev), decreasing = T)
  ev.top <- sorted.ev[1:5] %>%
    as_tibble() %>%
    mutate(variable=rownames(corr.matrix)[order(-abs(ev))][1:5]) %>%
    select(variable, value)
  
  cat("\nTop 5 contributors to multicollinearity at the condition idx:\n")
  
  for (row in 1:5) {
    cat("\t", ev.top[row,]$variable, ":", signif(ev.top[row,]$value, 3), "\n")
  }
  
  cat("\n------------------------\n")
}
```

<div class="language-plaintext highlighter-rouge"><pre class="highlight code-print"><code class="highlight">
## Eigenvalue: 0.00135666182336175 
## Condition Index: 65.3042479620253 
## 
## Top 5 contributors to multicollinearity at the condition idx:
##   avg_prcp : 0.814 
##   tmax : 0.459 
##   tmin : 0.356 
##   dem : 0.0108 
##   Winter_NDVI : 0.00657 
## 
## ------------------------
## Eigenvalue: 1.97064586234132e-15 
## Condition Index: 54184234.4382929 
## 
## Top 5 contributors to multicollinearity at the condition idx:
##   shrub_scrub : 0.508 
##   evergreen_forest : 0.5 
##   herbaceous : 0.473 
##   cultivated_crops : 0.318 
##   deciduous_forest : 0.207 
## 
## ------------------------
</code></pre></div><br>

```
corrplot::corrplot(corr.matrix)
```

<img class="sdm-eda-img" src="{{ site.baseurl }}/assets/plots/sdm-eda-15.png">

## Feature Engineering

As described previously, feature engineering is an essential step in predictive modeling, as it can help to reduce possible problems due to multi-collinearity or spatial autocorrelation. In addition, it can help reduce overfitting, as well as provide new insights unavailable through basic pre-processing of the data.
Setup

First, create a space to output any “engineered” rasters:

```
output.dir <- "artifacts/feature_engineering"
if (!dir.exists(output.dir)) dir.create(output.dir)
```

### Hierarchical Categories for Land Cover

<div>
    <p>
        Recall the <a href="https://www.mrlc.gov/data/legends/national-land-cover-database-class-legend-and-description" target="_blank">hierarchical categories for the land cover data</a>:
    </p>
    <ul>
        <li>Water:
            <ul>
                <li>11: Open Water</li>
                <li>12: Perennial Ice/Snow</li>
            </ul>
        </li>
        <li>Developed
            <ul>
                <li>21: Developed, Open Space</li>
                <li>22: Developed, Low Intensity</li>
                <li>23: Developed, Medium Intensity</li>
                <li>24: Developed, High Intensity</li>
            </ul>
        </li>
        <li>Barren
            <ul>
                <li>31: Barren Land (Rock/Sand/Clay)</li>
            </ul>
        </li>
        <li>Forest
            <ul>
                <li>41: Deciduous Forest</li>
                <li>42: Evergreen Forest</li>
                <li>43: Mixed Forest</li>
            </ul>
        </li>
        <li>Shrubland
            <ul>
                <li>51: Dwarf Shrub</li>
                <li>52: Shurb/Scrub</li>
            </ul>
        </li>
        <li>Herbaceous
            <ul>
                <li>71: Grassland/Herbaceous</li>
                <li>72: Sedge/Herbaceous</li>
                <li>73: Lichens</li>
                <li>74: Moss</li>
            </ul>
        </li>
        <li>Planted/Cultivated
            <ul>
                <li>81: Pasture/Hay</li>
                <li>82: Cultivated Crops</li>
            </ul>
        </li>
        <li>Wetlands
            <ul>
                <li>90: Woody Wetlands</li>
                <li>95: Emergent Herbaceous Wetlands</li>
            </ul>
        </li>
    </ul>
</div>

Using these categories, feature reduction can be accomplished using a heuristic methodology. Other redundant rasters, such as the Waterbody and Urban Imperviousness rasters, can also be combined with the respective related land cover rasters to further reduce multicollinearity between rasters.

```
create.parent.category.rasters <- function(input.raster, 
                                           state, 
                                           output.dir,
                                           layer.name="NLCD_Land") {
  # Define the category mappings
  categories <- list(
    water = c("Open Water"),
    ice.snow = c("Perennial Ice/Snow"),
    developed = c("Developed, Open Space",
                  "Developed, Low Intensity",
                  "Developed, Medium Intensity",  
                  "Developed, High Intensity"),
    barren = c("Barren Land"),
    forest = c("Mixed Forest", 
               "Deciduous Forest", 
               "Evergreen Forest"),
    shrubland = c("Shrub/Scrub", "Dwarf Shrub"),
    herbaceous = c("Herbaceous", "Grassland/Herbaceous",
                   "Sedge/Herbaceous", "Lichens", "Moss"),
    planted.cultivated = c("Cultivated Crops", "Hay/Pasture"),
    wetlands = c("Woody Wetlands", "Emergent Herbaceous Wetlands")
  )
  
  # Iterate through each category to create and save the binary raster
  for (cat.name in names(categories)) {
    # Define the output file path
    output.file <- file.path(output.dir, paste0(state, "_", cat.name, ".tif"))
    
    if (!file.exists(output.file)) {
      cat("Generating raster for", cat.name, "in", state, "\n")
      if (cat.name != "developed") {
        # Create a binary raster for the current category
        out.raster <- terra::lapp(input.raster[[layer.name]], 
                                  fun = function(x) {
                                    case_when(
                                      is.na(x) ~ NA,
                                      x %in% categories[[cat.name]] ~ 1.,
                                      T ~ 0)
                                  })
        names(out.raster) <- cat.name
        if (cat.name == "water") {
          # Combine with waterbody raster layer
          wb.raster <- input.raster[[paste0("waterbody_", state)]]
          out.raster <- (out.raster + wb.raster) / 2
          names(out.raster) <- "water"
        }
      } else {
        # Set developed raster to be a value, scale ranging from 0.25-1 by 0.25
        out.raster <- terra::lapp(input.raster[[layer.name]], 
                                  fun = function(x) {
                                    case_when(
                                      is.na(x) ~ NA,
                                      x == "Developed, Open Space" ~ 0.25,
                                      x == "Developed, Low Intensity" ~ 0.5,
                                      x == "Developed, Medium Intensity" ~ 0.75,
                                      x == "Developed, High Intensity" ~ 1.,
                                      T ~ 0)
                                  })
        # Combine with urban imperviousness
        ui.raster <- input.raster[[paste0("urban_imperviousness_", state)]]
        ui.min <- terra::minmax(ui.raster)[[1]]
        ui.max <- terra::minmax(ui.raster)[[2]]
        ui.raster.scaled <- (ui.raster - ui.min) / (ui.max - ui.min)
        out.raster <- (out.raster + ui.raster.scaled) / 2
        names(out.raster) <- "developed"
      }
      
      # Save the output raster to the specified output path
      writeRaster(out.raster, output.file, overwrite=T)
    }
  }
}

for (state in states) {
  create.parent.category.rasters(r.list[[state]], state, output.dir)
}
```

### Aggregate Seasonal NDVI Rasters (Min/Max)

NDVI, or Normalized Difference Vegetation Index, is a common measure of vegetation health. For each state, there are four separate NDVI rasters representing different seasons. Aggregating these rasters provides two essential insights: the minimum and maximum NDVI values. The minimum NDVI may correspond to the period when vegetation is least vigorous (e.g., winter or drought), while the maximum may indicate the period of peak vegetation growth (e.g., spring or summer). By scaling these values between 0 and 1, and then storing them separately as minimum and maximum NDVI rasters, a standardized measure of vegetation dynamics throughout the year for each state is achieved.

```
# Convert 4 separate NDVI rasters to a single raster
for (state in states) {
  output.file.min <- file.path(output.dir, paste0(state, "_min_NDVI.tif"))
  output.file.max <- file.path(output.dir, paste0(state, "_max_NDVI.tif"))
  if (!file.exists(output.file.min) | !file.exists(output.file.max)) {
    r.ndvi <- r.list[[state]][[names(r.list[[state]]) %like% "NDVI"]]
    # Parse seasons
    for (s in names(r.ndvi)) {
      r <- r.ndvi[[s]]
      .min <- terra::minmax(r)[[1]]
      .max <- terra::minmax(r)[[2]]
      # Scale to be from 0 to 1
      r.ndvi[[s]] <- (r - .min) / (.max - .min)
    }
    # Get min/max scaled NDVI values
    r.max <- max(r.ndvi)
    names(r.max) <- "max.ndvi"
    r.min <- min(r.ndvi)
    names(r.min) <- "min.ndvi"
    writeRaster(r.max, output.file.max, overwrite=T)
    writeRaster(r.min, output.file.min, overwrite=T)
  }
}
```

### Combine Temperature Rasters (Difference)

In this section, the difference between the maximum and minimum temperatures for each state is calculated. This temperature difference, often referred to as diurnal temperature range, can provide insights into the variability of daily temperatures. A high difference may suggest drastic temperature changes from day to night, which can impact ecosystems and human health. By creating a scaled raster of these differences (ranging from 0 to 1), the temperature variability can be compared across different states.

```
# Get the difference between the min and max temperatures as a raster
for (state in states) {
  output.file <- file.path(output.dir, paste0(state, "_tdiff.tif"))
  if (!file.exists(output.file)) {
    # Get difference
    r.tdiff <- r.list[[state]][[paste0("tmax_", state)]] - 
      r.list[[state]][[paste0("tmin_", state)]]
    # Get min/max differences
    .min <- terra::minmax(r.tdiff)[[1]]
    .max <- terra::minmax(r.tdiff)[[2]]
    # Scale to be from 0 to 1
    r.tdiff <- (r.tdiff - .min) / (.max - .min)
    names(r.tdiff) <- "temp.diff"
    # Write raster
    writeRaster(r.tdiff, output.file, overwrite=T)
  }
}
```

### Spatial Filters

These are components derived from the spatial coordinates, which can capture and account for spatial structures in the data. E.g., a polynomial term based on latitude and longitude could account for some of the spatial trend.

```
for (state in states) {
  output.file.lon <- file.path(output.dir, paste0(state, "_lon.tif"))
  output.file.lat <- file.path(output.dir, paste0(state, "_lat.tif"))
  output.file.lon2 <- file.path(output.dir, paste0(state, "_lon_poly.tif"))
  output.file.lat2 <- file.path(output.dir, paste0(state, "_lat_poly.tif"))
  output.file.lli <- file.path(output.dir, paste0(state, "_lon_lat_interaction.tif"))
  if (!all(file.exists(c(output.file.lon, output.file.lat, 
                         output.file.lon2, output.file.lat2,
                         output.file.lli)))) {
    # Get raster as df
    r.df <- terra::as.data.frame(r.list[[state]], xy=T, cells=T) %>% 
      rename(lon=x, lat=y) %>%
      select(lon, lat, cell) %>%
      st_as_sf(coords = c("lon", "lat"), crs = st_crs(r.list[[state]])) %>%
      st_transform(crs=4326) %>%
      cbind(st_coordinates(.)) %>%
      rename(lon="X", lat="Y")
    
    # Initialize empty raster templates
    r.lat <- ext(r.list[[state]]) %>% 
            rast(res=res(r.list[[state]]), crs=crs(r.list[[state]]))
    r.lon <- ext(r.list[[state]]) %>% 
            rast(res=res(r.list[[state]]), crs=crs(r.list[[state]]))
    # Fill with lat/lon values
    r.lat[r.df$cell] <- r.df$lat
    names(r.lat) <- "lat"
    r.lon[r.df$cell] <- r.df$lon
    names(r.lon) <- "lon"
    # Get polynomial and interaction terms
    r.lat2 <- r.lat^2 
    names(r.lat2) <- "lat.sqrd"
    r.lon2 <- r.lon^2 
    names(r.lon2) <- "lon.sqrd"
    r.lon.lat <- r.lon * r.lat
    names(r.lon.lat) <- "lon.lat.interaction"
    
    # Write outputs
    writeRaster(r.lat, output.file.lat, overwrite=T)
    writeRaster(r.lon, output.file.lon, overwrite=T)
    writeRaster(r.lat2, output.file.lat2, overwrite=T)
    writeRaster(r.lon2, output.file.lon2, overwrite=T)
    writeRaster(r.lon.lat, output.file.lli, overwrite=T)
  }
}
```

### De-trend DEM

Digital Elevation Models (DEM) give us a detailed overview of the terrain elevation. However, sometimes, larger spatial trends (e.g., the general increase in elevation from the coast to the mountains) can overshadow more localized features or variations. By de-trending the DEM, the underlying larger spatial trend is removed, highlighting only the local elevation differences. In the given section, a linear model is fitted to the DEM using latitude and longitude (and their polynomial terms) to capture the broader spatial trends. The residuals from this model represent the local variations in elevation, devoid of the overarching trend. These de-trended values are then scaled between 0 and 1 to create a standardized representation of the localized elevation changes for each state.

```
for (state in states) {
  output.file <- file.path(output.dir, paste0(state, "_detrend_dem.tif"))
  if (!file.exists(output.file)) {
    # Get raster as df
    r.df <- terra::as.data.frame(r.list[[state]], xy=T, cells=T) %>% 
      rename(lon=x, lat=y) %>%
      select(lon, lat, cell, !!sym(paste0("dem_", state))) %>%
      st_as_sf(coords = c("lon", "lat"), crs = st_crs(r.list[[state]])) %>%
      st_transform(crs=4326) %>%
      cbind(st_coordinates(.)) %>%
      rename(lon="X", lat="Y", dem=paste0("dem_", state)) %>%
      # Convert to data.table
      setDT()
    
    # Fit model based on location, with dem as response
    fit <- lm(dem ~ lat * lon + I(lat^2) + I(lon^2), 
              data=r.df[!is.na(dem)])
    # Get residuals from the model as "de-trended" dem values
    r.df[!is.na(dem), dem.detrended := residuals(fit)]
    # Scale de-trended values
    r.df[, dem.detrended := (dem.detrended - min(dem.detrended, na.rm=T)) / 
           (max(dem.detrended, na.rm=T) - min(dem.detrended, na.rm=T))]
    # Initialize empty raster template
    r.dem <- ext(r.list[[state]]) %>% 
            rast(res=res(r.list[[state]]), crs=crs(r.list[[state]]))
    # Fill with de-trended dem values
    r.dem[r.df$cell] <- r.df$dem.detrended
    names(r.dem) <- "dem.detrended"
    
    # Write outputs
    writeRaster(r.dem, output.file, overwrite=T)
  }
}
```

## Feature Selection

Feature selection is the process of determining the most relevant input variables to use in modeling or analysis. By selecting the most informative features, one can improve the model’s accuracy, reduce overfitting, reduce issues due to multi-collinearity and spatial autocorrelation, and decrease the computational cost of training.
Load all Original Features and Feature Engineering Outputs

In this section, the final dataset that will be used for feature selection is prepared. This involves:

- Loading the pseudo-absence and observation datasets.
- Fetching the feature-engineered raster data (geospatial data layers) and associating these rasters with the observational data.
- Loading the original raster layers, as these should still be included during the feature selection process.
- Converting the original land cover feature into binary variables.
- Saving the final aggregated dataset.

```
final.output.dir <- "artifacts/final_data"
if (!dir.exists(final.output.dir)) dir.create(final.output.dir)
final.fpath <- file.path(final.output.dir, "final_data.rds")

if (!file.exists(final.fpath)) {
  # Combine observation/pseudo-absence data
  df <- c(list.files(file.path("data", "final_pseudo_absence"), full.names = T),
          list.files(file.path("data", "final"), full.names = T)) %>%
    purrr::map_df(readRDS)
  
  # Get feature engineered raster data
  output.dir <- "artifacts/feature_engineered_final"
  if (!dir.exists(output.dir)) dir.create(output.dir)
  get.fe <- all(case_when(is_empty(list.files(output.dir)) ~ T,
                          list.files(output.dir) < length(states) ~ T,
                          T ~ F))
  fe.r.list <- list()
  if (get.fe) {
    for (state in states) {
      fe.r.list[[state]] <- list.files(
        "artifacts/feature_engineering", full.names = T)[
          grepl(
            paste0("^", state, "_"), 
            list.files("artifacts/feature_engineering", full.names = F)
          )
        ] %>% rast()
      # Save as multi-layer rasters by state
      writeRaster(fe.r.list[[state]], 
                  file.path(output.dir, paste0(state, ".tif")),
                  overwrite=T)
    }
  } else {
    for (state in states) {
      fe.r.list[[state]] <- rast(file.path(output.dir, paste0(state, ".tif")))
    }
  }
  
  # Add empty fields in dataframe for each new raster 
  purrr::map(fe.r.list, names) %>% 
    reduce(c) %>% 
    unique() %>% 
    sort() %>%
    purrr::walk(function(n) {
      if (!(n %in% names(df))) df[[n]] <<- 0
    })
  
  # Update point crs to match fe rasters
  df <- st_transform(df, st_crs(fe.r.list[[1]])) 
  df.list <- list()
  # From state multi-layer rasters, extract values from each point in df
  for (st in states) {
    cat(sprintf("Extracting points to values for %s...\n", st))
    r <- fe.r.list[[st]]
    df.list[[st]] <- df %>% filter(state == st)
    # Extract raster values
    for (r.name in names(r)) {
      cat("\tExtracting", r.name, "values for", st, "\n")
      x <- terra::extract(r[[r.name]], df.list[[st]])[[r.name]]
      df.list[[st]] <- mutate(df.list[[st]], !!r.name := x)
    }
  }
  # Combine list
  df <- do.call("rbind", df.list)
  rm(df.list)
  gc()
  # Update crs back
  df <- st_transform(df, 4326)
  # Remove rownames
  rownames(df) <- NULL
  
  # Convert land cover to binary variables
  binary.lc <- caret::dummyVars(~NLCD_Land, data=df, sep="") %>% 
    predict(., df) %>%
    as.data.frame()
  
  names(binary.lc) <- gsub("NLCD_Land", "", names(binary.lc)) %>% 
    gsub("\\/| |,", "_", .) %>% 
    gsub("__", "_", .) %>% 
    tolower()
  
  # Remove duplicates from feature engineering
  binary.lc <- binary.lc %>%
    select(-c("unclassified", "perennial_snow_ice", "barren_land",
              "shrub_scrub", "herbaceous"))
  # Combine with dataframe, remove land cover categorical var
  df <- df %>% select(-NLCD_Land) %>% cbind(binary.lc)
  saveRDS(df, final.fpath)
} else {
  df <- readRDS(final.fpath)
}

# Convert to data.table
df %>% setDT()

# View output
df %>% as_tibble()
```

<div data-pagedtable="false" pagedtable-page="0" class="pagedtable-wrapper">
  <script src="{{ site.baseurl }}/assets/data/sdm-2-eda-fs.js"></script>
  <script data-pagedtable-source data-global="sdm2FS" type="application/json"></script>
</div>

### Correlation

```
In this section, correlation between the predictors and the dependent variable (observations) for each species are examined. The correlation measure assists in identifying which predictors have a linear relationship with the response variable and can, therefore, impact model predictions.

obs.cor <- purrr::map(sort(unique(df$common.name)), function(spec) {
  corr.matrix <- cor(dplyr::select(df[common.name == spec], 
                            -c("state", "common.name", "geometry")))
  obs.cor <- corr.matrix[which(rownames(corr.matrix) == "observations"),] %>%
    as.data.frame()
  obs.cor$variable <- rownames(obs.cor)
  obs.cor %>%
    filter(variable != "observations") %>%
    rename(!!spec := `.`) %>%
    arrange(-abs(!!sym(spec))) %>%
    mutate(!!spec := paste0(variable, ": ", signif(!!sym(spec), 2))) %>%
    dplyr::select(!!sym(spec)) 
}) %>% 
  do.call("cbind", .)

obs.cor %>% as_tibble()
```

<div data-pagedtable="false" pagedtable-page="0" class="pagedtable-wrapper">
<script data-pagedtable-source="" type="application/json">
{"columns":[{"label":["Belted Kingfisher"],"name":[1],"type":["chr"],"align":["left"]},{"label":["Cedar Waxwing"],"name":[2],"type":["chr"],"align":["left"]},{"label":["Downy Woodpecker"],"name":[3],"type":["chr"],"align":["left"]},{"label":["Ruddy Duck"],"name":[4],"type":["chr"],"align":["left"]},{"label":["Sanderling"],"name":[5],"type":["chr"],"align":["left"]},{"label":["Sandhill Crane"],"name":[6],"type":["chr"],"align":["left"]},{"label":["Sharp-shinned Hawk"],"name":[7],"type":["chr"],"align":["left"]},{"label":["Wild Turkey"],"name":[8],"type":["chr"],"align":["left"]}],"data":[{"1":"urban_imperviousness: 0.082","2":"urban_imperviousness: 0.073","3":"urban_imperviousness: 0.059","4":"open_water: 0.076","5":"coastline: 0.2","6":"woody_wetlands: 0.054","7":"mixed_forest: 0.044","8":"dem.detrended: -0.073"},{"1":"dem.detrended: -0.071","2":"developed: 0.071","3":"developed: 0.058","4":"water: 0.071","5":"urban_imperviousness: 0.17","6":"wetlands: 0.044","7":"tmin: 0.016","8":"forest: -0.036"},{"1":"developed: 0.069","2":"dem: -0.055","3":"dem.detrended: -0.037","4":"waterbody: 0.044","5":"max.ndvi: -0.13","6":"waterbody: 0.043","7":"Fall_NDVI: 0.016","8":"developed: 0.035"},{"1":"water: 0.066","2":"dem.detrended: -0.044","3":"open_water: 0.034","4":"dem.detrended: -0.04","5":"Fall_NDVI: -0.13","6":"water: 0.033","7":"forest: 0.016","8":"avg_prcp: -0.035"},{"1":"open_water: 0.061","2":"evergreen_forest: -0.038","3":"developed_open_space: 0.034","4":"forest: -0.038","5":"Winter_NDVI: -0.13","6":"lat: -0.02","7":"Winter_NDVI: 0.015","8":"tmax: -0.035"},{"1":"waterbody: 0.049","2":"open_water: 0.035","3":"shrubland: -0.033","4":"evergreen_forest: -0.031","5":"Summer_NDVI: -0.12","6":"lat.sqrd: -0.019","7":"min.ndvi: 0.014","8":"tmin: -0.034"},{"1":"emergent_herbaceous_wetlands: 0.049","2":"shrubland: -0.035","3":"water: 0.031","4":"hay_pasture: 0.028","5":"tmax: 0.11","6":"tmin: -0.018","7":"avg_prcp: 0.014","8":"hay_pasture: 0.029"},{"1":"forest: -0.043","2":"developed_medium_intensity: 0.034","3":"dem: -0.03","4":"max.ndvi: -0.026","5":"developed: 0.11","6":"Winter_NDVI: -0.018","7":"max.ndvi: 0.014","8":"urban_imperviousness: 0.028"},{"1":"evergreen_forest: -0.038","2":"water: 0.032","3":"evergreen_forest: -0.029","4":"urban_imperviousness: 0.025","5":"avg_prcp: 0.1","6":"forest: -0.018","7":"temp.diff: 0.013","8":"developed_high_intensity: 0.028"},{"1":"Winter_NDVI: -0.037","2":"forest: -0.03","3":"lon.sqrd: -0.025","4":"Fall_NDVI: -0.025","5":"tmin: 0.088","6":"avg_prcp: -0.018","7":"dem: -0.013","8":"Winter_NDVI: -0.025"},{"1":"shrubland: -0.035","2":"Summer_NDVI: 0.029","3":"lon: 0.025","4":"shrubland: -0.024","5":"dem: -0.084","6":"tmax: -0.016","7":"lon: 0.013","8":"deciduous_forest: -0.024"},{"1":"wetlands: 0.031","2":"developed_low_intensity: 0.027","3":"emergent_herbaceous_wetlands: 0.025","4":"developed: 0.022","5":"min.ndvi: -0.082","6":"dem: 0.015","7":"Summer_NDVI: 0.012","8":"open_water: 0.024"},{"1":"max.ndvi: -0.027","2":"developed_open_space: 0.027","3":"developed_medium_intensity: 0.025","4":"developed_medium_intensity: 0.019","5":"Spring_NDVI: -0.079","6":"dem.detrended: -0.014","7":"lat: 0.012","8":"lon: -0.024"},{"1":"developed_medium_intensity: 0.026","2":"emergent_herbaceous_wetlands: 0.023","3":"cultivated_crops: -0.024","4":"Winter_NDVI: -0.019","5":"lon.sqrd: 0.079","6":"cultivated_crops: 0.014","7":"lat.sqrd: 0.012","8":"lon.sqrd: 0.023"},{"1":"developed_low_intensity: 0.024","2":"developed_high_intensity: 0.022","3":"developed_low_intensity: 0.023","4":"lon: -0.018","5":"lon.lat.interaction: -0.077","6":"evergreen_forest: -0.014","7":"lon.sqrd: -0.011","8":"lon.lat.interaction: -0.022"},{"1":"dem: -0.024","2":"herbaceous: -0.021","3":"Summer_NDVI: 0.023","4":"lon.sqrd: 0.017","5":"barren: 0.075","6":"emergent_herbaceous_wetlands: 0.014","7":"tmax: 0.011","8":"temp.diff: -0.022"},{"1":"developed_high_intensity: 0.022","2":"waterbody: 0.021","3":"lon.lat.interaction: 0.02","4":"temp.diff: -0.017","5":"developed_low_intensity: 0.074","6":"shrubland: -0.013","7":"herbaceous: -0.011","8":"max.ndvi: -0.022"},{"1":"min.ndvi: -0.019","2":"tmin: 0.016","3":"waterbody: 0.02","4":"herbaceous: -0.015","5":"lon: -0.074","6":"max.ndvi: -0.012","7":"shrubland: -0.01","8":"developed_medium_intensity: 0.019"},{"1":"Spring_NDVI: -0.018","2":"cultivated_crops: -0.016","3":"Fall_NDVI: 0.018","4":"mixed_forest: -0.015","5":"emergent_herbaceous_wetlands: 0.067","6":"planted.cultivated: 0.012","7":"cultivated_crops: -0.0079","8":"min.ndvi: -0.016"},{"1":"cultivated_crops: -0.018","2":"tmax: 0.015","3":"herbaceous: -0.018","4":"Summer_NDVI: -0.015","5":"forest: -0.055","6":"Summer_NDVI: -0.011","7":"open_water: 0.0076","8":"developed_open_space: 0.015"},{"1":"coastline: 0.015","2":"Spring_NDVI: 0.015","3":"forest: -0.017","4":"lon.lat.interaction: -0.015","5":"wetlands: 0.048","6":"herbaceous: -0.0099","7":"water: 0.0073","8":"planted.cultivated: 0.015"},{"1":"Fall_NDVI: -0.013","2":"hay_pasture: 0.015","3":"tmin: 0.015","4":"tmin: -0.014","5":"lat.sqrd: 0.039","6":"min.ndvi: -0.0093","7":"dem.detrended: -0.0061","8":"mixed_forest: -0.015"},{"1":"temp.diff: -0.013","2":"max.ndvi: 0.015","3":"planted.cultivated: -0.014","4":"avg_prcp: -0.014","5":"lat: 0.037","6":"Spring_NDVI: -0.0086","7":"lon.lat.interaction: 0.0059","8":"evergreen_forest: -0.014"},{"1":"developed_open_space: 0.013","2":"avg_prcp: 0.015","3":"developed_high_intensity: 0.013","4":"dem: -0.013","5":"deciduous_forest: -0.037","6":"temp.diff: -0.0084","7":"deciduous_forest: -0.0055","8":"Fall_NDVI: -0.014"},{"1":"Summer_NDVI: -0.012","2":"min.ndvi: 0.013","3":"avg_prcp: 0.011","4":"wetlands: 0.013","5":"planted.cultivated: -0.034","6":"lon.lat.interaction: 0.0071","7":"planted.cultivated: -0.0055","8":"herbaceous: 0.013"},{"1":"mixed_forest: -0.012","2":"wetlands: 0.012","3":"max.ndvi: 0.011","4":"cultivated_crops: -0.013","5":"dem.detrended: -0.032","6":"deciduous_forest: -0.0069","7":"urban_imperviousness: 0.0053","8":"water: 0.012"},{"1":"deciduous_forest: -0.012","2":"Fall_NDVI: 0.011","3":"Winter_NDVI: -0.0091","4":"woody_wetlands: 0.013","5":"waterbody: 0.031","6":"mixed_forest: -0.0055","7":"Spring_NDVI: 0.0052","8":"emergent_herbaceous_wetlands: 0.011"},{"1":"planted.cultivated: -0.011","2":"deciduous_forest: 0.0063","3":"woody_wetlands: -0.0089","4":"deciduous_forest: -0.012","5":"cultivated_crops: -0.026","6":"coastline: -0.0037","7":"waterbody: 0.0049","8":"shrubland: -0.011"},{"1":"herbaceous: -0.0071","2":"lat.sqrd: 0.006","3":"tmax: 0.0086","4":"developed_open_space: -0.011","5":"temp.diff: 0.025","6":"developed_open_space: -0.0036","7":"emergent_herbaceous_wetlands: 0.004","8":"lat: 0.011"},{"1":"avg_prcp: -0.0058","2":"lat: 0.0052","3":"hay_pasture: 0.0075","4":"tmax: -0.011","5":"mixed_forest: -0.024","6":"developed: -0.0031","7":"developed: 0.0034","8":"cultivated_crops: -0.011"},{"1":"tmax: -0.0053","2":"lon: 0.0052","3":"temp.diff: 0.0066","4":"developed_low_intensity: 0.01","5":"developed_medium_intensity: 0.024","6":"developed_low_intensity: -0.0026","7":"woody_wetlands: -0.0026","8":"lat.sqrd: 0.0095"},{"1":"lon.sqrd: -0.0044","2":"temp.diff: 0.0046","3":"wetlands: 0.0064","4":"planted.cultivated: 0.0088","5":"evergreen_forest: -0.022","6":"urban_imperviousness: -0.0025","7":"coastline: 0.0026","8":"Summer_NDVI: -0.0042"},{"1":"barren: -0.0038","2":"Winter_NDVI: -0.0042","3":"deciduous_forest: 0.0055","4":"min.ndvi: -0.0083","5":"water: 0.022","6":"Fall_NDVI: -0.0024","7":"evergreen_forest: -0.002","8":"Spring_NDVI: -0.0039"},{"1":"tmin: -0.0038","2":"lon.sqrd: -0.0041","3":"Spring_NDVI: -0.0054","4":"emergent_herbaceous_wetlands: 0.0059","5":"hay_pasture: -0.021","6":"open_water: 0.0024","7":"developed_medium_intensity: 0.0016","8":"woody_wetlands: -0.0034"},{"1":"hay_pasture: 0.0035","2":"barren: -0.0038","3":"barren: -0.0051","4":"lat: 0.0055","5":"developed_open_space: -0.015","6":"developed_high_intensity: 0.0019","7":"barren: -0.001","8":"wetlands: 0.003"},{"1":"lon: 0.0034","2":"mixed_forest: -0.0037","3":"mixed_forest: 0.0049","4":"lat.sqrd: 0.0052","5":"herbaceous: -0.0089","6":"barren: -0.0019","7":"developed_low_intensity: 0.00087","8":"barren: -0.0027"},{"1":"lon.lat.interaction: 0.0033","2":"coastline: 0.003","3":"coastline: -0.0042","4":"barren: -0.0047","5":"shrubland: -0.0076","6":"developed_medium_intensity: -0.0015","7":"hay_pasture: 0.00063","8":"developed_low_intensity: 0.0014"},{"1":"woody_wetlands: -0.0019","2":"woody_wetlands: -0.0029","3":"lat: -0.0032","4":"Spring_NDVI: -0.0037","5":"woody_wetlands: -0.0052","6":"hay_pasture: 0.0012","7":"wetlands: 0.00054","8":"coastline: -0.0014"},{"1":"lat.sqrd: -0.0016","2":"lon.lat.interaction: 0.0012","3":"lat.sqrd: -0.0032","4":"coastline: 0.0019","5":"open_water: 0.0013","6":"lon: -0.0011","7":"developed_high_intensity: -0.00026","8":"waterbody: 0.00066"},{"1":"lat: -0.0012","2":"planted.cultivated: 0.00072","3":"min.ndvi: 0.0026","4":"developed_high_intensity: 8e-04","5":"developed_high_intensity: 0.00016","6":"lon.sqrd: -0.00057","7":"developed_open_space: -3.3e-08","8":"dem: -0.00012"},{"1":"ice.snow: NA","2":"ice.snow: NA","3":"ice.snow: NA","4":"ice.snow: NA","5":"ice.snow: NA","6":"ice.snow: NA","7":"ice.snow: NA","8":"ice.snow: NA"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
</script>
</div>

### Train/Test Split

Ensuring the quality and robustness of a model requires that it is evaluated on unseen data. Similarly, it is important to avoid “data leakage”, which occurs when information from the test set inadvertently influences the training process. Reasons for data leakage can range from pre-processing steps being applied to the entire dataset instead of separately to train and test sets, to more subtle causes like target leakage where future information inadvertently gets used in the past due to data preparation mistakes.

In this section, the data is split into a training set and a test set. The training set is used for training the model, while the test set is reserved to evaluate the model’s performance. The data is stratified based on latitude, longitude, species, state, and absence to ensure that the distribution of these variables is consistent between the train and test sets. This stratification helps maintain representative samples and mitigates potential biases.

```
# Set seed for splitting and modeling
set.seed(19)

stratified.split.idx <- function(df, p=0.7, lat.lon.bins=25) {
  # Cut along lat/lon values to create grids (lat.bin & lon.bin)
  # lat.lon.bins is the number of divisions you want
  df$lat.bin <- cut(df$lat, breaks=lat.lon.bins, labels = F)
  df$lon.bin <- cut(df$lon, breaks=lat.lon.bins, labels = F)
  df$absence <- df$observations == 0
  
  # Create a new variable combining the stratification variables
  df %>%
    # stratify on lat/lon bins, species, state, and absence
    mutate(strata = paste(lat.bin, lon.bin, common.name, state, absence)) %>%
    pull(strata) %>%
    # Create the data partitions
    createDataPartition(., p = p, list = F) %>%
    suppressWarnings()
}

prepare.data <- function(df, p=.7, lat.lon.bins=25) {
  train.index <- stratified.split.idx(df, p=p, lat.lon.bins = lat.lon.bins)
  df.train <- df[train.index, ]
  df.test <- df[-train.index, ]
  
  list(train = df.train, 
       test = df.test,
       index = train.index)
}

# Get train/test indices
train.test <- prepare.data(df, .7)

# Split; Remove non-predictive variables
df.train <- df[train.test$index,] 
df.train[, `:=` (geometry=NULL)]
df.test <- df[-train.test$index,]
df.test[, `:=` (geometry=NULL)]

train.x <- df.train %>% dplyr::select(-observations)
train.y <- df.train$observations

test.x <- df.test %>% dplyr::select(-observations)
test.y <- df.test$observations
```

For computational efficiency, models are cached. This section provides a function to retrieve the model from cache if it exists or save it to the cache if it doesn’t. This approach speeds up the modeling process, especially when iterating and fine-tuning various models, by avoiding retraining on the same dataset.

```
# Get model from cache if it has been saved before
get.model <- function(model, file.name, model.path) {
  f.path <- file.path(model.path, file.name)
  if (!dir.exists(model.path)) {
    dir.create(model.path)
  }
  # Model cache check
  if (ifelse(file.exists(f.path), T, F)) {
    model <- readRDS(f.path)
  } else {
    saveRDS(model, f.path)
  }
  model
}
```

### LASSO Models - Importance from Coefficients

LASSO (Least Absolute Shrinkage and Selection Operator) regression is a type of linear regression that uses shrinkage. Here, data values are shrunk towards a central point, like the mean. The lasso procedure encourages simple, sparse models (i.e., models with fewer parameters).

In this section, a LASSO regression model is fitted for each species in each state. The regularization strength is controlled by the lambda parameter, where a value of zero results in regular linear regression and increasing values result in more regularization. By fitting LASSO models and examining the coefficients, we can determine the importance of each variable and variable interaction. Variables/interactions with non-zero coefficients are considered important, while those with coefficients shrunk to zero are considered non-informative.

```
species <- sort(unique(df.train$common.name))
if (!dir.exists("artifacts/models")) dir.create("artifacts/models")

# Define the control for the train function
ctrl <- trainControl(method = "cv", number = 5)

lasso.list <- purrr::map(species, function(spec) {
  spec.state.fit <- purrr::map(states, function(st) {
    cat("Fitting LASSO model for", spec, "in", st, "\n")
    spec.df <- df.train[common.name == spec & state == st][
      , `:=` (state=NULL, common.name = NULL)]
    
    # Remove any columns where all values are the same
    .remove <- names(which(sapply(spec.df, function(x) length(unique(x)) <= 1)))
    .remove <- .remove[.remove != "observations"]
    if (!is_empty(.remove)) {
      spec.df <- spec.df %>% dplyr::select(-.remove)
    }
    
    # Scale data
    # Identify fields to center/scale
    to.scale <- sapply(spec.df, function(x) is.numeric(x) && 
                         length(unique(x)) > 2)
    to.scale$observations <- F
    to.scale <- names(spec.df) %in% names(which(unlist(to.scale)))
    
    # Define the pre-processing steps (use the training data to avoid data leakage)
    # Apply centering and scaling to the non-binary fields and non-target
    preproc <- preProcess(spec.df[, ..to.scale], 
                          method = c("center", "scale"))
    
    # Center/Scale the training data
    df.train.scaled <- bind_cols(spec.df[, !(..to.scale)],
                                 predict(preproc, spec.df[, ..to.scale]))
    
    # LASSO (L1); Elastic Net, where alpha = 1
    fit.l1 <- get.model(
      train(observations ~ (.)^2, 
            data = df.train.scaled, 
            method = "glmnet",
            trControl = ctrl,
            tuneGrid = expand.grid(alpha = 1, 
                                   lambda = seq(0, 1, by = 0.1)),
            metric = "RMSE"),
      file.name=paste0(tolower(gsub(" ", "_", spec)), "_", st, "_regr_l1.rds"),
      model.path="artifacts/models/lasso_fs")
    
    gc()
    fit.l1
  })
  names(spec.state.fit) <- states
  spec.state.fit
})
names(lasso.list) <- species

# Get variable importance for LASSO models, based on coefficients
spec.state <- expand.grid(common.name=species, state=states, stringsAsFactors=F)
lasso.imp <- purrr::map_df(1:nrow(spec.state), function(i) {
  spec <- spec.state[i,]$common.name
  st <- spec.state[i,]$state
  fit <- lasso.list[[spec]][[st]]
  coef.df <- coef(fit$finalModel, s = fit$bestTune$lambda) %>%
    as.matrix() %>%
    as.data.frame()
  # Remove the intercept
  coef.df <- coef.df[-1, , drop = F]
  
  # Create a data frame of variable importance
  var.importance <- tibble(
    common.name = spec,
    state = st,
    variable = rownames(coef.df),
    importance = abs(coef.df[,1])
  ) %>%
    # Rank variables by importance
    arrange(state, common.name, -importance, variable) %>%
    # Only look at variables where imp. > 0
    filter(importance > 0)
})

lasso.imp
```

<div data-pagedtable="false" pagedtable-page="0" class="pagedtable-wrapper">
  <script src="{{ site.baseurl }}/assets/data/sdm-2-eda-lasso.js"></script>
  <script data-pagedtable-source data-global="sdm2Lasso" type="application/json"></script>
</div>

### Decision Tree Models - Variable Importance

Decision trees are a class of machine learning models that recursively split the dataset into subsets based on the value of input variables. The goal is to make the data within each subset as homogeneous as possible regarding the target variable. The importance of a variable in a decision tree model can be assessed by the frequency and depth at which it is used to split the data. Variables used frequently and closer to the root of the tree are typically considered more important.

The primary advantage of using a tree-based model as opposed to a model like LASSO is that the relationships do not need to be linear, and interactions are naturally accounted for due to the hierarchical nature of the model itself.

In the code below, Decision Tree models are trained for each species in each state. The variable importance is extracted from the tree structure itself, which offers an intuitive understanding of the relationships and hierarchies within the data.

```
# Decision Tree Variable Importance
tree.list <- purrr::map(species, function(spec) {
  spec.state.fit <- purrr::map(states, function(st) {
    cat("Fitting Decision Tree model for", spec, "in", st, "\n")
    spec.df <- df.train[common.name == spec & state == st][
      , `:=` (state=NULL, common.name = NULL)]
    
    # Remove any columns where all values are the same
    .remove <- names(which(sapply(spec.df, function(x) length(unique(x)) <= 1)))
    .remove <- .remove[.remove != "observations"]
    if (!is_empty(.remove)) {
      spec.df <- spec.df %>% dplyr::select(-.remove)
    }

    # Fit the Decision Tree model
    fit.tree <- get.model(
      rpart::rpart(observations ~ ., data=spec.df, method="anova",
                   control=rpart::rpart.control(cp=0.001)),
      file.name=paste0(tolower(gsub(" ", "_", spec)), "_", st, "_tree.rds"),
      model.path="artifacts/models/trees_fs")
    
    gc()
    fit.tree
  })
  names(spec.state.fit) <- states
  spec.state.fit
})
names(tree.list) <- species
```

```
# Get variable importance for decision tree models
tree.imp <- purrr::map_df(1:nrow(spec.state), function(i) {
  spec <- spec.state[i,]$common.name
  st <- spec.state[i,]$state
  fit <- tree.list[[spec]][[st]]
  vi <- fit$variable.importance
  # Create a data frame of variable importance
  var.importance <- tibble(
    common.name = spec,
    state = st,
    variable = names(vi),
    importance = (vi - min(vi)) / (max(vi) - min(vi))
  ) %>%
    # Rank variables by importance
    arrange(state, common.name, -importance, variable)
})

tree.imp
```

<div data-pagedtable="false" pagedtable-page="0" class="pagedtable-wrapper">
<script data-pagedtable-source="" type="application/json">
{"columns":[{"label":["common.name"],"name":[1],"type":["chr"],"align":["left"]},{"label":["state"],"name":[2],"type":["chr"],"align":["left"]},{"label":["variable"],"name":[3],"type":["chr"],"align":["left"]},{"label":["importance"],"name":[4],"type":["dbl"],"align":["right"]}],"data":[{"1":"Belted Kingfisher","2":"CO","3":"water","4":"1.0000000000"},{"1":"Belted Kingfisher","2":"CO","3":"dem.detrended","4":"0.8545118262"},{"1":"Belted Kingfisher","2":"CO","3":"temp.diff","4":"0.8476030006"},{"1":"Belted Kingfisher","2":"CO","3":"waterbody","4":"0.7946465968"},{"1":"Belted Kingfisher","2":"CO","3":"dem","4":"0.5451922439"},{"1":"Belted Kingfisher","2":"CO","3":"Fall_NDVI","4":"0.4012280915"},{"1":"Belted Kingfisher","2":"CO","3":"urban_imperviousness","4":"0.3221800637"},{"1":"Belted Kingfisher","2":"CO","3":"developed","4":"0.2831951900"},{"1":"Belted Kingfisher","2":"CO","3":"Winter_NDVI","4":"0.2698033769"},{"1":"Belted Kingfisher","2":"CO","3":"emergent_herbaceous_wetlands","4":"0.2697761852"},{"1":"Belted Kingfisher","2":"CO","3":"wetlands","4":"0.2697761852"},{"1":"Belted Kingfisher","2":"CO","3":"lon","4":"0.2520857179"},{"1":"Belted Kingfisher","2":"CO","3":"avg_prcp","4":"0.1826944382"},{"1":"Belted Kingfisher","2":"CO","3":"tmax","4":"0.1504288087"},{"1":"Belted Kingfisher","2":"CO","3":"Summer_NDVI","4":"0.1389604607"},{"1":"Belted Kingfisher","2":"CO","3":"lon.sqrd","4":"0.1289302404"},{"1":"Belted Kingfisher","2":"CO","3":"lon.lat.interaction","4":"0.0994649773"},{"1":"Belted Kingfisher","2":"CO","3":"open_water","4":"0.0784023257"},{"1":"Belted Kingfisher","2":"CO","3":"tmin","4":"0.0677655580"},{"1":"Belted Kingfisher","2":"CO","3":"lat","4":"0.0317966774"},{"1":"Belted Kingfisher","2":"CO","3":"lat.sqrd","4":"0.0317966774"},{"1":"Belted Kingfisher","2":"CO","3":"max.ndvi","4":"0.0219857423"},{"1":"Belted Kingfisher","2":"CO","3":"Spring_NDVI","4":"0.0114437997"},{"1":"Belted Kingfisher","2":"CO","3":"min.ndvi","4":"0.0031473708"},{"1":"Belted Kingfisher","2":"CO","3":"developed_high_intensity","4":"0.0000000000"},{"1":"Cedar Waxwing","2":"CO","3":"urban_imperviousness","4":"1.0000000000"},{"1":"Cedar Waxwing","2":"CO","3":"lat","4":"0.9575829604"},{"1":"Cedar Waxwing","2":"CO","3":"lon","4":"0.8501076749"},{"1":"Cedar Waxwing","2":"CO","3":"lat.sqrd","4":"0.7819617049"},{"1":"Cedar Waxwing","2":"CO","3":"dem","4":"0.7105478410"},{"1":"Cedar Waxwing","2":"CO","3":"Fall_NDVI","4":"0.6646716842"},{"1":"Cedar Waxwing","2":"CO","3":"developed","4":"0.6627924461"},{"1":"Cedar Waxwing","2":"CO","3":"lon.sqrd","4":"0.6159764330"},{"1":"Cedar Waxwing","2":"CO","3":"dem.detrended","4":"0.6011065416"},{"1":"Cedar Waxwing","2":"CO","3":"water","4":"0.5230809435"},{"1":"Cedar Waxwing","2":"CO","3":"max.ndvi","4":"0.5165668890"},{"1":"Cedar Waxwing","2":"CO","3":"lon.lat.interaction","4":"0.4815474451"},{"1":"Cedar Waxwing","2":"CO","3":"Spring_NDVI","4":"0.4074230713"},{"1":"Cedar Waxwing","2":"CO","3":"min.ndvi","4":"0.3502675256"},{"1":"Cedar Waxwing","2":"CO","3":"Summer_NDVI","4":"0.3383295444"},{"1":"Cedar Waxwing","2":"CO","3":"waterbody","4":"0.3266696734"},{"1":"Cedar Waxwing","2":"CO","3":"Winter_NDVI","4":"0.2855642115"},{"1":"Cedar Waxwing","2":"CO","3":"tmax","4":"0.2824664378"},{"1":"Cedar Waxwing","2":"CO","3":"avg_prcp","4":"0.2817523379"},{"1":"Cedar Waxwing","2":"CO","3":"shrubland","4":"0.2409818787"},{"1":"Cedar Waxwing","2":"CO","3":"tmin","4":"0.2363824007"},{"1":"Cedar Waxwing","2":"CO","3":"developed_low_intensity","4":"0.1657761454"},{"1":"Cedar Waxwing","2":"CO","3":"temp.diff","4":"0.1621125758"},{"1":"Cedar Waxwing","2":"CO","3":"open_water","4":"0.1335184102"},{"1":"Cedar Waxwing","2":"CO","3":"herbaceous","4":"0.1231602899"},{"1":"Cedar Waxwing","2":"CO","3":"emergent_herbaceous_wetlands","4":"0.0280599177"},{"1":"Cedar Waxwing","2":"CO","3":"developed_open_space","4":"0.0037935896"},{"1":"Cedar Waxwing","2":"CO","3":"developed_medium_intensity","4":"0.0019371521"},{"1":"Cedar Waxwing","2":"CO","3":"evergreen_forest","4":"0.0010139162"},{"1":"Cedar Waxwing","2":"CO","3":"developed_high_intensity","4":"0.0000000000"},{"1":"Downy Woodpecker","2":"CO","3":"dem","4":"1.0000000000"},{"1":"Downy Woodpecker","2":"CO","3":"dem.detrended","4":"0.8385421570"},{"1":"Downy Woodpecker","2":"CO","3":"water","4":"0.6755992568"},{"1":"Downy Woodpecker","2":"CO","3":"waterbody","4":"0.6458342049"},{"1":"Downy Woodpecker","2":"CO","3":"lon","4":"0.5736972240"},{"1":"Downy Woodpecker","2":"CO","3":"avg_prcp","4":"0.5716101416"},{"1":"Downy Woodpecker","2":"CO","3":"tmax","4":"0.5511130947"},{"1":"Downy Woodpecker","2":"CO","3":"lat","4":"0.5483662718"},{"1":"Downy Woodpecker","2":"CO","3":"lon.lat.interaction","4":"0.4954054524"},{"1":"Downy Woodpecker","2":"CO","3":"temp.diff","4":"0.4752087409"},{"1":"Downy Woodpecker","2":"CO","3":"lon.sqrd","4":"0.4376096039"},{"1":"Downy Woodpecker","2":"CO","3":"lat.sqrd","4":"0.3760914941"},{"1":"Downy Woodpecker","2":"CO","3":"Spring_NDVI","4":"0.3268294158"},{"1":"Downy Woodpecker","2":"CO","3":"tmin","4":"0.3010768497"},{"1":"Downy Woodpecker","2":"CO","3":"Summer_NDVI","4":"0.2732620180"},{"1":"Downy Woodpecker","2":"CO","3":"urban_imperviousness","4":"0.2544661917"},{"1":"Downy Woodpecker","2":"CO","3":"developed","4":"0.2359686026"},{"1":"Downy Woodpecker","2":"CO","3":"Fall_NDVI","4":"0.2219070383"},{"1":"Downy Woodpecker","2":"CO","3":"min.ndvi","4":"0.1754939841"},{"1":"Downy Woodpecker","2":"CO","3":"Winter_NDVI","4":"0.1523440801"},{"1":"Downy Woodpecker","2":"CO","3":"max.ndvi","4":"0.1335580113"},{"1":"Downy Woodpecker","2":"CO","3":"emergent_herbaceous_wetlands","4":"0.0951092239"},{"1":"Downy Woodpecker","2":"CO","3":"wetlands","4":"0.0766144515"},{"1":"Downy Woodpecker","2":"CO","3":"shrubland","4":"0.0527162180"},{"1":"Downy Woodpecker","2":"CO","3":"open_water","4":"0.0190355757"},{"1":"Downy Woodpecker","2":"CO","3":"deciduous_forest","4":"0.0100592431"},{"1":"Downy Woodpecker","2":"CO","3":"developed_open_space","4":"0.0085635813"},{"1":"Downy Woodpecker","2":"CO","3":"woody_wetlands","4":"0.0062747663"},{"1":"Downy Woodpecker","2":"CO","3":"developed_low_intensity","4":"0.0055411409"},{"1":"Downy Woodpecker","2":"CO","3":"developed_medium_intensity","4":"0.0029664693"},{"1":"Downy Woodpecker","2":"CO","3":"developed_high_intensity","4":"0.0000000000"},{"1":"Ruddy Duck","2":"CO","3":"water","4":"1.0000000000"},{"1":"Ruddy Duck","2":"CO","3":"waterbody","4":"0.6745622552"},{"1":"Ruddy Duck","2":"CO","3":"lon","4":"0.4753321906"},{"1":"Ruddy Duck","2":"CO","3":"dem.detrended","4":"0.4750941015"},{"1":"Ruddy Duck","2":"CO","3":"dem","4":"0.4412385640"},{"1":"Ruddy Duck","2":"CO","3":"max.ndvi","4":"0.4240337433"},{"1":"Ruddy Duck","2":"CO","3":"Summer_NDVI","4":"0.3615195679"},{"1":"Ruddy Duck","2":"CO","3":"lat","4":"0.3030882974"},{"1":"Ruddy Duck","2":"CO","3":"urban_imperviousness","4":"0.2855859820"},{"1":"Ruddy Duck","2":"CO","3":"developed","4":"0.2750348720"},{"1":"Ruddy Duck","2":"CO","3":"lon.sqrd","4":"0.2708945169"},{"1":"Ruddy Duck","2":"CO","3":"Winter_NDVI","4":"0.2605278060"},{"1":"Ruddy Duck","2":"CO","3":"Fall_NDVI","4":"0.2487934604"},{"1":"Ruddy Duck","2":"CO","3":"lat.sqrd","4":"0.2180903897"},{"1":"Ruddy Duck","2":"CO","3":"lon.lat.interaction","4":"0.1977943313"},{"1":"Ruddy Duck","2":"CO","3":"open_water","4":"0.1622558586"},{"1":"Ruddy Duck","2":"CO","3":"temp.diff","4":"0.1195762748"},{"1":"Ruddy Duck","2":"CO","3":"tmax","4":"0.0959956912"},{"1":"Ruddy Duck","2":"CO","3":"avg_prcp","4":"0.0823728594"},{"1":"Ruddy Duck","2":"CO","3":"tmin","4":"0.0297459157"},{"1":"Ruddy Duck","2":"CO","3":"min.ndvi","4":"0.0223953250"},{"1":"Ruddy Duck","2":"CO","3":"developed_low_intensity","4":"0.0158539723"},{"1":"Ruddy Duck","2":"CO","3":"wetlands","4":"0.0126343447"},{"1":"Ruddy Duck","2":"CO","3":"hay_pasture","4":"0.0070693549"},{"1":"Ruddy Duck","2":"CO","3":"developed_medium_intensity","4":"0.0000000000"},{"1":"Sanderling","2":"CO","3":"water","4":"1.0000000000"},{"1":"Sanderling","2":"CO","3":"waterbody","4":"0.8924349534"},{"1":"Sanderling","2":"CO","3":"max.ndvi","4":"0.2008918478"},{"1":"Sanderling","2":"CO","3":"open_water","4":"0.0882474452"},{"1":"Sanderling","2":"CO","3":"lon.lat.interaction","4":"0.0725406006"},{"1":"Sanderling","2":"CO","3":"lat","4":"0.0561696623"},{"1":"Sanderling","2":"CO","3":"Fall_NDVI","4":"0.0557786671"},{"1":"Sanderling","2":"CO","3":"lat.sqrd","4":"0.0496883437"},{"1":"Sanderling","2":"CO","3":"temp.diff","4":"0.0247240494"},{"1":"Sanderling","2":"CO","3":"tmax","4":"0.0247240494"},{"1":"Sanderling","2":"CO","3":"Spring_NDVI","4":"0.0122419022"},{"1":"Sanderling","2":"CO","3":"dem","4":"0.0012962637"},{"1":"Sanderling","2":"CO","3":"lon","4":"0.0012962637"},{"1":"Sanderling","2":"CO","3":"lon.sqrd","4":"0.0012962637"},{"1":"Sanderling","2":"CO","3":"avg_prcp","4":"0.0000000000"},{"1":"Sandhill Crane","2":"CO","3":"waterbody","4":"1.0000000000"},{"1":"Sandhill Crane","2":"CO","3":"water","4":"0.9932176443"},{"1":"Sandhill Crane","2":"CO","3":"woody_wetlands","4":"0.1686907425"},{"1":"Sandhill Crane","2":"CO","3":"lon.lat.interaction","4":"0.1659704445"},{"1":"Sandhill Crane","2":"CO","3":"lat","4":"0.1322290755"},{"1":"Sandhill Crane","2":"CO","3":"lat.sqrd","4":"0.1322290755"},{"1":"Sandhill Crane","2":"CO","3":"lon","4":"0.0804319553"},{"1":"Sandhill Crane","2":"CO","3":"lon.sqrd","4":"0.0804319553"},{"1":"Sandhill Crane","2":"CO","3":"tmin","4":"0.0569770881"},{"1":"Sandhill Crane","2":"CO","3":"Summer_NDVI","4":"0.0520303509"},{"1":"Sandhill Crane","2":"CO","3":"Fall_NDVI","4":"0.0217640563"},{"1":"Sandhill Crane","2":"CO","3":"avg_prcp","4":"0.0068309431"},{"1":"Sandhill Crane","2":"CO","3":"tmax","4":"0.0039214674"},{"1":"Sandhill Crane","2":"CO","3":"Spring_NDVI","4":"0.0028608883"},{"1":"Sandhill Crane","2":"CO","3":"Winter_NDVI","4":"0.0011652994"},{"1":"Sandhill Crane","2":"CO","3":"dem","4":"0.0011652994"},{"1":"Sandhill Crane","2":"CO","3":"dem.detrended","4":"0.0000000000"},{"1":"Sharp-shinned Hawk","2":"CO","3":"water","4":"1.0000000000"},{"1":"Sharp-shinned Hawk","2":"CO","3":"dem.detrended","4":"0.8333907422"},{"1":"Sharp-shinned Hawk","2":"CO","3":"waterbody","4":"0.7346013868"},{"1":"Sharp-shinned Hawk","2":"CO","3":"dem","4":"0.6965433304"},{"1":"Sharp-shinned Hawk","2":"CO","3":"lon","4":"0.5142104376"},{"1":"Sharp-shinned Hawk","2":"CO","3":"developed","4":"0.5027308556"},{"1":"Sharp-shinned Hawk","2":"CO","3":"urban_imperviousness","4":"0.4857823229"},{"1":"Sharp-shinned Hawk","2":"CO","3":"tmax","4":"0.4380800374"},{"1":"Sharp-shinned Hawk","2":"CO","3":"min.ndvi","4":"0.3962761873"},{"1":"Sharp-shinned Hawk","2":"CO","3":"Summer_NDVI","4":"0.3839120094"},{"1":"Sharp-shinned Hawk","2":"CO","3":"Spring_NDVI","4":"0.3806453001"},{"1":"Sharp-shinned Hawk","2":"CO","3":"lat","4":"0.3692269171"},{"1":"Sharp-shinned Hawk","2":"CO","3":"lon.sqrd","4":"0.3672184516"},{"1":"Sharp-shinned Hawk","2":"CO","3":"lon.lat.interaction","4":"0.3516245456"},{"1":"Sharp-shinned Hawk","2":"CO","3":"lat.sqrd","4":"0.3479490658"},{"1":"Sharp-shinned Hawk","2":"CO","3":"Fall_NDVI","4":"0.3388357496"},{"1":"Sharp-shinned Hawk","2":"CO","3":"avg_prcp","4":"0.3327998794"},{"1":"Sharp-shinned Hawk","2":"CO","3":"tmin","4":"0.3195752699"},{"1":"Sharp-shinned Hawk","2":"CO","3":"temp.diff","4":"0.3008659644"},{"1":"Sharp-shinned Hawk","2":"CO","3":"max.ndvi","4":"0.2676805018"},{"1":"Sharp-shinned Hawk","2":"CO","3":"Winter_NDVI","4":"0.1859190293"},{"1":"Sharp-shinned Hawk","2":"CO","3":"emergent_herbaceous_wetlands","4":"0.1813187288"},{"1":"Sharp-shinned Hawk","2":"CO","3":"open_water","4":"0.1522642484"},{"1":"Sharp-shinned Hawk","2":"CO","3":"wetlands","4":"0.0901066902"},{"1":"Sharp-shinned Hawk","2":"CO","3":"hay_pasture","4":"0.0395404734"},{"1":"Sharp-shinned Hawk","2":"CO","3":"shrubland","4":"0.0209174026"},{"1":"Sharp-shinned Hawk","2":"CO","3":"developed_open_space","4":"0.0063162771"},{"1":"Sharp-shinned Hawk","2":"CO","3":"forest","4":"0.0054562755"},{"1":"Sharp-shinned Hawk","2":"CO","3":"developed_low_intensity","4":"0.0053688356"},{"1":"Sharp-shinned Hawk","2":"CO","3":"developed_medium_intensity","4":"0.0037897663"},{"1":"Sharp-shinned Hawk","2":"CO","3":"developed_high_intensity","4":"0.0000000000"},{"1":"Wild Turkey","2":"CO","3":"Fall_NDVI","4":"1.0000000000"},{"1":"Wild Turkey","2":"CO","3":"max.ndvi","4":"0.8694374349"},{"1":"Wild Turkey","2":"CO","3":"dem.detrended","4":"0.4093010328"},{"1":"Wild Turkey","2":"CO","3":"Winter_NDVI","4":"0.3780587667"},{"1":"Wild Turkey","2":"CO","3":"lon","4":"0.3442523492"},{"1":"Wild Turkey","2":"CO","3":"lon.sqrd","4":"0.3442523492"},{"1":"Wild Turkey","2":"CO","3":"dem","4":"0.2331206011"},{"1":"Wild Turkey","2":"CO","3":"Summer_NDVI","4":"0.2219929974"},{"1":"Wild Turkey","2":"CO","3":"min.ndvi","4":"0.1982257652"},{"1":"Wild Turkey","2":"CO","3":"Spring_NDVI","4":"0.1526173978"},{"1":"Wild Turkey","2":"CO","3":"herbaceous","4":"0.1519922723"},{"1":"Wild Turkey","2":"CO","3":"temp.diff","4":"0.1429272543"},{"1":"Wild Turkey","2":"CO","3":"water","4":"0.1423808887"},{"1":"Wild Turkey","2":"CO","3":"developed_high_intensity","4":"0.0998946805"},{"1":"Wild Turkey","2":"CO","3":"avg_prcp","4":"0.0985603267"},{"1":"Wild Turkey","2":"CO","3":"tmax","4":"0.0928474269"},{"1":"Wild Turkey","2":"CO","3":"waterbody","4":"0.0786079251"},{"1":"Wild Turkey","2":"CO","3":"tmin","4":"0.0689261978"},{"1":"Wild Turkey","2":"CO","3":"forest","4":"0.0637592832"},{"1":"Wild Turkey","2":"CO","3":"evergreen_forest","4":"0.0617772499"},{"1":"Wild Turkey","2":"CO","3":"lon.lat.interaction","4":"0.0402370242"},{"1":"Wild Turkey","2":"CO","3":"developed","4":"0.0358065327"},{"1":"Wild Turkey","2":"CO","3":"shrubland","4":"0.0332250903"},{"1":"Wild Turkey","2":"CO","3":"lat","4":"0.0282340483"},{"1":"Wild Turkey","2":"CO","3":"urban_imperviousness","4":"0.0207948295"},{"1":"Wild Turkey","2":"CO","3":"lat.sqrd","4":"0.0140700629"},{"1":"Wild Turkey","2":"CO","3":"open_water","4":"0.0071715435"},{"1":"Wild Turkey","2":"CO","3":"woody_wetlands","4":"0.0000000000"},{"1":"Belted Kingfisher","2":"NC","3":"avg_prcp","4":"1.0000000000"},{"1":"Belted Kingfisher","2":"NC","3":"tmax","4":"0.9900927983"},{"1":"Belted Kingfisher","2":"NC","3":"lat","4":"0.5955567708"},{"1":"Belted Kingfisher","2":"NC","3":"lat.sqrd","4":"0.5417044151"},{"1":"Belted Kingfisher","2":"NC","3":"Winter_NDVI","4":"0.4481653456"},{"1":"Belted Kingfisher","2":"NC","3":"dem","4":"0.4142435754"},{"1":"Belted Kingfisher","2":"NC","3":"tmin","4":"0.4060997696"},{"1":"Belted Kingfisher","2":"NC","3":"lon.lat.interaction","4":"0.3832630445"},{"1":"Belted Kingfisher","2":"NC","3":"urban_imperviousness","4":"0.3654578126"},{"1":"Belted Kingfisher","2":"NC","3":"Fall_NDVI","4":"0.2800683021"},{"1":"Belted Kingfisher","2":"NC","3":"max.ndvi","4":"0.2767976653"},{"1":"Belted Kingfisher","2":"NC","3":"developed","4":"0.2521604912"},{"1":"Belted Kingfisher","2":"NC","3":"Spring_NDVI","4":"0.2033126675"},{"1":"Belted Kingfisher","2":"NC","3":"deciduous_forest","4":"0.2021470384"},{"1":"Belted Kingfisher","2":"NC","3":"lon","4":"0.1346714307"},{"1":"Belted Kingfisher","2":"NC","3":"lon.sqrd","4":"0.1205323605"},{"1":"Belted Kingfisher","2":"NC","3":"min.ndvi","4":"0.1187644294"},{"1":"Belted Kingfisher","2":"NC","3":"dem.detrended","4":"0.0953502146"},{"1":"Belted Kingfisher","2":"NC","3":"temp.diff","4":"0.0800527015"},{"1":"Belted Kingfisher","2":"NC","3":"Summer_NDVI","4":"0.0594568102"},{"1":"Belted Kingfisher","2":"NC","3":"hay_pasture","4":"0.0536310282"},{"1":"Belted Kingfisher","2":"NC","3":"water","4":"0.0395765145"},{"1":"Belted Kingfisher","2":"NC","3":"waterbody","4":"0.0394139799"},{"1":"Belted Kingfisher","2":"NC","3":"planted.cultivated","4":"0.0180594562"},{"1":"Belted Kingfisher","2":"NC","3":"developed_medium_intensity","4":"0.0172202137"},{"1":"Belted Kingfisher","2":"NC","3":"open_water","4":"0.0118140775"},{"1":"Belted Kingfisher","2":"NC","3":"developed_low_intensity","4":"0.0000000000"},{"1":"Cedar Waxwing","2":"NC","3":"tmax","4":"1.0000000000"},{"1":"Cedar Waxwing","2":"NC","3":"avg_prcp","4":"0.8692027483"},{"1":"Cedar Waxwing","2":"NC","3":"min.ndvi","4":"0.7479857365"},{"1":"Cedar Waxwing","2":"NC","3":"tmin","4":"0.7274344489"},{"1":"Cedar Waxwing","2":"NC","3":"lat","4":"0.6919517948"},{"1":"Cedar Waxwing","2":"NC","3":"lat.sqrd","4":"0.6919517948"},{"1":"Cedar Waxwing","2":"NC","3":"Spring_NDVI","4":"0.6635825655"},{"1":"Cedar Waxwing","2":"NC","3":"Winter_NDVI","4":"0.6171791210"},{"1":"Cedar Waxwing","2":"NC","3":"developed","4":"0.4700508494"},{"1":"Cedar Waxwing","2":"NC","3":"Fall_NDVI","4":"0.4237838871"},{"1":"Cedar Waxwing","2":"NC","3":"Summer_NDVI","4":"0.3589015436"},{"1":"Cedar Waxwing","2":"NC","3":"lon","4":"0.3163960242"},{"1":"Cedar Waxwing","2":"NC","3":"lon.sqrd","4":"0.2698902377"},{"1":"Cedar Waxwing","2":"NC","3":"max.ndvi","4":"0.2046090900"},{"1":"Cedar Waxwing","2":"NC","3":"dem.detrended","4":"0.1870020201"},{"1":"Cedar Waxwing","2":"NC","3":"urban_imperviousness","4":"0.1445527153"},{"1":"Cedar Waxwing","2":"NC","3":"dem","4":"0.0895717234"},{"1":"Cedar Waxwing","2":"NC","3":"developed_high_intensity","4":"0.0861890603"},{"1":"Cedar Waxwing","2":"NC","3":"lon.lat.interaction","4":"0.0634839451"},{"1":"Cedar Waxwing","2":"NC","3":"temp.diff","4":"0.0397296774"},{"1":"Cedar Waxwing","2":"NC","3":"emergent_herbaceous_wetlands","4":"0.0395059028"},{"1":"Cedar Waxwing","2":"NC","3":"coastline","4":"0.0356334762"},{"1":"Cedar Waxwing","2":"NC","3":"wetlands","4":"0.0142262898"},{"1":"Cedar Waxwing","2":"NC","3":"developed_open_space","4":"0.0129507456"},{"1":"Cedar Waxwing","2":"NC","3":"waterbody","4":"0.0126675810"},{"1":"Cedar Waxwing","2":"NC","3":"water","4":"0.0124937291"},{"1":"Cedar Waxwing","2":"NC","3":"developed_medium_intensity","4":"0.0000000000"},{"1":"Downy Woodpecker","2":"NC","3":"lat","4":"1.0000000000"},{"1":"Downy Woodpecker","2":"NC","3":"lat.sqrd","4":"0.9720565167"},{"1":"Downy Woodpecker","2":"NC","3":"urban_imperviousness","4":"0.8275239033"},{"1":"Downy Woodpecker","2":"NC","3":"dem","4":"0.7350117483"},{"1":"Downy Woodpecker","2":"NC","3":"lon.lat.interaction","4":"0.5720541156"},{"1":"Downy Woodpecker","2":"NC","3":"max.ndvi","4":"0.5019635413"},{"1":"Downy Woodpecker","2":"NC","3":"avg_prcp","4":"0.4924928374"},{"1":"Downy Woodpecker","2":"NC","3":"lon","4":"0.4441324287"},{"1":"Downy Woodpecker","2":"NC","3":"Winter_NDVI","4":"0.4193783935"},{"1":"Downy Woodpecker","2":"NC","3":"Fall_NDVI","4":"0.4027057001"},{"1":"Downy Woodpecker","2":"NC","3":"tmin","4":"0.3511270937"},{"1":"Downy Woodpecker","2":"NC","3":"Spring_NDVI","4":"0.3484566823"},{"1":"Downy Woodpecker","2":"NC","3":"lon.sqrd","4":"0.3032648211"},{"1":"Downy Woodpecker","2":"NC","3":"developed","4":"0.3001982586"},{"1":"Downy Woodpecker","2":"NC","3":"min.ndvi","4":"0.2752205576"},{"1":"Downy Woodpecker","2":"NC","3":"tmax","4":"0.2237159893"},{"1":"Downy Woodpecker","2":"NC","3":"dem.detrended","4":"0.1614255894"},{"1":"Downy Woodpecker","2":"NC","3":"Summer_NDVI","4":"0.1613861832"},{"1":"Downy Woodpecker","2":"NC","3":"emergent_herbaceous_wetlands","4":"0.1577429361"},{"1":"Downy Woodpecker","2":"NC","3":"waterbody","4":"0.1218234900"},{"1":"Downy Woodpecker","2":"NC","3":"water","4":"0.1188389204"},{"1":"Downy Woodpecker","2":"NC","3":"developed_open_space","4":"0.0650980559"},{"1":"Downy Woodpecker","2":"NC","3":"temp.diff","4":"0.0619954279"},{"1":"Downy Woodpecker","2":"NC","3":"developed_low_intensity","4":"0.0358853875"},{"1":"Downy Woodpecker","2":"NC","3":"developed_medium_intensity","4":"0.0357999374"},{"1":"Downy Woodpecker","2":"NC","3":"evergreen_forest","4":"0.0268193851"},{"1":"Downy Woodpecker","2":"NC","3":"deciduous_forest","4":"0.0063967535"},{"1":"Downy Woodpecker","2":"NC","3":"developed_high_intensity","4":"0.0055730878"},{"1":"Downy Woodpecker","2":"NC","3":"open_water","4":"0.0028156718"},{"1":"Downy Woodpecker","2":"NC","3":"wetlands","4":"0.0000000000"},{"1":"Downy Woodpecker","2":"NC","3":"woody_wetlands","4":"0.0000000000"},{"1":"Ruddy Duck","2":"NC","3":"Spring_NDVI","4":"1.0000000000"},{"1":"Ruddy Duck","2":"NC","3":"Fall_NDVI","4":"0.9204373676"},{"1":"Ruddy Duck","2":"NC","3":"Winter_NDVI","4":"0.8665517002"},{"1":"Ruddy Duck","2":"NC","3":"min.ndvi","4":"0.8041865783"},{"1":"Ruddy Duck","2":"NC","3":"max.ndvi","4":"0.6233550623"},{"1":"Ruddy Duck","2":"NC","3":"lon","4":"0.5824738659"},{"1":"Ruddy Duck","2":"NC","3":"dem","4":"0.5735902364"},{"1":"Ruddy Duck","2":"NC","3":"tmin","4":"0.5087970057"},{"1":"Ruddy Duck","2":"NC","3":"dem.detrended","4":"0.4180471169"},{"1":"Ruddy Duck","2":"NC","3":"lat","4":"0.2878992967"},{"1":"Ruddy Duck","2":"NC","3":"lon.sqrd","4":"0.2855463112"},{"1":"Ruddy Duck","2":"NC","3":"lon.lat.interaction","4":"0.0726024250"},{"1":"Ruddy Duck","2":"NC","3":"avg_prcp","4":"0.0253192261"},{"1":"Ruddy Duck","2":"NC","3":"Summer_NDVI","4":"0.0217573338"},{"1":"Ruddy Duck","2":"NC","3":"deciduous_forest","4":"0.0068276214"},{"1":"Ruddy Duck","2":"NC","3":"tmax","4":"0.0000000000"},{"1":"Sanderling","2":"NC","3":"dem","4":"1.0000000000"},{"1":"Sanderling","2":"NC","3":"dem.detrended","4":"0.4530310718"},{"1":"Sanderling","2":"NC","3":"lon.lat.interaction","4":"0.3741558829"},{"1":"Sanderling","2":"NC","3":"lat","4":"0.3708860894"},{"1":"Sanderling","2":"NC","3":"lat.sqrd","4":"0.3708860894"},{"1":"Sanderling","2":"NC","3":"max.ndvi","4":"0.2702967921"},{"1":"Sanderling","2":"NC","3":"Spring_NDVI","4":"0.2467318524"},{"1":"Sanderling","2":"NC","3":"lon","4":"0.2330645030"},{"1":"Sanderling","2":"NC","3":"min.ndvi","4":"0.2323215573"},{"1":"Sanderling","2":"NC","3":"Winter_NDVI","4":"0.2124026393"},{"1":"Sanderling","2":"NC","3":"Fall_NDVI","4":"0.1966487577"},{"1":"Sanderling","2":"NC","3":"Summer_NDVI","4":"0.1730838180"},{"1":"Sanderling","2":"NC","3":"avg_prcp","4":"0.0666547525"},{"1":"Sanderling","2":"NC","3":"woody_wetlands","4":"0.0075807194"},{"1":"Sanderling","2":"NC","3":"tmin","4":"0.0000000000"},{"1":"Sandhill Crane","2":"NC","3":"waterbody","4":"1.0000000000"},{"1":"Sandhill Crane","2":"NC","3":"avg_prcp","4":"0.2104934281"},{"1":"Sandhill Crane","2":"NC","3":"water","4":"0.1342184979"},{"1":"Sandhill Crane","2":"NC","3":"tmax","4":"0.0588702231"},{"1":"Sandhill Crane","2":"NC","3":"max.ndvi","4":"0.0439745982"},{"1":"Sandhill Crane","2":"NC","3":"Fall_NDVI","4":"0.0431347868"},{"1":"Sandhill Crane","2":"NC","3":"coastline","4":"0.0196200663"},{"1":"Sandhill Crane","2":"NC","3":"lon.lat.interaction","4":"0.0179404434"},{"1":"Sandhill Crane","2":"NC","3":"Spring_NDVI","4":"0.0087025175"},{"1":"Sandhill Crane","2":"NC","3":"Winter_NDVI","4":"0.0070228946"},{"1":"Sandhill Crane","2":"NC","3":"temp.diff","4":"0.0034594502"},{"1":"Sandhill Crane","2":"NC","3":"tmin","4":"0.0024855092"},{"1":"Sandhill Crane","2":"NC","3":"lat","4":"0.0005376274"},{"1":"Sandhill Crane","2":"NC","3":"lat.sqrd","4":"0.0005376274"},{"1":"Sandhill Crane","2":"NC","3":"open_water","4":"0.0000000000"},{"1":"Sharp-shinned Hawk","2":"NC","3":"lon","4":"1.0000000000"},{"1":"Sharp-shinned Hawk","2":"NC","3":"lon.sqrd","4":"0.9740540583"},{"1":"Sharp-shinned Hawk","2":"NC","3":"dem","4":"0.7689940294"},{"1":"Sharp-shinned Hawk","2":"NC","3":"dem.detrended","4":"0.6325724449"},{"1":"Sharp-shinned Hawk","2":"NC","3":"tmax","4":"0.5350567743"},{"1":"Sharp-shinned Hawk","2":"NC","3":"Winter_NDVI","4":"0.4577869377"},{"1":"Sharp-shinned Hawk","2":"NC","3":"temp.diff","4":"0.4425122509"},{"1":"Sharp-shinned Hawk","2":"NC","3":"Spring_NDVI","4":"0.4173709521"},{"1":"Sharp-shinned Hawk","2":"NC","3":"avg_prcp","4":"0.4120958942"},{"1":"Sharp-shinned Hawk","2":"NC","3":"tmin","4":"0.3738678813"},{"1":"Sharp-shinned Hawk","2":"NC","3":"min.ndvi","4":"0.3693675354"},{"1":"Sharp-shinned Hawk","2":"NC","3":"Fall_NDVI","4":"0.3253084840"},{"1":"Sharp-shinned Hawk","2":"NC","3":"max.ndvi","4":"0.2671369664"},{"1":"Sharp-shinned Hawk","2":"NC","3":"lat","4":"0.1860118298"},{"1":"Sharp-shinned Hawk","2":"NC","3":"lat.sqrd","4":"0.1860118298"},{"1":"Sharp-shinned Hawk","2":"NC","3":"developed_medium_intensity","4":"0.1344900431"},{"1":"Sharp-shinned Hawk","2":"NC","3":"urban_imperviousness","4":"0.0974421072"},{"1":"Sharp-shinned Hawk","2":"NC","3":"developed","4":"0.0935382373"},{"1":"Sharp-shinned Hawk","2":"NC","3":"Summer_NDVI","4":"0.0605633873"},{"1":"Sharp-shinned Hawk","2":"NC","3":"mixed_forest","4":"0.0242430663"},{"1":"Sharp-shinned Hawk","2":"NC","3":"lon.lat.interaction","4":"0.0211973638"},{"1":"Sharp-shinned Hawk","2":"NC","3":"coastline","4":"0.0075692182"},{"1":"Sharp-shinned Hawk","2":"NC","3":"waterbody","4":"0.0031041947"},{"1":"Sharp-shinned Hawk","2":"NC","3":"evergreen_forest","4":"0.0009974984"},{"1":"Sharp-shinned Hawk","2":"NC","3":"water","4":"0.0000000000"},{"1":"Wild Turkey","2":"NC","3":"tmax","4":"1.0000000000"},{"1":"Wild Turkey","2":"NC","3":"lon.lat.interaction","4":"0.4202815005"},{"1":"Wild Turkey","2":"NC","3":"avg_prcp","4":"0.3363239621"},{"1":"Wild Turkey","2":"NC","3":"dem","4":"0.2188951793"},{"1":"Wild Turkey","2":"NC","3":"lat","4":"0.1384551696"},{"1":"Wild Turkey","2":"NC","3":"lon","4":"0.1044790145"},{"1":"Wild Turkey","2":"NC","3":"Fall_NDVI","4":"0.0970534733"},{"1":"Wild Turkey","2":"NC","3":"Spring_NDVI","4":"0.0936821344"},{"1":"Wild Turkey","2":"NC","3":"Summer_NDVI","4":"0.0799186620"},{"1":"Wild Turkey","2":"NC","3":"developed_high_intensity","4":"0.0781518473"},{"1":"Wild Turkey","2":"NC","3":"temp.diff","4":"0.0634217210"},{"1":"Wild Turkey","2":"NC","3":"lat.sqrd","4":"0.0557619801"},{"1":"Wild Turkey","2":"NC","3":"developed","4":"0.0305571337"},{"1":"Wild Turkey","2":"NC","3":"tmin","4":"0.0221305143"},{"1":"Wild Turkey","2":"NC","3":"lon.sqrd","4":"0.0217858250"},{"1":"Wild Turkey","2":"NC","3":"max.ndvi","4":"0.0171693516"},{"1":"Wild Turkey","2":"NC","3":"Winter_NDVI","4":"0.0103960923"},{"1":"Wild Turkey","2":"NC","3":"min.ndvi","4":"0.0092442243"},{"1":"Wild Turkey","2":"NC","3":"urban_imperviousness","4":"0.0053955974"},{"1":"Wild Turkey","2":"NC","3":"developed_medium_intensity","4":"0.0033581990"},{"1":"Wild Turkey","2":"NC","3":"coastline","4":"0.0000000000"},{"1":"Belted Kingfisher","2":"OR","3":"dem","4":"1.0000000000"},{"1":"Belted Kingfisher","2":"OR","3":"lon.lat.interaction","4":"0.7135136900"},{"1":"Belted Kingfisher","2":"OR","3":"tmin","4":"0.6511305273"},{"1":"Belted Kingfisher","2":"OR","3":"temp.diff","4":"0.5173822565"},{"1":"Belted Kingfisher","2":"OR","3":"dem.detrended","4":"0.4939502013"},{"1":"Belted Kingfisher","2":"OR","3":"lon","4":"0.4652149262"},{"1":"Belted Kingfisher","2":"OR","3":"herbaceous","4":"0.4608245666"},{"1":"Belted Kingfisher","2":"OR","3":"lat","4":"0.4474459851"},{"1":"Belted Kingfisher","2":"OR","3":"tmax","4":"0.4101740567"},{"1":"Belted Kingfisher","2":"OR","3":"developed","4":"0.4078548325"},{"1":"Belted Kingfisher","2":"OR","3":"lat.sqrd","4":"0.3812718703"},{"1":"Belted Kingfisher","2":"OR","3":"avg_prcp","4":"0.3451914148"},{"1":"Belted Kingfisher","2":"OR","3":"Winter_NDVI","4":"0.3218915149"},{"1":"Belted Kingfisher","2":"OR","3":"min.ndvi","4":"0.3133487099"},{"1":"Belted Kingfisher","2":"OR","3":"Fall_NDVI","4":"0.3063446675"},{"1":"Belted Kingfisher","2":"OR","3":"max.ndvi","4":"0.2930945544"},{"1":"Belted Kingfisher","2":"OR","3":"Summer_NDVI","4":"0.2788955455"},{"1":"Belted Kingfisher","2":"OR","3":"lon.sqrd","4":"0.2775425283"},{"1":"Belted Kingfisher","2":"OR","3":"coastline","4":"0.2585445817"},{"1":"Belted Kingfisher","2":"OR","3":"Spring_NDVI","4":"0.1744778261"},{"1":"Belted Kingfisher","2":"OR","3":"urban_imperviousness","4":"0.1095784774"},{"1":"Belted Kingfisher","2":"OR","3":"waterbody","4":"0.0999770677"},{"1":"Belted Kingfisher","2":"OR","3":"water","4":"0.0989029144"},{"1":"Belted Kingfisher","2":"OR","3":"emergent_herbaceous_wetlands","4":"0.0369616004"},{"1":"Belted Kingfisher","2":"OR","3":"wetlands","4":"0.0223225959"},{"1":"Belted Kingfisher","2":"OR","3":"developed_open_space","4":"0.0063623159"},{"1":"Belted Kingfisher","2":"OR","3":"developed_low_intensity","4":"0.0012940304"},{"1":"Belted Kingfisher","2":"OR","3":"developed_medium_intensity","4":"0.0000000000"},{"1":"Cedar Waxwing","2":"OR","3":"Spring_NDVI","4":"1.0000000000"},{"1":"Cedar Waxwing","2":"OR","3":"Summer_NDVI","4":"0.6164387302"},{"1":"Cedar Waxwing","2":"OR","3":"max.ndvi","4":"0.6153948992"},{"1":"Cedar Waxwing","2":"OR","3":"developed","4":"0.5856734921"},{"1":"Cedar Waxwing","2":"OR","3":"Winter_NDVI","4":"0.5479319950"},{"1":"Cedar Waxwing","2":"OR","3":"water","4":"0.3370673113"},{"1":"Cedar Waxwing","2":"OR","3":"dem.detrended","4":"0.3075706725"},{"1":"Cedar Waxwing","2":"OR","3":"urban_imperviousness","4":"0.3007571144"},{"1":"Cedar Waxwing","2":"OR","3":"dem","4":"0.2790853032"},{"1":"Cedar Waxwing","2":"OR","3":"avg_prcp","4":"0.2452725442"},{"1":"Cedar Waxwing","2":"OR","3":"waterbody","4":"0.2022515164"},{"1":"Cedar Waxwing","2":"OR","3":"lon.lat.interaction","4":"0.1946561710"},{"1":"Cedar Waxwing","2":"OR","3":"lat","4":"0.1874800762"},{"1":"Cedar Waxwing","2":"OR","3":"lat.sqrd","4":"0.1659004298"},{"1":"Cedar Waxwing","2":"OR","3":"emergent_herbaceous_wetlands","4":"0.1523446162"},{"1":"Cedar Waxwing","2":"OR","3":"herbaceous","4":"0.1461168698"},{"1":"Cedar Waxwing","2":"OR","3":"tmax","4":"0.1381875015"},{"1":"Cedar Waxwing","2":"OR","3":"min.ndvi","4":"0.1355636241"},{"1":"Cedar Waxwing","2":"OR","3":"open_water","4":"0.1323618704"},{"1":"Cedar Waxwing","2":"OR","3":"tmin","4":"0.1013858639"},{"1":"Cedar Waxwing","2":"OR","3":"lon","4":"0.0655119690"},{"1":"Cedar Waxwing","2":"OR","3":"wetlands","4":"0.0628902318"},{"1":"Cedar Waxwing","2":"OR","3":"developed_medium_intensity","4":"0.0612662040"},{"1":"Cedar Waxwing","2":"OR","3":"temp.diff","4":"0.0602576637"},{"1":"Cedar Waxwing","2":"OR","3":"developed_low_intensity","4":"0.0561298803"},{"1":"Cedar Waxwing","2":"OR","3":"developed_open_space","4":"0.0561298803"},{"1":"Cedar Waxwing","2":"OR","3":"lon.sqrd","4":"0.0373842527"},{"1":"Cedar Waxwing","2":"OR","3":"coastline","4":"0.0000000000"},{"1":"Downy Woodpecker","2":"OR","3":"tmax","4":"1.0000000000"},{"1":"Downy Woodpecker","2":"OR","3":"lat","4":"0.9414470611"},{"1":"Downy Woodpecker","2":"OR","3":"lon","4":"0.8924486051"},{"1":"Downy Woodpecker","2":"OR","3":"urban_imperviousness","4":"0.8720701776"},{"1":"Downy Woodpecker","2":"OR","3":"tmin","4":"0.8543455667"},{"1":"Downy Woodpecker","2":"OR","3":"avg_prcp","4":"0.8334453355"},{"1":"Downy Woodpecker","2":"OR","3":"dem","4":"0.7689199092"},{"1":"Downy Woodpecker","2":"OR","3":"dem.detrended","4":"0.7223609070"},{"1":"Downy Woodpecker","2":"OR","3":"Spring_NDVI","4":"0.5467188833"},{"1":"Downy Woodpecker","2":"OR","3":"developed","4":"0.5241067441"},{"1":"Downy Woodpecker","2":"OR","3":"Winter_NDVI","4":"0.4342581373"},{"1":"Downy Woodpecker","2":"OR","3":"Summer_NDVI","4":"0.4319641349"},{"1":"Downy Woodpecker","2":"OR","3":"min.ndvi","4":"0.3523755423"},{"1":"Downy Woodpecker","2":"OR","3":"lon.lat.interaction","4":"0.2299682895"},{"1":"Downy Woodpecker","2":"OR","3":"max.ndvi","4":"0.2204017573"},{"1":"Downy Woodpecker","2":"OR","3":"lat.sqrd","4":"0.2037118811"},{"1":"Downy Woodpecker","2":"OR","3":"temp.diff","4":"0.1931012322"},{"1":"Downy Woodpecker","2":"OR","3":"water","4":"0.1786801451"},{"1":"Downy Woodpecker","2":"OR","3":"open_water","4":"0.1255551694"},{"1":"Downy Woodpecker","2":"OR","3":"lon.sqrd","4":"0.1124076390"},{"1":"Downy Woodpecker","2":"OR","3":"herbaceous","4":"0.1023900148"},{"1":"Downy Woodpecker","2":"OR","3":"developed_medium_intensity","4":"0.0756607025"},{"1":"Downy Woodpecker","2":"OR","3":"developed_low_intensity","4":"0.0487500257"},{"1":"Downy Woodpecker","2":"OR","3":"waterbody","4":"0.0443552792"},{"1":"Downy Woodpecker","2":"OR","3":"Fall_NDVI","4":"0.0251254568"},{"1":"Downy Woodpecker","2":"OR","3":"coastline","4":"0.0044597128"},{"1":"Downy Woodpecker","2":"OR","3":"emergent_herbaceous_wetlands","4":"0.0041054287"},{"1":"Downy Woodpecker","2":"OR","3":"hay_pasture","4":"0.0000000000"},{"1":"Ruddy Duck","2":"OR","3":"urban_imperviousness","4":"1.0000000000"},{"1":"Ruddy Duck","2":"OR","3":"developed","4":"0.5686106841"},{"1":"Ruddy Duck","2":"OR","3":"Winter_NDVI","4":"0.4442952883"},{"1":"Ruddy Duck","2":"OR","3":"Summer_NDVI","4":"0.4285798796"},{"1":"Ruddy Duck","2":"OR","3":"max.ndvi","4":"0.4172321144"},{"1":"Ruddy Duck","2":"OR","3":"min.ndvi","4":"0.3648609320"},{"1":"Ruddy Duck","2":"OR","3":"tmin","4":"0.3413724834"},{"1":"Ruddy Duck","2":"OR","3":"Spring_NDVI","4":"0.3037544626"},{"1":"Ruddy Duck","2":"OR","3":"lat","4":"0.1310815898"},{"1":"Ruddy Duck","2":"OR","3":"lat.sqrd","4":"0.1310815898"},{"1":"Ruddy Duck","2":"OR","3":"dem","4":"0.0550111720"},{"1":"Ruddy Duck","2":"OR","3":"developed_low_intensity","4":"0.0369642038"},{"1":"Ruddy Duck","2":"OR","3":"water","4":"0.0369642038"},{"1":"Ruddy Duck","2":"OR","3":"waterbody","4":"0.0369642038"},{"1":"Ruddy Duck","2":"OR","3":"lon.lat.interaction","4":"0.0086637184"},{"1":"Ruddy Duck","2":"OR","3":"developed_medium_intensity","4":"0.0034044228"},{"1":"Ruddy Duck","2":"OR","3":"developed_high_intensity","4":"0.0000000000"},{"1":"Sanderling","2":"OR","3":"coastline","4":"1.0000000000"},{"1":"Sanderling","2":"OR","3":"min.ndvi","4":"0.9834684301"},{"1":"Sanderling","2":"OR","3":"max.ndvi","4":"0.7792357087"},{"1":"Sanderling","2":"OR","3":"Winter_NDVI","4":"0.7671347109"},{"1":"Sanderling","2":"OR","3":"Spring_NDVI","4":"0.6993589393"},{"1":"Sanderling","2":"OR","3":"Summer_NDVI","4":"0.6263309098"},{"1":"Sanderling","2":"OR","3":"Fall_NDVI","4":"0.6233102436"},{"1":"Sanderling","2":"OR","3":"lon.lat.interaction","4":"0.4086675110"},{"1":"Sanderling","2":"OR","3":"lon","4":"0.3575731863"},{"1":"Sanderling","2":"OR","3":"lon.sqrd","4":"0.3575731863"},{"1":"Sanderling","2":"OR","3":"dem.detrended","4":"0.2699603143"},{"1":"Sanderling","2":"OR","3":"lat","4":"0.2699603143"},{"1":"Sanderling","2":"OR","3":"lat.sqrd","4":"0.2699603143"},{"1":"Sanderling","2":"OR","3":"dem","4":"0.2557247930"},{"1":"Sanderling","2":"OR","3":"temp.diff","4":"0.0021829575"},{"1":"Sanderling","2":"OR","3":"tmax","4":"0.0021829575"},{"1":"Sanderling","2":"OR","3":"tmin","4":"0.0021829575"},{"1":"Sanderling","2":"OR","3":"avg_prcp","4":"0.0000000000"},{"1":"Sandhill Crane","2":"OR","3":"dem","4":"1.0000000000"},{"1":"Sandhill Crane","2":"OR","3":"tmin","4":"0.4572272073"},{"1":"Sandhill Crane","2":"OR","3":"lat","4":"0.4014621049"},{"1":"Sandhill Crane","2":"OR","3":"dem.detrended","4":"0.3036785358"},{"1":"Sandhill Crane","2":"OR","3":"lon.lat.interaction","4":"0.2951310739"},{"1":"Sandhill Crane","2":"OR","3":"avg_prcp","4":"0.2399568434"},{"1":"Sandhill Crane","2":"OR","3":"lon","4":"0.2210271818"},{"1":"Sandhill Crane","2":"OR","3":"urban_imperviousness","4":"0.2210271818"},{"1":"Sandhill Crane","2":"OR","3":"tmax","4":"0.1840152333"},{"1":"Sandhill Crane","2":"OR","3":"lat.sqrd","4":"0.1553369880"},{"1":"Sandhill Crane","2":"OR","3":"emergent_herbaceous_wetlands","4":"0.0466278207"},{"1":"Sandhill Crane","2":"OR","3":"water","4":"0.0021759688"},{"1":"Sandhill Crane","2":"OR","3":"developed","4":"0.0000000000"},{"1":"Sharp-shinned Hawk","2":"OR","3":"dem","4":"1.0000000000"},{"1":"Sharp-shinned Hawk","2":"OR","3":"dem.detrended","4":"0.7992871995"},{"1":"Sharp-shinned Hawk","2":"OR","3":"waterbody","4":"0.4056020482"},{"1":"Sharp-shinned Hawk","2":"OR","3":"lon.lat.interaction","4":"0.2669567254"},{"1":"Sharp-shinned Hawk","2":"OR","3":"water","4":"0.2668753998"},{"1":"Sharp-shinned Hawk","2":"OR","3":"developed","4":"0.2457149214"},{"1":"Sharp-shinned Hawk","2":"OR","3":"urban_imperviousness","4":"0.2453609518"},{"1":"Sharp-shinned Hawk","2":"OR","3":"max.ndvi","4":"0.2075077629"},{"1":"Sharp-shinned Hawk","2":"OR","3":"lat","4":"0.1855630142"},{"1":"Sharp-shinned Hawk","2":"OR","3":"Fall_NDVI","4":"0.1821267625"},{"1":"Sharp-shinned Hawk","2":"OR","3":"Summer_NDVI","4":"0.1696669497"},{"1":"Sharp-shinned Hawk","2":"OR","3":"lat.sqrd","4":"0.1452583600"},{"1":"Sharp-shinned Hawk","2":"OR","3":"min.ndvi","4":"0.1326632239"},{"1":"Sharp-shinned Hawk","2":"OR","3":"lon","4":"0.1321724664"},{"1":"Sharp-shinned Hawk","2":"OR","3":"Spring_NDVI","4":"0.1016633430"},{"1":"Sharp-shinned Hawk","2":"OR","3":"open_water","4":"0.0923866886"},{"1":"Sharp-shinned Hawk","2":"OR","3":"lon.sqrd","4":"0.0834140260"},{"1":"Sharp-shinned Hawk","2":"OR","3":"temp.diff","4":"0.0750408622"},{"1":"Sharp-shinned Hawk","2":"OR","3":"tmax","4":"0.0447877440"},{"1":"Sharp-shinned Hawk","2":"OR","3":"avg_prcp","4":"0.0360172745"},{"1":"Sharp-shinned Hawk","2":"OR","3":"herbaceous","4":"0.0205018548"},{"1":"Sharp-shinned Hawk","2":"OR","3":"evergreen_forest","4":"0.0123066482"},{"1":"Sharp-shinned Hawk","2":"OR","3":"developed_open_space","4":"0.0105521195"},{"1":"Sharp-shinned Hawk","2":"OR","3":"hay_pasture","4":"0.0095181263"},{"1":"Sharp-shinned Hawk","2":"OR","3":"planted.cultivated","4":"0.0095181263"},{"1":"Sharp-shinned Hawk","2":"OR","3":"developed_low_intensity","4":"0.0064897628"},{"1":"Sharp-shinned Hawk","2":"OR","3":"developed_medium_intensity","4":"0.0027075687"},{"1":"Sharp-shinned Hawk","2":"OR","3":"forest","4":"0.0014478410"},{"1":"Sharp-shinned Hawk","2":"OR","3":"tmin","4":"0.0007239205"},{"1":"Sharp-shinned Hawk","2":"OR","3":"Winter_NDVI","4":"0.0000000000"},{"1":"Wild Turkey","2":"OR","3":"tmax","4":"1.0000000000"},{"1":"Wild Turkey","2":"OR","3":"temp.diff","4":"0.7888669111"},{"1":"Wild Turkey","2":"OR","3":"avg_prcp","4":"0.7294564907"},{"1":"Wild Turkey","2":"OR","3":"lat","4":"0.6163529120"},{"1":"Wild Turkey","2":"OR","3":"tmin","4":"0.5983447553"},{"1":"Wild Turkey","2":"OR","3":"lat.sqrd","4":"0.5798323733"},{"1":"Wild Turkey","2":"OR","3":"Winter_NDVI","4":"0.5517444440"},{"1":"Wild Turkey","2":"OR","3":"Fall_NDVI","4":"0.5194155447"},{"1":"Wild Turkey","2":"OR","3":"min.ndvi","4":"0.3999314458"},{"1":"Wild Turkey","2":"OR","3":"dem","4":"0.3858568138"},{"1":"Wild Turkey","2":"OR","3":"lon","4":"0.2991719689"},{"1":"Wild Turkey","2":"OR","3":"dem.detrended","4":"0.2580743733"},{"1":"Wild Turkey","2":"OR","3":"Spring_NDVI","4":"0.2566058777"},{"1":"Wild Turkey","2":"OR","3":"lon.sqrd","4":"0.2439830335"},{"1":"Wild Turkey","2":"OR","3":"water","4":"0.0750517086"},{"1":"Wild Turkey","2":"OR","3":"lon.lat.interaction","4":"0.0613636354"},{"1":"Wild Turkey","2":"OR","3":"max.ndvi","4":"0.0456741246"},{"1":"Wild Turkey","2":"OR","3":"waterbody","4":"0.0448859091"},{"1":"Wild Turkey","2":"OR","3":"developed","4":"0.0419888868"},{"1":"Wild Turkey","2":"OR","3":"Summer_NDVI","4":"0.0367842782"},{"1":"Wild Turkey","2":"OR","3":"urban_imperviousness","4":"0.0115785470"},{"1":"Wild Turkey","2":"OR","3":"forest","4":"0.0022753070"},{"1":"Wild Turkey","2":"OR","3":"open_water","4":"0.0014010842"},{"1":"Wild Turkey","2":"OR","3":"barren","4":"0.0000000000"},{"1":"Belted Kingfisher","2":"VT","3":"tmin","4":"1.0000000000"},{"1":"Belted Kingfisher","2":"VT","3":"urban_imperviousness","4":"0.8930344029"},{"1":"Belted Kingfisher","2":"VT","3":"temp.diff","4":"0.8659338690"},{"1":"Belted Kingfisher","2":"VT","3":"developed","4":"0.7731502112"},{"1":"Belted Kingfisher","2":"VT","3":"lon","4":"0.7160923521"},{"1":"Belted Kingfisher","2":"VT","3":"lon.sqrd","4":"0.7160923521"},{"1":"Belted Kingfisher","2":"VT","3":"avg_prcp","4":"0.5776997692"},{"1":"Belted Kingfisher","2":"VT","3":"dem","4":"0.4611335639"},{"1":"Belted Kingfisher","2":"VT","3":"tmax","4":"0.4112538866"},{"1":"Belted Kingfisher","2":"VT","3":"lat","4":"0.3339166776"},{"1":"Belted Kingfisher","2":"VT","3":"lat.sqrd","4":"0.3339166776"},{"1":"Belted Kingfisher","2":"VT","3":"lon.lat.interaction","4":"0.2152920781"},{"1":"Belted Kingfisher","2":"VT","3":"Spring_NDVI","4":"0.2048251671"},{"1":"Belted Kingfisher","2":"VT","3":"max.ndvi","4":"0.2036154214"},{"1":"Belted Kingfisher","2":"VT","3":"Fall_NDVI","4":"0.2035882186"},{"1":"Belted Kingfisher","2":"VT","3":"min.ndvi","4":"0.2030763946"},{"1":"Belted Kingfisher","2":"VT","3":"Winter_NDVI","4":"0.0783913438"},{"1":"Belted Kingfisher","2":"VT","3":"Summer_NDVI","4":"0.0265098572"},{"1":"Belted Kingfisher","2":"VT","3":"waterbody","4":"0.0187058240"},{"1":"Belted Kingfisher","2":"VT","3":"dem.detrended","4":"0.0175313754"},{"1":"Belted Kingfisher","2":"VT","3":"emergent_herbaceous_wetlands","4":"0.0156026790"},{"1":"Belted Kingfisher","2":"VT","3":"water","4":"0.0036965130"},{"1":"Belted Kingfisher","2":"VT","3":"herbaceous","4":"0.0000000000"},{"1":"Cedar Waxwing","2":"VT","3":"lon","4":"1.0000000000"},{"1":"Cedar Waxwing","2":"VT","3":"dem","4":"0.7973696964"},{"1":"Cedar Waxwing","2":"VT","3":"lon.sqrd","4":"0.7531549056"},{"1":"Cedar Waxwing","2":"VT","3":"Spring_NDVI","4":"0.5855109007"},{"1":"Cedar Waxwing","2":"VT","3":"Winter_NDVI","4":"0.5647742936"},{"1":"Cedar Waxwing","2":"VT","3":"min.ndvi","4":"0.5413023309"},{"1":"Cedar Waxwing","2":"VT","3":"avg_prcp","4":"0.4730540191"},{"1":"Cedar Waxwing","2":"VT","3":"temp.diff","4":"0.4487769623"},{"1":"Cedar Waxwing","2":"VT","3":"Summer_NDVI","4":"0.3826943526"},{"1":"Cedar Waxwing","2":"VT","3":"tmin","4":"0.2865894877"},{"1":"Cedar Waxwing","2":"VT","3":"max.ndvi","4":"0.2822429681"},{"1":"Cedar Waxwing","2":"VT","3":"Fall_NDVI","4":"0.2768346721"},{"1":"Cedar Waxwing","2":"VT","3":"tmax","4":"0.2441940337"},{"1":"Cedar Waxwing","2":"VT","3":"dem.detrended","4":"0.1907982944"},{"1":"Cedar Waxwing","2":"VT","3":"lon.lat.interaction","4":"0.1799635362"},{"1":"Cedar Waxwing","2":"VT","3":"waterbody","4":"0.1238827613"},{"1":"Cedar Waxwing","2":"VT","3":"developed_low_intensity","4":"0.0958934367"},{"1":"Cedar Waxwing","2":"VT","3":"lat","4":"0.0690363940"},{"1":"Cedar Waxwing","2":"VT","3":"water","4":"0.0557120776"},{"1":"Cedar Waxwing","2":"VT","3":"emergent_herbaceous_wetlands","4":"0.0512904478"},{"1":"Cedar Waxwing","2":"VT","3":"developed","4":"0.0276894579"},{"1":"Cedar Waxwing","2":"VT","3":"urban_imperviousness","4":"0.0239917861"},{"1":"Cedar Waxwing","2":"VT","3":"forest","4":"0.0160943370"},{"1":"Cedar Waxwing","2":"VT","3":"cultivated_crops","4":"0.0000000000"},{"1":"Downy Woodpecker","2":"VT","3":"dem","4":"1.0000000000"},{"1":"Downy Woodpecker","2":"VT","3":"Winter_NDVI","4":"0.7557493356"},{"1":"Downy Woodpecker","2":"VT","3":"Fall_NDVI","4":"0.7443769655"},{"1":"Downy Woodpecker","2":"VT","3":"max.ndvi","4":"0.7340272585"},{"1":"Downy Woodpecker","2":"VT","3":"tmax","4":"0.7044135992"},{"1":"Downy Woodpecker","2":"VT","3":"temp.diff","4":"0.6857273676"},{"1":"Downy Woodpecker","2":"VT","3":"avg_prcp","4":"0.6832801543"},{"1":"Downy Woodpecker","2":"VT","3":"tmin","4":"0.6337457942"},{"1":"Downy Woodpecker","2":"VT","3":"lon","4":"0.6179708548"},{"1":"Downy Woodpecker","2":"VT","3":"dem.detrended","4":"0.6175059035"},{"1":"Downy Woodpecker","2":"VT","3":"lon.sqrd","4":"0.5402907294"},{"1":"Downy Woodpecker","2":"VT","3":"Summer_NDVI","4":"0.3989386236"},{"1":"Downy Woodpecker","2":"VT","3":"urban_imperviousness","4":"0.2895786211"},{"1":"Downy Woodpecker","2":"VT","3":"lat","4":"0.2676125534"},{"1":"Downy Woodpecker","2":"VT","3":"water","4":"0.2323331454"},{"1":"Downy Woodpecker","2":"VT","3":"lat.sqrd","4":"0.1899324280"},{"1":"Downy Woodpecker","2":"VT","3":"developed_high_intensity","4":"0.1164183567"},{"1":"Downy Woodpecker","2":"VT","3":"lon.lat.interaction","4":"0.1010282362"},{"1":"Downy Woodpecker","2":"VT","3":"developed","4":"0.0884527721"},{"1":"Downy Woodpecker","2":"VT","3":"waterbody","4":"0.0735238680"},{"1":"Downy Woodpecker","2":"VT","3":"min.ndvi","4":"0.0619516288"},{"1":"Downy Woodpecker","2":"VT","3":"herbaceous","4":"0.0600599599"},{"1":"Downy Woodpecker","2":"VT","3":"Spring_NDVI","4":"0.0484708903"},{"1":"Downy Woodpecker","2":"VT","3":"developed_open_space","4":"0.0156853465"},{"1":"Downy Woodpecker","2":"VT","3":"developed_low_intensity","4":"0.0096863017"},{"1":"Downy Woodpecker","2":"VT","3":"open_water","4":"0.0018522644"},{"1":"Downy Woodpecker","2":"VT","3":"developed_medium_intensity","4":"0.0000000000"},{"1":"Ruddy Duck","2":"VT","3":"water","4":"1.0000000000"},{"1":"Ruddy Duck","2":"VT","3":"waterbody","4":"0.9975372711"},{"1":"Ruddy Duck","2":"VT","3":"max.ndvi","4":"0.3327821422"},{"1":"Ruddy Duck","2":"VT","3":"Summer_NDVI","4":"0.2367679928"},{"1":"Ruddy Duck","2":"VT","3":"dem","4":"0.2278378307"},{"1":"Ruddy Duck","2":"VT","3":"Fall_NDVI","4":"0.2029515602"},{"1":"Ruddy Duck","2":"VT","3":"lon.lat.interaction","4":"0.1958168400"},{"1":"Ruddy Duck","2":"VT","3":"avg_prcp","4":"0.1823489467"},{"1":"Ruddy Duck","2":"VT","3":"dem.detrended","4":"0.1665076825"},{"1":"Ruddy Duck","2":"VT","3":"lat","4":"0.1647373813"},{"1":"Ruddy Duck","2":"VT","3":"lat.sqrd","4":"0.1561312642"},{"1":"Ruddy Duck","2":"VT","3":"min.ndvi","4":"0.1172155885"},{"1":"Ruddy Duck","2":"VT","3":"lon","4":"0.0998497648"},{"1":"Ruddy Duck","2":"VT","3":"lon.sqrd","4":"0.0154275157"},{"1":"Ruddy Duck","2":"VT","3":"temp.diff","4":"0.0128897541"},{"1":"Ruddy Duck","2":"VT","3":"Winter_NDVI","4":"0.0062693014"},{"1":"Ruddy Duck","2":"VT","3":"Spring_NDVI","4":"0.0052764692"},{"1":"Ruddy Duck","2":"VT","3":"open_water","4":"0.0000000000"},{"1":"Sanderling","2":"VT","3":"Spring_NDVI","4":"1.0000000000"},{"1":"Sanderling","2":"VT","3":"Fall_NDVI","4":"0.4332540361"},{"1":"Sanderling","2":"VT","3":"max.ndvi","4":"0.4291567545"},{"1":"Sanderling","2":"VT","3":"min.ndvi","4":"0.4277909940"},{"1":"Sanderling","2":"VT","3":"Summer_NDVI","4":"0.2833729820"},{"1":"Sanderling","2":"VT","3":"dem","4":"0.2724468978"},{"1":"Sanderling","2":"VT","3":"Winter_NDVI","4":"0.0000000000"},{"1":"Sandhill Crane","2":"VT","3":"temp.diff","4":"1.0000000000"},{"1":"Sandhill Crane","2":"VT","3":"urban_imperviousness","4":"0.5795669962"},{"1":"Sandhill Crane","2":"VT","3":"water","4":"0.5384384527"},{"1":"Sandhill Crane","2":"VT","3":"waterbody","4":"0.5006780728"},{"1":"Sandhill Crane","2":"VT","3":"lon.lat.interaction","4":"0.4534210468"},{"1":"Sandhill Crane","2":"VT","3":"lat","4":"0.4244407253"},{"1":"Sandhill Crane","2":"VT","3":"lat.sqrd","4":"0.4217163171"},{"1":"Sandhill Crane","2":"VT","3":"Fall_NDVI","4":"0.3631758856"},{"1":"Sandhill Crane","2":"VT","3":"max.ndvi","4":"0.3631758856"},{"1":"Sandhill Crane","2":"VT","3":"developed","4":"0.2797761185"},{"1":"Sandhill Crane","2":"VT","3":"dem.detrended","4":"0.1802732528"},{"1":"Sandhill Crane","2":"VT","3":"Spring_NDVI","4":"0.1754423530"},{"1":"Sandhill Crane","2":"VT","3":"min.ndvi","4":"0.1754423530"},{"1":"Sandhill Crane","2":"VT","3":"Summer_NDVI","4":"0.1305268416"},{"1":"Sandhill Crane","2":"VT","3":"Winter_NDVI","4":"0.0834450028"},{"1":"Sandhill Crane","2":"VT","3":"dem","4":"0.0231231618"},{"1":"Sandhill Crane","2":"VT","3":"tmin","4":"0.0039029548"},{"1":"Sandhill Crane","2":"VT","3":"lon","4":"0.0029528480"},{"1":"Sandhill Crane","2":"VT","3":"avg_prcp","4":"0.0026127936"},{"1":"Sandhill Crane","2":"VT","3":"lon.sqrd","4":"0.0000000000"},{"1":"Sharp-shinned Hawk","2":"VT","3":"Summer_NDVI","4":"1.0000000000"},{"1":"Sharp-shinned Hawk","2":"VT","3":"tmax","4":"0.8078641027"},{"1":"Sharp-shinned Hawk","2":"VT","3":"Winter_NDVI","4":"0.7271529028"},{"1":"Sharp-shinned Hawk","2":"VT","3":"waterbody","4":"0.7061891659"},{"1":"Sharp-shinned Hawk","2":"VT","3":"water","4":"0.6983617633"},{"1":"Sharp-shinned Hawk","2":"VT","3":"lat","4":"0.6182690849"},{"1":"Sharp-shinned Hawk","2":"VT","3":"avg_prcp","4":"0.2576599569"},{"1":"Sharp-shinned Hawk","2":"VT","3":"dem","4":"0.2482900117"},{"1":"Sharp-shinned Hawk","2":"VT","3":"dem.detrended","4":"0.2481208563"},{"1":"Sharp-shinned Hawk","2":"VT","3":"lon","4":"0.1858817991"},{"1":"Sharp-shinned Hawk","2":"VT","3":"Spring_NDVI","4":"0.1855929402"},{"1":"Sharp-shinned Hawk","2":"VT","3":"min.ndvi","4":"0.1844801022"},{"1":"Sharp-shinned Hawk","2":"VT","3":"lon.lat.interaction","4":"0.1438018689"},{"1":"Sharp-shinned Hawk","2":"VT","3":"Fall_NDVI","4":"0.1360506606"},{"1":"Sharp-shinned Hawk","2":"VT","3":"temp.diff","4":"0.0895641313"},{"1":"Sharp-shinned Hawk","2":"VT","3":"tmin","4":"0.0773644711"},{"1":"Sharp-shinned Hawk","2":"VT","3":"lon.sqrd","4":"0.0720290027"},{"1":"Sharp-shinned Hawk","2":"VT","3":"lat.sqrd","4":"0.0402689210"},{"1":"Sharp-shinned Hawk","2":"VT","3":"developed","4":"0.0148988555"},{"1":"Sharp-shinned Hawk","2":"VT","3":"open_water","4":"0.0103794882"},{"1":"Sharp-shinned Hawk","2":"VT","3":"wetlands","4":"0.0063880546"},{"1":"Sharp-shinned Hawk","2":"VT","3":"woody_wetlands","4":"0.0063880546"},{"1":"Sharp-shinned Hawk","2":"VT","3":"emergent_herbaceous_wetlands","4":"0.0005043652"},{"1":"Sharp-shinned Hawk","2":"VT","3":"urban_imperviousness","4":"0.0000000000"},{"1":"Wild Turkey","2":"VT","3":"lat","4":"1.0000000000"},{"1":"Wild Turkey","2":"VT","3":"lat.sqrd","4":"0.9806686806"},{"1":"Wild Turkey","2":"VT","3":"dem","4":"0.9321054892"},{"1":"Wild Turkey","2":"VT","3":"Fall_NDVI","4":"0.5479245846"},{"1":"Wild Turkey","2":"VT","3":"dem.detrended","4":"0.4795097420"},{"1":"Wild Turkey","2":"VT","3":"temp.diff","4":"0.4429678007"},{"1":"Wild Turkey","2":"VT","3":"lon","4":"0.4340015893"},{"1":"Wild Turkey","2":"VT","3":"Spring_NDVI","4":"0.4144041885"},{"1":"Wild Turkey","2":"VT","3":"Winter_NDVI","4":"0.3583332364"},{"1":"Wild Turkey","2":"VT","3":"tmax","4":"0.3120291567"},{"1":"Wild Turkey","2":"VT","3":"avg_prcp","4":"0.2595023048"},{"1":"Wild Turkey","2":"VT","3":"tmin","4":"0.2060927143"},{"1":"Wild Turkey","2":"VT","3":"max.ndvi","4":"0.1809206765"},{"1":"Wild Turkey","2":"VT","3":"lon.lat.interaction","4":"0.1756924013"},{"1":"Wild Turkey","2":"VT","3":"forest","4":"0.1548397013"},{"1":"Wild Turkey","2":"VT","3":"min.ndvi","4":"0.1480134964"},{"1":"Wild Turkey","2":"VT","3":"Summer_NDVI","4":"0.1304238552"},{"1":"Wild Turkey","2":"VT","3":"open_water","4":"0.1238162330"},{"1":"Wild Turkey","2":"VT","3":"water","4":"0.1229245708"},{"1":"Wild Turkey","2":"VT","3":"lon.sqrd","4":"0.1067784411"},{"1":"Wild Turkey","2":"VT","3":"urban_imperviousness","4":"0.0299649433"},{"1":"Wild Turkey","2":"VT","3":"cultivated_crops","4":"0.0130296673"},{"1":"Wild Turkey","2":"VT","3":"developed","4":"0.0024898771"},{"1":"Wild Turkey","2":"VT","3":"developed_high_intensity","4":"0.0000000000"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
</script>
</div>

## Conclusion

Predictive modeling is iterative by nature, and thus conclusively selecting features to be included in the modeling process is impractical. The purpose of this exploratory analysis is just that - exploratory. The insights and results of this part of the study will be utilized in the iterative modeling process, and adjusted as needed.