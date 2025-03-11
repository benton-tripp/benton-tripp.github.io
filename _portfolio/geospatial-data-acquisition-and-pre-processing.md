---
layout: page
order: 6
title: Geospatial Data Acquisition and Pre-Processing
permalink: /portfolio/geospatial-data-acquisition-and-pre-processing
---

## Problem

Effective geospatial analysis requires the integration of diverse datasets—from satellite imagery and environmental rasters to observational records and official boundaries—into a unified framework. In my projects, I faced the challenge of acquiring raw data from multiple sources and transforming it into standardized, high-quality datasets ready for analysis. Whether preparing bird observation data for species distribution models or generating consistent spatial tiles for urban growth monitoring, the core objective was to ensure that all data layers share a common coordinate system, resolution, and extent.

This work addresses not only the technical difficulties inherent in managing heterogeneous datasets but also the need for reproducibility and scalability in geospatial projects. By establishing rigorous data acquisition and pre-processing workflows, I have created the foundation necessary to support advanced spatial analyses and modeling efforts across varied study areas.

## Analysis Procedures

My approach to geospatial data management begins with sourcing data from reliable, diverse repositories—ranging from official state boundaries and remote sensing imagery to crowd-sourced observation records and weather data. The acquired datasets are first cleaned and filtered to remove errors and inconsistencies, ensuring that only high-quality data moves forward in the workflow.

Next, I standardize the data by reprojecting all layers into a common coordinate reference system, resampling raster data to uniform resolutions, and masking or cropping them using consistent spatial boundaries. For instance, state boundary data are segmented into individual units for precise masking, while multispectral imagery and environmental layers are mosaicked and aligned to common extents. These procedures guarantee that disparate datasets can be accurately overlaid and analyzed together, paving the way for subsequent modeling and spatial analysis.

## Project Examples

**North Carolina Woodpecker Distribution Modeling**  

*(April, 2023; see [ArcGIS Python Workflow for Woodpecker Distribution in North Carolina](/portfolio/arcgis-python-workflow-for-woodpecker-distribution))* 

One specific example of complex data acquisition in this project is the comprehensive workflow developed for processing the FeederWatch dataset. The process begins by verifying whether a local CSV containing species codes already exists; if not, an Excel file is downloaded from a public Google Drive link, the “Species Codes” sheet is extracted, and key columns are renamed (e.g., `SPECIES_CODE` becomes `species_code`). Simultaneously, raw observation data are standardized by converting all field names to lowercase and ensuring that expected columns are present (assigning NaN where necessary). The dataset is then filtered to include only valid records, with additional criteria applied to exclude observations marked with plus codes and to retain only the targeted bird species. This integrated approach guarantees that both species metadata and observation records are consistent and ready for modeling.

**Species Distribution Modeling Benchmarks**  

*(December, 2023; see [Species Distribution Model Benchmark Study](/portfolio/species-distribution-model-benchmark-study))*  

A notable example in the SDM Benchmark Study involves the aggregation of weather raster data. The workflow downloads annual weather datasets at a 4km resolution for 2017–2019 and 30-year monthly normals at an 800m resolution. For each weather variable—such as precipitation, maximum temperature, and minimum temperature—the annual rasters are aggregated using appropriate statistical functions (e.g., maximum for `tmax`, minimum for `tmin`, mean for precipitation). These aggregated datasets are then reprojected to a consistent coordinate system and resampled to harmonize spatial resolution. Finally, a weighted average method is applied to combine the 4km and 800m data, where the weights reflect the relative importance and resolution differences between the two datasets. This detailed and multi-step process ensures that the environmental variables used in model benchmarking are accurately represented and spatially aligned.

**Utilizing the Clay Foundation Model and Sentinel-2 Imagery for Urban Growth Monitoring in Johnston County, North Carolina**  

*(December, 2024; see [Utilizing the Clay Foundation Model and Sentinel-2 Imagery for Urban Growth Monitoring in Johnston County, North Carolina](/portfolio/utilizing-the-clay-foundation-model))*  

In this project, one complex pre-processing example involves generating a uniform grid of 600×600-meter tiles based on the official county boundary, which is reprojected to UTM Zone 17N (EPSG:32617). Each tile is used as a spatial unit to crop high-resolution (10m) Sentinel-2 imagery, where the images are normalized using band-specific means and standard deviations to prepare them for deep learning analysis. In parallel, urban imperviousness data from the NLCD—originally at a 30-meter resolution—is reprojected and resampled to 200 meters using bilinear interpolation. The alignment and integration of these disparate datasets allow for the creation of consistent datacubes for each tile, which are then fed into the pretrained Clay Foundation Model to extract spatial embeddings that capture both spectral and contextual information essential for urban growth analysis.

<img src="{{ site.baseurl }}/assets/images/data-management-img-3.png" style="width:80%; max-width: 700px; min-width: 300px;">

*Figure 1: A 600×600-meter tile grid over Johnston County, North Carolina, used for systematically retrieving and processing Sentinel-2 imagery for urban growth analysis. Each tile defines a spatial unit for consistent data acquisition and subsequent deep learning integration.*

## Results

The outcome of these data management workflows is a suite of robust, standardized datasets that seamlessly integrate into my geospatial analyses. In one project, processed bird observation data were combined with environmental rasters to produce a consistent input for species distribution models. In another, a high-resolution, tile-based dataset was created from Sentinel-2 imagery and urban imperviousness data, enabling detailed monitoring of urban growth patterns.

These harmonized datasets not only enhance the accuracy and reliability of subsequent analyses but also streamline the entire modeling process. By ensuring that each data layer is consistent in terms of CRS, resolution, and spatial extent, the pre-processing workflows reduce the likelihood of errors and allow for more efficient and reproducible geospatial analysis.

## Reflection

This experience has shown the critical importance of rigorous data acquisition and pre-processing in geospatial projects. Through the systematic standardization and integration of different data sources, I have developed workflows that significantly improve the quality and reproducibility of spatial analyses. Managing complex datasets—from raw satellite imagery to intricate observational records—requires careful planning and execution, and the benefits are evident in the robustness of the final analytical outputs. By mastering these data management techniques, I have not only enhanced my technical capabilities but also positioned myself to tackle more complex spatial problems. This foundation in data acquisition and pre-processing is essential for advanced geospatial analysis and modeling, reinforcing the value of strong data management practices in both academic research and practical applications.

*Goals*: In the short term, I plan to expand my data validation scripts to automate error-checking across multiple repositories, while in the long term, I aim to incorporate machine learning-based data cleaning methods that can handle increasingly diverse and voluminous geospatial datasets.



