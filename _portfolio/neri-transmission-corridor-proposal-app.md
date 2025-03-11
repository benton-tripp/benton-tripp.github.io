---
layout: page
order: 7
title: NERI Transmission Corridor Proposal App
permalink: /portfolio/neri-transmission-corridor-proposal-app
---

*December, 2024*

## Problem

The primary objective of this project was to develop a GIS-based workflow for efficiently managing and analyzing proposed transmission line corridors within the New River Gorge National Park and Preserve (NERI). Planners required a robust system capable of evaluating environmental and cultural impacts of proposed infrastructure projects. The solution needed to integrate diverse spatial and non-spatial datasets, streamline geospatial analysis, and offer clear visual and tabular outputs to support informed decision-making while ensuring compliance with conservation priorities.

## Analysis Procedures

A PostgreSQL geospatial database was established, consolidating spatial data (e.g., vertebrate observations, culturally significant buildings, roads, trails, existing utilities) and associated non-spatial information such as sampling events and species metadata. Data processing leveraged ArcGIS Pro, including spatial joins, buffer creation, and attribute calculations. ArcGIS Server hosted customized geoprocessing services to evaluate intersections between proposed transmission line corridors and selected environmental or cultural resources, classifying each proposed segment as "Accepted" or "Rejected."

**Backend Processing and Analysis**

The backend workflow utilized Python scripting within ArcGIS Model Builder to automate critical geospatial analyses. Three distinct models supported this process:

- **Delete Proposed Lines:** Removes all existing proposed lines and associated dependent datasets to provide a clean environment for new analyses.
- **Reset Buffers and Proposal Statuses:** Clears buffer analyses and resets proposal statuses, facilitating iterative adjustments without restarting the entire workflow.
- **Process Proposal Status:** Accepts user-defined parameters such as buffer distance, species habitats, and building contribution statuses to conduct spatial intersections and update line acceptance statuses accordingly.

**Frontend Interface and Interactivity**

A responsive web application was developed using JavaScript and the ArcGIS JavaScript SDK. This frontend interface provided:

- **Interactive Map Layers:** Allowing users to visualize and selectively display environmental and infrastructural data.
- **Dynamic Proposal Editing:** Integrated editor widgets enable users to create, modify, or remove proposed transmission line segments interactively.
- **Real-time Buffer Analysis:** Interactive inputs enabled users to run analyses directly from the interface, dynamically updating results.
- **Data Exploration Tables:** Interactive data tables featuring dropdown filters and sorting capabilities were included, enabling users to easily navigate and explore detailed records for species observations and building intersections.
- **Map Exporting Functionality:** Users could generate and export customized maps in PDF or PNG formats for documentation purposes, further supporting reporting needs.

<img src="{{ site.baseurl }}/assets/images/web-services-img-1.png" style="width: 80%; max-width: 600px; min-width: 300px;">

*Figure 1: Interactive web-based geospatial interface demonstrating analysis capabilities, including selection of geoprocessing parameters such as buffer radius, species habitat type, and building contribution status for proposed transmission lines within NERI.*

<img src="{{ site.baseurl }}/assets/images/web-services-img-2.png" style="width: 80%; max-width: 600px; min-width: 300px;">

*Figure 2: Exported map illustrating the results of a geoprocessing analysis displaying proposed transmission lines, buffered areas indicating acceptance or rejection based on intersections with aquatic, terrestrial, and combined species observations within the New River Gorge National Park and Preserve (NERI).*

## Results

The resulting web-based application provides interactive mapping capabilities, enabling users to:

- Create, edit, and delete proposed transmission lines directly within the map interface.
- Conduct buffer analyses with adjustable distances, species habitat, and building contribution criteria.
- Clearly visualize analysis outcomes with color-coded buffers (green for Accepted, dark-red for Rejected).
- View detailed interactive tables of species observations and buildings within corridors, enhanced with sortable columns and dropdown filters for efficient data exploration.
- Export customizable map outputs in PDF or PNG formats for reporting and documentation purposes.

These functionalities significantly streamline the analysis process, supporting iterative refinement and facilitating comprehensive environmental and cultural assessments.

## Reflection

Developing this application underscored the importance of a tightly integrated database and GIS architecture, particularly the dynamic interaction between spatial analysis and user-friendly interfaces. Challenges included optimizing geoprocessing workflows for responsiveness, effectively managing complex intersecting datasets, and ensuring clarity in user-driven analysis outcomes. Future enhancements might include incorporating real-time data feeds or expanding the analysis to other infrastructure types, thus further extending the application's utility and adaptability.
