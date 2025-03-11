---
layout: page
order: 5
title: Managing Data in an Enterprise GeoDatabase
permalink: /portfolio/managing-data-in-an-enterprise-geodatabase
---

*December, 2024*

## Problem

Planning and conservation efforts within the New River Gorge National Park and Preserve (NERI) demanded a centralized and efficient system for handling both spatial and non-spatial data related to proposed transmission line corridors. This system needed to integrate datasets such as species observations, building attributes, park boundaries, and utility infrastructure into a single repository. By deploying an enterprise geodatabase in PostgreSQL, users could simultaneously ensure data integrity, maintain consistent relationships across multiple tables, and streamline interactions between ArcGIS Pro and a custom web application.

For additional details on the broader application (e.g., web services and geoprocessing workflows), see the NERI Transmission Corridor Proposal App.

## Analysis Procedures

**Database Schema Design**

- **Table Creation and Constraints:**
  - Key tables—such as `tbl_data`, `tbl_locations`, `tbl_events`, and `lu_species`—were defined using SQL CREATE TABLE statements with primary keys, foreign keys, and check constraints. This enforced referential integrity and prevented orphaned records.
  - Spatial feature classes (e.g., `fc_proposed_lines`, `fc_neri_buildings`) were managed via ArcGIS Pro, ensuring geometry fields and projection information remained consistent.

- **Views for Data Normalization and Filtering:**
  - `v_all_observations`: Merged raw observation data (`tbl_data`), location coordinates, event information, and species details into a unified record set. Data cleansing logic—such as mapping habitat strings to standardized values—was embedded in SQL CASE statements.
  - `v_transmission_lines`: Filtered the base utilities feature class to include only power transmission lines, simplifying subsequent geoprocessing tasks and map rendering.

**Data Loading and Maintenance**

- **Bulk Inserts and Updates:**
  - CSV files (e.g., `tbl_data.txt`, `tbl_locations.txt`) were loaded into PostgreSQL using SQL COPY statements. This facilitated large-scale data ingestion and minimized manual errors.
  - Foreign key constraints (e.g., `DataHaveSpecies`) ensured that each observation record referenced valid species IDs, preventing mismatched or duplicate entries.

- **Spatial Data Integration with ArcGIS Pro:**
  - Feature classes derived from the SQL views (e.g., `fc_all_observations`) allowed immediate visualization and editing within ArcGIS.
  - ArcGIS tools were used to project, dissolve, or clip data as needed, then push updates back into the PostgreSQL database—preserving a single source of truth for all spatial assets.

**Enterprise-Wide Access and Geoprocessing**

- **ArcGIS Server Services:**
  - Custom geoprocessing models (e.g., `Process Proposal Status`) interacted with the PostgreSQL data via REST endpoints. The models updated or reset records (such as buffer distances or proposal statuses) directly in the enterprise geodatabase.

- **Web Application Integration:**
  - A JavaScript front-end allowed end-users to propose or delete lines, which triggered updates to `fc_proposed_lines` and related tables. This real-time exchange of data was possible because of the well-structured schema and REST-based communication layer.

<img src="{{ site.baseurl }}/assets/images/NERI-uml.png" style="width:80%; max-width:1000px; min-width:300px;">

*Figure 1: UML Diagram portraying data structure, types, and relationships.*

<table style="border-collapse: collapse; width: 100%; font-family: Arial, sans-serif; font-size: 14px; color: #333;">
  <thead>
    <tr style="background-color: #aaa; text-align: left;">
      <!-- Keep the header in the regular font -->
      <th style="border: 1px solid #ccc; padding: 8px; white-space: nowrap; min-width: 150px;">Table / View</th>
      <th style="border: 1px solid #ccc; padding: 8px; white-space: nowrap; min-width: 150px;">Spatial / Non-Spatial</th>
      <th style="border: 1px solid #ccc; padding: 8px;">Details</th>
    </tr>
  </thead>
  <tbody>
    <!-- Row 1 -->
    <tr style="background-color: #ffffff;">
      <td style="border: 1px solid #ccc; padding: 8px; font-family: 'Courier New', Courier, monospace;">fc_all_observations</td>
      <td style="border: 1px solid #ccc; padding: 8px;">Spatial</td>
      <td style="border: 1px solid #ccc; padding: 8px;">Derived in ArcGIS Pro from the <span style="font-family: 'Courier New', Courier, monospace;">v_all_observations</span> view</td>
    </tr>
    <!-- Row 2 -->
    <tr style="background-color: #fafafa;">
      <td style="border: 1px solid #ccc; padding: 8px; font-family: 'Courier New', Courier, monospace;">fc_intersecting_buildings</td>
      <td style="border: 1px solid #ccc; padding: 8px;">Spatial</td>
      <td style="border: 1px solid #ccc; padding: 8px;">Derived in ArcGIS Pro from the buildings with user-selected contribution status filters that intersect with the proposed corridors.</td>
    </tr>
    <!-- Row 3 -->
    <tr style="background-color: #ffffff;">
      <td style="border: 1px solid #ccc; padding: 8px; font-family: 'Courier New', Courier, monospace;">fc_intersecting_species</td>
      <td style="border: 1px solid #ccc; padding: 8px;">Spatial</td>
      <td style="border: 1px solid #ccc; padding: 8px;">Derived in ArcGIS Pro from the species observation point data with user-selected habitat filters that intersect with the proposed corridors.</td>
    </tr>
    <!-- Row 4 -->
    <tr style="background-color: #fafafa;">
      <td style="border: 1px solid #ccc; padding: 8px; font-family: 'Courier New', Courier, monospace;">fc_neri_buildings</td>
      <td style="border: 1px solid #ccc; padding: 8px;">Spatial</td>
      <td style="border: 1px solid #ccc; padding: 8px;">NERI base data of buildings added to the database via ArcGIS Pro.</td>
    </tr>
    <!-- Row 5 -->
    <tr style="background-color: #ffffff;">
      <td style="border: 1px solid #ccc; padding: 8px; font-family: 'Courier New', Courier, monospace;">fc_park_boundary</td>
      <td style="border: 1px solid #ccc; padding: 8px;">Spatial</td>
      <td style="border: 1px solid #ccc; padding: 8px;">NERI base data of the park boundary added to the database via ArcGIS Pro.</td>
    </tr>
    <!-- Row 6 -->
    <tr style="background-color: #fafafa;">
      <td style="border: 1px solid #ccc; padding: 8px; font-family: 'Courier New', Courier, monospace;">fc_proposal_boundary</td>
      <td style="border: 1px solid #ccc; padding: 8px;">Spatial</td>
      <td style="border: 1px solid #ccc; padding: 8px;">Derived in ArcGIS Pro as the bounding box of all of the spatial data used in this project; Used as a reference for a “Proposal Area”.</td>
    </tr>
    <!-- Row 7 -->
    <tr style="background-color: #ffffff;">
      <td style="border: 1px solid #ccc; padding: 8px; font-family: 'Courier New', Courier, monospace;">fc_proposed_lines</td>
      <td style="border: 1px solid #ccc; padding: 8px;">Spatial</td>
      <td style="border: 1px solid #ccc; padding: 8px;">User-defined feature class representing proposed transmission lines; created in the web application and created/edited/deleted from the database interactively and through geoprocessing services.</td>
    </tr>
    <!-- Row 8 -->
    <tr style="background-color: #fafafa;">
      <td style="border: 1px solid #ccc; padding: 8px; font-family: 'Courier New', Courier, monospace;">fc_proposed_lines_buffer</td>
      <td style="border: 1px solid #ccc; padding: 8px;">Spatial</td>
      <td style="border: 1px solid #ccc; padding: 8px;">User-defined feature class representing the corridor around the proposed transmission lines; created/edited/deleted from the database using geoprocessing services.</td>
    </tr>
    <!-- Row 9 -->
    <tr style="background-color: #ffffff;">
      <td style="border: 1px solid #ccc; padding: 8px; font-family: 'Courier New', Courier, monospace;">fc_railroads</td>
      <td style="border: 1px solid #ccc; padding: 8px;">Spatial</td>
      <td style="border: 1px solid #ccc; padding: 8px;">NERI base data of the railroads in the region, added to the database via ArcGIS Pro.</td>
    </tr>
    <!-- Row 10 -->
    <tr style="background-color: #fafafa;">
      <td style="border: 1px solid #ccc; padding: 8px; font-family: 'Courier New', Courier, monospace;">fc_roads_and_trails</td>
      <td style="border: 1px solid #ccc; padding: 8px;">Spatial</td>
      <td style="border: 1px solid #ccc; padding: 8px;">NERI base data of roads and trails in the region, added to the database via ArcGIS Pro.</td>
    </tr>
    <!-- Row 11 -->
    <tr style="background-color: #ffffff;">
      <td style="border: 1px solid #ccc; padding: 8px; font-family: 'Courier New', Courier, monospace;">fc_streams</td>
      <td style="border: 1px solid #ccc; padding: 8px;">Spatial</td>
      <td style="border: 1px solid #ccc; padding: 8px;">NERI base data of streams in the region, added to the database via ArcGIS Pro.</td>
    </tr>
    <!-- Row 12 -->
    <tr style="background-color: #fafafa;">
      <td style="border: 1px solid #ccc; padding: 8px; font-family: 'Courier New', Courier, monospace;">fc_utilities</td>
      <td style="border: 1px solid #ccc; padding: 8px;">Spatial</td>
      <td style="border: 1px solid #ccc; padding: 8px;">NERI base data of utilities in the region, added to the database via ArcGIS Pro.</td>
    </tr>
    <!-- Row 13 -->
    <tr style="background-color: #ffffff;">
      <td style="border: 1px solid #ccc; padding: 8px; font-family: 'Courier New', Courier, monospace;">fc_water_bodies</td>
      <td style="border: 1px solid #ccc; padding: 8px;">Spatial</td>
      <td style="border: 1px solid #ccc; padding: 8px;">NERI base data of water bodies in the region, added to the database via ArcGIS Pro.</td>
    </tr>
    <!-- Row 14 -->
    <tr style="background-color: #fafafa;">
      <td style="border: 1px solid #ccc; padding: 8px; font-family: 'Courier New', Courier, monospace;">tbl_data</td>
      <td style="border: 1px solid #ccc; padding: 8px;">Non-Spatial</td>
      <td style="border: 1px solid #ccc; padding: 8px;">Contains raw observation data including species, locations, and sampling events; loaded from CSV.</td>
    </tr>
    <!-- Row 15 -->
    <tr style="background-color: #ffffff;">
      <td style="border: 1px solid #ccc; padding: 8px; font-family: 'Courier New', Courier, monospace;">tbl_events</td>
      <td style="border: 1px solid #ccc; padding: 8px;">Non-Spatial</td>
      <td style="border: 1px solid #ccc; padding: 8px;">Stores metadata about sampling events, such as the year of occurrence; loaded from CSV.</td>
    </tr>
    <!-- Row 16 -->
    <tr style="background-color: #fafafa;">
      <td style="border: 1px solid #ccc; padding: 8px; font-family: 'Courier New', Courier, monospace;">tbl_locations</td>
      <td style="border: 1px solid #ccc; padding: 8px;">Non-Spatial</td>
      <td style="border: 1px solid #ccc; padding: 8px;">Describes sampling locations with spatial coordinates and habitat details; loaded from CSV.</td>
    </tr>
    <!-- Row 17 -->
    <tr style="background-color: #ffffff;">
      <td style="border: 1px solid #ccc; padding: 8px; font-family: 'Courier New', Courier, monospace;">v_all_observations</td>
      <td style="border: 1px solid #ccc; padding: 8px;">Non-Spatial</td>
      <td style="border: 1px solid #ccc; padding: 8px;">Contains species metadata, including Latin and common names, with associated categories; loaded from CSV.</td>
    </tr>
    <!-- Row 18 -->
    <tr style="background-color: #fafafa;">
      <td style="border: 1px solid #ccc; padding: 8px; font-family: 'Courier New', Courier, monospace;">v_transmission_lines</td>
      <td style="border: 1px solid #ccc; padding: 8px;">Non-Spatial</td>
      <td style="border: 1px solid #ccc; padding: 8px;">A view filtering the utilities layer to include only power transmission lines; created via SQL.</td>
    </tr>
  </tbody>
</table>

*Table 1: Summary table of the enterprise geodatabase feature classes, tables, and views used in the NERI Transmission Corridor system. Items marked “Spatial” are feature classes managed in ArcGIS, while “Non-Spatial” items are standard relational tables or SQL views.*

## Results

By centralizing all data in a PostgreSQL enterprise geodatabase, the application users gained:

- **Data Consistency and Integrity:** Referenced foreign keys and consistent naming conventions (e.g., `SpeciesID`, `LocationID`) enforced strict data relationships. Views such as `v_all_observations` provided a unified lens for verifying species records without duplicating data.
- **Efficient Spatial-Non-Spatial Integration:** ArcGIS Pro users could seamlessly edit and visualize both tabular attributes (e.g., building contribution statuses) and geometry-based features (e.g., proposed line buffers) in the same environment.
- **Scalable Geoprocessing Workflows:** Automated models and scripts could rapidly ingest new CSV files, generate buffer polygons, or reset proposal statuses. The database design facilitated smooth concurrency—multiple teams could run analyses or add data without conflicts.

## Reflection

This data management solution highlights the power of an enterprise geodatabase in uniting diverse geospatial and tabular information. Key lessons include:

- **Schema-Driven Reliability:** Defining robust table structures, relationships, and constraints from the outset ensured stable, high-quality data.
- **SQL Views for Data Normalization:** Strategically placed logic in views (e.g., converting habitat strings to standardized labels) streamlined subsequent GIS processing, saving time in repeated tasks.
- **Interoperability:** The combination of ArcGIS Pro, PostgreSQL, and a custom web application exemplified a hybrid approach—capitalizing on proprietary GIS strengths while retaining the flexibility of open-source database technologies.

Overall, employing an enterprise geodatabase architecture provided a cohesive, scalable foundation for managing the dynamic and iterative data requirements of a complex environmental planning project. 
