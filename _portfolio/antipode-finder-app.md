---
layout: page
order: 8
title: Antipode Finder App
permalink: /portfolio/antipode-finder-app
---

*July, 2024*

<i>See <a href="https://antipode-finder.com" target="_blank">https://antipode-finder.com</a></i>

## Problem

The Antipode Finder project addresses the need for a fast, user-friendly web service that computes and visualizes the antipodal point for any given location. The application is built with Flask and leverages geospatial computations to determine the point exactly opposite a user’s chosen location on Earth. With a focus on responsive design and interactive mapping, the project serves as an accessible tool for users to explore global spatial relationships.

Developed as part of a Masters in Geospatial Information Science and Technology, the web service integrates multiple functionalities. It not only calculates antipodal coordinates using straightforward mathematical transformations but also offers interactive maps—featuring both Web Mercator and Azimuthal Equal Area projections—to present results visually.

## Analysis Procedure

The backend is developed in Python using Flask, where multiple route handlers manage distinct sections of the site (e.g., Home, About, Contact, Antipodal Cities, and Antipode of USA). The core logic reads a compressed CSV containing city antipode data and converts it into JSON-like structures for use within the templates. Additionally, a suite of RESTful endpoints (e.g., `/save_location`) allows the application to receive and process user location data via AJAX requests.

On the client side, extensive use of JavaScript (with jQuery, Leaflet, and D3.js) creates a dynamic user experience. The scripts enable real-time geolocation, interactive map drawing, and smooth transitions between different map projections. The implementation emphasizes performance through techniques such as debouncing and responsive design adjustments. The inclusion of meta tags, canonical URLs, sitemaps, and proper handling of static files further demonstrate the robust web services framework underlying the project.

## Results

The project successfully delivers a fully functional web service that enables users to quickly find and visualize the antipodal point for any location. The integration of Flask with client-side libraries results in a seamless experience where users can input an address or use geolocation to generate immediate visual feedback via interactive maps. Both the Web Mercator and Azimuthal Equal Area projections are provided, allowing users to compare spatial representations and better understand map distortions.

From a web services perspective, the project demonstrates reliable route handling, dynamic content generation, and effective performance optimizations. The service also produces a sitemap, manages static assets, and supports essential files (such as `robots.txt` and `ads.txt`) for smooth search engine indexing and compliance with modern web standards. This comprehensive approach ensures that the application is both user-centric and robust in its technical implementation.

<img src="{{ site.baseurl }}/assets/images/web-services-img-3.png" style="width: 80%; max-width: 1000px; min-width: 300px;">

*Figure 1: Web Mercator projection displaying a selected point near North Carolina, USA (blue marker), with its corresponding antipode indicated in red near Australia.*

<img src="{{ site.baseurl }}/assets/images/web-services-img-4.png" style="width: 80%; max-width: 1000px; min-width: 300px;">

*Figure 2: Azimuthal Equal Area projection illustrating the selected point in North Carolina (blue) and its antipodal location in the Indian Ocean (red). This projection accurately represents areas, allowing direct spatial comparison between antipodal points.*

<img src="{{ site.baseurl }}/assets/images/web-services-img-5.png" style="width: 80%; max-width: 1000px; min-width: 300px;">

*Figure 3: World map using the Mollweide projection, displaying global landmasses (gray) overlaid with their corresponding antipodal geography (purple). This visualization highlights the relationship and symmetry between geographic locations and their antipodes.*

## Reflection

Working on the Antipode Finder has been instrumental in honing my skills in building scalable web services within a geospatial context. The project provided hands-on experience with Flask, RESTful APIs, and dynamic client-side interactions using JavaScript and D3.js. Integrating multiple data sources and ensuring compatibility across different map projections challenged me to consider both user experience and performance optimization. Reflecting on the development process, I gained valuable insights into designing web applications that effectively merge backend logic with rich, interactive frontends. The iterative testing and debugging of route management, AJAX data handling, and map rendering deepened my understanding of web services best practices. This project not only expanded my technical repertoire but also reinforced the importance of accessibility in modern web design.

