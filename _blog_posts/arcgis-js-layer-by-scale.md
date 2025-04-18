---
layout: post
title: "Efficient Scale-Dependent Layer Switching with ArcGIS JavaScript"
permalink: /blog_posts/arcgis-js-layer-by-scale
---

Switching between ArcGIS layers based on the map scale is a common technique to ensure efficient performance in web mapping applications. This post illustrates how you can dynamically alternate between a feature layer and a map image layer as users zoom in and out, maintaining popup functionality through spatial queries.

## Problem

Rendering complex polygon layers at various scales can significantly impact application performance. Feature layers provide detailed interactivity at large scales (zoomed-in views) but can be resource-intensive for large geographic areas. Map image layers render quickly at small scales (zoomed-out views) but typically lack direct interactive popup capabilities.

The solution is to intelligently switch between these layer types based on the map scale and retain popup interactions by querying the underlying data.

## Implementation

We use the ArcGIS JavaScript SDK to achieve scale-dependent layer visibility:

### Define Scale Threshold

Set a specific scale (`breakScale`) at which layers switch visibility:

```javascript
const breakScale = 50000; // Layers switch visibility at this scale threshold
```

### Feature Layer Setup

Configure the feature layer for detailed zoomed-in interactions, popups, and rendering:

```javascript
const featureLayer = new FeatureLayer({
    url: layerPath,            // URL to specific feature layer within the map service
    minScale: breakScale,      // Minimum scale at which this layer becomes visible (zoomed-in views)
    maxScale: 500,             // Maximum zoom-in scale for visibility
    visible: true,             // Initially set visible (controlled by scale)
    outFields: ["*"],          // All fields retrieved; can be reduced for performance
    renderer: renderer,        // Defines polygon styling (optional)
    popupTemplate              // Defines popup content and format (optional)
});
```

### Map Image Layer Setup

The map image layer ensures quick rendering for zoomed-out views:

```javascript
const imageLayer = new MapImageLayer({
    url: serviceRoot,          // Root URL of the map service
    minScale: 0,               // Visible from all scales down to specified max scale
    maxScale: breakScale,      // Maximum scale at which this layer is visible (zoomed-out views)
    sublayers: [{
        id: 9,                 // Sublayer ID within the map service
        visible: true,
        renderer: renderer
    }]
});
```

Both layers are added simultaneously, but only one is visible at any given scale.

## Retaining Popup Functionality

Popups at smaller scales (map image layer visibility) are achieved by executing spatial queries against the feature layer:

```javascript
view.popup.autoOpenEnabled = false;
view.on("click", function (event) {
    const query = featureLayer.createQuery();
    query.geometry = event.mapPoint;
    query.spatialRelationship = "intersects";
    query.outFields = ["*"];
    query.returnGeometry = true;

    featureLayer.queryFeatures(query)
        .then(function (featureSet) {
            if (featureSet.features.length > 0) {
                const selectedFeature = featureSet.features[0];
                view.popup.open({
                    features: [selectedFeature],
                    location: event.mapPoint
                });
            }
        })
        .catch(err => console.error("Query failed: ", err));
});
```

## Full Example

Below is the full implementation of a web map with layers that dynamically render depending on the scale. *This example uses data provided by [Wake County, NC Flood Zones](https://maps.wake.gov/arcgis/rest/services/Boundaries/ZipCodes/MapServer/0).*

```
<!DOCTYPE html>
<html lang="en">
<!-- ArcGIS JavaScript SDK Scale Dependent Layers Example -->

<head>
    <meta charset="utf-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Wake County Flood Plains</title>
    <link rel="stylesheet" href="https://js.arcgis.com/4.31/esri/themes/light/main.css">
    <script src="https://js.arcgis.com/4.31/"></script>
    <style>
        html,
        body {
            background-color: #333;
            padding: 0;
            margin: 0;
            width: 100%;
            min-width: 320px;
        }

        #viewDiv {
            height: 650px;
            width: 90%;
            min-width: 310px;
            margin: 10px auto;
            display: flex;
            justify-content: center;
            border: 1px solid #000;
            border-radius: 5px;
        }
    </style>
</head>

<body>
    <div id="viewDiv">
        <!-- MAP HERE -->
    </div>
</body>

<script>
    const serviceRoot = "https://maps.wake.gov/arcgis/rest/services/EnerGov/portal/MapServer";
    const layerPath = serviceRoot + "/9";
    const breakScale = 50000;
    const outFields = ["*"];

    require([
        "esri/Map", 
        "esri/views/MapView",
        "esri/layers/MapImageLayer",
        "esri/layers/FeatureLayer",
        "esri/widgets/Legend",
        "esri/widgets/Home",
        "esri/widgets/ScaleBar",
        "esri/widgets/Expand",
        "esri/widgets/Compass"
    ], (Map, MapView, MapImageLayer, FeatureLayer, Legend, Home, ScaleBar, Expand, Compass) => {

        // Define the map and view
        const map = new Map({
            basemap: "gray-vector"
        });

        const view = new MapView({
            container: "viewDiv",
            map,
            center: [-78.640885,35.781924],
            scale: 270000
        });

        // Define popup template
        const popupTemplate = {
            title: "Floodplain",
            content: [{
                type: "fields",
                fieldInfos: [
                { fieldName: "FLOODARID",    label: "Flood Area ID" },
                { fieldName: "FLOODZONE",    label: "FEMA Flood Zone" },
                { fieldName: "ZONE_IMAPS",   label: "Map Flood Zone" },
                { fieldName: "VERDATE",      label: "Verified Date",    format: { dateFormat: "short-date" } }
                ]
            }]
        };

        // Define the styling renderer
        const renderer = {
            type: "simple",
            symbol: {
                type: "simple-fill",
                color: [0, 134, 242, 0.9],
                outline: { color: "#000", width: 0.1 }
            }
        };

        // Feature layer (large‑scale, or zoomed-in) 
        featureLayer = new FeatureLayer({
            title: "Wake County Flood Plains – Feature Layer",
            url: layerPath,
            minScale: breakScale,
            maxScale: 500,
            visible: true,
            outFields: outFields,
            renderer: renderer,
            popupTemplate
        });

        // Map‑image layer (small‑scale, or zoomed-out) 
        const imageLayer = new MapImageLayer({
        url: serviceRoot,
        title: "Wake County Flood Plains – Image Layer",
        minScale:0,
        maxScale: breakScale, 
        sublayers: [    
            {
                id: 9,      
                visible: true,
                renderer: renderer,
                title: null
            }
        ]
        });

        // Add the layers to the map
        map.addMany([imageLayer, featureLayer]);

        // Add widgets
        const legend = new Legend({
            view: view,
            style: "classic"
        });

        const legendExpand = new Expand({
            view: view,
            content: legend,
            expandTooltip: "Expand Legend",
            expanded: true
        });
        view.ui.add(legendExpand, "bottom-left");

        const scaleBar = new ScaleBar({
            view: view,
            unit: "dual"
        });
        view.ui.add(scaleBar, "bottom-right");

        const compass = new Compass({
            view: view,
            visible: true
        });
        view.ui.add(compass, "bottom-right");

        const homeWidget = new Home({
            view: view
        });
        view.ui.add(homeWidget, "top-left");

        // Show pop-up (queries the feature layer, but works for both the image and feature layer)
        view.popup.autoOpenEnabled = false; // Disable default
        view.on("click", function (event) {

            const query = featureLayer.createQuery();
            query.geometry = event.mapPoint;
            query.spatialRelationship = "intersects";
            query.outFields = outFields;
            query.returnGeometry = true;

            featureLayer.queryFeatures(query)
                .then(function (featureSet) {
                    if (featureSet.features.length > 0) {
                        selectedFeature = featureSet.features[0];
                        view.popup.open({
                            features: [selectedFeature],
                            location: event.mapPoint
                        });
                    }
                })
                .catch(err => {
                    console.error("Query failed: ", err);
                });
        });

    });
</script>

</html>
```

Interact with the embedded map below to see the implementation of the code above in action:


<div>
    <iframe style="width:100%; height:690px;"
        src="/assets/iframe_src/example-arcgis-js-layer-by-scale.html" 
        width="100%" 
        height="700px" 
        style="border:0.5px solid #bbb; border-radius:3px;">
    </iframe>
</div>

## Summary

Implementing scale-dependent layers significantly enhances application responsiveness, combining fast rendering at smaller scales with detailed interactivity at larger scales. Using spatial queries ensures that popups can still be used within the `MapImageLayer` as well as the `FeatureLayer`.
