---
layout: page
order: 1
title: Utilizing the Clay Foundation Model and Sentinel-2 Imagery for Urban Growth Monitoring in Johnston County, North Carolina
permalink: /portfolio/utilizing-the-clay-foundation-model
---

*December, 2024*

Github Page: <a href="https://github.com/benton-tripp/urban-clay" target="_blank">https://github.com/benton-tripp/urban-clay</a>

## Problem

Many regions worldwide were experiencing rapid urban expansion, prompting a need for timely tracking of impervious surfaces—critical indicators of urban growth. Traditional datasets such as the U.S. National Land Cover Database (NLCD) often updated only every five years, creating data gaps in fast-changing areas. To address this issue, I developed a near-real-time framework that combined Sentinel-2 satellite imagery with the Clay Foundation Model—a pretrained deep learning model for Earth observation.

I chose Johnston County, North Carolina, as a use case because of its substantial population growth and diverse land-use patterns. By training neural networks on Clay embeddings of Sentinel-2 spectral imagery, I aimed to generate frequent, high-resolution imperviousness estimates that would help fill critical data gaps for planning authorities and researchers.

## Analysis Procedures

I focused on Johnston County, which covered about 2,050 km² in eastern North Carolina. After reprojecting the county boundary into UTM Zone 17N (EPSG:32617), I created a 600×600 m tile grid over the entire region. Each tile was then paired with two core datasets:

1. Sentinel-2 Imagery (four spectral bands at 10 m resolution), where I selected images with minimal cloud cover for multiple dates.  
2. NLCD Imperviousness data, which I clipped and resampled to 200 m resolution, producing a 3×3 target matrix per tile.

Figure 1 shows an example of the Sentinel-2 tile grid for April 2016, with each tile representing a consistent spatial unit of analysis.

<img src="{{site.baseurl}}/assets/images/tripp_b_image3.png" style="width: 80%; height:auto; max-width: 500px; min-width: 300px;">

<i style="width: 80%; max-width: 500px;">
Figure 1: Sentinel-2 RGB imagery of Johnston County in April 2016, overlaid with a 600x600-meter tile grid used and the county boundary.
</i>

### Embeddings Using the Clay Foundation Model

I leveraged the Clay Foundation Model to transform each tile’s multi-temporal, multi-spectral Sentinel-2 data into a 768-dimensional embedding vector. First, I compiled the tile’s spectral-spatial data, temporal embeddings, spatial embeddings, wave embeddings, and a scalar for ground sample distance, as summarized below: 

<script src="https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-mml-chtml.js"></script>
<div style="margin:5px; padding: 5px 5px 5px 15px; border: 1px solid #ababab; max-width: 800px; min-width: 300px;">
    <div>
        <span>Spectral-Spatial Data (pixels):</span>
        <ul style="margin-left:30px;">
            <li>4 bands (Blue, Green, Red, NIR)</li>
            <li>Each band normalized by known mean and standard deviation</li>
            <li>Spatial resolution: 10 m</li>
            <li>Tile size: 600 m × 600 m</li>
            <li>Resulting pixel matrix: 60 × 60</li>
            <li>(Overall spectral-spatial input: 4 × 60 × 60)</li>
        </ul>
        <span>Temporal Embeddings:</span>
        <ul style="margin-left:30px;">
            <li>Hour of day (\( h \)), normalized as:</li>
            <li>\(\sin(h \cdot 2\pi/24), \cos(h \cdot 2\pi/24)\)</li>
            <li>Week of year (\( w \)), normalized as:</li>
            <li>\(\sin(w \cdot 2\pi/52), \cos(w \cdot 2\pi/52)\)</li>
        </ul>
        <span>Spatial Embeddings:</span>
        <ul style="margin-left:30px;">
            <li>Tile centroid latitude (\( \text{lat} \)) and longitude (\( \text{lon} \)) in
                radians</li>
            <li>Normalized as:</li>
            <li>\(\sin(\text{lat}), \cos(\text{lat})\)</li>
            <li>\(\sin(\text{lon}), \cos(\text{lon})\)</li>
        </ul>
        <span>Ground Sample Distance (GSD):</span>
        <ul style="margin-left:30px;">
            <li>10 meters</li>
        </ul>
        <span><i>* Sine/cosine transformations ensure cyclical patterns are properly
                represented</i></span>
    </div>
</div>

I then fed this combined input into the Clay Foundation Model to extract high-level features describing local spectral signatures. Figure 2 illustrates how each tile’s inputs were constructed and passed through the model, ultimately producing a 768-dimensional embedding vector for each tile and date. 

<img src="{{site.baseurl}}/assets/images/tripp_b_image2.png" style="width: 80%; height:auto; max-width: 900px; min-width: 300px;">

<i style="width: 80%; max-width: 900px;">
Figure 2: Workflow illustrating the creation of datacubes for each date and tile in the study. Each datacube includes spectral-spatial tensors, temporal embeddings, spatial embeddings, wave embeddings, and ground sample distance (GSD). These datacubes serve as input to the Clay Foundation Model, which outputs a 768-dimensional vector of embedded data for each date and tile.
</i>

### Model Training

Once I had the 768-dimensional Clay embeddings, I trained two types of fully connected neural networks (simple vs. deep) to predict imperviousness. The target variable was the 3×3 imperviousness matrix for each tile, flattened into a single nine-value output vector. I used mean squared error (MSE) as the main loss function and performed hyperparameter tuning (adjusting dropout, weight decay, and layer sizes) to manage overfitting. I also split the data by space (cluster-based partitioning) and time (pre-2017 vs. post-2017) to ensure that model evaluation was robust and geographically diverse.

## Results

1. Small Feature Challenges: I found it challenging to accurately estimate smaller impervious features—like roads—because their narrow footprints were diluted within each 600×600 m tile.  
2. Seasonal Influences: Winter imagery revealed higher imperviousness estimates (less vegetation cover to obscure surfaces), whereas in spring/fall, vegetation tended to lower the estimates. Figure 3 showcases this seasonal difference in predicted imperviousness for April 2023 (left) versus November 2023 (right).
3. Model Performance:
    - Deep Neural Networks (DNNs) captured broad spatial patterns effectively, showing lower root mean squared error (RMSE) overall.  
    - Simple Neural Networks (SNNs) occasionally excelled in mean absolute error (MAE), suggesting they were better at minimizing smaller, consistent errors.
    - On average, from 2016 to 2023, my models estimated about a 2.44% yearly increase in imperviousness, which aligned reasonably well with the documented rapid growth in Johnston County.

<img src="{{site.baseurl}}/assets/images/tripp_b_image1.png" style="width: 80%; height:auto; max-width: 900px; min-width: 300px;">

<i style="width: 80%; height:auto; max-width: 900px;">
Figure 3: Model-estimated urban imperviousness for Johnston County. Left - April 2023 estimates (spring). Right - November 2023 estimates (winter).
</i>

## Reflection

By combining the Clay Foundation Model with frequent, high-resolution Sentinel-2 imagery, I demonstrated that near-real-time imperviousness monitoring was possible in a fast-growing region. DNNs proved highly effective in capturing complex spatial variability, whereas certain SNN configurations handled smaller, uniform errors more efficiently. Seasonal variations highlighted the importance of carefully selecting imagery and potentially incorporating vegetation indices to improve my embeddings’ ability to account for leaf-on/leaf-off conditions. The difficulties in detecting smaller impervious features hinted at opportunities for alternative architectures—such as convolutional networks—or more detailed sampling strategies. Nevertheless, this proof-of-concept showed that foundation models could effectively supplement traditional land cover products, especially in areas where timely urban development data were critical for planning and management.
