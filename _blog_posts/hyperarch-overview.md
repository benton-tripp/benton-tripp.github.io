---
layout: post
title: Forecast Reconciliation in Python
permalink: /blog_posts/hyperarch-overview
---

Reconcile hierarchal forecasts into coherent forecasts using Python.

## Introduction
In my experience, it is somewhat unusual to work with real time-series data that does not have underlying "levels". For example, consider the data that might be produced via transactions at a grocery store. Transactions can be described at an individual product level, a shopper level, or at a store level. The products purchased might also be categorized into a specific type of product, and those categories might in turn fall into broader categories. For a business owner, this complex hierarchal structure can complicate making accurate and unbiased forecasts for their business. A data-driven person will most likely produce forecasts for each of the different levels, but they don't always add up. This is where reconciliation becomes necessary.

The majority of my explanations throughout this post are based on the book [*Forecasting: Principles and Practice*](https://otexts.com/fpp3/) by Rob J. Hyndman and George Athanasopoulos - an excellent resource for forecasting in general and completely free online. I highly recommend that you spend some time reading Hyndman's [more detailed explanation](https://otexts.com/fpp3/hierarchical.html) of reconciling hierarchal and grouped forecasts at some point. It is not my intention here to replace this book, although I still intend to give a comprehensive explanation of forecasting reconciliation. My intention in writing this blog-post is to provide whomever is reading this with the background needed to develop their own forecast reconciliation code. I have found that the [existing tools]() I have used are not (in my opinion) at a stage in their development where they can be considered reliable in a production environement, or cannot be generalized to every scenario (I know this from personal experience). I have also found that staring at a mathematical formula or applying some pre-defined function that somebody else created are not the best ways to gain a firm understanding of what is actually happening. Hopefully by sharing some short explanations of how one might go about writing their own Python code, I can also provide some additional insight that might not be attained otherwise.

## Hierarchal Time Series
To begin, I will give a brief introduction to the differences between hierarchal and grouped time series (no, they are not quite the same thing). A time series is considered *hierarchal*
when lower levels withing the hierarchy only fall under one domain. For example consider the following hierarchy tree:

```
from treelib import Tree

tree = Tree()
tree.create_node("Total", "total")  # root node
tree.create_node("Cat 1", "cat1", parent="total")
tree.create_node("Cat 2", "cat2", parent="total")
tree.create_node("SubCategory 1", "sub1", parent="cat1")
tree.create_node("SubCategory 2", "sub2", parent="cat1")
tree.create_node("SubCategory 3", "sub3", parent="cat2")
tree.create_node("SubCategory 4", "sub4", parent="cat2")
tree.create_node("SubCategory 5", "sub5", parent="cat2")
tree.show()
```
<div class="language-plaintext highlighter-rouge"><pre class="highlight code-print"><code class="highlight">
Total
├── Cat 1
│   ├── SubCategory 1
│   └── SubCategory 2
└── Cat 2
    ├── SubCategory 3
    ├── SubCategory 4
    └── SubCategory 5
</code></pre></div><br>

`Sub-Category 1` falls exclusively under `Category 1`, `Sub-Category 4` falls exlusively under `Sub-Category 2`, etc. This *exclusivity* is what defines a hierarchal time series. It is also important to note that mathematically, each of the sub-categories should add up to the category above them, and `Total` is the total sum of each of the bottom-level categories. The mathematical expression looks like this:

*T* = *C1* + *C2* = *SC1* + *SC2* + *SC3* + *SC4* + *SC5* <br>
*C1* = *SC1* + *SC2* <br>
*C2* = *SC3* + *SC4* + *SC5* <br>

On the other hand, a *grouped* time series is when that exclusivity between sub-domains does not exist. For example, consider the same hierarchy tree but with only three sub-categories spread across each of the two primary categories:

```
tree = Tree()
tree.create_node("Total", "total")  # root node
tree.create_node("Category 1", "cat1", parent="total")
tree.create_node("Category 2", "cat2", parent="total")
tree.create_node("SubCategory 1", "sub1", parent="cat1")
tree.create_node("SubCategory 2", "sub2", parent="cat1")
tree.create_node("SubCategory 3", "sub3", parent="cat1")
tree.create_node("SubCategory 1", "sub4", parent="cat2")
tree.create_node("SubCategory 2", "sub5", parent="cat2")
tree.create_node("SubCategory 3", "sub6", parent="cat2")
tree.show()
```
<div class="language-plaintext highlighter-rouge"><pre class="highlight code-print"><code class="highlight">
Total
├── Category 1
│   ├── SubCategory 1
│   ├── SubCategory 2
│   └── SubCategory 3
└── Category 2
    ├── SubCategory 1
    ├── SubCategory 2
    └── SubCategory 3
</code></pre></div><br>

The same hierarchy could also be described by the following tree:

```
tree = Tree()
tree.create_node("Total", "total")  # root node
tree.create_node("SubCategory 1", "sub1", parent="total")
tree.create_node("SubCategory 2", "sub2", parent="total")
tree.create_node("SubCategory 3", "sub3", parent="total")
tree.create_node("Category 1", "cat1", parent="sub1")
tree.create_node("Category 2", "cat2", parent="sub1")
tree.create_node("Category 1", "cat3", parent="sub2")
tree.create_node("Category 2", "cat4", parent="sub2")
tree.create_node("Category 1", "cat5", parent="sub3")
tree.create_node("Category 2", "cat6", parent="sub3")
tree.show()
```
<div class="language-plaintext highlighter-rouge"><pre class="highlight code-print"><code class="highlight">
Total
├── SubCategory 1
│   ├── Category 1
│   └── Category 2
├── SubCategory 2
│   ├── Category 1
│   └── Category 2
└── SubCategory 3
    ├── Category 1
    └── Category 2
</code></pre></div><br>

This complexity of multiple arangements of the hierarchy groups means that the original formulas used are no longer valid. Hyndman [describes](https://otexts.com/fpp3/hts.html#grouped-time-series) this structural concept as *"not naturally disaggregat[ing] in a unique hierarchical manner."* 

The difference between hierarchal and grouped time series is important to understand, because the summing matrix (explained briefly in the next section) is dependent on the structure of the time series.

## Building-Blocks of Coherent Forecasts
Coherent - or *reconciled* - forecasts are constructed from a few key components:

**Base Forecasts:** <br>
Forecasts at each hierarchal level represented as an *m* x *n* matrix (*m* rows and *n* columns), where the columns of the matrix represent the hierarchal levels of the data and the rows represent each time period of the forecast horizon.

**Summing matrix:** <br>
The summing matrix describes the hierarchal/grouped structure of the data. For hierarchal summing matrices, the number of columns matches the number of unique categories on the bottom level. The number of rows is determined by the total number of unique categories across all levels. The values of the summing matrix are binary values that represent which bottom-level (column) categories map to each of the hierarchies across all levels. Consider the example shared previously of hierarchal data with two categories and five sub-categories. The summing matrix would look like the following:

```
import pandas as pd

pd.DataFrame(
    data = {
        'Sub-Category 1': [1,1,0,1,0,0,0,0],
        'Sub-Category 2': [1,1,0,0,1,0,0,0],
        'Sub-Category 3': [1,0,1,0,0,1,0,0],
        'Sub-Category 4': [1,0,1,0,0,0,1,0],
        'Sub-Category 5': [1,0,1,0,0,0,0,1]
        },
    index = [
        'Total', 'Category 1', 'Category 2', 
        'Cat 1 -> Sub-Cat 1', 
        'Cat 1 -> Sub-Cat 2',
        'Cat 2 -> Sub-Cat 3',
        'Cat 2 -> Sub-Cat 4',
        'Cat 2 -> Sub-Cat 5'
        ])
```

<div style="min-width: 320px; overflow-x: auto; padding-bottom:10px;">
<table class="hts-dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>Sub-Category 1</th>
      <th>Sub-Category 2</th>
      <th>Sub-Category 3</th>
      <th>Sub-Category 4</th>
      <th>Sub-Category 5</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>Total</th>
      <td>1</td>
      <td>1</td>
      <td>1</td>
      <td>1</td>
      <td>1</td>
    </tr>
    <tr>
      <th>Category 1</th>
      <td>1</td>
      <td>1</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>Category 2</th>
      <td>0</td>
      <td>0</td>
      <td>1</td>
      <td>1</td>
      <td>1</td>
    </tr>
    <tr>
      <th>Cat 1 -&gt; Sub-Cat 1</th>
      <td>1</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>Cat 1 -&gt; Sub-Cat 2</th>
      <td>0</td>
      <td>1</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>Cat 2 -&gt; Sub-Cat 3</th>
      <td>0</td>
      <td>0</td>
      <td>1</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>Cat 2 -&gt; Sub-Cat 4</th>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>1</td>
      <td>0</td>
    </tr>
    <tr>
      <th>Cat 2 -&gt; Sub-Cat 5</th>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>1</td>
    </tr>
  </tbody>
</table>
</div>

Because grouped time series data does not disaggregate in a "unique manner," there is not one unique bottom level. This means that the summing matrix needs additional columns as well as additional rows. The grouped data from earlier in the post would look like this:

```
pd.DataFrame(
    data = {
        'Cat 1 / Sub-Cat 1': [1,1,0,1,0,0,1,0,0,0,0,0], 
        'Cat 1 / Sub-Cat 2': [1,1,0,0,1,0,0,1,0,0,0,0],
        'Cat 1 / Sub-Cat 3': [1,1,0,0,0,1,0,0,1,0,0,0],
        'Cat 2 / Sub-Cat 1': [1,0,1,1,0,0,0,0,0,1,0,0],
        'Cat 2 / Sub-Cat 2': [1,0,1,0,1,0,0,0,0,0,1,0],
        'Cat 2 / Sub-Cat 3': [1,0,1,0,0,1,0,0,0,0,0,1]
        },
    index = [
        'Total', 
        'Category 1', 
        'Category 2', 
        'Sub-Category 1',
        'Sub-Category 2',
        'Sub-Category 3',
        'Cat 1 / Sub-Cat 1', 
        'Cat 1 / Sub-Cat 2',
        'Cat 1 / Sub-Cat 3',
        'Cat 2 / Sub-Cat 1',
        'Cat 2 / Sub-Cat 2',
        'Cat 2 / Sub-Cat 3'
        ])
```

<div style="min-width: 320px; overflow-x: auto; padding-bottom:10px;">
<table class="hts-dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>Cat 1 / Sub-Cat 1</th>
      <th>Cat 1 / Sub-Cat 2</th>
      <th>Cat 1 / Sub-Cat 3</th>
      <th>Cat 2 / Sub-Cat 1</th>
      <th>Cat 2 / Sub-Cat 2</th>
      <th>Cat 2 / Sub-Cat 3</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>Total</th>
      <td>1</td>
      <td>1</td>
      <td>1</td>
      <td>1</td>
      <td>1</td>
      <td>1</td>
    </tr>
    <tr>
      <th>Category 1</th>
      <td>1</td>
      <td>1</td>
      <td>1</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>Category 2</th>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>1</td>
      <td>1</td>
      <td>1</td>
    </tr>
    <tr>
      <th>Sub-Category 1</th>
      <td>1</td>
      <td>0</td>
      <td>0</td>
      <td>1</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>Sub-Category 2</th>
      <td>0</td>
      <td>1</td>
      <td>0</td>
      <td>0</td>
      <td>1</td>
      <td>0</td>
    </tr>
    <tr>
      <th>Sub-Category 3</th>
      <td>0</td>
      <td>0</td>
      <td>1</td>
      <td>0</td>
      <td>0</td>
      <td>1</td>
    </tr>
    <tr>
      <th>Cat 1 / Sub-Cat 1</th>
      <td>1</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>Cat 1 / Sub-Cat 2</th>
      <td>0</td>
      <td>1</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>Cat 1 / Sub-Cat 3</th>
      <td>0</td>
      <td>0</td>
      <td>1</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>Cat 2 / Sub-Cat 1</th>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>1</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>Cat 2 / Sub-Cat 2</th>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>1</td>
      <td>0</td>
    </tr>
    <tr>
      <th>Cat 2 / Sub-Cat 3</th>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>1</td>
    </tr>
  </tbody>
</table>
</div>

**Mapping Matrix:** <br>
 The key componenent of forecast reconciliation is the mapping matrix. This matrix varies depending on the reconciliation method used, but the principle remains the same. When you multiply the mapping matrix with the summing and base forecast matrices (for either hierarchal or grouped time series), the result is a set of coherent forecasts. The challenge then is finding an "optimal" mapping matrix that can be used to reconcile the forecasts with the least amount of variance.

## Reconciliation Methods
Many different reconciliation methods exist, and some are better than others. I will not go into detail on the majority of them, but I will still list a few. Reconciliation methods fall under two categories: [single-level](https://otexts.com/fpp3/single-level.html) and [minimum trace](https://otexts.com/fpp3/reconciliation.html#the-mint-optimal-reconciliation-approach) methods. 

Some single-level methods include:
- Bottom-Up
- Top-Down
- Proportions of historical averages (method exists for proportions of historical actuals as well as historical forecasts)
- Middle-Out 

Some minimum trace methods include:
- Ordinary least squares (OLS)
- Weighted least squares (WLS) with variance scaling
- WLS with structrual scaling

For a more rigorous explanation of the different minimum trace methods - along with a few additional approaches - feel free to check out [this](https://robjhyndman.com/papers/mint.pdf) paper.

## Available Reconciliation Tools
As mentioned previously, there are a couple of existing tools/libraries that can apply the majority of reconciliation methods. The aforementioned Rob Hyndman helped to develop an R package called [`fabletools`](https://fabletools.tidyverts.org/reference/reconcile.html) that in my opinion is one of the best forecasting packages, and the best tool out there for forecast reconciliation. An effort is also underway to develop a Python package called [`scikit-hts`](https://scikit-hts.readthedocs.io/en/latest/readme.html). However, it is still pretty buggy and doesn't have all of the capabilities that exist in `fabletools`.

## Forecast Reconciliation in Python
I will not be sharing an explanation for coding all of the different methods available, nor will I be strictly showing all of the typical parts of machine learning (i.e. I won't be splitting data into training/test sets, generating forecasts, looking at accuracy metrics, etc.). My intention is to demonstrate how one might go about programming their own reconciliation algorithm, so I am going to be working with a sample dataset with pre-existing forecasts. <br>

Below are plots of each of the hierarchal levels in the sample dataset:

```
import pandas as pd
import numpy as np
import matplotlib.pylab as plt
plt.style.use('bmh')

df = pd.read_csv('../data/sample_hierarchal_data.csv').set_index('date')
df.index = pd.to_datetime(df.index).date 

def plot_hier(dataset, level='top'):
    data = dataset.copy()
    fig, ax = plt.subplots(figsize=(16, 5))
    data['pred'] = ~data.actual
    data.loc[max(data.loc[(data.pred == False)].index), 'pred'] = True
    if level == 'top':
        l1, = ax.plot(data.loc[data.actual==True]['total'], lw=1.5, c='blue')
        l1.set_label(f'totals - actual')
        l2, = ax.plot(data.loc[data.pred==True]['total'], lw=1.5, c='darkred')
        l2.set_label('totals - forecast')
        ax.set_title('total')
        ax.legend(bbox_to_anchor=(1, 1), loc='upper left', shadow=True, prop={'size': 11})
    elif level == 'middle':
        l1, = ax.plot(data.loc[data.actual==True]['category_1'], lw=1.5, c='green')
        l1.set_label('cat 1 - actual')
        l2, = ax.plot(data.loc[data.actual==True]['category_2'], lw=1.5, c='purple')
        l2.set_label('cat 2 - actual')
        l3, = ax.plot(data.loc[data.pred==True]['category_1'], lw=1.5, c='orange')
        l3.set_label('cat 1 - forecast')
        l4, = ax.plot(data.loc[data.pred==True]['category_2'], lw=1.5, c='cyan')
        l4.set_label('cat 2 - forecast')
        ax.set_title('categories')
        ax.legend(bbox_to_anchor=(1, 1), loc='upper left', shadow=True, prop={'size': 11})
    elif level == 'bottom':
        l1, = ax.plot(data.loc[data.actual==True]['subcategory_1'], lw=1.5, c='lime')
        l1.set_label('subcat 1 - actual')
        l2, = ax.plot(data.loc[data.actual==True]['subcategory_2'], lw=1.5, c='navy')
        l2.set_label('subcat 2 - actual')
        l3, = ax.plot(data.loc[data.actual==True]['subcategory_3'], lw=1.5, c='rebeccapurple')
        l3.set_label('subcat 3 - actual')
        l4, = ax.plot(data.loc[data.actual==True]['subcategory_4'], lw=1.5, c='tomato')
        l4.set_label('subcat 4 - actual')
        l5, = ax.plot(data.loc[data.actual==True]['subcategory_5'], lw=1.5, c='darkgoldenrod')
        l5.set_label('subcat 5 - actual')

        l6, = ax.plot(data.loc[data.pred==True]['subcategory_1'], lw=1.5, c='darkorange')
        l6.set_label('subcat 1 - forecast')
        l7, = ax.plot(data.loc[data.pred==True]['subcategory_2'], lw=1.5, c='blueviolet')
        l7.set_label('subcat 2 - forecast')
        l8, = ax.plot(data.loc[data.pred==True]['subcategory_3'], lw=1.5, c='gold')
        l8.set_label('subcat 3 - forecast')
        l9, = ax.plot(data.loc[data.pred==True]['subcategory_4'], lw=1.5, c='darkcyan')
        l9.set_label('subcat 4 - forecast')
        l10, = ax.plot(data.loc[data.pred==True]['subcategory_5'], lw=1.5, c='yellow')
        l10.set_label('subcat 5 - forecast')
        ax.set_title('subcategories')
        ax.legend(bbox_to_anchor=(1, 1), loc='upper left', shadow=True, prop={'size': 11})
    else:
        raise ValueError('level must be "top", "middle", or "bottom"')

plot_hier(df, level='top')
plot_hier(df, level='middle')
plot_hier(df, level='bottom')
```

<img src="{{ site.baseurl }}/assets/plots/hts-plt-1.png" 
    style="background-color: #ddd; width:98%; margin:5px; min-width: 320px; max-width: 1000px; height: auto; border: 1px solid #aaa;">
<img src="{{ site.baseurl }}/assets/plots/hts-plt-2.png" 
    style="background-color: #ddd; width:98%; margin:5px; min-width: 320px; max-width: 1000px; height: auto; border: 1px solid #aaa;">
<img src="{{ site.baseurl }}/assets/plots/hts-plt-3.png" 
    style="background-color: #ddd; width:98%; margin:5px; min-width: 320px; max-width: 1000px; height: auto; border: 1px solid #aaa;">


The first steps in reconciling these forecasts are to define the Base Forecast and Summing Matrices. In order to do this, we need to define the bottom/middle/top levels of the hierarchical structure. I have found that having the data in the following format makes this process pretty simple to do:

```
import pandas as pd

# read data
df = pd.read_csv('../data/hierarchical_data.csv')
display(df)
```

<div style="max-height:600px; overflow-y: auto; min-width: 320px; overflow-x: auto; padding-bottom:10px; padding-right:5px;">
<table class="hts-dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>date</th>
      <th>parent</th>
      <th>child</th>
      <th>forecast</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>1/1/2017</td>
      <td>category_1</td>
      <td>subcategory_1</td>
      <td>54.309109</td>
    </tr>
    <tr>
      <th>1</th>
      <td>2/1/2017</td>
      <td>category_1</td>
      <td>subcategory_1</td>
      <td>57.023606</td>
    </tr>
    <tr>
      <th>2</th>
      <td>3/1/2017</td>
      <td>category_1</td>
      <td>subcategory_1</td>
      <td>57.296805</td>
    </tr>
    <tr>
      <th>3</th>
      <td>4/1/2017</td>
      <td>category_1</td>
      <td>subcategory_1</td>
      <td>59.326200</td>
    </tr>
    <tr>
      <th>4</th>
      <td>5/1/2017</td>
      <td>category_1</td>
      <td>subcategory_1</td>
      <td>61.666806</td>
    </tr>
    <tr>
      <th>5</th>
      <td>6/1/2017</td>
      <td>category_1</td>
      <td>subcategory_1</td>
      <td>63.596201</td>
    </tr>
    <tr>
      <th>6</th>
      <td>1/1/2017</td>
      <td>category_1</td>
      <td>subcategory_2</td>
      <td>51.471773</td>
    </tr>
    <tr>
      <th>7</th>
      <td>2/1/2017</td>
      <td>category_1</td>
      <td>subcategory_2</td>
      <td>52.272320</td>
    </tr>
    <tr>
      <th>8</th>
      <td>3/1/2017</td>
      <td>category_1</td>
      <td>subcategory_2</td>
      <td>53.810966</td>
    </tr>
    <tr>
      <th>9</th>
      <td>4/1/2017</td>
      <td>category_1</td>
      <td>subcategory_2</td>
      <td>55.726676</td>
    </tr>
    <tr>
      <th>10</th>
      <td>5/1/2017</td>
      <td>category_1</td>
      <td>subcategory_2</td>
      <td>58.749791</td>
    </tr>
    <tr>
      <th>11</th>
      <td>6/1/2017</td>
      <td>category_1</td>
      <td>subcategory_2</td>
      <td>61.995237</td>
    </tr>
    <tr>
      <th>12</th>
      <td>1/1/2017</td>
      <td>category_1</td>
      <td>subcategory_3</td>
      <td>49.628702</td>
    </tr>
    <tr>
      <th>13</th>
      <td>2/1/2017</td>
      <td>category_1</td>
      <td>subcategory_3</td>
      <td>49.317622</td>
    </tr>
    <tr>
      <th>14</th>
      <td>3/1/2017</td>
      <td>category_1</td>
      <td>subcategory_3</td>
      <td>50.058703</td>
    </tr>
    <tr>
      <th>15</th>
      <td>4/1/2017</td>
      <td>category_1</td>
      <td>subcategory_3</td>
      <td>51.134935</td>
    </tr>
    <tr>
      <th>16</th>
      <td>5/1/2017</td>
      <td>category_1</td>
      <td>subcategory_3</td>
      <td>55.290929</td>
    </tr>
    <tr>
      <th>17</th>
      <td>6/1/2017</td>
      <td>category_1</td>
      <td>subcategory_3</td>
      <td>62.709843</td>
    </tr>
    <tr>
      <th>18</th>
      <td>1/1/2017</td>
      <td>category_2</td>
      <td>subcategory_4</td>
      <td>55.511010</td>
    </tr>
    <tr>
      <th>19</th>
      <td>2/1/2017</td>
      <td>category_2</td>
      <td>subcategory_4</td>
      <td>56.083683</td>
    </tr>
    <tr>
      <th>20</th>
      <td>3/1/2017</td>
      <td>category_2</td>
      <td>subcategory_4</td>
      <td>56.218949</td>
    </tr>
    <tr>
      <th>21</th>
      <td>4/1/2017</td>
      <td>category_2</td>
      <td>subcategory_4</td>
      <td>56.143276</td>
    </tr>
    <tr>
      <th>22</th>
      <td>5/1/2017</td>
      <td>category_2</td>
      <td>subcategory_4</td>
      <td>56.608268</td>
    </tr>
    <tr>
      <th>23</th>
      <td>6/1/2017</td>
      <td>category_2</td>
      <td>subcategory_4</td>
      <td>57.260983</td>
    </tr>
    <tr>
      <th>24</th>
      <td>1/1/2017</td>
      <td>category_2</td>
      <td>subcategory_5</td>
      <td>53.306174</td>
    </tr>
    <tr>
      <th>25</th>
      <td>2/1/2017</td>
      <td>category_2</td>
      <td>subcategory_5</td>
      <td>53.071540</td>
    </tr>
    <tr>
      <th>26</th>
      <td>3/1/2017</td>
      <td>category_2</td>
      <td>subcategory_5</td>
      <td>55.659420</td>
    </tr>
    <tr>
      <th>27</th>
      <td>4/1/2017</td>
      <td>category_2</td>
      <td>subcategory_5</td>
      <td>56.842727</td>
    </tr>
    <tr>
      <th>28</th>
      <td>5/1/2017</td>
      <td>category_2</td>
      <td>subcategory_5</td>
      <td>56.126760</td>
    </tr>
    <tr>
      <th>29</th>
      <td>6/1/2017</td>
      <td>category_2</td>
      <td>subcategory_5</td>
      <td>53.772583</td>
    </tr>
    <tr>
      <th>30</th>
      <td>1/1/2017</td>
      <td>total</td>
      <td>category_1</td>
      <td>171.273741</td>
    </tr>
    <tr>
      <th>31</th>
      <td>2/1/2017</td>
      <td>total</td>
      <td>category_1</td>
      <td>172.425305</td>
    </tr>
    <tr>
      <th>32</th>
      <td>3/1/2017</td>
      <td>total</td>
      <td>category_1</td>
      <td>175.925120</td>
    </tr>
    <tr>
      <th>33</th>
      <td>4/1/2017</td>
      <td>total</td>
      <td>category_1</td>
      <td>184.916548</td>
    </tr>
    <tr>
      <th>34</th>
      <td>5/1/2017</td>
      <td>total</td>
      <td>category_1</td>
      <td>206.228804</td>
    </tr>
    <tr>
      <th>35</th>
      <td>6/1/2017</td>
      <td>total</td>
      <td>category_1</td>
      <td>234.987138</td>
    </tr>
    <tr>
      <th>36</th>
      <td>1/1/2017</td>
      <td>total</td>
      <td>category_2</td>
      <td>114.194868</td>
    </tr>
    <tr>
      <th>37</th>
      <td>2/1/2017</td>
      <td>total</td>
      <td>category_2</td>
      <td>114.802626</td>
    </tr>
    <tr>
      <th>38</th>
      <td>3/1/2017</td>
      <td>total</td>
      <td>category_2</td>
      <td>120.535837</td>
    </tr>
    <tr>
      <th>39</th>
      <td>4/1/2017</td>
      <td>total</td>
      <td>category_2</td>
      <td>122.949692</td>
    </tr>
    <tr>
      <th>40</th>
      <td>5/1/2017</td>
      <td>total</td>
      <td>category_2</td>
      <td>122.393983</td>
    </tr>
    <tr>
      <th>41</th>
      <td>6/1/2017</td>
      <td>total</td>
      <td>category_2</td>
      <td>118.695822</td>
    </tr>
    <tr>
      <th>42</th>
      <td>1/1/2017</td>
      <td>NaN</td>
      <td>total</td>
      <td>275.827205</td>
    </tr>
    <tr>
      <th>43</th>
      <td>2/1/2017</td>
      <td>NaN</td>
      <td>total</td>
      <td>278.333631</td>
    </tr>
    <tr>
      <th>44</th>
      <td>3/1/2017</td>
      <td>NaN</td>
      <td>total</td>
      <td>286.530407</td>
    </tr>
    <tr>
      <th>45</th>
      <td>4/1/2017</td>
      <td>NaN</td>
      <td>total</td>
      <td>294.169931</td>
    </tr>
    <tr>
      <th>46</th>
      <td>5/1/2017</td>
      <td>NaN</td>
      <td>total</td>
      <td>304.508837</td>
    </tr>
    <tr>
      <th>47</th>
      <td>6/1/2017</td>
      <td>NaN</td>
      <td>total</td>
      <td>315.620680</td>
    </tr>
  </tbody>
</table>
</div>

Notice I have included `date`, `parent`, `child`, and `forecast` as columns in the dataset. Of course, however you choose to structure your data is completely up to you. But the methods I use to gather all of the information needed to reconcile the forecasts require that the data is formatted this way. <br>

Next, we can go ahead and define the different hierarchies within this dataset, and restructure the data into the correct base forecast matrix structure:

```
###-----------------------------------------------------------------------------
### pre-process so that it is formatted correctly for the base forecast matrix
### define hierarchy structure
###-----------------------------------------------------------------------------

# create new column combining parent & child
df['parent -> child'] = df.apply(
    lambda x: x.child if (x.parent is np.nan) else f'{x.parent}{" -> "}{x.child}', axis=1
    )

# get lists of hierarchy values (`parent`, `child`, and `parent -> child`)
parents = df.loc[~df.parent.isna()].parent.unique().tolist()
children = df.child.unique().tolist()
hier = df['parent -> child'].unique().tolist()

# define bottom level (any child values that don't exist in the parent column)
btm = [v for v in hier if (len(v.split(' -> ')) > 1) and (v.split(' -> ')[1] not in parents)]

# pivot data, set column order to be the same order as `labs` (this is an important step)
df = df.pivot(index='date', columns='parent -> child', values='forecast')[hier]
df.columns.names = [None]

# define base forecast matrix
bf = df.values

print('Restructured data:')
display(df)
print(f'Bottom level: {btm}\n')
print(f'Labels: {hier}\n')
print(f'Base forecast matrix: \n{bf}')
```

<div style="min-width:320px; overflow-x:auto; padding-bottom: 10px;">
<table class="hts-dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>category_1 -&gt; subcategory_1</th>
      <th>category_1 -&gt; subcategory_2</th>
      <th>category_1 -&gt; subcategory_3</th>
      <th>category_2 -&gt; subcategory_4</th>
      <th>category_2 -&gt; subcategory_5</th>
      <th>total -&gt; category_1</th>
      <th>total -&gt; category_2</th>
      <th>total</th>
    </tr>
    <tr>
      <th>date</th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>1/1/2017</th>
      <td>54.309109</td>
      <td>51.471773</td>
      <td>49.628702</td>
      <td>55.511010</td>
      <td>53.306174</td>
      <td>171.273741</td>
      <td>114.194868</td>
      <td>275.827205</td>
    </tr>
    <tr>
      <th>2/1/2017</th>
      <td>57.023606</td>
      <td>52.272320</td>
      <td>49.317622</td>
      <td>56.083683</td>
      <td>53.071540</td>
      <td>172.425305</td>
      <td>114.802626</td>
      <td>278.333631</td>
    </tr>
    <tr>
      <th>3/1/2017</th>
      <td>57.296805</td>
      <td>53.810966</td>
      <td>50.058703</td>
      <td>56.218949</td>
      <td>55.659420</td>
      <td>175.925120</td>
      <td>120.535837</td>
      <td>286.530407</td>
    </tr>
    <tr>
      <th>4/1/2017</th>
      <td>59.326200</td>
      <td>55.726676</td>
      <td>51.134935</td>
      <td>56.143276</td>
      <td>56.842727</td>
      <td>184.916548</td>
      <td>122.949692</td>
      <td>294.169931</td>
    </tr>
    <tr>
      <th>5/1/2017</th>
      <td>61.666806</td>
      <td>58.749791</td>
      <td>55.290929</td>
      <td>56.608268</td>
      <td>56.126760</td>
      <td>206.228804</td>
      <td>122.393983</td>
      <td>304.508837</td>
    </tr>
    <tr>
      <th>6/1/2017</th>
      <td>63.596201</td>
      <td>61.995237</td>
      <td>62.709843</td>
      <td>57.260983</td>
      <td>53.772583</td>
      <td>234.987138</td>
      <td>118.695822</td>
      <td>315.620680</td>
    </tr>
  </tbody>
</table>
</div>
<div class="language-plaintext highlighter-rouge"><pre class="highlight code-print"><code class="highlight">
Bottom level: [
    'category_1 -> subcategory_1', 
    'category_1 -> subcategory_2', 
    'category_1 -> subcategory_3', 
    'category_2 -> subcategory_4', 
    'category_2 -> subcategory_5'
]
</code></pre></div>
<div class="language-plaintext highlighter-rouge"><pre class="highlight code-print"><code class="highlight">
Labels: [
    'category_1 -> subcategory_1', 
    'category_1 -> subcategory_2', 
    'category_1 -> subcategory_3', 
    'category_2 -> subcategory_4', 
    'category_2 -> subcategory_5', 
    'total -> category_1', 
    'total -> category_2', 'total'
]
</code></pre></div>
<div class="language-plaintext highlighter-rouge"><pre class="highlight code-print"><code class="highlight">
Base forecast matrix: 
[[ 54.30910919  51.47177337  49.62870215  55.51101017  53.30617447
  171.2737405  114.1948681  275.8272049 ]
 [ 57.02360626  52.27232002  49.3176221   56.083683    53.0715397
  172.4253051  114.8026261  278.3336314 ]
 [ 57.29680483  53.81096637  50.05870335  56.21894853  55.65941955
  175.92512    120.5358369  286.530407  ]
 [ 59.32619958  55.72667601  51.13493525  56.14327624  56.84272682
  184.9165477  122.9496924  294.169931  ]
 [ 61.66680586  58.74979071  55.29092923  56.60826839  56.12676037
  206.2288039  122.3939829  304.508837  ]
 [ 63.59620071  61.99523706  62.70984315  57.26098344  53.7725831
  234.9871381  118.695822   315.6206805 ]]
</code></pre></div><br>

From the restructured dataset, defining the summing matrix can easily be accomplished using the following code:

```
import numpy as np

#define Summing Matrix

def middle_sm_vals(col):
    """Set the 'middle' hierarchy values to 1 in summing matrix"""
    tmp = col.reset_index()
    tmp.loc[tmp['index'].str.contains(col.name.split(' -> ')[0]), col.name] = 1
    return tmp.set_index('index')[col.name]

# stack each different level matrices into one summing matrix according to label order (in this case, top is on the bottom)
sm = np.vstack((
            # bottom level
            np.identity(len(btm)),

            # middle hierarchies
            pd.DataFrame(
                index=[h for h in hier if h not in btm and ' -> ' in h], 
                columns=btm, data=0
                ).apply(lambda x: middle_sm_vals(x)).values, 
    
            # top level all 1's
            np.ones(len(btm))))

print(sm)
```
<div class="language-plaintext highlighter-rouge"><pre class="highlight code-print"><code class="highlight">
[[1. 0. 0. 0. 0.]
 [0. 1. 0. 0. 0.]
 [0. 0. 1. 0. 0.]
 [0. 0. 0. 1. 0.]
 [0. 0. 0. 0. 1.]
 [1. 1. 1. 0. 0.]
 [0. 0. 0. 1. 1.]
 [1. 1. 1. 1. 1.]]
</code></pre></div>

Finally, we can define our mapping matrix. For this post, I will be demonstrating how to reconcile the forecasts using the Ordinary Least Squares method. The other methods can easily be implemented without too many major adjustments to the code now that you have the base forecast and summing matrices defined. If you would like to use any of the other methods, you can refer to [this](https://otexts.com/fpp3/reconciliation.html#the-mint-optimal-reconciliation-approach) section of Hyndman's book for the different optimal reconciliation approaches.

Using the OLS reconciliation method, we will use the following formula:

Let the summing matrix `sm` = *S*, the base forecast matrix `bf` = *F*, the forecast horizon = *h*, and the reconciled forecast output = *R*.

The mapping matrix *M* = (*S<sup>T</sup>S*)*<sup>-1</sup>S<sup>T</sup>*

Thus, *R* = *SMF<sub>h</sub>*

```
# define M, multiply, multiply by S
m = np.dot(np.dot(sm, np.linalg.inv(np.dot(np.transpose(sm), sm))), np.transpose(sm))

# reconcile by taking inner product between M and F for every h
rec = np.array([np.dot(m, bf[x, :]) for x in range(bf.shape[0])])

# output to dataframe
df_rec = pd.DataFrame(rec, index=df.index, columns=df.columns)

display(df_rec)
```

<div style="min-width: 320px; overflow-x:auto; padding-bottom:10px;">
<table class="hts-dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>category_1 -&gt; subcategory_1</th>
      <th>category_1 -&gt; subcategory_2</th>
      <th>category_1 -&gt; subcategory_3</th>
      <th>category_2 -&gt; subcategory_4</th>
      <th>category_2 -&gt; subcategory_5</th>
      <th>total -&gt; category_1</th>
      <th>total -&gt; category_2</th>
      <th>total</th>
    </tr>
    <tr>
      <th>date</th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>1/1/2017</th>
      <td>57.873479</td>
      <td>55.036143</td>
      <td>53.193072</td>
      <td>56.768012</td>
      <td>54.563177</td>
      <td>166.102694</td>
      <td>111.331189</td>
      <td>277.433882</td>
    </tr>
    <tr>
      <th>2/1/2017</th>
      <td>60.108384</td>
      <td>55.357098</td>
      <td>52.402400</td>
      <td>57.475269</td>
      <td>54.463126</td>
      <td>167.867882</td>
      <td>111.938395</td>
      <td>279.806277</td>
    </tr>
    <tr>
      <th>3/1/2017</th>
      <td>60.639390</td>
      <td>57.153552</td>
      <td>53.401289</td>
      <td>58.642004</td>
      <td>58.082475</td>
      <td>171.194231</td>
      <td>116.724479</td>
      <td>287.918710</td>
    </tr>
    <tr>
      <th>4/1/2017</th>
      <td>63.419464</td>
      <td>59.819940</td>
      <td>55.228199</td>
      <td>58.679279</td>
      <td>59.378730</td>
      <td>178.467603</td>
      <td>118.058009</td>
      <td>296.525612</td>
    </tr>
    <tr>
      <th>5/1/2017</th>
      <td>67.924990</td>
      <td>65.007974</td>
      <td>61.549113</td>
      <td>57.998405</td>
      <td>57.516897</td>
      <td>194.482077</td>
      <td>115.515303</td>
      <td>309.997380</td>
    </tr>
    <tr>
      <th>6/1/2017</th>
      <td>72.801796</td>
      <td>71.200833</td>
      <td>71.915439</td>
      <td>56.527244</td>
      <td>53.038843</td>
      <td>215.918068</td>
      <td>109.566087</td>
      <td>325.484155</td>
    </tr>
  </tbody>
</table>
</div>

To double check, you can double check that everything adds up coherently:

*T* = *C1* + *C2* = *SC1* + *SC2* + *SC3* + *SC4* + *SC5* <br>
*C1* = *SC1* + *SC2* <br>
*C2* = *SC3* + *SC4* + *SC5* <br>

And that's it! You have successfully reconciled the forecasts.

## References
[1] Rob J. Hyndman and George Athanasopoulos, *Forecasting: Principles and Practice 3rd ed.* (2022), [https://otexts.com/fpp3/](https://otexts.com/fpp3/) <br>
[2] Shanika L. Wickramasuriya, George Athanasopoulos, and Rob Hyndman; *Optimal Forecast Reconciliation for Hierarchical and Grouped Time Series Through Trace Minimization* (2019); [https://doi.org/10.1080/01621459.2018.1448825](https://doi.org/10.1080/01621459.2018.1448825)<br>
[3] Mitchell O'Hara-Wild, Rob Hyndman, Earo Wang; *fabletools* (2022); [https://fabletools.tidyverts.org/reference/reconcile.html](https://fabletools.tidyverts.org/reference/reconcile.html)<br>
[4] Carlo Mazzaferro, *scikit-hts* (2019), [https://scikit-hts.readthedocs.io/en/latest/readme.html](https://scikit-hts.readthedocs.io/en/latest/readme.html)