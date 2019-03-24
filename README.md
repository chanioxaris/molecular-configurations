## Overview

## Metric Distance

### Discrete Frechet Distance (DFD)
The [Fréchet distance](https://en.wikipedia.org/wiki/Fr%C3%A9chet_distance) is a measure of similarity between curves that takes into account the location and ordering of the points along the curves. It is named after [Maurice Fréchet](https://en.wikipedia.org/wiki/Maurice_Fr%C3%A9chet).

Imagine a man traversing a finite curved path while walking his dog on a leash, with the dog traversing a separate one. Assume that the dog varies her speed to keep as much slack in her leash as possible: the Fréchet distance between the two curves is the length of the shortest leash sufficient for both to traverse their separate paths. Note that the definition is symmetric with respect to the two curves—the Frechet distance would be the same if the dog was walking her owner.

![Frechet Distance](https://github.com/chanioxaris/MolecularConfigurations-Clustering/blob/master/img/frechet_distance.jpg)

## Evaluation (Silhouette)
Silhouette refers to a method of interpretation and validation of consistency within clusters of data. The technique provides a succinct graphical representation of how well each object lies within its cluster.
The silhouette value is a measure of how similar an object is to its own cluster (cohesion) compared to other clusters (separation). The silhouette ranges from −1 to +1, where a high value indicates that the object is well matched to its own cluster and poorly matched to neighboring clusters. If most objects have a high value, then the clustering configuration is appropriate. If many points have a low or negative value, then the clustering configuration may have too many or too few clusters.


![Silhouette](https://github.com/chanioxaris/MolecularConfigurations-Clustering/blob/master/img/silhouette.jpg)

## Data files

### Input file 
The format of input text file, described by the following structure:
```
numConform: <integer>
N: <integer>
x11 y11 z11
x12 y12 z12
...
x1N y1N z1N
x21 y21 z21
...
x<numConform><N> y<numConform><N> z<numConform><N>
```
where ```numConform``` is the total molecular configurations, ```N``` the number of points in each molecular configuration

### Output file 
The format of output text file, described by the following structure:
```
k: <integer>
s: <real in [-1,1]>
ID11 ID12 ID13 ID14
...
IDk1 IDk2 IDk3
```
where ```k``` the number of clusters, ```s``` the silhouette value  

## Compile

`./makefile`

## Usage

`./proteins –i [input file] -frechet (optional)`
