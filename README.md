## Overview

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
where ```k``` the number of clusters created, ```s``` the silhouette value  

## Compile

`./makefile`

## Usage

`./proteins –i [input file] -frechet (optional)`
