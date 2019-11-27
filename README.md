# s3RNA: Subcellular swine small RNA data

This is the R pacakge for analysis of the granulosa small RNA sequencing dataset.

## Installation
s3RNA is a package for the R statisical computing language, and was built and tested on versions >= 3.5, and certain dependencies did not work during testing on 3.3. To install s3RNA, type either of these commands into the console:

``` r
devtools::install_github("derektoms/s3RNA")
```

...or...

``` r
source("https://install-github.me/derektoms/s3RNA")
```

No packages imported by s3RNA, although it suggests rnaseqGene and DESeq2 for analysis and visualization.

## Run the Shiny app

Two functions are provided with this package, each of which loads a specific dataset:
```r
miRNA()
```
and 
```r
snoRNA()
```
Each of these returns a list of two data frames containing unprocessed read counts (```raw.count```) and column data (```col.dat```). Analysis has been performed using the DESeq2 package, although additional approaches are possible. An example use would be to load the data, define an analysis model, compute significant differences and check the distribution of p-values. This is shown below:
```r
library(s3RNA)
library(ggplot)

sno_mat <- snoRNA()
sno <- DESeqDataSetFromMatrix(countData = sno_mat$raw.count, colData = sno_mat$col.dat, design=~size*subcell+batch)

sno.d <- DESeq(sno)

## Full linear model, two factors
size <- results(dds1, contrast=c("size","small","large"))
local <- results(dds1, contrast=c("subcell","cytosol","nucleus"))
sno.loc <- local[which(local$padj<0.1),]

## Check distribution of p-values
ggplot(as(local, "data.frame"), aes(x = pvalue)) 
+ geom_histogram(binwidth = 0.01, fill = "darkslategray", boundary = 0)
```

---

# License
This package is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License, version 3, as published by the Free Software Foundation.

This program is distributed in the hope that it will be useful, but without any warranty; without even the implied warranty of merchantability or fitness for a particular purpose. See the GNU General Public License for more details.

A copy of the GNU General Public License, version 3, is available at https://www.r-project.org/Licenses/GPL-3

