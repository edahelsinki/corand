# CORAND webdemo

## Prerequisites 

First, you should have [R](https://www.r-project.org/) and the following R libraries and their prerequisites installed:

* [shiny](https://cran.r-project.org/package=shiny)
* [DT](https://cran.r-project.org/package=DT)

To run the web demo you should first clone this repository and then execute the following
command from this directory: `Rscript --vanilla run.R`
Then you can just point your browser to <http://127.0.0.1:7000> and everything should work smoothly.

Alternatively, you can install [RStudio](https://www.rstudio.com/), set the current working directory
([`setwd`](https://stat.ethz.ch/R-manual/R-devel/library/base/html/getwd.html)) to the directory where
run.R resides, and then issuing command `runApp("webdemo",port=7000)`. 


## Walk-through of the UI

The system is experimental and meant to demonstrate the computational concepts. It is not therefore
meant to be used as a "production data analysis system", even though it works for that purpose as well.
The user interface (UI) is still quite rough, as you will notice.

We will now walk through the UI.

### Left-hand-side items

#### select rows

A central idea in the UI is that at each moment you can select a subset of rows and columns, which form the `current` selection. Here you can select the rows by painting the data items in the right hand side top scatterplot. If the selection
mode is `add` painting of rows add rows to the selection, if the selection mode is `remove` the painting of rows 
removes rows from the current selection. `clear current selection` empties the current selection of rows
and `reverse current selection` selects all unselected rows and unselects the selected rows.

#### select columns

As `select rows`, but for columns. The only difference is that there is a `toggle factors` button that removes
the factors from the lower right hand side scatterplot, if the data contains factors.

#### recompute projection

This option recomputes the projection used to make the right hand side scatterplots. The projections 
can either be updated by making PCA on the current selection (`PCA project on selection`) or by using
the tiling projection as described in [Puolamäki et al. (2018)](https://arxiv.org/abs/1805.07725) 
(`tiling projection`).

#### parallel coordinates

This interaction controls the parallel coordinate plot on the right hand side. You can adjust the number of columns
to show as well as the ordering (either via the first projection coordinate or the difference of the current selection
with respect to all data). Here as well you can add or remove columns from the current selection.

#### tiling

This interaction shows the status of the current pair of hypothesis. See [Puolamäki et al. (2018)](https://arxiv.org/abs/1805.07725) for a discussion.

#### tiles

This interaction contains one button (`add current tile`) by which you can save the current tile, after which it will
appear to the table at the bottom of the right hand side.

#### datafile

Here you can choose the datafile to use. The demo comes with some pre-installed data files, but you can easily
add your own data files by adding data frames, saved as rds files, to the data directory.

Here you can also dump the current datafile as well as state of the system to a rds file.

### Right-hand-side items

#### rows

Here you see a scatterplot to the directions found either by PCA or the tiling projection showing the data points.
You can make selections of rows from this plot by painting the appropriate points.

#### columns

Here you can see a scatterplot to the directions found either by PCA or the tiling projection showing the attributes.
You can make selections of rows from this plot by painting the appropriate points.
If the data attributes happen to have zero mean and unit variance the xy-coordinates of an attribute are simply
the attributes weights in the first and second projection direction, respectively.

#### parallel coordinate plot

Here you can see selected attributes and how their value compares for the full data vs. selection. 
You can select columns here also by painting with a mouse.

#### tile table

Here you can see the table of tiles stored by the system. The `current` tile is the current selection and `ALL` 
tile is a special tile containing all rows and columns. The system automatically generates a tile that correspond
to each of the factors. From the table, you can see the name of the tile, and the number of rows and columns in the tile.
You can also see the similarity of the tile to the current selection (recall, precision, jaccard). 

Columns
H1 and H2 show if the tile is included in the corresponding hypothesis pair. The contents of H1 and H2 affect how
the tiling projection is computed.

By default, H1 contains only the special tile `ALL` and H2 is empty. At his limit the tiling projection corresponds
to the PCA of the correlation matrix of the data.

If both of the values of H1 and H2 for a tile are false, then the tile is ignored in computing the tiling projection.

## Walk-through of the German example

Here we provide a walk-through of the experiment described in Puolamäki et al. (2019) or Appendix B of 
[Puolamäki et al. (2018)](https://arxiv.org/abs/1805.07725). 

After starting the web interface, change the dataset to `socio_economic_germany_cut.rds` by
choosing it in the `datafile` drop down menu. Then hit the `tiling projection` button, after
which you should reproduce Figure 9.

You can select the yellow points (eastern rural districts) either by mouse or by following clever trick: 
at the tile table at the bottom right, click with your mouse the number of rows (`sizeR`) of the following
tiles: `TypeUrban`, `RegionNorth`, `RegionSouth`, and `RegionWest`. Then click `reverse current selection` 
under **select rows**. You should now have in your current selection the districts that are both from
`RegionEast` and of `TypeRural`. (If you make a mistace you can always clear your selection by hitting
`clear current selection` 
under **select rows**.)

We will then select the rows using the parallel coordinate plot. With the rural eastern districts selected,
push `order by X1` under **parallel coordinates**. Select, e.g., the top 8 columns from the parallel coordinate plot
(from Unemploy. to Income). You notice that the distinctive features of the rural eastern district include, e.g.,
high unemplyment, low income, and small number of children.

With the rows and columns selected, add the current tile to the memory by clicking `add current tile` under **tiles**.
The tile will appear to the tile table under name `user01`. Add the tile to the background information by clicking
both H1 and H2 for tile `user01` from false to true, and then hit `tiling projection` under **recompute projection**. 
You should have now produced something similar to Figure 8. If you click `clear current selection` under both
**select rows** and **select columns** to clear selections and then click the number of rows (`sizeR`) of `TypeUrban`
in tile table you should highlight all urban locations. (Notice that you may not reproduce Figure 8 exactly, because 
the detailed projection depends on the columns chosen!)




