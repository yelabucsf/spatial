## this function plots the results of an xyzeq experiment
## requires the imager library

## input
## should be 5 columns: gene, x, y, cell_id, count
library(imager);
library(data.table);
library(sp);
install.packages("rgeos")
library(raster)
install.packages("geos-config")
library(rgeos);
##library(reshape);
library(Matrix);
library(latticeExtra);

## output
make_grid <- function(x, cell_diameter) {
  ext <- as(extent(x), "SpatialPolygons")
  ext <- as(extent(x)+cell_diameter, "SpatialPolygons");
  projection(ext) <- projection(x)

  # generate array of hexagon centers
  g <- spsample(ext,type = "hexagonal", offset = c(0,-0.5), cellsize=cell_diameter);
##                offset = c(0.5, 0.5))
  # convert center points to hexagons
  g <- HexPoints2SpatialPolygons(g)

  # clean up feature IDs
  row.names(g) <- as.character(1:length(g))
  ##browser();
  return(g)
}

input <- "/Users/Christa.Caggiano/Documents/UCSF_year1/Ye-rotation2/XYZ_20180319_cell_2_report/read_human_mouse/spatial_output.csv";

mat <- fread(input);
mat_m = mat[, c("V1", "V2", "V4")]
mat_h = mat[, c("V1", "V2", "V3")]
colnames(mat_h) <- c("x", "y","value")
colnames(mat_m) <- c("x", "y","value")

# mat$value <- mat$mouse;

# @TODO remove log transform 
mat.agg <- aggregate(mat$value, by=list(mat$x,mat$y), FUN=function(x) {log(sum(x+1))});
colnames(mat.agg) <- c("y","x","value")

mat.h.agg <- aggregate(mat_h$value, by=list(mat_h$x,mat_h$y), FUN=function(x) {log(sum(x+1))});
mat.m.agg <- aggregate(mat_m$value, by=list(mat_m$x,mat_m$y), FUN=function(x) {log(sum(x+1))});
colnames(mat.h.agg) <- c("x","y","value")
colnames(mat.m.agg) <- c("x","y","value")

m = merge(mat.h.agg, mat.m.agg, by=c("y", "x"), all=T )
m[is.na(m)] = 0

m$value = m$value.x/(m$value.y+m$value.x)

m[which(!is.finite(m))] = 0
mat.agg = m[, c("y", "x", "value")]
## let's offset x by 0.5 for even numbered ys
##mat$x[which(mat$y %% 2 == 1)] <- mat$x[which(mat$y %% 2 == 1)] + 1;

img <- as.cimg(mat.agg);

## let's make the grid based on what we know about the dimensions
x.n <- 14;
y.n <- 28;
x.n <- 18;
y.n <- 44;

## if the diameter is 1 on the x axis
## the radius is 1/sqrt(3)
## the extra bit is r/2 = 1/(2*sqrt(3))
## we should then scale everything in the y direction by 3/(2*sqrt(3))

##mat.coords <- expand.grid(1:x.n, 0:(3*(y.n)/(2*sqrt(3))));
##mat.coords$Var1[which(mat.coords$Var2 %% 2 == 1)] <- mat.coords$Var1[which(mat.coords$Var2 %% 2 == 1)]-1

mat.coords <- rbind(c(1,1),c(x.n,1),c(x.n,3*(y.n)/(2*sqrt(3))),c(1,3*(y.n)/(2*sqrt(3))));
##mat.coords.offset <- rbind(c(0,0),c(x.n,0),c(x.n,y.n),c(0,y.n));

img_area <- SpatialPolygons(list(Polygons(list(Polygon(mat.coords)), "x")));
img_grid <- make_grid(img_area, cell_diameter=1);

mat.agg.fill <- expand.grid(1:x.n, 1:y.n);
colnames(mat.agg.fill) <- c("x","y");
mat.agg.fill$value <- 0;
matched <- match(paste(mat.agg$y,mat.agg$x,sep="."),paste(mat.agg.fill$y,mat.agg.fill$x,sep="."));
mat.agg.fill$value[matched] <- mat.agg$value;

img_poly <- SpatialPolygonsDataFrame(img_grid, data.frame(mat.agg.fill$value), match.ID=F);

##plot(xlim=c(0,14),ylim=c(0,28));
spplot(img_poly);
##plot(img_grid, border="orange", add=T);

plot(mat$V4, mat$V3, xlab="human", ylab="mouse")


##hex_points <- spsampspsample(img_area, type = "hexagonal", cellsize = 0.5)

##pdf("test.pdf"); plot(img); dev.off();
## this function plots the results of an xyzeq experiment
## requires the imager library

## input
## should be 5 columns: gene, x, y, cell_id, count
library(imager);
library(data.table);
library(sp);
library(raster)
library(rgeos);
##library(reshape);
library(Matrix);
library(latticeExtra);

## output
make_grid <- function(x, cell_diameter) {
  ext <- as(extent(x), "SpatialPolygons")
  ext <- as(extent(x)+cell_diameter, "SpatialPolygons");
  projection(ext) <- projection(x)
  
  # generate array of hexagon centers
  g <- spsample(ext,type = "hexagonal", offset = c(0,-0.5), cellsize=cell_diameter);
  ##                offset = c(0.5, 0.5))
  # convert center points to hexagons
  g <- HexPoints2SpatialPolygons(g)
  
  # clean up feature IDs
  row.names(g) <- as.character(1:length(g))
  ##browser();
  return(g)
}


# make_grid <- function(x, n, cell_diameter, cell_area, clip = FALSE) {
#   if (missing(cell_diameter)) {
#     if (missing(cell_area)) {
#       stop("Must provide cell_diameter or cell_area")
#     } else {
#       cell_diameter <- sqrt(2 * cell_area / sqrt(3))
#     }
#   }
#   ##ext <- as(extent(x) + cell_diameter, "SpatialPolygons")
#   ext <- as(extent(x), "SpatialPolygons")
#   projection(ext) <- projection(x)
#   # generate array of hexagon centers
#   g <- spsample(ext, n=392, type = "hexagonal", cellsize = cell_diameter,
#                   offset = c(0,0));
# ##                offset = c(0.5, 0.5))
#   # convert center points to hexagons
#   g <- HexPoints2SpatialPolygons(g, dx = cell_diameter)
#   # clip to boundary of study area
#   if (clip) {
#     g <- gIntersection(g, x, byid = TRUE)
#   } else {
#     g <- g[x, ]
#   }
#   # clean up feature IDs
#   row.names(g) <- as.character(1:length(g))
#   ##browser();
#   return(g)
# }

##input <- "~/Box\ Sync/XYZeq\ Project/XYZ_20180122_456_2/count.txt";
##input <- "~/Box\ Sync/XYZeq\ Project/XYZ_20180122_456_2/spatial_output.csv";
input <- "/Users/Christa.Caggiano/Documents/UCSF_year1/Ye-rotation2/XYZ_20180328_pan12_report/human_mouse_gene_count/spatial_output.csv";
require(data.table)
mat <- fread(input);

##mat$value <- mat$mouse;

mat.agg <- aggregate(mat$V1, by=list(mat$V2,mat$V3), FUN=function(x) {log(sum(x+1))});
colnames(mat.agg) <- c("y","x","value")


mat.mouse.agg <- aggregate(mat$mouse, by=list(mat$x,mat$y), FUN=function(x) {sum(x)});
mat.human.agg <- aggregate(mat$human, by=list(mat$x,mat$y), FUN=function(x) {sum(x)});
mat.agg <- as.data.frame(cbind(mat.mouse.agg$Group.1, mat.mouse.agg$Group.2, value=mat.mouse.agg$x/(mat.mouse.agg$x+mat.human.agg$x)))
colnames(mat.agg) <- c("y","x","value")

## let's offset x by 0.5 for even numbered ys
##mat$x[which(mat$y %% 2 == 1)] <- mat$x[which(mat$y %% 2 == 1)] + 1;

img <- as.cimg(mat.agg);

## let's make the grid based on what we know about the dimensions
x.n <- 14;
y.n <- 27;
# x.n <- 18;
# y.n <- 44;

## if the diameter is 1 on the x axis
## the radius is 1/sqrt(3)
## the extra bit is r/2 = 1/(2*sqrt(3))
## we should then scale everything in the y direction by 3/(2*sqrt(3))

##mat.coords <- expand.grid(1:x.n, 0:(3*(y.n)/(2*sqrt(3))));
##mat.coords$Var1[which(mat.coords$Var2 %% 2 == 1)] <- mat.coords$Var1[which(mat.coords$Var2 %% 2 == 1)]-1

mat.coords <- rbind(c(1,1),c(x.n,1),c(x.n,3*(y.n)/(2*sqrt(3))),c(1,3*(y.n)/(2*sqrt(3))));
##mat.coords.offset <- rbind(c(0,0),c(x.n,0),c(x.n,y.n),c(0,y.n));

img_area <- SpatialPolygons(list(Polygons(list(Polygon(mat.coords)), "x")));
img_grid <- make_grid(img_area, cell_diameter=1);

mat.agg.fill <- expand.grid(1:x.n, 1:y.n);
colnames(mat.agg.fill) <- c("x","y");
mat.agg.fill$value <- 0;
matched <- match(paste(mat.agg$y,mat.agg$x,sep="."),paste(mat.agg.fill$y,mat.agg.fill$x,sep="."));
mat.agg.fill$value[matched] <- mat.agg$value;

img_poly <- SpatialPolygonsDataFrame(img_grid, data.frame(mat.agg$value), match.ID=F);

##plot(xlim=c(0,14),ylim=c(0,28));
spplot(img_poly);
##plot(img_grid, border="orange", add=T);

##hex_points <- spsampspsample(img_area, type = "hexagonal", cellsize = 0.5)

##pdf("test.pdf"); plot(img); dev.off();
