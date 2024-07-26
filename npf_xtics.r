#load the ncdf4 package 
library(ncdf4)

#open a netCDF file 
mydata <- nc_open('/home/coliewo/Desktop/analysis/combined_test1.nc')

#extract date and diameter
date <- ncvar_get(mydata,"time") 
diameter <- ncvar_get(mydata,"diameter")

neg_particles <- ncvar_get(mydata,"neg_particles")

library(raster) # package for raster manipulation
library(rgdal) # package for geospatial analysis
library(ggplot2) # package for plotting

plot(neg_particles)
