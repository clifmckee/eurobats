##########################
# Code for determining which species distributions overlap
# R version 3.5.2 'Eggshell Igloo'
# Rtools version 3.5
# Works: 2019.02.22
##########################

###################
### Import data ###
###################

# Necessary packages
library(sp) # Version 1.3-1
library(rgeos) # Version 0.3-28
library(maps) # Version 3.3.0
library(mapdata) # Version 2.3.0
library(maptools) # Version 0.9-3
library(geosphere) # Version 1.5-7
library(ggplot2) # Version 3.0.0
library(ggthemes) # Version 4.0.1
library(reshape2) # Version 1.4.3
library(phangorn) # Version 2.4.0
library(viridisLite) # Version 0.3.0
library(viridis) # Version 0.5.1

# Set working directory
# setwd('C:/Users/Clif/Dropbox/Eurobats')

# Print out session info
sessionInfo()

# Read in the shape file
# Available from IUCN at http://www.iucnredlist.org/technical-documents/spatial-data/
# readShapeSpatial() is in the 'sp' package
# Shape files for the distributions of every terrestrial mammal in the IUCN database
# mammterr = readShapeSpatial('~/TERRESTRIAL_MAMMALS/TERRESTRIAL_MAMMALS.shp')
load(file='./Data/mammterr.RData')
summary(mammterr)
class(mammterr)

# Create a list of species names to filter mammterr
allnames = c('Eptesicus nilssonii',
            'Eptesicus serotinus',
            'Miniopterus schreibersii',
            'Myotis bechsteinii',
            'Myotis blythii',
            'Myotis capaccinii',
            'Myotis dasycneme',
            'Myotis emarginatus',
            'Myotis daubentonii',
            'Myotis myotis',
            'Myotis mystacinus',
            'Myotis nattereri',
            'Nyctalus noctula',
            'Pipistrellus nathusii',
            'Pipistrellus pipistrellus',
            'Pipistrellus pygmaeus',
            'Plecotus auritus',
            'Rhinolophus blasii',
            'Rhinolophus euryale',
            'Rhinolophus ferrumequinum',
            'Rhinolophus mehelyi') # The species I care about
# Filter mammterr down to just binomial names
all.binomial = mammterr$binomial

# Which rows of the data are the species I care about
keep = list()
for(i in 1:length(allnames)){
  keep[[i]] = which(all.binomial == allnames[i])
}

x = keep[[1]]

for(i in 2:length(allnames)){
  x = c(x, keep[[i]])
}

keep = x
myspecies.distr = mammterr[keep,]

# Helpful to output the 'keep' numbers to Excel so you can manually inspect polygons and
# see which species is assigned to each shapefile
write.csv(x=keep, file='./Results/keep_species.csv')

################################
### Checking data for errors ###
################################

# See Bat_polygons.R for code to resolve self-intersections and other errors
source(file='./Code/Bat_polygons.R')

# Set species order
sp.order = c('Ept.nil', 'Ept.ser', 'Min.sch', 'Myo.bec', 'Myo.bly', 'Myo.cap', 'Myo.das',
            'Myo.dau', 'Myo.ema', 'Myo.myo', 'Myo.mys', 'Myo.nat', 'Nyc.noc', 'Pip.nat',
            'Pip.pip', 'Pip.pyg', 'Ple.aur', 'Rhi.bla', 'Rhi.eur', 'Rhi.fer', 'Rhi.meh')

# Read in centroids
centroids = read.csv('./Results/Bat_species_centroids.csv', header=T)
centroids = centroids[,c(1, 2, 4, 3)]
centroids$X = sp.order; colnames(centroids) = c('species', 'polygon', 'lat', 'lon')

# Calculate line distance between centroids
index = expand.grid(1:length(sp.order), 1:length(sp.order))
linedistance = matrix(NA, length(sp.order), length(sp.order), dimnames=list(sp.order, sp.order))
for(i in 1:dim(index)[1]){
  linedistance[index[i,1], index[i,2]]=distm(c(centroids$lat[index[i,1]], centroids$lon[index[i,1]]),
                                             c(centroids$lat[index[i,2]], centroids$lon[index[i,2]]))
}
write.csv(linedistance,'./Results/linedistance.csv')

# Make UPGMA tree from linedistance matrix
upgma.tree = upgma(linedistance)
plot(upgma.tree)
# Save as Nexus format
write.nexus(upgma.tree, file='./Results/UPGMA_distance.nex')

##################################
### Species distribution plots ###
##################################

# Plot all species distributions as transparent layers
greentrans = rgb(0, 100, 0, 50, maxColorValue=255)
dev.off()
png('./Results/species_map.png', height=5, width=7, units='in', res=300)
map('world', xlim=c(-15,155), ylim=c(15,70), col='grey')
plot(myspecies.distr, add=T, col=greentrans)
dev.off()

# Plot species distributions separately
png('./Results/sep_species_map.png', height=11, width=8.5, units='in', res=300)
par(mfrow=c(7, 3), mar=c(0, 0, 0, 0))
for(i in 1:21){
  map('world', xlim=c(-15,155), ylim=c(15,70), col='grey')
  plot(get(paste('polygon.sp', i, sep='')), col=greentrans, add=T)
  text(x=105, y=65, labels=allnames[i])
}
dev.off()

######################
### Sampling sites ###
######################

# Read in sampling site locations
sites = read.csv('./Data/Eurobats_map_points_filtered.csv')
sites$Prevalence = sites$Positives/sites$Samples

# Make world map
world = fortify(map_data('world'))

png('./Results/sample_map.png', height=10, width=13.333, units='in', res=300)
ggplot() +
  geom_map(data=world, map=world, aes(map_id=region), fill=NA, colour='grey') +
  xlim(0, 30) +
  ylim(40, 55) +
  theme_map(base_size=20) +
  geom_count(data=sites, aes(x=long, y=lat, size=Samples), alpha=0.6, colour='black') +
  scale_size(range=c(1, 12), breaks=c(1,10,50,100,300), name='Sample sizes') +
  annotate('text', x=9, y=53, label = 'Netherlands', size=6) +
  annotate('segment', x=7.3, xend=5, y=53, yend=52) +
  annotate('text', x=8, y=51, label = 'Belgium', size=6) +
  annotate('segment', x=6.8, xend=4.8, y=51, yend=50.7) +
  annotate('text', x=25, y=46, label = 'Romania', size=6) +
  annotate('text', x=19, y=47, label = 'Hungary', size=6) +
  theme(legend.title=element_text(colour='black', size=24, face='bold'))
dev.off()

png('./Results/sample_map_prevalence.png', height=10, width=13.333, units='in', res=300)
ggplot() +
  geom_map(data=world, map=world, aes(map_id=region), fill=NA, colour='grey') +
  xlim(0, 30) +
  ylim(40, 55) +
  theme_map(base_size=20) +
  geom_count(data=sites, aes(x=long, y=lat, size=Samples, colour=Prevalence)) +
  scale_colour_viridis(alpha=0.6, direction=-1, option='D') +
  scale_size(range=c(1, 12), breaks=c(1,10,50,100,300), name='Sample sizes') +
  annotate('text', x=9, y=53, label = 'Netherlands', size=6) +
  annotate('segment', x=7.3, xend=5, y=53, yend=52) +
  annotate('text', x=8, y=51, label = 'Belgium', size=6) +
  annotate('segment', x=6.8, xend=4.8, y=51, yend=50.7) +
  annotate('text', x=25, y=46, label = 'Romania', size=6) +
  annotate('text', x=19, y=47, label = 'Hungary', size=6) +
  theme(legend.title=element_text(colour='black', size=24, face='bold'))
dev.off()

#####################################
### Overlap and area calculations ###
#####################################

# gArea(), gLength(), and gIntersection() use the 'rgeos' package
# over() uses the 'sp' package

# Calculate area of polygon (test)
plot(polygon.sp1, col='red')
plot(polygon.sp2, col='blue',add=T)
plot(polygon.sp8, col='green',add=T)
gArea(polygon.sp21)
gLength(polygon.sp21)

# Calculate area of union between two polygons (test)
gArea(gUnion(polygon.sp1, polygon.sp2)) # Eptesicus nilssonii and Eptesicus serotinus

# Gives line segments that overlap (test)
# All NAs means no overlaps
over(polygon.sp1, polygon.sp2) # Eptesicus nilssonii and Eptesicus serotinus

# If they overlap, give a 1, if not 0 (test)
ifelse(sum(over(polygon.sp1, polygon.sp2), na.rm=TRUE)>0, 1, 0) # Eptesicus nilssonii and Eptesicus serotinus

##############################
### Overlap and area loops ###
##############################

### Range loop
rangematrix = matrix(NA, length(allnames), 2); rownames(rangematrix) = allnames
for(i in 1:dim(rangematrix)[1]){
  rangematrix[i] = gArea(get(paste('polygon.sp', i, sep='')))
  rangematrix[i, 2] = log(rangematrix[i, 1], 10)
}
colnames(rangematrix) = c('range.size', 'log10.range.size')
write.csv(rangematrix, './Results/range_size.csv')

### Fragmentation loop
fragmatrix = matrix(NA, length(allnames), 2); rownames(fragmatrix) = allnames
for(i in 1:dim(fragmatrix)[1]){
  fragmatrix[i] = gArea(get(paste('polygon.sp', i, sep='')))/gLength(get(paste('polygon.sp', i, sep='')))
  fragmatrix[i, 2] = log(fragmatrix[i, 1], 10)
}
colnames(fragmatrix) = c('range.frag', 'log10.range.frag')
write.csv(fragmatrix, './Results/range_fragmentation.csv')

### Overlap loop
# Set up an empty matrix
overlapmatrix = matrix(NA, length(allnames), length(allnames), dimnames=list(allnames, allnames))

# Create an index for the loop
index = expand.grid(1:length(allnames), 1:length(allnames))

# Loops through index and calculates whether two shapefiles have lines that overlap
for(i in 1:dim(index)[1]){
if(exists(paste('polygon.sp', index[i, 1], sep='')) & exists(paste('polygon.sp', index[i, 2], sep=''))){
  overlapmatrix[index[i, 1], index[i, 2]] = ifelse(sum(over(get(paste('polygon.sp', index[i, 1], sep='')),
                                                            get(paste('polygon.sp', index[i,2], sep=''))),
                                                       na.rm=TRUE)>0, 1, 0)
  print(overlapmatrix[index[i, 1], index[i, 2]]) # Print as loop goes just to make sure it is working
    }
}

# Output matrix is sent to Excel sheet
write.csv(overlapmatrix, file='./Results/Bat_species_overlap.csv')

### Area difference loop
# Start with matrix of 0s and 1s from overlap matrix
areamatrix = overlapmatrix

# Loops through overlap matrix and calculates absolute difference in total area between shapefiles
for(i in 1:dim(index)[1]){
  areamatrix[index[i, 1], index[i, 2]] =
    abs(
      gArea(get(paste('polygon.sp', index[i, 1], sep='')))-
      gArea(get(paste('polygon.sp', index[i, 2], sep='')))
    )
  print(areamatrix[index[i, 1], index[i, 2]]) # Print as loop goes just to make sure it is working
}

# Output matrix is sent to Excel sheet
write.csv(areamatrix, file='./Results/Bat_species_diff_area.csv')

### Intersection area loop
# Start with matrix of 0s and 1s from overlap matrix
sectmatrix = overlapmatrix

# Loops through overlap matrix and calculates area of intersection between overlapping shapefiles
for(i in 1:dim(index)[1]){
  if(exists(paste('polygon.sp', index[i, 1], sep='')) & exists(paste('polygon.sp', index[i, 2], sep=''))
     & overlapmatrix[index[i, 1], index[i, 2]]!=0){
    sectmatrix[index[i, 1], index[i, 2]] =
      gArea(
        gIntersection(get(paste('polygon.sp', index[i, 1], sep='')),
                      get(paste('polygon.sp', index[i, 2], sep=''))))
    print(sectmatrix[index[i, 1], index[i, 2]]) # Print as loop goes just to make sure it is working
  }
}

# Output matrix is sent to Excel sheet
write.csv(sectmatrix, file='./Results/Bat_species_overlap_area.csv')

### Percent overlap loop
# Start with matrix of 0s and 1s from overlap matrix
pctmatrix = overlapmatrix

# Loops through overlap matrix and calculates percent area of intersection between overlapping shapefiles
for(i in 1:dim(index)[1]){
  if(exists(paste('polygon.sp', index[i, 1], sep='')) & exists(paste('polygon.sp', index[i, 2], sep=''))
     & overlapmatrix[index[i, 1], index[i, 2]]!=0){
    pctmatrix[index[i, 1], index[i, 2]] =
      gArea(
        gIntersection(get(paste('polygon.sp', index[i, 1], sep='')),
                      get(paste('polygon.sp', index[i, 2], sep=''))))/
      gArea(
        gUnion(get(paste('polygon.sp', index[i, 1], sep='')),
               get(paste('polygon.sp', index[i, 2], sep=''))))
    print(pctmatrix[index[i, 1], index[i, 2]]) # Print as loop goes just to make sure it is working
  }
}

# Output matrix is sent to Excel sheet
write.csv(pctmatrix, file='./Results/Bat_species_overlap_pct.csv')

###################
### End of code ###
###################