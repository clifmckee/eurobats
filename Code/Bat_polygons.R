##########################
# Code for combining shapefiles for each species
# R version 3.5.2 'Eggshell Igloo'
# Rtools version 3.5
# Works: 2019.02.22
##########################

# Necessary packages
# gUnaryUnion() and gSimplify() are in the 'rgeos' package
library(rgeos) # Version 0.3-28

# Set working directory
# setwd('C:/Users/Clif/Dropbox/Eurobats')

# Print out session info
sessionInfo()

##########################
# Instructions:
# 1. Start out by using the gUnaryUnion() function on all species polygons.
#     This will collapse all polygons attributed to a species to just a single polygon.
# 2. Run loop A below to check all shapefiles for self-intersections.
#     These create holes in the polygons that interfere with calculating intersections.
# 3. Examine which shapefiles produce self-intersections.
# 4. Manually change the shapefiles with self-intersections with the gSimplify() function
#     Use tol=0.0001 which is the setting for the rebuffering
##########################

polygon.sp1 = gUnaryUnion(mammterr[mammterr$binomial == allnames[1],]) # Eptesicus nilssonii
polygon.sp2 = gUnaryUnion(gSimplify(mammterr[mammterr$binomial == allnames[2],], tol=0.0001)) # Eptesicus serotinus
polygon.sp3 = gUnaryUnion(mammterr[mammterr$binomial == allnames[3],]) # Miniopterus schreibersii
polygon.sp4 = gUnaryUnion(mammterr[mammterr$binomial == allnames[4],]) # Myotis bechsteinii
polygon.sp5 = gUnaryUnion(mammterr[mammterr$binomial == allnames[5],]) # Myotis blythii
polygon.sp6 = gUnaryUnion(mammterr[mammterr$binomial == allnames[6],]) # Myotis capaccinii
polygon.sp7 = gUnaryUnion(mammterr[mammterr$binomial == allnames[7],]) # Myotis dasycneme
polygon.sp8 = gUnaryUnion(mammterr[mammterr$binomial == allnames[8],]) # Myotis daubentonii
polygon.sp9 = gUnaryUnion(mammterr[mammterr$binomial == allnames[9],]) # Myotis emarginatus
polygon.sp10 = gUnaryUnion(mammterr[mammterr$binomial == allnames[10],]) # Myotis myotis
polygon.sp11 = gUnaryUnion(mammterr[mammterr$binomial == allnames[11],]) # Myotis mystacinus
polygon.sp12 = gUnaryUnion(mammterr[mammterr$binomial == allnames[12],]) # Myotis nattereri
polygon.sp13 = gUnaryUnion(mammterr[mammterr$binomial == allnames[13],]) # Nyctalus noctula
polygon.sp14 = gUnaryUnion(mammterr[mammterr$binomial == allnames[14],]) # Pipistrellus nathusii
polygon.sp15 = gUnaryUnion(mammterr[mammterr$binomial == allnames[15],]) # Pipistrellus pipistrellus
polygon.sp16 = gUnaryUnion(mammterr[mammterr$binomial == allnames[16],]) # Pipistrellus pygmaeus
polygon.sp17 = gUnaryUnion(mammterr[mammterr$binomial == allnames[17],]) # Plecotus auritus
polygon.sp18 = gUnaryUnion(mammterr[mammterr$binomial == allnames[18],]) # Rhinolophus blasii
polygon.sp19 = gUnaryUnion(mammterr[mammterr$binomial == allnames[19],]) # Rhinolophus euryale
polygon.sp20 = gUnaryUnion(mammterr[mammterr$binomial == allnames[20],]) # Rhinolophus ferrumequinum
polygon.sp21 = gUnaryUnion(mammterr[mammterr$binomial == allnames[21],]) # Rhinolophus mehelyi

# Example code if tolerance needs to be reset
# gSimplify(mammterr[:], tol=0.0001)

#########################################
### Error checking and centroid loops ###
#########################################

# A. Loop for checking all species shapefiles for self-intersections
for (i in 1:length(allnames)){
  out = gIsValid(get(paste('polygon.sp', i, sep='')), reason=T, byid=T)
  print(i)
  print(out)
}

centroids = array(0, dim=c(21,3))

# B. Loop for calculating the centroid of all species shapefiles
for(i in 1:length(allnames)){
  out = coordinates(get(paste('polygon.sp', i, sep='')))
  print(i)
  print(out)
  centroids[i,1:3] = c(paste('polygon.sp', i, sep=''), out[1], out[2])
}
write.csv(x=centroids, file='./Results/Bat_species_centroids.csv')

###################
### End of code ###
###################
