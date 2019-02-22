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
setwd('C:/Users/Clif/Dropbox/Eurobats')

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

polygon.sp1 = gUnaryUnion(mammterr[19587:19594,])# Eptesicus nilssonii
polygon.sp2 = gUnaryUnion(gSimplify(mammterr[35043:35093,], tol=0.0001))# Eptesicus serotinus
polygon.sp3 = gUnaryUnion(mammterr[25553:25568,])# Miniopterus schreibersii
polygon.sp4 = gUnaryUnion(mammterr[23718:23730,])# Myotis bechsteinii
polygon.sp5 = gUnaryUnion(mammterr[41313:41392,])# Myotis blythii
polygon.sp6 = gUnaryUnion(mammterr[18481:18518,])# Myotis capaccinii
polygon.sp7 = gUnaryUnion(mammterr[17746:17752,])# Myotis dasycneme
polygon.sp8 = gUnaryUnion(mammterr[13913:14114,])# Myotis daubentonii
polygon.sp9 = gUnaryUnion(mammterr[38593:38608,])# Myotis emarginatus
polygon.sp10 = gUnaryUnion(mammterr[28254:28258,])# Myotis myotis
polygon.sp11 = gUnaryUnion(mammterr[29195:29256,])# Myotis mystacinus
polygon.sp12 = gUnaryUnion(mammterr[33124:33148,])# Myotis nattereri
polygon.sp13 = gUnaryUnion(mammterr[24898:24918,])# Nyctalus noctula
polygon.sp14 = gUnaryUnion(mammterr[37135:37157,])# Pipistrellus nathusii
polygon.sp15 = gUnaryUnion(mammterr[30014:30050,])# Pipistrellus pipistrellus
polygon.sp16 = gUnaryUnion(mammterr[32948:32973,])# Pipistrellus pygmaeus
polygon.sp17 = gUnaryUnion(mammterr[2442:2449,])#Plecotus auritus
polygon.sp18 = gUnaryUnion(mammterr[989:1005,])#Rhinolophus blasii
polygon.sp19 = gUnaryUnion(mammterr[4613:4626,])#Rhinolophus euryale
polygon.sp20 = gUnaryUnion(mammterr[5629:5661,])#Rhinolophus ferrumequinum
polygon.sp21 = gUnaryUnion(mammterr[4327:4395,])#Rhinolophus mehelyi

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