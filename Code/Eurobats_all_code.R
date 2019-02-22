##########################
# Required code for Eurobats analyses
# R version 3.5.2 'Eggshell Igloo'
# Rtools version 3.5
# Works: 2019.02.22
##########################

# Clear environment
rm(list=ls())

# Load package
library(rsq) # Version 1.0.1
library(lattice) # Version 0.20-35
library(DAAG) # Version 1.22
library(MASS) # Version 7.3-50
library(boot) # Version 1.3-20
library(survey) # Version 3.33-2
library(grid) # Version 3.5.1
library(Matrix) # Version 1.2-14
library(survival) # Version 2.42-3
library(mitools) # Version 2.3
library(relaimpo) # Version 2.2-3
library(Rcpp) # Version 0.12.18
library(ape) # Version 5.1
library(permute) # Version 0.9-4
library(vegan) # Version 2.5-2
library(paco) # Version 0.3.2
library(reshape2) # Version 1.4.3
library(igraph) # Version 1.2.2
library(nlme) # Version 3.1-137
library(picante) # Version 1.7
library(HiveR) # Version 0.3.42
library(cluster) # Version 2.0.7-1
library(gridExtra) # Version 2.3
library(MuMIn) # Version 1.42.1
library(ggplot2) # Version 3.0.0
library(cowplot) # Version 0.9.3
library(devtools) # Version 1.13.6
install_github('gastonstat/arcdiagram', force=F)
library(arcdiagram) # Version 0.1.11
library(viridisLite) # Version 0.3.0
library(viridis) # Version 0.5.1
library(coda) # Version 0.19-1

# Print session info
sessionInfo()

# Set working directory
# setwd('C:/Users/Clif/Dropbox/Eurobats')

# Set seed
set.seed(20180902)

##############################
### Phylogenetic diversity ###
##############################

# Read in Bartonella tree and host association table
bart = read.nexus('./Data/2_Eurobats_strict_IBD-SAMP-OUT.trees.txt')
bart_samp = read.nexus('./Data/2_Eurobats_strict_IBD-SAMP.trees.txt')
samp100 = seq(1, 9000, 90)
bart_samp100 = bart_samp[samp100]
hpmat = read.csv('./Data/hpmat_new.csv', row.names=1)
# Relabel columns with correct tree tip labels
colnames(hpmat) = bart$tip.label
# Calculate Faith's PD
hppd = pd(hpmat, bart)
# Calculate Faith's PD with tree uncertainty
faith_samp100 = as.data.frame(array(NA, c(100, 21)))
for(i in 1:100){
  colnames(hpmat) = bart_samp100[[i]]$tip.label
  hppd_samp100 = pd(hpmat, bart_samp100[[i]])
  faith_samp100[i,] = hppd_samp100$PD
}
colnames(faith_samp100) = rownames(hpmat)
interval = as.data.frame(array(NA, c(21, 3)))
for(i in 1:21){
  interval[i,] = quantile(faith_samp100[,i], probs=c(.025, .5, .975))
}
hppd_summary = data.frame(mean=apply(faith_samp100, 2, mean),
                          lower=interval[,1],
                          median=interval[,2],
                          upper=interval[,3])

# Enter additional data for the number of Bartonella OTUs, ectoparasites, and reviewed studies
hppd_summary$samples = hppd$SR
hppd_summary$OTUs = c(1, 5, 9, 4, 9, 1, 4, 8, 5, 5, 1, 2, 3, 3, 2, 2, 1, 2, 4, 6, 4)
hppd_summary$ectos = c(5, 11, 7, 8, 12, 6, 11, 15, 8, 14, 12, 11, 9, 7, 9, 3, 13, 4, 6, 10, 4)
hppd_summary$studies = c(20, 34, 45, 36, 49, 29, 24, 81, 17, 85, 23, 52, 24, 15, 48, 13, 42, 14, 23, 38, 12)
write.csv(hppd_summary, './Results/hppd.csv')

# Test correlations between phylogenetic diversity and covariates
# Test for number of samples
cor.test(hppd_summary$mean, log(hppd_summary$samples))
plot(log(hppd_summary$samples), hppd_summary$mean)
# Test for number of OTUs
cor.test(hppd_summary$mean, hppd_summary$OTUs)
plot(hppd_summary$OTUs, hppd_summary$mean)
# Test for number of ectoparasites
cor.test(hppd_summary$mean, hppd_summary$ectos)
plot(hppd_summary$ectos, hppd_summary$mean)
# Test for number of ectoparasites vs. number of reviewed studies
cor.test(hppd_summary$ectos, log(hppd_summary$studies))
plot(log(hppd_summary$studies), hppd_summary$ectos)
# Create linear models to explain phylogenetic diversity
test1 = lm(mean~log(samples)+ectos, hppd_summary); summary(test1)
test2 = lm(OTUs~log(samples)+ectos, hppd_summary); summary(test2)

####################################
### Ectoparasite sampling effort ###
####################################

# Read in ectoparasite sampling effort data
samp = read.csv('./Data/TableS9_new.csv', header=T)
# Test correlation between sampling and number of ectoparasites per species
cor.test(samp$log.count, samp$records.minus)
plot(samp$records.minus, samp$log.count)

########################################
### Community detection using igraph ###
########################################

# Import edge and weight list
data = as.data.frame(read.csv('./Data/Eurobats_host-vector_matrix - vect_adj_new.csv'))
# data = data[which(data$V3>=1),]
data = data[, c(1, 2, 4)]

# Create igraph network
data.net = graph.data.frame(data,directed=F)
# Assign colors based on vertex type: yellow, ectoparasites; red, Bartonella; green, bats
V(data.net)$color = c(rep('yellow', 17), rep('red', 20), rep('green', 21))

# igraph node labels
sp.order = c('Ept.nil', 'Ept.ser', 'Min.sch', 'Myo.bec', 'Myo.bly', 'Myo.cap', 'Myo.das',
            'Myo.dau', 'Myo.ema', 'Myo.myo', 'Myo.mys', 'Myo.nat', 'Nyc.noc', 'Pip.nat',
            'Pip.pip', 'Pip.pyg', 'Ple.aur', 'Rhi.bla', 'Rhi.eur', 'Rhi.fer', 'Rhi.meh')

# View igraph network
plot.igraph(data.net,axes=F,
            edge.arrow.size=0.2,
            vertex.size=15,
            vertex.label.cex=0.6,
            vertex.label.color='black',
            frame=F,layout=layout.gem)

# Community detection of igraph network
clusbet = cluster_edge_betweenness(data.net, weights=E(data.net)$weight)
clusgreed = cluster_fast_greedy(data.net, weights=E(data.net)$weight)
clusinf = cluster_infomap(data.net, e.weights=E(data.net)$weight, nb.trials=10)
cluslab = cluster_label_prop(data.net, weights=E(data.net)$weight)
cluseig = cluster_leading_eigen(data.net, weights=E(data.net)$weight)
cluslouv = cluster_louvain(data.net, weights=E(data.net)$weight)
clusopt = cluster_optimal(data.net, weights=E(data.net)$weight)
clusspin = cluster_spinglass(data.net, weights=E(data.net)$weight, spins=25)
cluswalk = cluster_walktrap(data.net, weights=E(data.net)$weight, steps=4)
# Combine communities into data frame
communities = as.data.frame(cbind(membership(clusbet), membership(clusgreed), membership(clusinf),
                                 membership(cluslab), membership(cluseig), membership(cluslouv),
                                 membership(clusopt), membership(clusspin), membership(cluswalk)))
colnames(communities)=c('betweenness', 'greedy', 'infomap',
                        'label', 'eigen', 'louvain',
                        'optimal', 'spinglass', 'walktrap')
# Calculate network modularity for selected algorithms
graph.modularity = c(modularity(clusinf), modularity(cluslouv), modularity(clusopt), modularity(clusspin))
write.csv(communities, './Results/network_identified_communities.csv')

# Plot Louvain communities
dev.off()
png('./Results/bart_comm.png', height=8, width=8, units='in', res=300)
plot(cluslouv, data.net,
     edge.width=25*E(data.net)$weight,
     vertex.size=10,
     vertex.label.cex=0.5,
     vertex.label.color='black',
     layout=layout.gem)
dev.off()

## Create community graph, edge weights are the number of edges
cg = contract.vertices(data.net, membership(clusinf))
cg2 = simplify(cg, remove.loops=TRUE)

for(i in 1:length(V(cg2)$name)){
  V(cg2)$name[[i]]=i
}

## Plot the community graph
plot(cg2, vertex.label=V(cg2)$name,
     edge.width=10*E(cg2)$weight, layout=layout.circle)
adj.cg2 = data.frame(get.edgelist(cg2))
adj.cg2$weight = E(cg2)$weight

# Extra calculations from igraph network
# Calculate degree and plot histogram
graph.degree = degree(graph=data.net, v=V(data.net), mode='total')
hist(graph.degree, xlab='Node degree', main='', breaks=25, las=1, freq=F)
# Calculate weighted degree and plot histogram
graph.wdegree = strength(graph=data.net, vids=V(data.net), mode='total', loops=T, weights=E(data.net)$weight)
hist(graph.wdegree, xlab='Node weighted degree', main='', breaks=25, las=1, freq=F)
# Calculate PageRank and plot histogram
graph.page = page.rank(graph=data.net, algo='prpack', vids=V(data.net),
                      directed=F, damping=0.85, weights=E(data.net)$weight)
hist(graph.page$vector ,xlab='Node PageRank', main='', breaks=25, las=1, freq=F)
# Make summary data frame
graph.summary = data.frame(cbind(graph.degree, graph.wdegree, graph.page$vector))
colnames(graph.summary) = c('degree','wd', 'page')
# Select nodes with highest weighted degree
graph.summary[which(graph.summary$wd>=quantile(graph.summary$wd, .75)),]

# Make Hiveplot for entire matrix
dev.off()
png('./Results/net.png', width=11, height=8.5, units='in', res=300)
net = edge2HPD(data)
net1 = mineHPD(net, option='rad <- tot.edge.count')
net2 = mineHPD(net1, option='axis <- source.man.sink')

# Edit node order and color externally and read in
net.data1 = read.csv('./Data/net2_nodes_new.csv', header=T)
net2$nodes$radius = as.numeric(net.data1$radius)
net2$nodes$size = 2
net2$nodes$color = as.character(net.data1$color)

# Edit edge color externally and read in
net.data2 = read.csv('./Data/net2_edges_new.csv',header=T)
net2$edges$weight = 25*net2$edges$weight
net2$edges$color = numeric(dim(net.data2)[1])
net2$edges$color = as.character(net.data2$color)

net2$axis.cols = 'white'
# Annotations for sample types
plotHive(net2,axLabs = c('bartonella', 'hosts', 'ectoparasites'),ch=0,bkgnd='white',
         axLab.gpar = gpar(col = 'black',fontsize = 20,lwd = 3))
# Annotations for communities
grid.text('Communities', x = .9, y = .95, default.units = 'npc', gp = gpar(fontsize = 24, col = 'black'))
grid.text('Min/Myo', x = .9, y = .915, default.units = 'npc', gp = gpar(fontsize = 20, col = '#e7298a'))
grid.text('Rhi', x = .9, y = .88, default.units = 'npc', gp = gpar(fontsize = 20, col = '#7570b3'))
grid.text('VespA', x = .9, y = .845, default.units = 'npc', gp = gpar(fontsize = 20, col = '#1b9e77'))
grid.text('VespB', x = .9, y = .81, default.units = 'npc', gp = gpar(fontsize = 20, col = '#d95f02'))
grid.text('VespC', x = .9, y = .775, default.units = 'npc', gp = gpar(fontsize = 20, col = '#66a61e'))
grid.text('VespD', x = .9, y = .74, default.units = 'npc', gp = gpar(fontsize = 20, col = '#e6ab02'))
grid.text('VespE', x = .9, y = .705, default.units = 'npc', gp = gpar(fontsize = 20, col = '#a6761d'))
# Annotations for abbreviated species names
grid.text(c('OTU31', 'OTU1', 'OTU2', 'OTU42', 'OTU11',
            'OTU16', 'OTU4', 'OTU25', 'OTU5', 'OTU3',
            'OTU24', 'OTU32', 'OTU22', 'OTU19', 'OTU26',
            'OTU33', 'OTU21', 'OTU27', 'EsSk', 'OTU20'),
          x=.465, y=seq(.523, .523+.0193*19, .0193),
          default.units = 'npc', gp = gpar(fontsize = 12, col = 'black'))
grid.text(c('Myo.bly', 'Myo.myo', 'Min.sch', 'Myo.cap', 'Rhi.fer', 'Rhi.eur', 'Rhi.meh',
            'Rhi.bla', 'Myo.dau', 'Myo.nat', 'Myo.bec', 'Myo.das', 'Myo.mys', 'Nyc.noc',
            'Pip.pip', 'Pip.nat', 'Pip.pyg', 'Ept.ser', 'Ept.nil', 'Ple.aur', 'Myo.ema'),
          x=seq(.504, .504-.0166*20, -0.0166), y=seq(.451, .451-.0097*20, -.0097),
          rot=300, default.units = 'npc', gp = gpar(fontsize = 12, col = 'black'))
grid.text(c('Spi.myo', 'Nyb.sch', 'Pen.duf', 'Pen.con', 'Pht.bia', 'Nyb.kol',
            'Bas.nan', 'Bas.nat', 'Arg.ves', 'Cim.pip', 'Isc.var', 'Spi.and',
            'Ixo.ari', 'Mac.sp4', 'Spi.sp2', 'Spi.kol', 'Spi.ple'),
          x=seq(.53, .53+.0168*16, .0168), y=seq(.531, .531-.0095*16, -.0095),
          rot=60, default.units = 'npc', gp = gpar(fontsize = 12, col = 'black'))
dev.off()

###################################################
### Cophylogeny analyses using PACo and ParaFit ###
###################################################

# Input host tree
pruneH = c('Eptesicus_nilssoni', 'Eptesicus_serotinus', 'Miniopterus_schreibersii',
          'Myotis_bechsteinii', 'Myotis_blythii', 'Myotis_capaccinii',
          'Myotis_dasycneme', 'Myotis_daubentonii', 'Myotis_emarginatus',
          'Myotis_myotis', 'Myotis_mystacinus', 'Myotis_nattereri',
          'Nyctalus_noctula', 'Pipistrellus_nathusii', 'Pipistrellus_pipistrellus',
          'Pipistrellus_pygmaeus', 'Plecotus_auritus', 'Rhinolophus_blasii',
          'Rhinolophus_euryale', 'Rhinolophus_ferrumequinum', 'Rhinolophus_mehelyi')
TreeH = read.nexus('./Data/Shi_et_al_2015_Evolution.nex')
# Prune tree to specific species
TreeH.prune = drop.tip(TreeH, TreeH$tip.label[-match(pruneH, TreeH$tip.label)])
plot(TreeH.prune)

# Input parasite tree
pruneP = c('Min.sch_B44655', 'Min.sch_B44692', 'Myo.bly_B44657', 'Rhi.eur_B44540', 'Rhi.meh_Pht.bia_77CJ1334221',
          'Rhi.fer_B44671', 'Rhi.eur_B44537', 'Myo.bly_B44620', 'Myo.ema_B44544', 'Myo.das_Cim.pip_264BB',
          'Ept.ser_Spi.sp1_26S83', 'Myo.dau_KF003129', 'Myo.bly_B44622', 'Rhi.fer_B44613', 'Myo.ema_B44735',
          'Nyc.noc_AJ871615', 'Pip.pyg_B44734', 'Myo.dau_Nyb.kol_2BF78', 'Myo.ema_B44617', 'Myo.bly_B44721')
TreeP = read.nexus('./Data/2_Eurobats_strict_IBD-SAMP-OUT.trees.txt')
# Prune tree to specific species
TreeP.prune = drop.tip(TreeP, TreeP$tip.label[-match(pruneP, TreeP$tip.label)])
plot(TreeP.prune)
# Import parasite species for uncertainty
TreeP_samp = read.nexus('./Data/2_Eurobats_strict_IBD-SAMP.trees.txt')
samp100 = seq(1, 9000, 90)
TreeP_samp100 = TreeP_samp[samp100]
# Prune trees to specific species for cophylogeny
TreeP_samp100.prune = NULL
for(i in 1:length(samp100)){
  TreeP_samp100.prune[[i]] = drop.tip(TreeP_samp100[[i]],
                                    TreeP_samp100[[i]]$tip.label[-match(pruneP,
                                                                   TreeP_samp100[[i]]$tip.label)])
}
plot(TreeP_samp100.prune[[1]])
# Prune trees to specific species for tip-association tests
Tree3 = read.tree('./Data/3_Eurobats_mafft_NJ.nex')
TreeP_samp100.pruneTA = NULL
for(i in 1:100){
  TreeP_samp100.pruneTA[[i]] = drop.tip(TreeP_samp100[[i]],
                                        TreeP_samp100[[i]]$tip.label[-match(Tree3$tip.label,
                                                                                    TreeP_samp100[[i]]$tip.label)])
}
write.nexus(TreeP_samp100.pruneTA, file='./Data/Eurobats_MCC_samp100.nex', translate=F)

# Compute patristic distance of host tree
H.order = c('Ept.nil', 'Ept.ser', 'Min.sch', 'Myo.bec', 'Myo.bly', 'Myo.cap', 'Myo.das',
           'Myo.dau', 'Myo.ema', 'Myo.myo', 'Myo.mys', 'Myo.nat', 'Nyc.noc', 'Pip.nat',
           'Pip.pip', 'Pip.pyg', 'Ple.aur', 'Rhi.bla', 'Rhi.eur', 'Rhi.fer', 'Rhi.meh')
host.D = cophenetic (TreeH.prune)
# Adjust distances to be between 0 and 1
host.D = host.D/max(host.D)
# Reorder rows and columns
host.D = host.D[pruneH, pruneH]
rownames(host.D) = H.order; colnames(host.D) = H.order
# Rename tips
TreeH.prune$tip.label = c('Pip.pyg', 'Pip.pip', 'Nyc.noc', 'Pip.nat','Ept.nil', 'Ept.ser', 'Ple.aur',
                        'Myo.das', 'Myo.cap', 'Myo.nat', 'Myo.bly', 'Myo.myo', 'Myo.bec', 'Myo.dau',
                        'Myo.ema', 'Myo.mys', 'Min.sch', 'Rhi.fer', 'Rhi.eur', 'Rhi.meh', 'Rhi.bla')
plot(TreeH.prune)

# Compute patristic distance of parasite tree
P.order = c('OTU1', 'OTU2', 'OTU3', 'OTU4', 'OTU5',
            'OTU11', 'OTU16', 'OTU19', 'OTU20', 'OTU21',
            'EsSk', 'OTU22', 'OTU24', 'OTU25', 'OTU26',
            'OTU27', 'OTU31', 'OTU32', 'OTU33', 'OTU42')
para.D = cophenetic(TreeP.prune)
# Adjust distances to be between 0 and 1
para.D = para.D/max(para.D)
# Reorder rows and columns
para.D = para.D[pruneP, pruneP]
rownames(para.D) = P.order; colnames(para.D) = P.order
# Rename tips
TreeP.prune$tip.label = c('OTU42', 'OTU33', 'OTU31', 'OTU32', 'OTU3', 
                        'OTU2', 'OTU1', 'OTU4', 'OTU5', 'OTU16',
                        'OTU11', 'OTU26', 'OTU27', 'EsSk', 'OTU21',
                        'OTU20', 'OTU19', 'OTU22', 'OTU25', 'OTU24')
plot(TreeP.prune)
# Compute patristic distance of parasite trees for uncertainty
para.D_samp100 = NULL
for(i in 1:length(samp100)){
  para.D_samp100[[i]] = cophenetic(TreeP_samp100.prune[[i]])
}
# Adjust distances to be between 0 and 1
for(i in 1:length(samp100)){
  para.D_samp100[[i]] = para.D_samp100[[i]]/max(para.D_samp100[[i]])
}
# Reorder rows and columns
for(i in 1:length(samp100)){
  para.D_samp100[[i]] = para.D_samp100[[i]][pruneP, pruneP]
}

# Input host association matrix
HP = as.matrix(read.csv('./Data/host_associations_OTUs_new.csv', row.names=1, header=T)) 
rownames(HP) = H.order
# Reorder rows and columns
host.D = host.D[rownames(HP), rownames(HP)]
para.D = para.D[colnames(HP), colnames(HP)]
# Input host association matrix for uncertainty
HP_samp100 = as.matrix(read.csv('./Data/host_associations_OTUs_new_uncertainty.csv',
                                row.names=1, header=T)) 
rownames(HP_samp100) = H.order
# Reorder rows and columns
host.D = host.D[rownames(HP_samp100), rownames(HP_samp100)]
for(i in 1:length(samp100)){
  para.D_samp100[[i]] = para.D_samp100[[i]][colnames(HP_samp100), 
                                            colnames(HP_samp100)]
}

# Run PACo analysis
D = prepare_paco_data(host.D, para.D, HP)
D = add_pcoord(D, correction='cailliez')
D = PACo(D, nperm=1000, seed=20180902, method='r0', symmetric=F, proc.warnings=T, shuffled=T)
D = paco_links(D)
# Plot null versus observed
hist(D$shuffled, xlim=c(21, 24))
abline(v=D$gof$ss, col='red')
# Goodness of fit and p-value
D$gof
# Output residuals and jackknife values for individual links
D.resid = as.data.frame(residuals_paco(D$proc, type='interaction')); colnames(D.resid) = 'resid'
D.resid$mean = as.matrix(D$jackknife$mean)
D.resid$upper = as.matrix(D$jackknife$upper)
# Run PACo analysis for uncertainty
gof_samp100 = as.data.frame(array(NA, c(length(samp100), 2)))
for(i in 1:length(samp100)){
  D_samp100 = prepare_paco_data(host.D, para.D_samp100[[i]], HP_samp100)
  D_samp100 = add_pcoord(D_samp100, correction='cailliez')
  D_samp100 = PACo(D_samp100, nperm=1000, seed=20180902, method='r0',
                   symmetric=F, proc.warnings=T, shuffled=T)
  gof_samp100[i, 1] = D_samp100$gof$ss
  gof_samp100[i, 2] = D_samp100$gof$p
}
colnames(gof_samp100) = c('ss', 'p')
sum(gof_samp100$p<0.05)

# ParaFit calculations
PF.out = parafit(host.D, para.D, HP, nperm=999, test.links=T, correction='cailliez', silent=F)
# Goodness of fit and p-value
PF.out$ParaFitGlobal; PF.out$p.global
# Output individual link tests
PF.test = as.data.frame(PF.out$link.table)
PF.test = PF.test[order(PF.test$Parasite),]
# ParaFit calculations for uncertainty
ParaFit_samp100 = as.data.frame(array(NA, c(length(samp100), 2)))
for(i in 1:length(samp100)){
  PF.out_samp100 = parafit(host.D, para.D_samp100[[i]], HP, nperm=999,
                   test.links=F, correction='cailliez', silent=F)

  ParaFit_samp100[i, 1] = PF.out_samp100$ParaFitGlobal
  ParaFit_samp100[i, 2] = PF.out_samp100$p.global
}
colnames(ParaFit_samp100) = c('ParaFitGlobal', 'p.global')
sum(ParaFit_samp100$p.global<0.05)

# Output test significance results
D.resid = cbind(D.resid,PF.test)
D.resid$col = as.character('grey')
D.resid[which(D.resid$p.F1<0.1),]$col = 'yellow'
D.resid[which(D.resid$p.F1<0.05),]$col = 'green'
D.resid[which(D.resid$p.F1<0.01),]$col = 'green3'
D.resid[which(D.resid$p.F1<0.005),]$col = 'green4'
col = c('#00000000','#00000033','#000000CC','#000000FF')
D.resid$sig = as.character('#00000033')
D.resid[which((D.resid$mean^2)<mean(D.resid$resid^2)),]$sig = col[3]
D.resid[which((D.resid$upper^2)<mean(D.resid$resid^2)),]$sig = col[4]
write.csv(D.resid, './Results/PACo_resid.csv')
PACo.cols = read.csv('./Results/PACo_colors_new.csv', row.names=1, header=T)

# Barplot of residuals
dev.off()
png('./Results/bar_resid.png',width=11,height=8.5,units='in',res=300)
par(mar=c(5, 5, 1, 1))
PACo.cols2 = PACo.cols[order(PACo.cols$order),]
bar.resid = barplot(PACo.cols2$resid,las=2, names=rownames(PACo.cols2), cex.lab=1.5, cex.axis=1.5,
                   cex.names=0.6, ylim=c(0,1), ylab='PACo residuals')
abline(h=mean(PACo.cols$resid), lty=2, lwd=3, col='red')
dev.off()

# Barplot of jackknifed contributions
png('./Results/bar_jack.png', width=11, height=8.5, units='in', res=300)
par(mar=c(5, 5, 1, 1))
PACo.cols2 = PACo.cols[order(PACo.cols$order),]
bar.jack = barplot(as.vector(PACo.cols2$mean), las=2, names=rownames(PACo.cols2),
                  col=as.character(PACo.cols2$b.col), cex.axis=1.5, cex.lab=1.5,
                  cex.names=0.6, ylim=c(0,1), ylab='Squared residuals')
arrows(bar.jack, PACo.cols2[,7], bar.jack, PACo.cols2[,6], angle=90, code=1, length=0.05)
abline(h=mean(PACo.cols2$resid^2), lty=2, lwd=3, col='red')
MyCols = c('#e7298a', '#7570b3', '#1b9e77', '#d95f02',
           '#66a61e', '#e6ab02', '#a6761d', '#bebebe')
MyLabs = c('Min/Myo', 'Rhi', 'VespA', 'VespB',
           'VespC', 'VespD', 'VespE', 'not in same community')
Nfact = 8
Nrows = 2
Ncols = ceiling(Nfact/Nrows)
MyOrder = as.vector(matrix(1:(Nrows*Ncols), nrow=Nrows, ncol=Ncols, byrow=F))
legend('top',
       legend=MyLabs[MyOrder],
       col=MyCols[MyOrder],
       pch=15, cex=1.2, bty='n', ncol=Ncols, title='Communities')
dev.off()

# Reorder links by parasite
PACo.cols3 = PACo.cols[order(PACo.cols$Parasite),]

# Host ordination matrices
png('./Results/ordination.png', width=11, height=8.5, units='in', res=300)
par(mar=c(5, 5, 1, 1))
HostX = D$proc$X
ParY = D$proc$Yrot
# Plotting host and parasite ordinations
ordination = plot(HostX, pch=2, col='black',
                 xlab='PACo axis 1', ylab='PACo axis 2', las=1, cex.lab=1.5, cex.axis=1.5)
arrows(ParY[,1], ParY[,2], HostX[,1], HostX[,2], length=0.12, angle=15, xpd=F,
       col=as.character(PACo.cols3$c.col), lwd=25*PACo.cols3$weight)
points(HostX, pch=2, col='black')
points(ParY, pch=1, col='black')
# Unique() removes duplicated points
HostX = unique(D$proc$X)
ParY = unique(D$proc$Y)
host.groups = c('Vesp', 'Min', 'Myo', 'Myo', 'Vesp', 'Rhi', 'Myo',
                'Myo', 'Myo', 'Rhi', 'Rhi', 'Rhi', 'Myo', 'Myo',
                'Vesp', 'Vesp', 'Vesp', 'Vesp', 'Vesp', 'Myo', 'Myo')
ordihull(ord=HostX, groups=host.groups, display='sites', col='black')
legend(0.1, 0.43,
       legend=c('Min/Myo', 'Rhi', 'VespA', 'VespB',
                'VespC', 'VespD', 'VespE', 'not in same community'),
       col=c('#e7298a', '#7570b3', '#1b9e77', '#d95f02',
             '#66a61e', '#e6ab02', '#a6761d', '#bebebe'),
       lty=1, lwd=3, cex=1.2, bty='n', title='Communities')
legend(0.34, 0.43,
       legend=c('0.005', '0.05', '0.5'),
       lty=1, lwd=25*c(.01, 0.1, 0.5), cex=1.2, bty='n', title='Edge weight')
legend(0.49, 0.43,
       legend=c('not supported', 'supported', 'highly supported'),
       lty=1, lwd=3, col=c('#00000020', '#00000080', '#000000'), cex=1.2, bty='n', title='Link support')
text(0, 0.32, 'Vespertilioninae', col='black', cex=1.5)
text(-0.1, -0.2, 'Myotinae', col='black', cex=1.5)
text(0.1, 0.13, 'Miniopteridae', col='black', cex=1.5)
text(0.62, -0.03, 'Rhinolophidae', col='black', cex=1.5)
dev.off()

# Second association matrix
association = read.csv('./Data/associations_cophylo_new.csv', header=T)

# Cophyloplot
cophyloplot(TreeH.prune, TreeP.prune, assoc=association[,1:2], length.line=0, space=50,
            col=as.character(association[,3]), gap=3, rotate=T,
            lwd=25*association[,4], font=1)
# Must output PDF file manually

###########################
### Regression analyses ###
###########################

# Set host species order
sp.order = c('Ept.nil', 'Ept.ser', 'Min.sch', 'Myo.bec', 'Myo.bly', 'Myo.cap','Myo.das',
            'Myo.dau', 'Myo.ema', 'Myo.myo', 'Myo.mys', 'Myo.nat', 'Nyc.noc', 'Pip.nat',
            'Pip.pip', 'Pip.pyg', 'Ple.aur', 'Rhi.bla', 'Rhi.eur', 'Rhi.fer', 'Rhi.meh')

# Read in Bartonella data
bart.data = read.csv('./Data/Eurobats_host-vector_matrix - bart_count_new.csv',
                     header=T, row.names=1)

# Read in ectoparasite data
ecto.data = read.csv('./Data/Eurobats_host-vector_matrix - ecto_count_new.csv',
                     header=T, row.names=1)

# Filter data frames to remove rows without data
bart.data = bart.data[which(rownames(bart.data) %in% rownames(ecto.data)),]
bart.data = bart.data[, which(apply(bart.data, 2, sum)!=0)]
ecto.data = ecto.data[which(rownames(ecto.data) %in% rownames(bart.data)),]
ecto.data = ecto.data[, which(apply(ecto.data, 2, sum)!=0)]

# Choose one of the following for Bartonella community dissimilarity calculation and comment out others:
# bart.dist = as.matrix(vegdist(bart.data, 'binomial', binary=F))
# bart.dist = as.matrix(vegdist(bart.data, 'cao', binary=F))
bart.dist = 1-cor(t(bart.data), method='spearman')
# Remove lower triangle (redundant because symmetrical)
bart.dist[lower.tri(bart.dist)] = NA
# Melt data into long format
bart = melt(bart.dist)
# Filter out values for self dissimilarity and missing values (from lower triangle)
bart = bart[which(bart$Var1!=bart$Var2 & bart$value!='NA'),]
# Assign nicknames to values
bart$nick = paste(bart[,1], bart[,2], sep='_')

# Choose one of the following for ectoparasite community dissimilarity calculation and comment out others:
# ecto.dist = as.matrix(vegdist(ecto.data, 'binomial', binary=F))
# ecto.dist = as.matrix(vegdist(ecto.data, 'cao', binary=F))
ecto.dist = 1-cor(t(ecto.data), method='spearman')
# Remove lower triangle (redundant because symmetrical)
ecto.dist[lower.tri(ecto.dist)] = NA
# Melt data into long format
ecto = melt(ecto.dist)
# Filter out values for self dissimilarity and missing values (from lower triangle)
ecto = ecto[which(ecto$Var1!=ecto$Var2 & ecto$value!='NA'),]
# Assign nicknames to values
ecto$nick = paste(ecto[,1], ecto[,2], sep='_')

# Read in host species tree
shi = read.nexus('./Data/Shi_et_al_2015_Evolution.nex')
shi.species = c('Eptesicus_nilssoni', 'Eptesicus_serotinus', 'Miniopterus_schreibersii',
               'Myotis_bechsteinii', 'Myotis_blythii', 'Myotis_capaccinii',
               'Myotis_dasycneme', 'Myotis_daubentonii', 'Myotis_emarginatus',
               'Myotis_myotis', 'Myotis_mystacinus', 'Myotis_nattereri',
               'Nyctalus_noctula', 'Pipistrellus_nathusii', 'Pipistrellus_pipistrellus',
               'Pipistrellus_pygmaeus', 'Plecotus_auritus', 'Rhinolophus_blasii',
               'Rhinolophus_euryale', 'Rhinolophus_ferrumequinum', 'Rhinolophus_mehelyi')
# Filter to the chosen host species
prune.shi = drop.tip(shi, shi$tip.label[-match(shi.species, shi$tip.label)])
plot(prune.shi)
write.nexus(prune.shi, file='./Results/prune_shi.nex', translate=T)
# Calculate tree distance and adjust branch lengths to be proportion of maximum branch length
dist.shi = cophenetic.phylo(prune.shi)
dist.shi = dist.shi/max(dist.shi)
dist.shi = dist.shi[shi.species, shi.species]
rownames(dist.shi) = sp.order; colnames(dist.shi) = sp.order
# Assign host.mean
host.mean = dist.shi
write.csv(host.mean, './Results/host_mean.csv')

# Filter host.mean to the species that are in Bartonella community data frame
host.mean = host.mean[which(rownames(host.mean) %in% rownames(bart.data)),]
host.mean = host.mean[, which(colnames(host.mean) %in% rownames(bart.data))]
# Remove lower triangle (redundant because symmetrical)
host.mean[lower.tri(host.mean)] = NA
# Melt data into long format
host = melt(host.mean)
# Filter out values for self dissimilarity and missing values (from lower triangle)
host = host[which(host$Var1!=host$Var2 & host$value!='NA'),]
# Assign nicknames to values
host$nick = paste(host[,1], host[,2], sep='_')

# Read in bat species geographic overlap data
over.data = as.matrix(read.csv('./Data/Bat_species_overlap_pct.csv'), header=T)
rownames(over.data) = sp.order; colnames(over.data) = sp.order

# Filter over.data to the species that are in Bartonella community data frame
over.data = over.data[which(rownames(over.data) %in% rownames(bart.data)),]
over.data = over.data[, which(colnames(over.data) %in% rownames(bart.data))]
# Remove lower triangle (redundant because symmetrical)
over.data[lower.tri(over.data)] = NA
# Melt data into long format
over = melt(over.data)
# Filter out values for self dissimilarity and missing values (from lower triangle)
over = over[which(over$Var1!=over$Var2 & over$value!='NA'),]
# Assign nicknames to values
over$nick = paste(over[,1],over[,2],sep='_')

# Read in bat species co-roosting data
roost.data = read.csv('./Data/table_batroosting.csv', header=T, row.names=1)
# Reorder rows and columns
roost.data = roost.data[sp.order, sp.order]
rownames(roost.data) = sp.order; colnames(roost.data) = sp.order

# Filter roost.data to the species that are in Bartonella community data frame
roost.data = roost.data[which(rownames(roost.data) %in% rownames(bart.data)),]
roost.data = roost.data[, which(colnames(roost.data) %in% rownames(bart.data))]
# Remove lower triangle (redundant because symmetrical)
roost.data[lower.tri(roost.data)] = NA
# Melt data into long format
roost = melt(as.matrix(roost.data))
# Filter out values for self dissimilarity and missing values (from lower triangle)
roost = roost[which(roost$Var1!=roost$Var2 & roost$value!='NA'),]
# Assign nicknames to values
roost$nick = paste(roost[,1],roost[,2],sep='_')

# Create data frame
df = as.data.frame(cbind(
  as.numeric(as.character(bart$value)),
  as.numeric(as.character(ecto$value)),
  as.numeric(as.character(host$value)),
  as.numeric(as.character(over$value))))
# Rescale data to standard normal
df = data.frame(apply(df, 2, scale))
df = cbind(df, roost$value)
colnames(df) = c('bart', 'ecto', 'host', 'over', 'roost')
# Examine correlations
pairs(df)

# Test model for Bartonella community dissimilarity
reg1 = lm(bart~host+over+roost, data=df,  na.action=na.fail)
# Examine model fit (residual vs. fitted, Q-Q, etc.)
par(mfrow=c(2,2)); plot(reg1); par(mfrow=c(1,1))
# Histogram of model residuals
hist(reg1$resid, breaks=20)
# Shapiro-Wilk test for normality of residuals
shapiro.test(reg1$residuals)
# Model R2 and adjusted R2
R2 = rsq(reg1, adj=F); R2
R2adj = rsq(reg1, adj=T); R2adj
# Model summary, confidence intervals, etc.
summary(reg1)
confint(reg1)
summary(aov(reg1))
# Model partial R2 and adjusted R2
pR2 = rsq.partial(reg1, adj=F); pR2
pR2adj = rsq.partial(reg1, adj=T); pR2adj
# Model cross-validation
cv.lm(bart~host+over+roost, data=df, m=10)
# Model covariate relative importance and bootstrapping
calc.relimp(reg1, type='lmg', rela=T)
boot1 = boot.relimp(reg1, b=1000, type='lmg', rank=T, diff=T, rela=T)
booteval.relimp(boot1)
plot(booteval.relimp(boot1, sort=T))
# Model selection by AICc
dreg1 = dredge(reg1); dreg1
mreg1 = model.avg(dreg1, subset=delta<=2)
summary(mreg1)
confint(mreg1)

# Test model for ectoparasite community dissimilarity
reg2 = lm(ecto~host+over+factor(roost), data=df, na.action=na.fail)
# Examine model fit (residual vs. fitted, Q-Q, etc.)
par(mfrow=c(2,2)); plot(reg2); par(mfrow=c(1,1))
# Histogram of model residuals
hist(reg2$resid, breaks=20)
# Shapiro-Wilk test for normality of residuals
shapiro.test(reg2$residuals)
# Model R2 and adjusted R2
R2 = rsq(reg2, adj=F); R2
R2adj = rsq(reg2, adj=T); R2adj
# Model summary, confidence intervals, etc.
summary(reg2)
confint(reg2)
summary(aov(reg2))
# Model partial R2 and adjusted R2
pR2 = rsq.partial(reg2, adj=F); pR2
pR2adj = rsq.partial(reg2, adj=T); pR2adj
# Model cross-validation
cv.lm(ecto~host+over+roost, data=df, m=10)
# Model covariate relative importance and bootstrapping
calc.relimp(reg2, type='lmg', rela=T)
boot2 = boot.relimp(reg2, b=1000, type='lmg', rank=T, diff=T, rela=T)
booteval.relimp(boot2)
plot(booteval.relimp(boot2, sort=T))
# Model selection by AICc
dreg2 = dredge(reg2); dreg2
mreg2 = model.avg(dreg2, subset=delta<=2)
summary(mreg2)
confint(mreg2)

# Read in bat host species Bartonella host-switching rates
rates = as.data.frame(read.csv('./Data/3_Eurobats_all_bat_rates_logspprates.csv'))
# Test correlation between median rates and least sampled host species
cor.test(rates$median, log(rates$least))
plot(log(rates$least), rates$median)
# Test correlation between median rates and Bayes factors
cor.test(rates$median, 2*log(rates$BF))
plot(rates$median, 2*log(rates$BF))
# Filter data to make sure same species are in each data vector
ectosub = ecto[which(ecto$nick %in% rates$nick),]
ectosub = ectosub[order(ectosub$nick),]
hostsub = host[which(host$nick %in% rates$nick),]
hostsub = hostsub[order(hostsub$nick),]
oversub = over[which(over$nick %in% rates$nick),]
oversub = oversub[order(oversub$nick),]
roostsub = roost[which(roost$nick %in% rates$nick),]
roostsub = roostsub[order(roostsub$nick),]
rates$ecto = ectosub$value
rates$host = hostsub$value
rates$over = oversub$value
rates$roost = roostsub$value
# Make data frame
df2 = data.frame(mean=rates$mean,
                 median=rates$median,
                 least=log(rates$least),
                 ecto=rates$ecto,
                 host=rates$host,
                 over=rates$over)
# Rescale data to standard normal
df2 = data.frame(apply(df2, 2, scale))
df2$roost = rates$roost
# Examine correlations
pairs(df2)

# Test correlation between median and mean rates
cor.test(df2$mean, df2$median)
# Test model for Bartonella host-switching rates
reg3 = lm(median~least+host, data=df2, na.action=na.fail)
# Examine model fit (residual vs. fitted, Q-Q, etc.)
par(mfrow=c(2,2)); plot(reg3); par(mfrow=c(1,1))
# Histogram of model residuals
hist(reg3$resid, breaks=20)
# Shapiro-Wilk test for normality of residuals
shapiro.test(reg3$residuals)
# Model R2 and adjusted R2
R2 = rsq(reg3, adj=F); R2
R2adj = rsq(reg3, adj=T); R2adj
# Model summary, confidence intervals, etc.
summary(reg3)
confint(reg3)
summary(aov(reg3))
# Model partial R2 and adjusted R2
pR2 = rsq.partial(reg3, adj=F); pR2
pR2adj = rsq.partial(reg3, adj=T); pR2adj
# Model cross-validation
cv.lm(median~least+host, data=df2, m=10)
# Model covariate relative importance and bootstrapping
calc.relimp(reg3, type='lmg', rela=T)
boot = boot.relimp(reg3, b=1000, type='lmg', rank=T, diff=T, rela=T)
booteval.relimp(boot)
plot(booteval.relimp(boot))
# Model selection by AICc
dreg3 = dredge(reg3); dreg3
mreg3 = model.avg(dreg3, subset=delta<=2)
summary(mreg3)
confint(mreg3)

# Plots for linear regression summary
# Panel A: Bartonella dissimilarity
A1 = ggplot(data=df, aes(x=host, y=bart)) +
  geom_abline(intercept=reg1$coef[1], slope=reg1$coef[2],
              col='#E69F00', size=2) +
  geom_point(size=3, shape=1, alpha=0.4) +
  xlim(-4, 2) + xlab('Host phylogenetic distance') +
  ylim(-3, 2) + ylab('Bartonella dissimilarity') +
  theme_cowplot(font_size=12)
A2 = ggplot(data=df, aes(x=over, y=bart)) +
  geom_abline(intercept=reg1$coef[1], slope=reg1$coef[3],
              col='#56B4E9', size=2) +
  geom_point(size=3, shape=1, alpha=0.4) +
  xlim(-2, 4) + xlab('Host geographic range overlap') +
  ylim(-3, 2) + ylab('') +
  theme_cowplot(font_size=12)
A3 = ggplot(data=df, aes(x=as.factor(roost), y=bart)) +
  geom_violin() +
  geom_jitter(width=0.25, size=3, shape=1, alpha=0.4) +
  scale_x_discrete(labels=c('no sharing', 'sharing'), name='Host roost sharing') +
  ylim(-3, 2) + ylab('') +
  theme_cowplot(font_size=12)
# Combine panel A plots
A = plot_grid(A1, A2, A3, ncol=3)
# Panel B: Ectoparasite dissimilarity
B1 = ggplot(data=df, aes(x=host, y=ecto)) +
  geom_abline(intercept=reg1$coef[1], slope=reg2$coef[2],
              col='#E69F00', size=2) +
  geom_point(size=3, shape=1, alpha=0.4) +
  xlim(-4, 2) + xlab('Host phylogenetic distance') +
  ylim(-2, 3) + ylab('Ectoparasite dissimilarity') +
  theme_cowplot(font_size=12)
B2 = ggplot(data=df, aes(x=over, y=ecto)) +
  geom_abline(intercept=reg1$coef[1], slope=reg2$coef[3],
              col='#56B4E9', size=2) +
  geom_point(size=3, shape=1, alpha=0.4) +
  xlim(-2, 4) + xlab('Host geographic range overlap') +
  ylim(-2, 3) + ylab('') +
  theme_cowplot(font_size=12)
B3 = ggplot(data=df, aes(x=as.factor(roost), y=ecto)) +
  geom_violin() +
  geom_jitter(width=0.25, size=3, shape=1, alpha=0.4) +
  scale_x_discrete(labels=c('no sharing', 'sharing'), name='Host roost sharing') +
  ylim(-2, 3) + ylab('') +
  theme_cowplot(font_size=12)
# Combine panel B plots
B = plot_grid(B1, B2, B3, ncol=3)
# Panel C: Bartonella host-switching rates
C1 = ggplot(data=df2, aes(x=host, y=median)) +
  geom_abline(intercept=reg3$coef[1], slope=reg3$coef[3],
              col='#E69F00', size=2) +
  geom_point(size=3, shape=1, alpha=0.4) +
  xlim(-3, 2) + xlab('Host phylogenetic distance') +
  ylim(-2, 3) + ylab('Bartonella host-switching rate') +
  theme_cowplot(font_size=12)
C2 = ggplot(data=df2, aes(x=least, y=median)) +
  geom_abline(intercept=reg3$coef[1], slope=reg3$coef[2],
              col='#009E73', size=2) +
  geom_point(size=3, shape=1, alpha=0.4) +
  xlim(-2, 2) + xlab('Least sampled bat species') +
  ylim(-2, 3) + ylab('') +
  theme_cowplot(font_size=12)
# Combine panel C plots
C = plot_grid(C1, C2, ncol=2)
# Combine all panels
dev.off()
png('./Results/lin_reg.png', width=8.5, height=11, units='in', res=300)
plot_grid(A, B, C, labels=c('A', 'B', 'C'), label_size=24, nrow=3)
dev.off()

png('./Results/lin_reg_C.png', width=10, height=5, units='in', res=300)
par(mfrow=c(1, 2), mar=c(5, 5, 1, 1))
plot(df2$host, df2$median, las=1, pch=19, cex.axis=1.2, cex.lab=1.2,
     xlab='Host phylogenetic distance', ylab='Bartonella host-switching rate')
abline(a=reg3$coef[1], b=reg3$coef[3], col='#E69F0080', lwd=3)
plot(df2$least, df2$median, las=1, pch=19, cex.axis=1.2, cex.lab=1.2,
     xlab='Least sampled bat species', ylab='')
abline(a=reg3$coef[1], b=reg3$coef[2], col='#009E7380', lwd=3)
dev.off()

# Make transition rates network
# Import edge and weight list for median rates
rates = as.data.frame(read.csv('./Data/Eurobats_median_rates.csv'), header=T)
rates.net = graph.data.frame(rates, directed=F)
rates.edges = get.edgelist(rates.net)

# Arcplot for median rates
png('./Results/Eurobats_median_rates.png', width=8.5, height=11, units='in', res=300)
v.cols = c('#56b4e9', '#0072b2', '#a65628', '#e41a1c')
col=c('#00000000', '#00000033', '#000000AA', '#000000FF')
arcplot(rates.edges, show.nodes=T, horizontal=F, lwd.arcs=5,
        col.arcs=c(col[2], col[3], col[2], col[3], col[2], col[2], col[2], col[2], col[2], col[3],
                   col[2], col[2], col[2], col[3], col[3], col[4], col[2], col[2], col[2], col[2],
                   col[2], col[3], col[2], col[2], col[2], col[2], col[2], col[2], col[3], col[3]),
        ordering=c(21, 15, 12, 20, 9, 17, 6, 14, 5, 4, 3, 11, 10, 13, 8, 16, 19, 7, 18, 2, 1),
        col.labels='black', cex.labels=1.1,
        col.nodes=c(v.cols[1], v.cols[1], v.cols[2], v.cols[2], v.cols[2], v.cols[2], v.cols[1], 
                    v.cols[1], v.cols[3], v.cols[2], v.cols[2], v.cols[4], v.cols[2], v.cols[2],
                    v.cols[4], v.cols[1], v.cols[2], v.cols[1], v.cols[1], v.cols[4], v.cols[4]),
        lwd.nodes=10)
legend('topright',legend=c('<1', '1 to 1.5', '>1.5'), cex=1.4,
       lty=1, lwd=5,
       col=col[2:4],
       bty='n',
       title='Median rates')
legend('bottomright', legend=c('Vespertilioninae (Vespertilionidae)',
                            'Myotinae (Vespertilionidae)',
                            'Miniopteridae', 'Rhinolophidae'),
       cex=1.4, pch=19,
       col=v.cols[1:4],
       bty='n',
       title='Bat taxonomy')
dev.off()

# Import edge and weight list for BF
BF = as.data.frame(read.csv('./Data/Eurobats_Bayes_factors.csv'))
BF.net = graph.data.frame(BF, directed=F)
BF.edges = get.edgelist(BF.net)

# Arcplot for Bayes factors
png('./Results/Eurobats_Bayes_factors.png', width=8.5, height=11, units='in', res=300)
arcplot(BF.edges, show.nodes=T, horizontal=F, lwd.arcs=5,
        col.arcs=c(col[4], col[3], col[2], col[3], col[2], col[3], col[2], col[2], col[2], col[3],
                   col[3], col[2], col[2], col[3], col[4], col[4], col[3], col[2], col[2], col[2],
                   col[2], col[3], col[2], col[3], col[2], col[3], col[4], col[4], col[4], col[4]),
        ordering=c(21, 15, 12, 20, 9, 17, 6, 14, 5, 4, 3, 11, 10, 13, 8, 16, 19, 7, 18, 2, 1),
        col.labels='black', cex.labels=1.1,
        col.nodes=c(v.cols[1], v.cols[1], v.cols[2], v.cols[2], v.cols[2], v.cols[2], v.cols[1],
                    v.cols[1], v.cols[3], v.cols[2], v.cols[2], v.cols[4], v.cols[2], v.cols[2],
                    v.cols[4], v.cols[1], v.cols[2], v.cols[1], v.cols[1], v.cols[4], v.cols[4]),
        lwd.nodes=10)
legend('topright', legend=c('<6', '6 to 10', '>10'), cex=1.4,
       lty=1, lwd=5,
       col=col[2:4],
       bty='n',
       title='Bayes factors (2 ln K)')
dev.off()

###################
### End of code ###
###################