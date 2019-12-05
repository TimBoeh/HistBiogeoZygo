## Boehnert et al. 2019
## Historical Biogeography of South American Zygophyllaceae
## Published in Frontiers of Biogeography

# ######## install.packages() ############
# ::::: installation only need ones ::::::
install.packages("optimx", dependencies=TRUE, repos="http://cran.rstudio.com")
install.packages("snow")
install.packages("phylobase", dependencies=TRUE, repos="http://cran.rstudio.com")
install.packages("vegan")
install.packages("devtools")
install.packages('rexpokit')
install.packages('cladoRcpp')
install.packages("glue")
install.packages("Rcpp")
install.packages("roxygen2")
install.packages("fields")
install.packages("readPNG")
install.packages("png")
install.packages("colorspace")
install.packages("strap")

## Install the updated BioGeoBEARS from GitHub:
library(devtools)
devtools::install_github(repo = "nmatzke/BioGeoBEARS", INSTALL_opts = "--byte-compile", dependencies=T)

## Install Phyloch from fmichonneau GitHub repo; this is the same package as from Christophs homepage:
devtools::install_github(repo = "fmichonneau/phyloch", INSTALL_opts = "--byte-compile", dependencies=T)
#install.packages(pkgs="http://www.christophheibl.de/phyloch_1.5-5.zip", repos=NULL, type="source")

# ############ library() #############
# :::::::::: load packages  ::::::::::
library(BioGeoBEARS)
library(cladoRcpp)
library(rexpokit)
library(parallel)
library(optimx)
library(ape)
library(strap)    #for ploting geologic time scale
library(phyloch)
library(png)

# ############ setwd() ####################
# :::::  set Working Directory  :::::::::::
setwd("C:/YOUR/PASS/TO/WORKING/DIRECTORY")     ## for Windows
setwd("/home/YOUR/PASS/TO/WORKING/DIRECTORY")  ## for Linux
getwd()
list.files()

## ######### PLOT BEAST TREEs #########
## Load two Trees with two different calibrations schemes
## read Best tree with 4 secondary calibrations
ZygoBeastTree <- read.beast("TREE-FILE.tre")
ZygoBeastTree$root.time <- ZygoBeastTree$height[1]    # without this line you get error: "tree$root.time is missing, check tree is time scaled."
margin <- branching.times(ZygoBeastTree)[1]-85.8      # set margin according to oldest node

## set colors for clade rect() function
FagColT <- rgb(0.11, 0.67, 0.57, alpha = 0.3); FagCol <- rgb(0.11, 0.67, 0.57)
LarrColT <- rgb(0.85, 0.71, 0.11, alpha = 0.3); LarrCol <- rgb(0.85, 0.71, 0.11)

## Plot the dated phylogeny (4 calibrations) with HPD Bars
pdf("Zygophyllales_Beast_Tree_Appendix_4cal.pdf", height = 11, width = 8)
geoscalePhylo(tree = ladderize(ZygoBeastTree, right = F), boxes = "Epoch", units = c("Period", "Epoch"),
              cex.age = 0.7, cex.ts = 0.75, cex.tip = 0.4, x.lim = c(-15,85), y.lim = c(4,113), width = 0, 
              quat.rm = T, erotate = 270, tick.scale = 5, label.offset = 0.5)
HPDbars(ZygoBeastTree, label = "height_95%_HPD", col = "#8b948b", lwd = 5)
par(new = T)
rect(70.5, 15.5, 84.5, 41.5, col = LarrColT, border = NA); segments(84.5, 15.6, 84.5, 41.4, col = LarrCol, lwd = 2)
rect(70.5, 87.5, 84.5, 93.5, col = FagColT, border = NA); segments(84.5, 87.6, 84.5, 93.4, col = FagCol, lwd = 2)
rect(70.5, 85.5, 84.5, 86.5, col = FagColT, border = NA); segments(84.5, 85.6, 84.5, 86.4, col = FagCol, lwd = 2)
text(85.5, 89.25, expression("New World  "*italic(Fagonia)), srt = 90, cex = 0.5); text(85.5, 29, "Larreoideae", srt = 90, cex = 0.55)
par(new = T)
geoscalePhylo(tree = ladderize(ZygoBeastTree, right = F), units = c("Period", "Epoch"),
              cex.age = 0.7, cex.ts = 0.75, cex.tip = 0.4, x.lim = c(-15, 85), y.lim = c(4, 113), width = 1, 
              quat.rm = T, erotate = 270, tick.scale = 5, label.offset = 0.5)
## plot circles for calibrated nodes
symbols(-0.1, 13.125,  circles = 16, inches = 0.07, bg = "#286e9c", fg = "#114363", add = T)
symbols(5.8, 24.35,  circles = 16, inches = 0.07, bg = "#286e9c", fg = "#114363", add = T)
symbols(17.7, 43.4,  circles = 16, inches = 0.07, bg = "#286e9c", fg = "#114363", add = T)
symbols(33.7, 60.75,  circles = 16, inches = 0.07, bg = "#286e9c", fg = "#114363", add = T)
par(xpd = NA)
rect(-20, -25, margin, 119.6, -30, col = "white", border=NA); # to cut the time scale
rect(margin, -7, margin, 0.5, col = "white", border=NA)       # to cut the time scale
segments(margin, -18, margin, -0.35)                          # to put a line at the end of time scale
dev.off()

## ::::::::::::::::::::::::::::::::::::::::::::
## read Best tree with 1 secondary calibrations
ZygoBeastTree2 <- read.beast("Zygo_CP_Master_2019-11-15.tre")
ZygoBeastTree2$root.time <- ZygoBeastTree2$height[1]    # without this line you get error: "tree$root.time is missing, check tree is time scaled."
margin <- branching.times(ZygoBeastTree2)[1]-85.8      # set margin according to .... FEDERICO FRAGEN

## set colors for clade rect() function
FagColT <- rgb(0.11, 0.67, 0.57, alpha = 0.3); FagCol <- rgb(0.11, 0.67, 0.57)
LarrColT <- rgb(0.85, 0.71, 0.11, alpha = 0.3); LarrCol <- rgb(0.85, 0.71, 0.11)

## Plot the dated phylogeny (1 calibrations) with HPD Bars
pdf("Zygophyllales_Beast_Tree_Appendix_1cal.pdf", height = 11, width = 8)
geoscalePhylo(tree = ladderize(ZygoBeastTree2, right = F), boxes = "Epoch", units = c("Period", "Epoch"),
              cex.age = 0.7, cex.ts = 0.75, cex.tip = 0.4, x.lim = c(-15,85), y.lim = c(4,113), width = 0, 
              quat.rm = T, erotate = 270, tick.scale = 5, label.offset = 0.5)
HPDbars(ZygoBeastTree2, label = "height_95%_HPD", col = "#8b948b", lwd = 5)
par(new = T)
rect(64.5, 15.5, 78.5, 41.5, col = LarrColT, border = NA); segments(78.5, 15.6, 78.5, 41.4, col = LarrCol, lwd = 2)
rect(64.5, 47.5, 78.5, 53.5, col = FagColT, border = NA); segments(78.5, 47.6, 78.5, 53.4, col = FagCol, lwd = 2)
rect(64.5, 45.5, 78.5, 46.5, col = FagColT, border = NA); segments(78.5, 45.6, 78.5, 46.4, col = FagCol, lwd = 2)
text(79.5, 49, expression("New World  "*italic(Fagonia)), srt = 90, cex = 0.5); text(79.5, 29, "Larreoideae", srt = 90, cex = 0.55)
par(new = T)
geoscalePhylo(tree = ladderize(ZygoBeastTree2, right = F), units = c("Period", "Epoch"),
              cex.age = 0.7, cex.ts = 0.75, cex.tip = 0.4, x.lim = c(-15, 85), y.lim = c(4, 113), width = 1, 
              quat.rm = T, erotate = 270, tick.scale = 5, label.offset = 0.5)
## calibration points
symbols(-0.1, 13.8,  circles = 16, inches = 0.07, bg = "#286e9c", fg = "#114363", add = T)
par(xpd = NA)
rect(-25, -25, margin, 119.6, -30, col = "white", border=NA); # to cut the time scale
segments(margin, -18, margin, -0.35)                          # to put a line at the end of time scale
dev.off()

## ########### BIOGEOBEARS GLOBAL #######################
## ::::::::::::::::::  GLOBAL -- ZYGOPHYLLACEAE
## Load Tree and prun for analysis of Zygophyllaceae GLOBAL
z <- read.nexus("Zygo_CP_Master_2019-06-25.tre")          ## read Beast Tree with "read.nexus"
z$tip.label                                               ## get list of tiplabels
plot(ladderize(z, right = F), cex = 0.6)                  ## plot tree just to see the tree

## ladderize tree and extract Zygophyllaceae == exclude Krameria
z1 <- extract.clade(z, node = getMRCA(z, c("Seetzenia_lanata", "Zygophyllum_rosowii")))
z1 <- ladderize(z1, right = F)
plot(z1); nodelabels(cex = 0.7)

## read beast tree again but this time with "read.beast"
zB <- read.beast("Zygo_CP_Master_2019-06-25.tre")
zB1 <- ladderize(extract.clade(zB, node = getMRCA(z, c("Seetzenia_lanata", "Zygophyllum_rosowii"))), right = F)
zB0 <- extract.clade2(zB, node = getMRCA(z, c("Seetzenia_lanata", "Zygophyllum_rosowii")))

## make a .txt file only with the tip labels --> add distribution info outside R by hand
write.table(z1$tip.label, file = "BioGeoBEARS_TipTable_Zygo.txt")

## write a tree from extracted read.nexus() TREE as BioGeoBears needs an extra tree file
write.tree(z1, "BioGeoBears_Lad_Tree_Zygo.tre")

## :::::::::::::::::::::::::::::::::::::::::::::
## Analysis with 'maxareas' & 'max_range_size' = 2

## get the colors for the number of states and maxareas
states_list_0based_index = rcpp_areas_list_to_states_list(areas = c("A","B","C","D","E"), maxareas = 2)  ## adjust 'maxareas' and 'max_range_size' 
colors_matrix = get_colors_for_numareas(5)
colors_list_for_states = mix_colors_for_states(colors_matrix, states_list_0based_index)
colors_list_for_states

# read a table with color lists
coltipsZygo <- read.table("BioGeoBEARS_color_Zygo.txt", row.names = 1)
coltips <- coltipsZygo[z1$tip.label,]

## set parameter for BioGeoBears and run analysis
BioGeoBEARS_model_object = define_BioGeoBEARS_model_object(minval_anagenesis = 1e-15, minval_cladogenesis = 1e-03, maxval = 5)
BioGeoBEARS_run_object = define_BioGeoBEARS_run(BioGeoBEARS_model_object = BioGeoBEARS_model_object, trfn="BioGeoBears_Lad_Tree_Zygo.tre", 
                                                geogfn = "BioGeoBEARS_distr_Zygo.txt", num_cores_to_use = 8, max_range_size = 2)
BioGeoBEARS_run_object$geogfn = "BioGeoBEARS_distr_Zygo.txt"
zygoAArC = bears_optim_run(BioGeoBEARS_run_object)

# get the states with the highest likelihood at each node, use only 
# columns 2 to 16 because the first column is an order number
nodestates <- get_ML_states_from_relprobs(zygoAArC$ML_marginal_prob_each_state_at_branch_top_AT_node[,2:16],
                                          statenames = c("A","B","C","D","E","AB","AC","AD","AE","BC","BD","BE","CD","CE","DE"))

# read png map to plot later in the figure 
world <- readPNG("BGB_Zygo_Globalmap_5.png")
Africa <- readPNG("PIE_Africa.png")
SAm <- readPNG("PIE_SAm.png")

zB1$root.time <- zB0$height[1]
z1$root.time <- zB0$height[1]
zB1$"height_95%_HPD_MIN" <- zB0$"height_95%_HPD_MIN"
zB1$"height_95%_HPD_MAX" <- zB0$"height_95%_HPD_MAX"

weis <- rgb(1,1,1, alpha = 0.6)

## Plot BioGeoBears Results with maxareas=2
pdf("BGBears_BeastTree_Zygophyllaceae_2Areas.pdf", height = 11, width = 8)
geoscalePhylo(tree = z1, boxes = "Epoch", units = c("Period", "Epoch"),
              cex.age = 0.7, cex.ts = 0.75, cex.tip = 0.4, x.lim = c(-15,79), y.lim = c(4,112.25), width = 0, 
              quat.rm = T, erotate = 270, arotate = 90, tick.scale = 5, label.offset = 6)
HPDbars(zB1, label = "height_95%_HPD", col = "#8b948b", lwd = 5)
par(new = T)
geoscalePhylo(tree = z1, units = c("Period", "Epoch"),
              cex.age = 0.7, cex.ts = 0.75, cex.tip = 0.4, x.lim = c(-15,79), y.lim = c(4,112.25), width = 1, 
              quat.rm = T, erotate = 270, arotate = 90, tick.scale = 5, label.offset = 6)
par(xpd = NA)
rect(-21, -25, -16.5, 119, col = "white", border = NA); rect(-21, 86, 35, 119, col = "white", border = NA); rect(-21, 60, 7, 90, col = "white", border = NA)
segments(-16.5, -18, -16.5, -0.35)                      # to put a line at the end of time scale
tiplabels(pch = 22, bg = as.character(coltips$V2), cex = 1.0, adj = c(1.5, 0.51), lwd = 0.4)
tiplabels(pch = 22, bg = as.character(coltips$V3), cex = 1.0, adj = c(2.5, 0.51), lwd = 0.4)
tiplabels(pch = 22, bg = as.character(coltips$V4), cex = 1.0, adj = c(3.5, 0.51), lwd = 0.4)
tiplabels(pch = 22, bg = as.character(coltips$V5), cex = 1.0, adj = c(4.5, 0.51), lwd = 0.4)
tiplabels(pch = 22, bg = as.character(coltips$V6), cex = 1.0, adj = c(5.5, 0.51), lwd = 0.4)
mtext(side = 3, line = -0.8, at = c(65.34, 66.34, 67.34, 68.34, 69.34), text = c("A", "B", "C", "D", "E"), cex = 0.575, font = 2)
mtext(side = 1, line = 4.8, at = 63, text = "Quat.", cex = 0.5)
legend(-15, 84, legend = c("A: North and Central America","B: South America","C: Africa and the Mediterranean","D: Asia","E: Australia",
                           "A+B","A+C","A+D","A+E","B+C","B+D","B+E","C+D","C+E","D+E")
       [c(1:7,10,13,14)], pch = 22, pt.bg = colors_list_for_states[c(1:7,10,13,14)], bty = "n", cex = 0.75, pt.cex = 1.5)
plot_BioGeoBEARS_results(zygoAArC, plotwhat = "pie", tipcex = 0, statecex = 0.4, plotsplits = F, plotlegend = F,
                         legend_cex = 0.5, titlecex = 0, 
                         cornercoords_loc = "/home/tim/R/x86_64-pc-linux-gnu-library/3.4/BioGeoBEARS/extdata/a_scripts/", 
                         skiptree = T, tipboxes_TF = F)
nodelabels(nodestates[115:227], adj = c(1.75, -0.4), frame = "none", cex = 0.65)
rasterImage(world, -18, 86, 40, 117); rasterImage(Africa, 11, 57.6, 12.8, 59.5); rasterImage(SAm, 11, 23.15, 12.8, 25.05)
points(-6.5, 106, pch = 16, col = weis, cex = 2); text(-6.5, 106, "A", font = 2, cex = 0.85); text(0.5, 96.5, "B", font = 2, cex = 0.85)
text(13, 102, "C", font = 2, cex = 0.85); text(24.5, 105.5, "D", font = 2, cex = 0.85); text(31.9, 94.5, "E", font = 2, cex = 0.85)
text(10.25, 59.25, "C", cex = 0.7); text(10.4, 25.25, "B", cex = 0.7)
dev.off()

## :::::::::::::::::::::::::::::::::::::::::::::
## Analysis with 'maxareas' & 'max_range_size' = 3

## get the colors for the number of states and maxareas
states_list_0based_index = rcpp_areas_list_to_states_list(areas = c("A","B","C","D","E"), maxareas = 3)  ## adjust 'maxareas' and 'max_range_size' 
colors_matrix = get_colors_for_numareas(5)
colors_list_for_states = mix_colors_for_states(colors_matrix, states_list_0based_index)
colors_list_for_states

# read a table with color lists
coltipsZygo <- read.table("BioGeoBEARS_color_Zygo.txt", row.names = 1)
coltips <- coltipsZygo[z1$tip.label,]

## set parameter for BioGeoBears and run analysis
BioGeoBEARS_model_object = define_BioGeoBEARS_model_object(minval_anagenesis = 1e-15, minval_cladogenesis = 1e-03, maxval = 5)
BioGeoBEARS_run_object = define_BioGeoBEARS_run(BioGeoBEARS_model_object = BioGeoBEARS_model_object, trfn="BioGeoBears_Lad_Tree_Zygo.tre", 
                                                geogfn = "BioGeoBEARS_distr_Zygo.txt", num_cores_to_use = 8, max_range_size = 3)
BioGeoBEARS_run_object$geogfn = "BioGeoBEARS_distr_Zygo.txt"
zygoAArC = bears_optim_run(BioGeoBEARS_run_object)

## get the states with the highest likelihood at each node, use only 
## columns 2 to 16 because the first column is an order number
nodestates3 <- get_ML_states_from_relprobs(zygoAArC$ML_marginal_prob_each_state_at_branch_top_AT_node[,2:26],
                                           statenames = c("A","B","C","D","E","AB","AC","AD","AE","BC","BD","BE","CD","CE","DE", 
                                                          "ABC", "ABD", "ABE", "ACD", "ACE", "ADE", "BCD", "BCE", "BDE", "CDE"))

## read png map to plot later in the figure 
world <- readPNG("BGB_Zygo_Globalmap_5.png")
Africa <- readPNG("PIE_Africa.png")
SAm <- readPNG("PIE_SAm.png")

zB1$root.time <- zB0$height[1]
z1$root.time <- zB0$height[1]
zB1$"height_95%_HPD_MIN" <- zB0$"height_95%_HPD_MIN"
zB1$"height_95%_HPD_MAX" <- zB0$"height_95%_HPD_MAX"

weis <- rgb(1,1,1, alpha = 0.6)

pdf("BGBears_BeastTree_Zygophyllaceae_3Areas.pdf", height = 11, width = 8)
geoscalePhylo(tree = z1, boxes = "Epoch", units = c("Period", "Epoch"),
              cex.age = 0.7, cex.ts = 0.75, cex.tip = 0.4, x.lim = c(-15,79), y.lim = c(4,112.25), width = 0, 
              quat.rm = T, erotate = 270, arotate = 90, tick.scale = 5, label.offset = 6)
HPDbars(zB1, label = "height_95%_HPD", col = "#8b948b", lwd = 5)
par(new = T)
geoscalePhylo(tree = z1, units = c("Period", "Epoch"),
              cex.age = 0.7, cex.ts = 0.75, cex.tip = 0.4, x.lim = c(-15,79), y.lim = c(4,112.25), width = 1, 
              quat.rm = T, erotate = 270, arotate = 90, tick.scale = 5, label.offset = 6)
par(xpd = NA)
rect(-21, -25, -16.5, 119, col = "white", border = NA); rect(-21, 86, 35, 119, col = "white", border = NA); rect(-21, 60, 7, 90, col = "white", border = NA)
segments(-16.5, -18, -16.5, -0.35)                      # to put a line at the end of time scale
tiplabels(pch = 22, bg = as.character(coltips$V2), cex = 1.0, adj = c(1.5, 0.51), lwd = 0.4)
tiplabels(pch = 22, bg = as.character(coltips$V3), cex = 1.0, adj = c(2.5, 0.51), lwd = 0.4)
tiplabels(pch = 22, bg = as.character(coltips$V4), cex = 1.0, adj = c(3.5, 0.51), lwd = 0.4)
tiplabels(pch = 22, bg = as.character(coltips$V5), cex = 1.0, adj = c(4.5, 0.51), lwd = 0.4)
tiplabels(pch = 22, bg = as.character(coltips$V6), cex = 1.0, adj = c(5.5, 0.51), lwd = 0.4)
mtext(side = 3, line = -0.8, at = c(65.34, 66.34, 67.34, 68.34, 69.34), text = c("A", "B", "C", "D", "E"), cex = 0.575, font = 2)
mtext(side = 1, line = 4.8, at = 63, text = "Quat.", cex = 0.5)
legend(-15, 84, legend = c("A: North and Central America","B: South America","C: Africa and the Mediterranean","D: Asia","E: Australia",
                     "A+B","A+C","A+D","A+E","B+C","B+D","B+E","C+D","C+E","D+E", "ABC", "ABD", "ABE", "ACD", "ACE", "ADE", "BCD", "B+C+E", "BDE", "CDE")
       [c(1:7,10,13,14,23)], pch = 22, pt.bg = colors_list_for_states[c(1:7,10,13,14,23)], bty = "n", cex = 0.75, pt.cex = 1.5)
plot_BioGeoBEARS_results(zygoAArC, plotwhat = "pie", tipcex = 0, statecex = 0.4, plotsplits = T, plotlegend = FALSE,
                         legend_cex = 0.5, titlecex = 0, splitcex = 0.3,
                         cornercoords_loc = "/home/tim/R/x86_64-pc-linux-gnu-library/3.4/BioGeoBEARS/extdata/a_scripts/", 
                         skiptree = T, tipboxes_TF = F)
nodelabels(nodestates3[115:227], adj = c(1.75, -0.4), frame = "none", cex = 0.65)
rasterImage(world, -18, 86, 40, 117)#; rasterImage(Africa, 11, 57.6, 12.8, 59.5); rasterImage(SAm, 11, 23.15, 12.8, 25.05)
points(-6.5, 106, pch = 16, col = weis, cex = 2); text(-6.5, 106, "A", font = 2, cex = 0.85); text(0.5, 96.5, "B", font = 2, cex = 0.85)
text(13, 102, "C", font = 2, cex = 0.85); text(24.5, 105.5, "D", font = 2, cex = 0.85); text(31.9, 94.5, "E", font = 2, cex = 0.85)
text(10.25, 59.25, "C", cex = 0.7); text(10.4, 25.25, "B", cex = 0.7)
dev.off()


## ########### BIOGEOBEARS AMERICA #########
## ::::::::::::::::::  SOUTH AMERICA -- LARREOIDEAE
## Load Tree and prun for analysis of LARREOIDEAE SOUTH AMERICA
l <- read.nexus("Zygo_CP_Master_2019-06-25.tre")          ## read Beast Tree with "read.nexus"
l$tip.label                                               ## get list of tiplabels
plot(ladderize(l, right = F), cex = 0.6)                  ## plot tree just to see outgroup ...

## ladderize tree and extract only Larreoideae
l1 <- extract.clade(l, node = getMRCA(l, c("Metharme_lanata", "Guaiacum_coulteri")))
l1 <- ladderize(l1, right = F)

## read beast tree again but this time with "read.beast"
lB <- read.beast("Zygo_CP_Master_2019-06-25.tre")
lB1 <- ladderize(extract.clade(lB, node = getMRCA(l, c("Metharme_lanata", "Guaiacum_coulteri"))), right = F)
lB0 <- extract.clade2(lB, node = getMRCA(l, c("Metharme_lanata", "Guaiacum_coulteri")))

## make a txt. file only with the tip labels --> add distribution info
write.table(l1$tip.label, file = "BioGeoBEARS_TipTable_Larr.txt")

## write a tree from extracted read.nexus() TREE
write.tree(l1, "BioGeoBears_Lad_Tree_Larr.tre")


## :::::::::::::::::::::::::::::::::::::::::::::
## Analysis with 'maxareas' & 'max_range_size' = 2

## get the colors for the number of states and maxareas
states_list_0based_index = rcpp_areas_list_to_states_list(areas=c("A","B","C","D","E"), maxareas = 2)
colors_matrix = get_colors_for_numareas(5)
colors_list_for_states = mix_colors_for_states(colors_matrix, states_list_0based_index)
colors_list_for_states    # Use the first five color codes the color "BioGeoBEARS_color_Larr.txt" file

## read a table with color lists
coltipsLarr <- read.table("BioGeoBEARS_color_Larr.txt", row.names = 1)
coltips <- coltipsLarr[l1$tip.label,]

## set parameter for BioGeoBears and run analysis
BioGeoBEARS_model_object = define_BioGeoBEARS_model_object(minval_anagenesis = 1e-15, minval_cladogenesis = 1e-03, maxval = 5)
BioGeoBEARS_run_object = define_BioGeoBEARS_run(BioGeoBEARS_model_object = BioGeoBEARS_model_object, trfn="BioGeoBears_Lad_Tree_Larr.tre", 
                                                geogfn = "BioGeoBEARS_distr_Larr.txt", num_cores_to_use = 8, max_range_size = 2)
BioGeoBEARS_run_object$geogfn = "BioGeoBEARS_distr_Larr.txt"
larrAArC = bears_optim_run(BioGeoBEARS_run_object)

## get the states with the highest likelihood at each node, use only 
## columns 2 to 16 because the first column is an order number
nodestatesLarr <- get_ML_states_from_relprobs(larrAArC$ML_marginal_prob_each_state_at_branch_top_AT_node[,2:16],
                                          statenames = c("A","B","C","D","E","AB","AC","AD","AE","BC","BD","BE","CD","CE","DE"))

## read png map to plot later in the figure 
Americas <- readPNG("BGB_Larreoideae_Americas_5_21.png")

lB1$root.time <- lB0$height[1]
l1$root.time <- lB0$height[1]
lB1$"height_95%_HPD_MIN" <- lB0$"height_95%_HPD_MIN"
lB1$"height_95%_HPD_MAX" <- lB0$"height_95%_HPD_MAX"

weis <- rgb(1,1,1, alpha = 0.6)

pdf("BGBears_BeastTree_Larreoideae_2Areas.pdf", height = 7, width = 10)
geoscalePhylo(tree = l1, boxes = "Epoch", units = c("Period", "Epoch"),
              cex.age = 0.7, cex.ts = 0.75, cex.tip = 0.6, x.lim = c(-15,47), y.lim = c(1.25,26.5), width = 0, 
              quat.rm = T, erotate = 270, tick.scale = 5, label.offset = 6)
rect(-20, 9, -1, 28, col = "white", border = NA)  # make the area behind the map white!
rasterImage(Americas, -23, -1, 6, 28)
HPDbars(lB1, label = "height_95%_HPD", col = "#8b948b", lwd = 6)
par(new = T)
geoscalePhylo(tree = l1, units = c("Period", "Epoch"),
              cex.age = 0.7, cex.ts = 0.75, cex.tip = 0.6, x.lim = c(-15,47), y.lim = c(1.25,26.5), width = 1.55, 
              quat.rm = T, erotate = 270, tick.scale = 5, label.offset = 6)
par(xpd = NA)
rect(-22, -8, -17, 10, col = "white", border = NA)
segments(-17, -3.85, -17, 0.245)                      # to put a line at the end of time scale
tiplabels(pch = 22, bg = as.character(coltips$V2), cex = 1.85, adj = c(1.5, 0.51), lwd = 1)
tiplabels(pch = 22, bg = as.character(coltips$V3), cex = 1.85, adj = c(2.5, 0.51), lwd = 1)
tiplabels(pch = 22, bg = as.character(coltips$V4), cex = 1.85, adj = c(3.5, 0.51), lwd = 1)
tiplabels(pch = 22, bg = as.character(coltips$V5), cex = 1.85, adj = c(4.5, 0.51), lwd = 1)
tiplabels(pch = 22, bg = as.character(coltips$V6), cex = 1.85, adj = c(5.5, 0.51), lwd = 1)
mtext(side = 3, line = -1, at = c(29.65, 30.675, 31.7, 32.75, 33.8), text = c("A", "B", "C", "D", "E"), cex = 0.8, font = 2)
mtext(side = 1, line = 2.85, at = 27.5, text = "Quat.", cex = 0.725)
legend(-6, 27.75, legend = c("A: Northern America","B: Central America","C: southeast S. America (SE-SA)","D: Peruvian Andes & coastal Desert","E: Atacama Desert",
                             "A+B","A+C","A+D","A+E","B+C","B+D","B+E","C+D","C+E","D+E")
       [c(1:7,10,13,14)], pch = 22, pt.bg = colors_list_for_states[c(1:7,10,13,14)], bty = "n", cex = 0.75, pt.cex = 1.5)
plot_BioGeoBEARS_results(larrAArC, plotwhat = "pie", tipcex = 0, statecex = 0.4, plotsplits = T, splitcex = 0.4, plotlegend = FALSE,
                         legend_cex = 0.5, titlecex = 0, 
                         cornercoords_loc = "/home/tim/R/x86_64-pc-linux-gnu-library/3.4/BioGeoBEARS/extdata/a_scripts/", 
                         skiptree = T, tipboxes_TF = F)
nodelabels(nodestatesLarr[27:51], adj = c(1.75, -0.7), frame = "none", cex = 0.7)
points(-16, 21.75, pch = 16, col = weis, cex = 3.5); text(-16, 21.75, "A", font = 2, cex = 1.25); text(-8.15, 18.35, "B", font = 2, cex = 1.25); 
text(-3.25, 6.5, "C", font = 2, cex = 1.25); text(-9, 11.65, "D", font = 2, cex = 1.25); text(-6.57, 8, "E", font = 2, cex = 1.25)
dev.off()


## :::::::::::::::::::::::::::::::::::::::::::::
## Analysis with 'maxareas' & 'max_range_size' = 3

## get the colors for the number of states and maxareas
states_list_0based_index = rcpp_areas_list_to_states_list(areas=c("A","B","C","D","E"), maxareas = 3)
colors_matrix = get_colors_for_numareas(5)
colors_list_for_states = mix_colors_for_states(colors_matrix, states_list_0based_index)
colors_list_for_states    # Use the first five color codes the color "BioGeoBEARS_color_Larr.txt" file

## read a table with color lists
coltipsLarr <- read.table("BioGeoBEARS_color_Larr.txt", row.names = 1)
coltips <- coltipsLarr[l1$tip.label,]

## set parameter for BioGeoBears and run analysis
BioGeoBEARS_model_object = define_BioGeoBEARS_model_object(minval_anagenesis = 1e-15, minval_cladogenesis = 1e-03, maxval = 5)
BioGeoBEARS_run_object = define_BioGeoBEARS_run(BioGeoBEARS_model_object = BioGeoBEARS_model_object, trfn="BioGeoBears_Lad_Tree_Larr.tre", 
                                                geogfn = "BioGeoBEARS_distr_Larr.txt", num_cores_to_use = 8, max_range_size = 3)
BioGeoBEARS_run_object$geogfn = "BioGeoBEARS_distr_Larr.txt"
larrAArC = bears_optim_run(BioGeoBEARS_run_object)

## get the states with the highest likelihood at each node, use only 
## columns 2 to 16 because the first column is an order number
nodestatesLarr3 <- get_ML_states_from_relprobs(larrAArC$ML_marginal_prob_each_state_at_branch_top_AT_node[,2:26],
                                               statenames = c("A","B","C","D","E","AB","AC","AD","AE","BC","BD","BE","CD","CE","DE",
                                                              "ABC", "ABD", "ABE", "ACD", "ACE", "ADE", "BCD", "BCE", "BDE", "CDE"))

## read png map to plot later in the figure 
Americas <- readPNG("BGB_Larreoideae_Americas_5_21.png")

lB1$root.time <- lB0$height[1]
l1$root.time <- lB0$height[1]
lB1$"height_95%_HPD_MIN" <- lB0$"height_95%_HPD_MIN"
lB1$"height_95%_HPD_MAX" <- lB0$"height_95%_HPD_MAX"

weis <- rgb(1,1,1, alpha = 0.6)

pdf("BGBears_BeastTree_Larreoideae_3Areas.pdf", height = 7, width = 10)
geoscalePhylo(tree = l1, boxes = "Epoch", units = c("Period", "Epoch"),
              cex.age = 0.7, cex.ts = 0.75, cex.tip = 0.6, x.lim = c(-15,47), y.lim = c(1.25,26.5), width = 0, 
              quat.rm = T, erotate = 270, tick.scale = 5, label.offset = 6)
rect(-20, 9, -1, 28, col = "white", border = NA)  # make the area behind the map white!
rasterImage(Americas, -23, -1, 6, 28)
HPDbars(lB1, label = "height_95%_HPD", col = "#8b948b", lwd = 6)
par(new = T)
geoscalePhylo(tree = l1, units = c("Period", "Epoch"),
              cex.age = 0.7, cex.ts = 0.75, cex.tip = 0.6, x.lim = c(-15,47), y.lim = c(1.25,26.5), width = 1.55, 
              quat.rm = T, erotate = 270, tick.scale = 5, label.offset = 6)
par(xpd = NA)
rect(-22, -8, -17, 10, col = "white", border = NA)
segments(-17, -3.85, -17, 0.245)                      # to put a line at the end of time scale
tiplabels(pch = 22, bg = as.character(coltips$V2), cex = 1.85, adj = c(1.5, 0.51), lwd = 1)
tiplabels(pch = 22, bg = as.character(coltips$V3), cex = 1.85, adj = c(2.5, 0.51), lwd = 1)
tiplabels(pch = 22, bg = as.character(coltips$V4), cex = 1.85, adj = c(3.5, 0.51), lwd = 1)
tiplabels(pch = 22, bg = as.character(coltips$V5), cex = 1.85, adj = c(4.5, 0.51), lwd = 1)
tiplabels(pch = 22, bg = as.character(coltips$V6), cex = 1.85, adj = c(5.5, 0.51), lwd = 1)
mtext(side = 3, line = -1, at = c(29.65, 30.675, 31.7, 32.75, 33.8), text = c("A", "B", "C", "D", "E"), cex = 0.8, font = 2)
mtext(side = 1, line = 2.85, at = 27.5, text = "Quat.", cex = 0.725)
legend(-6, 27.75, legend = c("A: Northern America","B: Central America","C: southeast S. America (SE-SA)","D: Peruvian Andes & coastal Desert","E: Atacama Desert",
                             "A+B","A+C","A+D","A+E","B+C","B+D","B+E","C+D","C+E","D+E", "ABC", "ABD", "ABE", "A+C+D", "ACE", "ADE", "BCD", "B+C+E", "BDE", "C+D+E")
       [c(1:5,6,10,13,14,19,23,25)], pch = 22, pt.bg = colors_list_for_states[c(1:5,6,10,13,14,19,23,25)], bty = "n", cex = 0.75, pt.cex = 1.5)
plot_BioGeoBEARS_results(larrAArC, plotwhat = "pie", tipcex = 0, statecex = 0.4, plotsplits = T, splitcex = 0.3, plotlegend = FALSE,
                         legend_cex = 0.5, titlecex = 0, 
                         cornercoords_loc = "/home/tim/R/x86_64-pc-linux-gnu-library/3.4/BioGeoBEARS/extdata/a_scripts/", 
                         skiptree = T, tipboxes_TF = F)
nodelabels(nodestatesLarr3[27:51], adj = c(1.75, -0.7), frame = "none", cex = 0.7)
points(-16, 21.75, pch = 16, col = weis, cex = 3.5); text(-16, 21.75, "A", font = 2, cex = 1.25); text(-8.15, 18.35, "B", font = 2, cex = 1.25); 
text(-3.25, 6.5, "C", font = 2, cex = 1.25); text(-9, 11.65, "D", font = 2, cex = 1.25); text(-6.57, 8, "E", font = 2, cex = 1.25)
dev.off()
