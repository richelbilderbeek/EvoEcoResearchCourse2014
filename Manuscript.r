# Manuscript, script to generate figures and tables for the Evolution & Ecology Research course 2014
# Copyright (C) 2014 Richel Bilderbeek
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.See the
# GNU General Public License for more details.
# You should have received a copy of the GNU General Public License
# along with this program.If not, see <http://www.gnu.org/licenses/>.
#
# File can be downloaded from https://github.com/richelbilderbeek/EvoEcoResearchCourse2014
#
# Research question: Can I predict where species are from a redox potential?
# Experiment: 
# * Couple shallow (0-5 cm) benthos abundances with redox measurement of 2 cm
# * Couple deep (5-20 cm) benthos abundances with redox measurement of 10 cm
# Answer: 
# * No 

rm(list = ls())
setwd("~/GitHubs/EvoEcoResearchCourse2014")
library(reshape2)
library(testit)

# Create benthos data as such:
#  dist_m depth_cm species_name
#  2350   2        Bathyporeia_spec.
#  2400   10       Nereis_diversicolor
# Removes species at unknown depths
CreateDataBenthos <- function()
{
	data_benthos <- read.table("benthos_species_diversity.csv",header=TRUE,sep=",")
	#Remove useless columns
	data_benthos <-data_benthos[ c("dist_m","soil_layer","species_name") ]
	# Only keep rows with 'D' or 'S' (for deep and shallow/top part of the core)
	data_benthos<-data_benthos[data_benthos$soil_layer != "U" & data_benthos$soil_layer != "?", ]
	# Drop the unneeded levels 
	data_benthos$soil_layer <- droplevels(data_benthos$soil_layer)
	#Of data_benthod, change the column 'soil_layer' to 'depth_cm': 
	# if soil_layer == "S" -> depth_cm = 2
	# if soil_layer == "D" -> depth_cm = 10
	names(data_benthos) <- c("dist_m","depth_cm","species_name")
	data_benthos$depth_cm <- as.factor(data_benthos$depth_cm)
	# Rename the levels
	levels(data_benthos$depth_cm) <- c(2,10)
	data_benthos
}

# Count the species occurring at two depths
# species_name      2  10
# Arenicola_marina  12  1
# Bathyporeia_spec. 2   0
GetSpeciesCountAtDepths <- function()
{
	species_to_depth <- dcast(
		CreateDataBenthos(),
		species_name ~ depth_cm,
		value.var = "species_name",
		fill = 0, 
		fun.aggregate=length # fun.aggregate may also be mean
	) 
	species_to_depth
}

# Obtain the species names occuring at least thrice in both of the two depths in the soil
# Will be this:
# [1] Carcinus_maenas     Cerastoderma_edule  Crassostrea_gigas   Hemigropsus_takanoi Hydrobia_ulvae     
# [6] Littorina_littorea  Mytilus_edulis      Nereis_diversicolor
# 8 Levels: Carcinus_maenas Cerastoderma_edule Crassostrea_gigas Hemigropsus_takanoi ... Nereis_diversicolor
GetSelectedSpecies <- function()
{
	species_to_depth <- GetSpeciesCountAtDepths()
	selected_species_to_depth <- subset(species_to_depth,
		  species_to_depth$"2" > 2 & species_to_depth$"10" > 2)
	selected_species_to_depth <- droplevels(selected_species_to_depth)
	# Collect all species' names
	species_list <- levels(selected_species_to_depth$species_name)
	selected_species <- selected_species_to_depth$species_name
	selected_species <- droplevels(selected_species)
	selected_species
}

GetDataBenthosSelected <- function()
{
  data_benthos <- CreateDataBenthos()
	selected_species <- GetSelectedSpecies()
  #data_benthos_selected <- data_benthos[ data_benthos$species_name %in% selected_species,]
	data_benthos_selected <- data_benthos[ data_benthos$species_name %in% GetSelectedSpecies(),]
  data_benthos_selected
}

# Prepare the redox data
# dist_m depth_cm replicate redox_calib
# 220    2        1         304.89200
# 220    2        2         302.63036
CreateDataRedox <- function()
{	
	data_redox <- read.table("Redox.csv",header=TRUE,sep="\t")
	data_redox <- subset(data_redox,depth_cm != 5)
	data_redox
	#Remove the replicate column
	data_redox <- data_redox[c("dist_m","depth_cm","redox_calib") ]
	# Take the average redox values, so that every distance and depth has a single redox value
	data_redox <- aggregate(data_redox,list(data_redox$depth_cm,data_redox$dist_m),mean)
	data_redox <- data_redox[c("dist_m","depth_cm","redox_calib") ]
	data_redox
}

CreateDataRedoxAll <- function()
{	
	data_redox <- read.table("Redox.csv",header=TRUE,sep="\t")
	#Remove the replicate column
	data_redox <- data_redox[c("dist_m","depth_cm","redox_calib") ]
	# Take the average redox values, so that every distance and depth has a single redox value
	data_redox <- aggregate(data_redox,list(data_redox$depth_cm,data_redox$dist_m),mean)
	data_redox <- data_redox[c("dist_m","depth_cm","redox_calib") ]
	data_redox
}

GetRedoxPerDistance <- function()
{	
	data_redox <- read.table("Redox.csv",header=TRUE,sep="\t")
	#Remove the replicate column
	data_redox <- data_redox[c("dist_m","depth_cm","redox_calib") ]
	# Take the average redox values, so that every distance and depth has a single redox value
	data_redox <- aggregate(data_redox,list(data_redox$depth_cm,data_redox$dist_m),mean)
	data_redox <- data_redox[c("dist_m","depth_cm","redox_calib") ]
	data_redox
}

# Get the 32 redox values used in this study
# [1] -987.68238 -857.44961 -834.07933 -811.22734 -750.16306 -742.48291 -699.70022 -616.30224 -577.85436 -528.05116
# ...
# [31]  350.83156  370.14974
GetRedoxValues <- function()
{
	redox_values <- subset(CreateDataRedox(),select="redox_calib")
	redox_values <- unique(redox_values)
	redox_values <- sort(redox_values$redox_calib)
	redox_values
}

# Obtain the distances of all redox sites
#  [1]  220  250  300  350  400  500 1000 1150 1400 1600 1800 1970 2250 2300 2350 2400
GetDistances <- function()
{
	distances <- subset(CreateDataRedox(),select="dist_m")
	distances <- unique(distances)
	distances <- sort(distances$dist_m)
	distances
}

# Get all species per redox potential
# redox_calib   species_name
# -78.26751 Hydrobia_ulvae
# -78.26751 Hydrobia_ulvae
GetDataCombined <- function()
{
	data_combined <- merge(CreateDataRedox(),GetDataBenthosSelected(),by=c("dist_m","depth_cm"),all=FALSE)
	data_combined$species_name <- droplevels(data_combined$species_name)
	data_combined <- data_combined[ c("redox_calib","species_name") ]
	data_combined
}

# Get the redox potentials of Hydrobia
# [1] -78.26751 -78.26751 -78.26751 -78.26751 -78.26751 -78.26751 -78.26751 -78.26751 -78.26751 -78.26751
GetRedoxesHydrobia <- function()
{
  redoxes_hydrobia <- subset(GetDataCombined(),species_name == "Hydrobia_ulvae")
  redoxes_hydrobia <- subset(redoxes_hydrobia,select=c("redox_calib"))
  redoxes_hydrobia
}

GetRedoxesNereis <- function()
{
  redoxes_nereis <- subset(GetDataCombined(),species_name == "Nereis_diversicolor")
  redoxes_nereis <- subset(redoxes_nereis,select=c("redox_calib"))
	redoxes_nereis
}

TallySpeciesPerRedox <- function()
{
	redox_to_species <- dcast(
		GetDataCombined(),
		redox_calib ~ species_name,
		value.var = "species_name",
		fill = 0, 
		fun.aggregate=length # fun.aggregate may also be mean
	) 
	redox_calib <- redox_to_species$redox_calib
	# Note that not all species are present anymore that often. This is due to that not all distances are redoxed
	redox_to_species
}

TallySelectedSpeciesPerRedox <- function()
{
	redox_to_species <- TallySpeciesPerRedox()
	redox_to_species <- redox_to_species[, colSums(redox_to_species) > 6]
	redox_to_species <- cbind(TallySpeciesPerRedox(),redox_to_species)
	redox_to_species <- redox_to_species[, c("redox_calib","Hydrobia_ulvae","Nereis_diversicolor")]
	redox_to_species
}

# Can the redox potentials at a certain distane be assumed to be linear at a certain depth?
CalcOrderednessPerDistance <- function()
{
	orderedness_per_distance <- data.frame(dist_m = numeric(),order = factor())
	for (i in GetDistances())
	{
    data <- subset(CreateDataRedoxAll(),dist_m == i)
    low <- subset(data,depth_cm == 2)$redox_calib
    mid <- subset(data,depth_cm == 5)$redox_calib
    high <- subset(data,depth_cm == 10)$redox_calib
	  if (low < mid && mid < high) 
	  { 
	  	orderedness_per_distance <- rbind(orderedness_per_distance,data.frame(dist_m = i,order = "LTH"))
	  } else if(low > mid && mid > high) 
	  { 
	  	orderedness_per_distance <- rbind(orderedness_per_distance,data.frame(dist_m = i, order = "HTL"))
	  } else  { 
	  	orderedness_per_distance <- rbind(orderedness_per_distance,data.frame(dist_m = i, order = "UNO"))
	  }
	}
	orderedness_per_distance
}

# Generate figure for species abundances for the range of redox potentials
#   in two vertically aligned plots
CreateFigureSpeciesAbundancesSeperate <- function()
{
	svg(filename="Figure_species_abundances_seperate.svg")
	par(mfrow=c(2,1))
	plot(
		Hydrobia_ulvae ~ redox_calib, 
		data = TallySpeciesPerRedox(), 
		t = "b", 
		pch = 19, 
		col = "black",
		#main = "Hydrobia ulvae abundance",
		xlab = "Redox potential (mV)",
		ylab = "# Hydrobia ulvae"
	)
	plot(
		Nereis_diversicolor ~ redox_calib, 
		data = TallySpeciesPerRedox(), 
		t = "b", 
		pch = 19, 
		col = "black",
		#main = "Nereis diversicolor abundance",
		xlab = "Redox potential (mV)",
		ylab = "# Nereis Diversicolor",
		ylim = c(0,4)
	)
	par(mfrow=c(1,1))
	dev.off()
}

# Generate figure for species abundances for the range of redox potentials
#   in the same plot
CreateFigureSpeciesAbundances <- function()
{
	svg(filename="Figure_species_abundances.svg")
	par(mar = c(5, 4, 4, 4) + 0.3)  # Leave space for z axis
	par(new = FALSE) 
	plot(Hydrobia_ulvae ~ redox_calib, data = TallySpeciesPerRedox(), pch=19, axes=FALSE, xlab="Redox potential (mV)", ylab="", 
	   type="b",col="blue"
		#main="Species abundandances\nfor different redox potentials"
	)
	axis(1, col="black",las=1) #'las=1' align labels horizontally
	axis(2, col="blue",las=1) #'las=1' align labels horizontally
	mtext("# Hydrobia",side=2,line=2.5,col="blue")
	box()
	par(new = TRUE) # Prevents R from clearing the area 
	plot(
		Nereis_diversicolor ~ redox_calib, axes = FALSE, bty = "n", xlab = "", ylab = "", data = TallySpeciesPerRedox(),
		type="b",
		col = "red",
		yaxt="n", 
		pch = 19
	)
	axis(
		side=4, 
		at=seq(0,3,1), #Otherwise, 0.5 would be shown as a tick mark
		col="red",
		las=1 #Align labels horizontally
	)
	mtext("# Nereis", side=4, line=3,col="red")
	legend("topright",
	  inset=0.05,title = "Species name", c("Hydrobia ulvae", "Nereis diversicolor"),horiz=TRUE,
	  fill=c("blue","red"), 
		cex = 0.75
	)
	dev.off()
}

# Redox potentials along the transect for the two depths
CreateFigureRedoxPerDistance <- function()
{
	svg(filename="Figure_redox_per_distance.svg")
	dist_to_redox <- subset(CreateDataRedox(), CreateDataRedox()$dist_m %in% GetDistances())
	y_min <- min(dist_to_redox$redox_calib) - 100
	y_max <- max(dist_to_redox$redox_calib) + 1000
	y_text <- y_max - 200
	plot(
		dist_to_redox$redox_calib ~ dist_to_redox$dist_m,
		col = as.factor(dist_to_redox$depth_cm), 
		pch = 19,
		#main="Redox potentials along the transect\nfor the two depths",
		xlab="Distance along transect (m)",
		ylab="Redox potential (mV)",
		ylim=c(y_min,y_max)
	)
	segments(1125,y_min,1125,y_max,col="black",lty=3)
	segments(1900,y_min,1900,y_max,col="black",lty=3)
	segments(2225,y_min,2225,y_max,col="black",lty=3)
	text((1125 + 150) / 2,y_text,"Salt\nmarsh")
	text((1900 + 1125) / 2,y_text,"Mud\nflat")
	text((2225 + 1900) / 2,y_text,"Mussel\nbed")
	text(-50 + ((2600 + 2225) / 2),y_text,"Mud\nflat")
	legend("bottomleft",
	  inset=0.05,title = "Sampling depth", c("2 cm", "10 cm"),horiz=TRUE,
	  fill=c(1,2,3), cex = 1.0
	)
	rm(dist_to_redox)
	rm(y_min)
	rm(y_max)
	rm(y_text)
	dev.off()
}

assert("CreateDataBenthos: 863 individuals were scored at known depths",length(CreateDataBenthos()$dist_m) == 863)
assert("GetSpeciesCountAtDepths: All 863 individuals must be seperated correctly at their depths",length(CreateDataBenthos()$dist_m) == sum(GetSpeciesCountAtDepths()$"2") + sum(GetSpeciesCountAtDepths()$"10"))
assert("GetSpeciesCountAtDepths: 20 species were scored at all depths",length(GetSpeciesCountAtDepths()$species_name) == 20)
assert("GetSelectedSpecies: 8 species were found at both depths, at each depth occurring at least thrice",length(GetSelectedSpecies()) == 8)
assert("GetDataBenthosSelected: 740 individuals of the 8 selected species were scored at known depths",length(GetDataBenthosSelected()$species_name) == 740)
assert("CreateDataRedox must only consider depths of 2 and 10 cm",length(subset(CreateDataRedox(),depth_cm != 2 & depth_cm != 10)$depth_cm) == 0)
assert("GetRedoxValues: 32 redox potentials are investigated",length(GetRedoxValues()) == 32)
assert("GetDistances: 16 distances are investigated",length(GetDistances()) == 16)
assert("GetDataCombined: needs to be 1176 Hydrobia",length(subset(GetDataCombined(),species_name == "Hydrobia_ulvae")$species_name) == 294)
assert("GetDataCombined: needs to be 36 Nereis diversicolor",length(subset(GetDataCombined(),species_name == "Nereis_diversicolor")$species_name) == 9)
assert("GetRedoxesHydrobia: must be 294 values",length(GetRedoxesHydrobia()$redox_calib) == 294)
assert("GetRedoxesNereis: must be 9 values",length(GetRedoxesNereis()$redox_calib) == 9)
assert("TallySpeciesPerRedox: must be 9 columns (redoxes and each species its frequency",length(TallySpeciesPerRedox()) == 9)
assert("TallySelectedSpeciesPerRedox: must be 3 columns (redoxes, Hydrobia and Nereis",length(TallySelectedSpeciesPerRedox()) == 3)
assert("CalcOrderednessPerDistance: must have 16 distances",length(CalcOrderednessPerDistance()$dist_m) == 16)

write.csv(GetSpeciesCountAtDepths(),file="table_species_count_at_depth.csv")
write.csv(GetSelectedSpecies(),file="table_selected_species.csv")
write.csv(CalcOrderednessPerDistance(),file="table_orderedness_per_distance.csv")
write.csv(TallySpeciesPerRedox(),file="table_redox_to_species.csv")
write.csv(TallySelectedSpeciesPerRedox(),file="table_redox_to_selected_species.csv")

# Do statistics per species
shapiro.test(GetRedoxesHydrobia()$redox_calib)
# Shapiro-Wilk normality test
# 
# data:  redoxes_hydrobia$redox_calib
# W = 0.7567, p-value < 2.2e-16

shapiro.test(GetRedoxesNereis()$redox_calib)
# Shapiro-Wilk normality test
# 
# data:  redoxes_nereis$redox_calib
# W = 0.8341, p-value = 0.04965

# Creates Figure_species_abundances_seperate.svg 
CreateFigureSpeciesAbundancesSeperate()

# Create Figure_species_abundances.svg 
CreateFigureSpeciesAbundances()

# Create Figure_redox_per_distance.svg
CreateFigureRedoxPerDistance()
