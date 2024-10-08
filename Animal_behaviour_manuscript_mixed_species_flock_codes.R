setwd("D:/Papers/MSF_Acoustics/Analysis_Animal_behaviour/Repository/")
library(tidyverse)
## Phenotypic similarity analysis###

a<- read_csv("Matrix_visual_final.csv")
b<- read_csv("Body_variables_Birdtree_ID.csv")

colnames(b)
b<- b[,c("Species", "Body mass")]


dat<- left_join(b,a , by="Species")
names(dat)[1]<- ""




min(colSums(dat[,3:ncol(dat)]))
max(colSums(dat[,3:ncol(dat)]))#checking min and max for flock sizes in the dat file


func3<-function (d){
  sd.obs3<-numeric()
  for (i in 3:ncol(dat)){
    new <-d[d[,i] == 1,] #subset dat3 and create a "new" matrix such that the ith column has all 1s
    sd.obs3[i]<-sd(new$'Body mass') # in new look at the column bs and get the sd in bs for that flock size and store it in the i th position of the empty vector sd.obs
  }
  return(sd.obs3)
}

sd.obs3<-func3(dat) #apply func3 on dat and calculate sds for all flocks
#sd.obs3
sd.obs.mat3<-cbind(colSums(dat[,c(3:ncol(dat))]), sd.obs3[-c(1:2)]) #cbind sds and flock sizes for all flocks in dat
nrow(sd.obs.mat3)#checking no. of rows in the newly bound matrix
class(sd.obs.mat3)
sd.obs.mat3<-as.data.frame(sd.obs.mat3)

################# use this for final ses: obtaining sds for categories based on flock size 
sd.obs.mat3$cat<-cut(sd.obs.mat3[,1], breaks=c(0,5,10,22), labels = c("Small", "Medium", "Large"))

#Number of each flock richness class
cat_number<- sd.obs.mat3 %>% 
  group_by(cat) %>% 
  summarise(Number= length(cat))
cat_number

obs.sml.overall.last<-with(sd.obs.mat3, tapply(X=sd.obs.mat3[,2], INDEX = cat, FUN = mean))
###############################################################

library(EcoSimR)
emp.mat3<-lapply(1:1000,function(i) data.frame()) # creating a 1000 empty matrices
my.matrices3<-replicate(1000, sim5(dat[,3:ncol(dat)]),simplify = FALSE) # simulating 1000 matrices by randomizing the dat matrix

#binding the trait data which is present in the first three cols of dat to each of the newly randomized matrices and storing it in empty matrices that we created couple of steps ago
for(j in 1:length(emp.mat3)){
  emp.mat3[[j]]<-cbind(dat[,c(1:2)],my.matrices3[[j]])
}

list.sd3<- lapply(emp.mat3, func3) # applying func3 on all the newly created matrices. 
sd.null3<-matrix(unlist(list.sd3), ncol=1000, byrow=FALSE)

sd.null3<-sd.null3[-c(1:2),]
sd.null3<-as.data.frame(sd.null3)


#These are steps to club flock sizes into S, M, L without averaging within richness (i.e. avg for flock size 2, flock size 3 and so on)
fl.sd.sml<-matrix(nrow=3, ncol = 1000)




for (k in 1:1000){
  fl.sd.temp<-matrix(nrow=180,ncol = 2 )
  fl.sd.temp<-cbind(colSums(dat[,-c(1:2)]),sd.null3[,k]) #for each column in sd.null3, bind flock sizes to it
  fl.sd.temp<-as.data.frame(fl.sd.temp)
  fl.sd.temp$cat<-cut(fl.sd.temp[,1], breaks=c(0,5,10,22), labels = c("Small", "Medium", "Large"))# for the temporary vector create a new column cat and assign S,M,L according to sizes 2,3,4=small, 5-9=medium, >10=large
  fl.sd.sml[,k]<-with(fl.sd.temp, tapply(X=fl.sd.temp[,2], INDEX = cat, FUN = mean))# get mean for every "cat"
}



mean(fl.sd.sml[1,])# mean for row one which is SMALL flocks
mean(fl.sd.sml[2,])# mean for row two which is MEDIUM flocks
mean(fl.sd.sml[3,])# mean for row three which is LARGE flocks

sd(fl.sd.sml[1,])# SD for row one which is SMALL flocks
sd(fl.sd.sml[2,])# SD for row two which is MEDIUM flocks
sd(fl.sd.sml[3,])# SD for row three which is LARGE flocks


#creating category-wise summary for Obs mean, Exp mean and SD in exp mean 
null.means<-c(mean(fl.sd.sml[1,]),mean(fl.sd.sml[2,]),mean(fl.sd.sml[3,]))
null.sds<-c(sd(fl.sd.sml[1,]),sd(fl.sd.sml[2,]),sd(fl.sd.sml[3,]))
obs.means<-obs.sml.overall.last

fin<-cbind(null.means,null.sds,obs.means)

#SESs based on "fin"
effectsizes<-numeric() #numeric ses to store the SESs
for (y in 1:3){ # for all rows, SES=(mean(obs)-mean(exp))/sd(exp)
  effectsizes[y]<- ((fin[y,3]-fin[y,1])/fin[y,2])
}

effectsizes


##################END############################



##Calculate PC scores of acoustic variables###

##Calculating PC scores

dat<- read_csv("Acoustic_pcf_final.csv")

n_distinct(dat$Species)

dat_pcf<- read_csv("Acoustic_pcf_final.csv")

dat_pcf<- dat_pcf[, -c(17,18)]#Removing columns not needed for the analysis


dat.pca<- prcomp(dat_pcf[, c(7:14)], center = T, scale. = T)
summary(dat.pca)
dat.pca$rotation
library(factoextra)
get_eigenvalue(dat.pca)


data<- cbind(dat_pcf, dat.pca$x[,1:2])

##Extracting the mean Pc1 scores for each species
data1<- data %>% 
  group_by(Species) %>% 
  summarise(PC1= mean(PC1))

View(data1)

dat.pca<- prcomp(dat_pcf[, c(7:14)], center = T, scale. = T)
summary(dat.pca)
dat.pca$rotation
library(factoextra)
get_eigenvalue(dat.pca)


data<- cbind(dat_pcf, dat.pca$x[,1:2])

##Extracting the mean Pc2 scores for each species
data2<- data %>% 
  group_by(Species) %>% 
  summarise(PC1= mean(PC1))
View(data2)
##Saving the file
##Each species has an associated mean PC score## The analysis remains the same except for the different categories associated with flock rishness and vocal activity categories##

##################END#############################



##Correlation between beak and PC1 scores of species###

## PC1 vs beak morphology

dat<- read_csv("PC1_species.csv")
beak<- read_csv("Body_variables_Birdtree_ID.csv")
dat<- left_join(dat, beak, by="Species")
colnames(dat)

library(ggpubr)

ggscatter(dat, x =  "PC1" , y ="Beak.Length_Culmen",
          add = "reg.line", # Add regression line
          conf.int = TRUE, # Add confidence interval
          add.params = list(color = "red",
                            fill = "lightgray")
)+
  stat_cor(method = "spearman")+
  ylab("Beak length (culmen) (mm)")+
  xlab("Mean PC1")

ggscatter(dat, x =  "PC1" , y ="Beak.Length_Nares",
          add = "reg.line", # Add regression line
          conf.int = TRUE, # Add confidence interval
          add.params = list(color = "red",
                            fill = "lightgray")
)+
  stat_cor(method = "spearman")+
  ylab("Beak length (nares) (mm)")+
  xlab("Mean PC1")


ggscatter(dat, x = "PC1", y = "Beak.Width",
          add = "reg.line", # Add regression line
          conf.int = TRUE, # Add confidence interval
          add.params = list(color = "green",
                            fill = "lightgray")
)+
  stat_cor(method = "spearman")+
  ylab("Beak width (mm))")+
  xlab("Mean PC1")


ggscatter(dat, x = "PC1", y = "Beak.Depth",
          add = "reg.line", # Add regression line
          conf.int = TRUE, # Add confidence interval
          add.params = list(color = "orange",
                            fill = "lightgray")
)+
  stat_cor(method = "spearman")+
  ylab("Beak delpth (mm))")+
  xlab("Mean PC1")



## Phylogenetic correction

library(ape)
library(phytools)
library(nlme)
library(visreg)
mytree<- read.nexus("output.nex")
#
mytree
class(mytree)

#Using consensus.edge function
trees<- consensus.edges(mytree, method = "least.squares", rooted= TRUE)

plot(trees, cex=0.7)



# Vector of species (tips) to be removed##because there are only 54 vocalising species
tips_to_remove <- c("Turdus_merula","Coracina_melanoptera", "Oriolus_chinensis", "Oriolus_oriolus", "Harpactes_fasciatus", "Pycnonotus_melanicterus", "Picus_squamatus", "Eumyias_thalassinus")

# Drop the specified tips from the tree
trees <- drop.tip(trees, tips_to_remove, trim.internal = TRUE, 
                  subtree = FALSE, root.edge = 0, 
                  rooted = is.rooted(trees), 
                  collapse.singles = TRUE, interactive = FALSE)

dat$`Species Bird tree`

#Converting the species ID column into row names
mydata<- as.data.frame(dat)
rownames(mydata) <- mydata$`Species Bird tree`
#Reordering the rows to match the tree
mydata <- mydata[match(trees$tip.label,rownames(mydata)),]
#

trees$tip.label
mydata$`Species Bird tree`



#Independant contrasts
trees<- multi2di(trees) ##You can resolve multifurcations with branches of zero length using
#multi2di(), e.g.:
x1<- pic(dat$PC1, trees)
x1

y1<- pic(dat$Beak.Length_Culmen, trees)
y1

z1<- pic(dat$Beak.Length_Nares, trees)
z1

w1<- pic(dat$Beak.Width, trees)
p1<- pic(dat$Beak.Depth, trees)
phy_cor<- data.frame(x1, y1, z1, w1, p1)
names(phy_cor)<- c("PC1", "Beak length: culmen", "Beak length: nares", "Beak width", "Beak Depth")
head(phy_cor)

#Spearman's correlation
ggscatter(phy_cor, x = "PC1", y = "Beak length: culmen",
          add = "reg.line", # Add regression line
          conf.int = TRUE, # Add confidence interval
          add.params = list(color = "red",
                            fill = "lightgray")
)+
  stat_cor(method = "spearman")+
  ylab("PIC (Beak length: culmen (mm)")+
  xlab("PIC (Mean PC1)")

#Spearman's correlation
ggscatter(phy_cor, x = "PC1", y = "Beak length: nares",
          add = "reg.line", # Add regression line
          conf.int = TRUE, # Add confidence interval
          add.params = list(color = "green",
                            fill = "lightgray")
)+
  stat_cor(method = "spearman")+
  ylab(" PIC (Beak length: nares (mm)")+
  xlab("PIC (Mean PC1)")

#Spearman's correlation
ggscatter(phy_cor, x = "PC1", y = "Beak width",
          add = "reg.line", # Add regression line
          conf.int = TRUE, # Add confidence interval
          add.params = list(color = "orange",
                            fill = "lightgray")
)+
  stat_cor(method = "spearman")+
  ylab("PIC (Beak width (mm)")+
  xlab("PIC (Mean PC1)")

#Spearman's correlation
ggscatter(phy_cor, x = "PC1", y = "Beak Depth",
          add = "reg.line", # Add regression line
          conf.int = TRUE, # Add confidence interval
          add.params = list(color = "grey",
                            fill = "lightgray")
)+
  stat_cor(method = "spearman")+
  ylab("PIC(Beak depth (mm)")+
  xlab("PIC (Mean PC1)")



#################END###################################

###Faith's PD analysis####

library(tidyverse)

##Calculating Faiths PD flocks the flock matrix and randomizing the flocks 
#to check if what is observed is withing expectation by chance
##Code adapted from Rohit Caravansary's work

## Calculate PD for all flocks
library(picante)

mytree<- read.nexus("output1.nex")
#
mytree # prints basic tree information
class(mytree)

#Using consensus.edge function
trees<- consensus.edges(mytree, method = "least.squares", rooted= T)
plot(trees, cex=0.7)
trees$tip.label

outgroup<- "Spilornis_cheela"

# Root the tree using the outgroup
rooted_tree <- root(trees, outgroup = outgroup, resolve.root = TRUE)

# Plot the rooted tree
plot(rooted_tree, main = "Rooted Tree", cex=0.7)


##Input matrix 
dat<- read_csv("Matrix_visual_final_FPD.csv")

dat<- as.data.frame(dat)


# Extract the first column as character elements
row_names <- as.character(dat[, 1])

# Set the row names of the matrix to be the first column elements
rownames(dat) <- row_names

# Display the modified matrix with updated row names
print(dat)

# # #Reordering the rows to match the tree
dat <- dat[match(rooted_tree$tip.label,rownames(dat)),]
dat<- na.omit(dat) #Removing the outgroup because it is not present in the matrix
data<-t(dat)
#write.csv(data, "data_for_pd1_test.csv")
data<- read_csv("data_for_pd1.csv")

data<- as.data.frame(data)
data<- data[-1,]
# Extract the first column to be used as row names
rownames(data) <- data[, 1]

# Remove the first column
data <- data[, -1]


# Print the modified matrix
print(data)




# Faith's PD
PD<-pd(data,rooted_tree,include.root=T) #Faith's PD using community matrix and pruned tree
PD #worked

###Alternate path to calculate the Faith's pd ses values
a<- ses.pd(data, rooted_tree,null.model = "independentswap", runs=999,
           iterations=1000)

# 
ses<- a$pd.obs.z
library(ggpubr)

# Create a smoothed density plot with interesting colors
plot1 <- ggplot(data.frame(value = ses), aes(x = ses)) +
  geom_density(fill = "skyblue", color = "darkblue", alpha = 0.7) +
  xlab("SES values across flocks") +
  ylab("Probabiliy density")+
  theme_classic2()+
  theme(
    axis.title.x = element_text(size = 20),  # Increase x-axis title font size
    axis.title.y = element_text(size = 20),
    axis.text.x = element_text(size = 14),  # Increase x-axis labels font size
    axis.text.y = element_text(size = 14)   # Increase y-axis labels font size# Increase y-axis title font size
  )
plot1

############################END###############################


####Blomberg's K######

# Calculating Blomber's k for the bodymass of 62 flock participants using rooted tree



treebm<- read.nexus("output.nex")
treebm

treesbm<- consensus.edges(treebm, method = "least.squares", rooted= T)
plot(treesbm, cex=0.7)


bodymass<- read_csv("Body_variables_Birdtree_ID.csv")
bodymass<- as.data.frame(bodymass)
rownames(bodymass)<- bodymass$`Species Bird tree`

##Matching the rows to match the tree
treesbm$tip.label


dat<- bodymass[match(treesbm$tip.label, rownames(bodymass)),]
rownames(dat)

treesbm2<- multi2di(treesbm)




# Calculate Blomberg's K
k_result <- phylosig(treesbm, dat$`Body mass`, method = "K", test=TRUE)
print(k_result)


# Calculating Blomber's k for the PC1 of 54 vocalising participants



pc1<- read_csv("PC1_species.csv")
pc1<- left_join(pc1,bodymass, by= "Species")

pc1<- pc1[,c(7,3)]
pc1<- as.data.frame(pc1)
rownames(pc1)<- pc1$`Species Bird tree`



tree2<- read.nexus("output.nex")
tree2

trees2<- consensus.edges(tree2, method = "least.squares", rooted= TRUE)
plot(trees2, cex=0.7)


# Vector of species (tips) to be removed##because there are only 54 vocalising species
tips_to_remove <- c("Turdus_merula","Coracina_melanoptera", "Oriolus_chinensis", "Oriolus_oriolus", "Harpactes_fasciatus", "Pycnonotus_melanicterus", "Picus_squamatus", "Eumyias_thalassinus")

# Drop the specified tips from the tree
treespc <- drop.tip(trees2, tips_to_remove, trim.internal = TRUE, 
                    subtree = FALSE, root.edge = 0, 
                    rooted = is.rooted(trees2), 
                    collapse.singles = TRUE, interactive = FALSE)




##Matching the rows to match the tree

dat<- pc1[match(treespc$tip.label, rownames(pc1)),]

rownames(dat)
treespc$tip.label


treespc1<- multi2di(treespc)



# Calculate Blomberg's K
k_result <- phylosig(treespc1, dat$PC1, method = "K", test=TRUE)
print(k_result)

##################END#####################################


##Phylogenetically corrected body masss against vocal activity##

c<- read_csv("correlation.csv")

##Plotting scattaer plot using spearmans method

ggscatter(c, x = "Body mass", y = "Acoustic",
          add = "reg.line",                                 # Add regression line
          conf.int = TRUE,                                  # Add confidence interval
          add.params = list(color = "red",
                            fill = "lightgray")
)+
  stat_cor(method = "spearman")+
  ylab("Vocal activity index")+
  xlab("Body mass (g)")


mytree<-  read.nexus("output2.nex")
# 
mytree               # prints basic tree information



class(mytree)

#Using consensus.edge function
trees<- consensus.edges(mytree, method = "least.squares", rooted= TRUE)
plot(trees)



trees$tip.label

mydata<- read_csv("correlation.csv")
#Converting the species ID column into row names
mydata<- as.data.frame(mydata)

rownames(mydata) <- mydata$`Species ID`

#Reordering the rows to match the tree
mydata <- mydata[match(trees$tip.label,rownames(mydata)),]

# 
mydat2 <- mydata[, c(3,4)]

library(phytools)
dotTree(trees, as.matrix(mydat2)[,c("Acoustic")], method = "plotTree") # plot trait "Acoustic" at tree tips 




dotTree(trees, as.matrix(mydat2)[,c("Body mass")], method= "plotTree") # plot trait "Acoustic" at tree tips 

#Independant contrasts

trees<- multi2di(trees) ##You can resolve multifurcations with branches of zero length using
#multi2di(), e.g.:

x1<- pic(mydata$`Body mass`, trees)
x1
y1<- pic(mydata$Acoustic, trees)
y1

phy_cor<- data.frame(x1, y1)
names(phy_cor)<- c("Body mass", "Acoustic")


#Spearman's correlation

ggscatter(phy_cor, x = "Body mass", y = "Acoustic",
          add = "reg.line",                                 # Add regression line
          conf.int = TRUE,                                  # Add confidence interval
          add.params = list(color = "red",
                            fill = "lightgray")
)+
  stat_cor(method = "spearman")+
  ylab("Vocal activity index")+
  xlab("Body mass (g)")


###########END###################################



