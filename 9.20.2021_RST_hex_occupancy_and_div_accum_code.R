#####--------------------------------------------------------------------------------------------

### Baltimore_street_tree_redline ###

#####--------------------------------------------------------------------------------------------

# Last updated: 20 Sept 2021
# Author: Karin Burghardt
# Contact: kburghar@umd.edu

setwd("/Volumes/GoogleDrive/My Drive/Projects/redlining_analysis/redline_balt_street")
date<-"2021_09_20"

# Load required libraries
library(readr)
library(car)
library(effects)
library(emmeans)
library(FD)
library(tidyverse)
library(ggplot2)
library('scales')
library(vegan)
library(reshape2)
library(dplyr)
library(lme4)
library(lmerTest)
library(stringr)
library(BiodiversityR)
library(iNEXT)
library(janitor)
library(cowplot)
library(MuMIn)

theme_set(theme_bw(base_family="Helvetica"))
options(scipen = 999, digits = 4)

#function to specify "not in"
`%nin%` <- Negate(`%in%`)

# HOLC custom colors
holc_col <- c("A"='#92BC6B', "B"='#92C7C9', "C"='#E7DC6B', "D"='#E47D67')

#import dataset- 95,119 trees (park trees still included)
st_tree <- read_csv("street_trees_Baltimore_w_HOLC_grades_HEX_2021-03-15.csv")

###This first bit is data wrangling to create CSV for composition analysis

#create neighbor,poly, hex column
st_tree$id_hex_poly<-paste(st_tree$holc_id.hex_id, st_tree$poly_id, sep=".")

#combine needed species column- must be done before character is changed to factor
st_tree$SPP <- str_replace_all(st_tree$SPP, 'Acer tataricum ginnala', 'Acer tataricum')
st_tree$SPP <- str_replace_all(st_tree$SPP, 'Populus nigra Italica', 'Populus nigra')
st_tree$SPP <- str_replace_all(st_tree$SPP, 'persimmon, Japanese', 'Diospyros kaki')
st_tree$SPP <- str_replace_all(st_tree$SPP, 'X Cupressocyparis leylandii', 'Cupressocyparis leylandii')
st_tree$SPP <- str_replace_all(st_tree$SPP, 'Gleditsia triacanthos inermis', 'Gleditsia triacanthos')
st_tree$SPPorig<-st_tree$SPP#column to preserve original species designation

#Replace all trees rated with condition as stump or dead so that they stay in dataset but get removed from SPP column
st_tree$CONDITION<-as.factor(st_tree$CONDITION)
Condition_inHex_percent <- tabyl(st_tree,CONDITION, show_na = FALSE)
st_tree$SPP[st_tree$CONDITION %in% c("Dead", "Stump","Stump w")] <- "Dead"

#name factors
st_tree$holc_id<-as.factor(st_tree$holc_id)
st_tree$holc_id.hex_id<-as.factor(st_tree$holc_id.hex_id)
st_tree$holc_grade<-as.factor(st_tree$holc_grade)
st_tree$poly_id<-as.factor(st_tree$poly_id)
st_tree$hex_id<-as.factor(st_tree$hex_id)
st_tree$id_hex_poly<-as.factor(st_tree$id_hex_poly)
st_tree$SPP<-as.factor(st_tree$SPP)
st_tree$COMMON<-as.factor(st_tree$COMMON)
st_tree$GENUS<-as.factor(st_tree$GENUS)
st_tree$SPACE_TYPE<-as.factor(st_tree$SPACE_TYPE)
st_tree$LOC_TYPE<-as.factor(st_tree$LOC_TYPE)
st_tree<-separate(data = st_tree, col = SPP, into = c("Genus", "species"), sep = "\\ ",remove = FALSE)
st_tree$Genus<-as.factor(st_tree$Genus)
st_tree$COLLECTOR<-as.factor(st_tree$COLLECTOR)
st_tree$CULTIVAR<-as.factor(st_tree$CULTIVAR)
summary(st_tree)
levels(st_tree$SPP)

#remove park trees= 5911 trees
st_tree2<-st_tree%>%
  filter(LOC_TYPE=="Street")

#check trees per neighborhood- smallest=B7 with 92 trees
numbertrees_per_holcid <- tabyl(st_tree2,holc_id, show_na = FALSE)

#create csv file for Meghan
write.csv(st_tree2,file=sprintf("output/st_tree_inHex%s.csv",date), row.names = TRUE)

#create list of trees per species to check
SPP_list_inHex_with_potential <- tabyl(st_tree2,SPP, show_na = FALSE)
print(SPP_list_inHex_with_potential)
write.csv(SPP_list_inHex_with_potential,file=sprintf("output/SPP_list_inHex_with_potential%s.csv",date), row.names = FALSE)


################Occupancy analysis and figures#############


#create list of SPP designations to make zero but keep as potential site in hex for abundance analysis as potential sites of trees
undesired <- c('Vacant Site', 'Vacant Potential', 'Stump','Vacant Site Not Suitable','NA',"Z Add 01"," ","Dead")

#create list of designations to keep as living trees but remove for diversity analysis
undesiredSPP <- c('unknown shrub','unknown tree','Ficus spp.','Fraxinus spp.','Hydrangea spp.','Ilex spp.','Ilex x','Juniperus spp.','Magnolia x','Photinia spp.','Populus spp.','Quercus spp.','Quercus x','Salix spp.','Ulmus spp.')


#create dummy abundance column for pivoting later so each potential site is hex is counted
st_tree2$abundance<-1

#create columns to quantify abundance based on size classes (S,M,L,empty)- need to do this to retain potenial sites when pivoting later

st_tree2<-st_tree2%>%
  mutate(small=ifelse(SPP %nin% undesired&DBH>=0&DBH<=5, 1, 0))

st_tree2<-st_tree2%>%
  mutate(large=ifelse(SPP %nin% undesired&DBH>=20, 1, 0))

st_tree2<-st_tree2%>%
  mutate(medium=ifelse(SPP %nin% undesired&DBH>5&DBH<20, 1, 0))

st_tree2<-st_tree2%>%
  mutate(empty=ifelse(SPP%in%undesired, 1, 0))

#####add sanity check column to check if we aren't doublecounting any trees and all are in one category only- result= all good!
st_tree2<-st_tree2%>%
  mutate(sum = rowSums(.[51:54]))

#####create a series of species x hex matrices for both occupancy analysis- note: this keeps in undesiredSPP because they are living trees!#####

#create holc_id_id_hex_poly x species matrix with all potential sites included
com.pot<-pivot_wider(st_tree2, id_cols=holc_id.hex_id, names_from = SPP, values_from = c(abundance),values_fn = list(abundance = sum),values_fill = 0)
com.pot<-column_to_rownames(com.pot, var = "holc_id.hex_id")

######create matrix with only trees than less than 5dbh for size-based occupancy analysis
com.small<-pivot_wider(st_tree2, id_cols=holc_id.hex_id, names_from = SPP, values_from = c(small),values_fn = list(small = sum),values_fill = 0)
com.small<-column_to_rownames(com.small, var = "holc_id.hex_id")

######create matrix with only trees than more than 5 dbh and less than 20 for size-based occupancy analysis
com.medium<-pivot_wider(st_tree2, id_cols=holc_id.hex_id, names_from = SPP, values_from = c(medium),values_fn = list(medium = sum),values_fill = 0)
com.medium<-column_to_rownames(com.medium, var = "holc_id.hex_id")

######create matrix with only trees than less than more than 20 for size-based occupancy analysis
com.large<-pivot_wider(st_tree2, id_cols=holc_id.hex_id, names_from = SPP, values_from = c(large),values_fn = list(large = sum),values_fill = 0)
com.large<-column_to_rownames(com.large, var = "holc_id.hex_id")

######create matrix of vacancies
com.empty<-pivot_wider(st_tree2, id_cols=holc_id.hex_id, names_from = SPP, values_from = c(empty),values_fn = list(empty = sum),values_fill = 0)
com.empty<-column_to_rownames(com.empty, var = "holc_id.hex_id")

#remove potential spots with no tree- only live trees left- use for abundance values for proportion analysis
com.live <- com.pot %>%
  select(-one_of(undesired))
summary(com.live)

# This matrix removes SPP designations that are unclear, repetitive or at genera level when species level ident are present for other trees.Some no longer exist with new hex size. Use: diversity analysis by hex level 
com <- com.live %>%
  select(-one_of(undesiredSPP))
summary(com)


# Create an environmental dataframe that includes total for occupancy analysis, add holc_grade factor back to enviro dataframe, create columns with totals for size classes and empties
env<-pivot_wider(st_tree2, id_cols=holc_id.hex_id, names_from = holc_grade)
env$holc_grade<-substr(env$holc_id.hex_id, 1, 1)
env<-separate(data = env, col = holc_id.hex_id, into = c("holc_id", "hex_id"), sep = "\\.",remove = FALSE)
env$holc_grade<-as.factor(env$holc_grade)
env$hex_id<-as.factor(env$hex_id)
env$holc_id<-as.factor(env$holc_id)
env<-column_to_rownames(env, var = "holc_id.hex_id")
env$site.totals.pot <- apply(com.pot,1,sum)
env$site.totals.live <- apply(com.live,1,sum)
env$site.totals <- apply(com,1,sum)
env$site.totals.small <- apply(com.small,1,sum)
env$site.totals.medium <- apply(com.medium,1,sum)
env$site.totals.large <- apply(com.large,1,sum)
env$site.totals.empty <- apply(com.empty,1,sum)
env$proportion_of_possible_sites_occupied <- env$site.totals.live/env$site.totals.pot
env$proportion_of_possible_sites_occupied_by_small <- env$site.totals.small/env$site.totals.pot
env$proportion_of_possible_sites_occupied_by_medium <- env$site.totals.medium/env$site.totals.pot
env$proportion_of_possible_sites_occupied_by_large <- env$site.totals.large/env$site.totals.pot
env$proportion_of_possible_sites_empty <- env$site.totals.empty/env$site.totals.pot
summary(env)



#######Create Occupancy figures

####### figures by totals not proportions#####

########create plot of potential tree spots per hex by grade- illustrates differences in # of potential locations for trees in D- that is why standardization by number of sites is needed####

tree_number_grade_hex_pot<-ggplot(env, aes(x=holc_grade, y=site.totals.pot,fill=holc_grade)) +
  geom_jitter(width = 0.2, colour="black",alpha=.2)+
  geom_boxplot(notch=TRUE,outlier.shape = NA,weight=5)+scale_fill_manual(values = holc_col)+
  labs(x = "HOLC Grade", y = "Potential tree sites/3.3 ha hexagon", colour = "HOLC grade", shape = "HOLC grade")

ggsave(file=sprintf("output/tree_number_grade_hex_pot%s.tiff",date), plot=tree_number_grade_hex_pot, width=6, height=4)

########create plot of living trees per hex by grade- illustrates differences in # of potential locations for trees in D- that is why standardization by number of sites is needed####

tree_number_grade_hex_live<-ggplot(env, aes(x=holc_grade, y=site.totals.live,fill=holc_grade)) +
  geom_jitter(width = 0.2, colour="black",alpha=.2)+
  geom_boxplot(notch=TRUE,outlier.shape = NA,weight=5)+scale_fill_manual(values = holc_col)+
  labs(x = "HOLC Grade", y = "Living street trees/3.3 ha hexagon", colour = "HOLC grade", shape = "HOLC grade")

ggsave(file=sprintf("output/tree_number_grade_hex_live%s.tiff",date), plot=tree_number_grade_hex_live, width=6, height=4)



#######Create multipanel boxplot#######

#create plot of empty tree spots per hex by grade

tree_number_grade_hex_empty<-ggplot(env, aes(x=holc_grade, y=proportion_of_possible_sites_empty,fill=holc_grade))+ geom_jitter(width = 0.2, colour="black",alpha=.2)+
  geom_boxplot(notch=TRUE,outlier.shape = NA,weight=5)+scale_fill_manual(values = holc_col)+
  labs(x = "HOLC Grade", y = "proportion of potential sites/3.3 ha hexagon", fill = "HOLC grade", shape = "HOLC grade")+ theme(legend.position = "none")

#create plot of proportion small trees per hex by grade

tree_number_grade_hex_small<-ggplot(env, aes(x=holc_grade, y=proportion_of_possible_sites_occupied_by_small,fill=holc_grade))+ geom_jitter(width = 0.2, colour="black",alpha=.2)+
  geom_boxplot(notch=TRUE,outlier.shape = NA,weight=5)+scale_fill_manual(values = holc_col)+
  labs(x = "HOLC Grade", y = "proportion of potential sites/3.3 ha hexagon", colour = "HOLC grade", shape = "HOLC grade")+ theme(legend.position = "none")

#create plot of proportion medium trees per hex by grade

tree_number_grade_hex_medium<-ggplot(env, aes(x=holc_grade, y=proportion_of_possible_sites_occupied_by_medium,fill=holc_grade))+ geom_jitter(width = 0.2, colour="black",alpha=.2)+
  geom_boxplot(notch=TRUE,outlier.shape = NA,weight=5)+scale_fill_manual(values = holc_col)+
  labs(x = "HOLC Grade", y = "proportion of potential sites/3.3 ha hexagon", colour = "HOLC grade", shape = "HOLC grade")+ theme(legend.position = "none")

#create plot of proportion large trees per hex by grade

tree_number_grade_hex_large<-ggplot(env, aes(x=holc_grade, y=proportion_of_possible_sites_occupied_by_large,fill=holc_grade))+ geom_jitter(width = 0.2, colour="black",alpha=.2)+
  geom_boxplot(notch=TRUE,outlier.shape = NA,weight=5)+scale_fill_manual(values = holc_col)+
  labs(x = "HOLC Grade", y = "proportion of potential sites/3.3 ha hexagon", colour = "HOLC grade", shape = "HOLC grade")+ theme(legend.position = "none")

#individual graphs

ggsave(file=sprintf("output/tree_number_grade_hex_empty%s.pdf",date), plot=tree_number_grade_hex_empty, width=3, height=4)
ggsave(file=sprintf("output/tree_number_grade_hex_small%s.pdf",date), plot=tree_number_grade_hex_small, width=3, height=4)
ggsave(file=sprintf("output/tree_number_grade_hex_medium%s.pdf",date), plot=tree_number_grade_hex_medium, width=3, height=4)
ggsave(file=sprintf("output/tree_number_grade_hex_large%s.pdf",date), plot=tree_number_grade_hex_large, width=3, height=4)


####create multipanel plot of proportion all categories of trees per hex by grade- an R update broke this code but individual graphs above can be combined

Fig3boxplot<-plot_grid(tree_number_grade_hex_empty, tree_number_grade_hex_large, tree_number_grade_hex_medium, tree_number_grade_hex_small, labels = c('I. No living tree', 'II. Large tree','III. Medium tree','IV. Small tree'), label_size = 12,ncol = 4, nrow = 1,hjust = 0.01, label_x = 0.24,vjust = -.2)+
  theme(plot.margin = unit(c(1,0,0,0), "lines"))

ggsave(file=sprintf("output/Fig3boxplot_proportion%s.pdf",date), plot=Fig3boxplot, width=7.5, height=3.5)


###### Occupancy models####

#Do the number of potential sites differ across holc_grades?
#check distribution of data- much more normal than I would think for count data but poisson better still!
ggplot(env,aes(site.totals.pot,fill = holc_grade, colour = holc_grade))+
  geom_density(alpha = 0.1) +
  xlim(0, 160)+ scale_colour_manual(values = holc_col)+scale_fill_manual(values = holc_col)

# glmer poisson
mod.pot<-glmer(site.totals.pot~holc_grade+(1|holc_id),data=env,family="poisson")
summary(mod.pot)
plot(mod.pot)
ranef(mod.pot)
r.squaredGLMM(mod.pot)
drop1(mod.pot, test="Chisq")

#Yes- potential sites per area (3.3 ha) varies across neighborhoods with D neighborhoods with the most potential locations (likely due to denser roads)
#So- need to use proportional occupancy analysis rather then living tree densities to account for different number of potential sites

####Occupancy analysis
###LARGE TREES

mod.large<-glmer(site.totals.large/site.totals.pot~holc_grade+(1|holc_id),weights=site.totals.pot,
           data=env,family="binomial")
summary(mod.large)
ranef(mod.large)
coef(mod.large)
summ(mod.large)
r.squaredGLMM(mod.large)
#this is the intercept only model I am using drop1() Chisq to compare to:
#mod.large.nul<-glmer(site.totals.large/site.totals.pot~1+(1|holc_id),weights=site.totals.pot,data=env,family="binomial")

mod.large.chi<-drop1(mod.large, test="Chisq")
print(mod.large.chi)

#If holc_grade significant than further do paired comparisons with emmeans: it is
emms.large<-emmeans(mod.large,~holc_grade,type = "response")
summary(emms.large)
emms.large.df = as.data.frame(emms.large)
pairs(emms.large,ratios = TRUE, type="response")
plot(emms.large,comparisons = TRUE) + theme_bw() + 
  labs(x = "Estimated marginal mean (Large tree in location- back-transformed)", y = "HOLC grade")

####(A&B;C;D) are the groups.

####Medium trees

mod.medium<-glmer(site.totals.medium/site.totals.pot~holc_grade+(1|holc_id),weights=site.totals.pot,
                 data=env,family="binomial")
summary(mod.medium)
ranef(mod.medium)
coef(mod.medium)
summ(mod.medium)
r.squaredGLMM(mod.medium)

#this is the intercept only model I am using drop1() Chisq to compare to:
#mod.medium.nul<-glmer(site.totals.medium/site.totals.pot~1+(1|holc_id),weights=site.totals.pot,data=env,family="binomial")

mod.medium.chi<-drop1(mod.medium, test="Chisq")
mod.medium.chi
#HOLC_GRADE_NOT_SIG- do not run comparisons

#### EMPTY LOCATIONS

mod.empty<-glmer(site.totals.empty/site.totals.pot~holc_grade+(1|holc_id),weights=site.totals.pot,
                  data=env,family="binomial")
summary(mod.empty)
ranef(mod.empty)
coef(mod.empty)
r.squaredGLMM(mod.empty)
#this is the intercept only model I am using drop1() Chisq to compare to:
#mod.empty.nul<-glmer(site.totals.empty/site.totals.pot~1+(1|holc_id),weights=site.totals.pot,data=env,family="binomial")

mod.empty.chi<-drop1(mod.empty, test="Chisq")
mod.empty.chi
#If holc_grade significant than further do paired comparisons with emmeans: Not different.

####SMALL

mod.small<-glmer(site.totals.small/site.totals.pot~holc_grade+(1|holc_id),weights=site.totals.pot,
                 data=env,family="binomial")
summary(mod.small)
ranef(mod.small)
coef(mod.small)
r.squaredGLMM(mod.small)

#this is the intercept only model I am using drop1() Chisq to compare to:
#mod.small.nul<-glmer(site.totals.small/site.totals.pot~1+(1|holc_id),weights=site.totals.pot,data=env,family="binomial")

mod.small.chi<-drop1(mod.small, test="Chisq")
mod.small.chi

#If holc_grade significant than further do paired comparisons with emmeans: It is!
emms.small<-emmeans(mod.small,~holc_grade,type = "response")
summary(emms.small)
emms.small.df = as.data.frame(emms.small)
pairs(emms.small,ratios = TRUE, type="response")
plot(emms.small,comparisons = TRUE) + theme_bw() + 
  labs(x = "Estimated marginal mean (small tree in location- back-transformed)", y = "HOLC grade")
pwpp(emms.small)

#B&D differ; no other pair-wise diffs


###### FIGURE 2: Species accumulation curves-extrapolation with I-Next package: 

st_tree_next<-st_tree2%>%
  filter(!SPP %in% undesired) %>% 
  filter(!is.na(SPP))%>% 
  filter(!SPP %in% undesiredSPP)

#create grade x species matrix for all trees
x<-pivot_wider(st_tree_next, id_cols=SPP, names_from = holc_grade, values_from = c(abundance),values_fn = list(abundance = sum),values_fill = 0)
x<-column_to_rownames(x, var = "SPP")

#ALLTREES FIGURE 1- rarefaction and estimates at all 3 q levels- takes 20 + mins to run

#estimate exact species diversity across holc grades standardized both by number of trees and coverage

size_rare_estD<-estimateD(x, datatype = "abundance", base = "size", level = NULL,conf = 0.95)
print(size_rare_estD)
write.csv(size_rare_estD,file=sprintf("output/size_rare_estD%s.csv",date), row.names = FALSE)

SC_rare_estD<-estimateD(x, "abundance", base="coverage", level=NULL, conf=0.95)
print(SC_rare_estD)
write.csv(SC_rare_estD,file=sprintf("output/SC_rare_estD%s.csv",date), row.names = FALSE)

#create rarefaction/extrapolation curves for figure

out.all<- iNEXT(x, q=c(0,1,2),datatype="abundance")
out.all$DataInfo # showing basic data information.
out.all$AsyEst # showing asymptotic diversity estimates.
out.all$iNextEst # showing diversity estimates with rarefied and extrapolated.

#summary rare/extrapolation figure

qlabels <- c("0" = "richness (q=0)", "1" = "Shannon's EFN (q=1)","2" = "Simpson's EFN (q=2)")

accum_alltree_all_q<-ggiNEXT(out.all, type=1, facet.var="order") + theme_bw(base_size=10)+ theme_bw(base_size=10)+ xlim(c(0,20000))+scale_colour_manual(values = holc_col,name="HOLC grade")+scale_fill_manual(values = holc_col, name="HOLC grade")+
  labs(x = "Number of individual trees sampled", y = "Tree species diversity", colour = "HOLC grade", shape = "HOLC grade", fill="HOLC grade")+
  theme(legend.position="bottom",legend.title=element_blank())+facet_wrap(~order, scales="free",labeller=labeller(order = qlabels))

#save summary figure

ggsave(file=sprintf("output/accum_alltree_all_div_pooled_indv%s.pdf",date), plot=accum_alltree_all_q, width=7, height=4)
ggsave(file=sprintf("output/accum_alltree_all_div_pooled_indv%s.tiff",date), plot=accum_alltree_all_q, width=7, height=4)


#summary rarefaction figure- simply ends same curves without extrapolation
out.all.rare<- iNEXT(x, q=c(0,1,2),datatype="abundance", endpoint=5238)

accum_alltree_all_q_rare<-ggiNEXT(out.all.rare, type=1, facet.var="order") +geom_line(size = .1, alpha=.2)+ theme_bw(base_size=10)+ theme_bw(base_size=10)+ xlim(c(0,6000))+scale_colour_manual(values = holc_col,name="HOLC grade")+scale_fill_manual(values = holc_col, name="HOLC grade")+
  labs(x = "Number of individual trees sampled", y = "Tree species diversity", colour = "HOLC grade", shape = "HOLC grade", fill="HOLC grade")+
  theme(legend.position="bottom",legend.title=element_blank())+facet_wrap(~order, scales="free",labeller=labeller(order = qlabels))

#save summary figure
ggsave(file=sprintf("output/accum_alltree_all_div_pooled_indv_rare%s.pdf",date), plot=accum_alltree_all_q_rare, width=7, height=4)
ggsave(file=sprintf("output/accum_alltree_all_div_pooled_indv_rare%s.tiff",date), plot=accum_alltree_all_q_rare, width=7, height=4)


#summary SC rarefaction figure
accum_alltree_all_q_rareSC<-ggiNEXT(out.all.rare, type=3, facet.var="order")+geom_line(size = .5, alpha=.9)+ theme_bw(base_size=10)+ theme_bw(base_size=10)+scale_colour_manual(values = holc_col,name="HOLC grade")+scale_fill_manual(values = holc_col, name="HOLC grade")+
  labs(x = "Sample coverage", y = "Tree species diversity", colour = "HOLC grade", shape = "HOLC grade", fill="HOLC grade")+
  theme(legend.position="bottom",legend.title=element_blank())+facet_wrap(~order, scales="free",labeller=labeller(order = qlabels))

#save summary figure
ggsave(file=sprintf("output/accum_alltree_all_div_pooled_indv_rareSC%s.pdf",date), plot=accum_alltree_all_q_rareSC, width=7, height=4)
ggsave(file=sprintf("output/accum_alltree_all_div_pooled_indv_rareSC%s.tiff",date), plot=accum_alltree_all_q_rareSC, width=7, height=4)








