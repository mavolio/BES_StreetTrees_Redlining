#####--------------------------------------------------------------------------------------------

### Baltimore_street_tree_redline ###

#####--------------------------------------------------------------------------------------------

# Last updated: 3 May 2021
# Author: Karin Burghardt
# Contact: kburghar@umd.edu

setwd("/Volumes/GoogleDrive/My Drive/Projects/redlining_analysis/redline_balt_street")
date<-"2021_05_01"

# Load required libraries
library(readr)
library(car)
library(effects)
library(emmeans)
library(FD)
library(nlme)
library(piecewiseSEM)
library(plotrix)
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
library("cowplot")
theme_set(theme_bw(base_family="Helvetica"))
options(scipen = 999, digits = 4)

#function to specify "not in"
`%nin%` <- Negate(`%in%`)

# HOLC custom colors
holc_col <- c("A"='#92BC6B', "B"='#92C7C9', "C"='#E7DC6B', "D"='#E47D67')

#st_tree <- read_csv("street_trees_Baltimore_w_HOLC_grades_2021-01-15.csv")
#st_tree <- read_csv("street_trees_Baltimore_w_HOLC_grades_HEX_2021-01-27.csv")
#st_tree <- read_csv("street_trees_Baltimore_w_HOLC_grades_HEX_2021-02-03.csv")
#st_tree <- read_csv("street_trees_Baltimore_w_HOLC_grades_HEX_2021-03-04.csv")
st_tree <- read_csv("street_trees_Baltimore_w_HOLC_grades_HEX_2021-03-15.csv")


###This first bit is data wrangling to create CSV for Meghan's analysis

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

#create list of species to check
SPP_list_inHex <- tabyl(st_tree2,SPP, show_na = FALSE)
print(SPP_list_inHex)
write.csv(SPP_list_inHex,file=sprintf("output/SPP_list_inHex%s.csv",date), row.names = FALSE)


################Occupancy analysis and figures#############


#create list of SPP designations to make zero but keep as potenial site in hex for abundance analysis as potential sites of trees
undesired <- c('Vacant Site', 'Vacant Potential', 'Stump','Vacant Site Not Suitable','NA',"Z Add 01"," ","Dead")

#create list of designations to keep as living trees but remove for diversity analysis
undesiredSPP <- c('unknown shrub','unknown tree','Acer spp.','Carya spp.','Cornus spp.','Cornus x','Ficus spp.','Fraxinus spp.','Hydrangea spp.','Ilex spp.','Ilex x','Juniperus spp.','Magnolia spp.','Magnolia x','Photinia spp.','Picea spp.','Populus spp.','Quercus spp.','Quercus x','Rhus spp.','Salix spp.','Ulmus spp.')

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

# This matrix removes SPP designations that are unclear, repetitive or at genera level.Some no longer exist with new hex size. Use: diversity analysis by hex level 
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

########create plot of potential tree spots per hex by grade- illustrates differences in # of potential locations for trees in D- that is why standardization by number of sites is needed####

tree_number_grade_hex_pot<-ggplot(env, aes(x=holc_grade, y=site.totals.pot,fill=holc_grade)) +
  geom_jitter(width = 0.2, colour="black",alpha=.2)+
  geom_boxplot(notch=TRUE)+scale_fill_manual(values = holc_col)+
  labs(x = "HOLC Grade", y = "Possible street tree locations per hexagon", colour = "HOLC grade", shape = "HOLC grade")

#ggsave(file=sprintf("output/tree_number_grade_hex_pot%s.pdf",date), plot=tree_number_grade_hex_pot, width=6, height=4)

########create plot of living trees per hex by grade- illustrates differences in # of potential locations for trees in D- that is why standardization by number of sites is needed####

tree_number_grade_hex_live<-ggplot(env, aes(x=holc_grade, y=site.totals.live,fill=holc_grade)) +
  geom_jitter(width = 0.2, colour="black",alpha=.2)+
  geom_boxplot(notch=TRUE)+scale_fill_manual(values = holc_col)+
  labs(x = "HOLC Grade", y = "Living street trees per hexagon", colour = "HOLC grade", shape = "HOLC grade")

#ggsave(file=sprintf("output/tree_number_grade_hex_pot%s.pdf",date), plot=tree_number_grade_hex_pot, width=6, height=4)


#####standardized plots####

#create density plot by proportion of empty sites available
density_st_tree_hex_occupied_empty<-ggplot(env,aes(proportion_of_possible_sites_empty,fill = holc_grade, colour = holc_grade))+
  geom_density(alpha = 0.1) +
  xlim(0, 1)+ scale_colour_manual(values = holc_col)+scale_fill_manual(values = holc_col)+ theme(legend.position = "none")

#create density plot by proportion of possible site occupied by small trees
density_st_tree_hex_occupied_small<-ggplot(env,aes(proportion_of_possible_sites_occupied_by_small,fill = holc_grade, colour = holc_grade))+
  geom_density(alpha = 0.1) +
  xlim(0, 1)+ scale_colour_manual(values = holc_col)+scale_fill_manual(values = holc_col)+ theme(legend.position = "none")

#create density plot by proportion of possible site occupied by medium trees
density_st_tree_hex_occupied_medium<-ggplot(env,aes(proportion_of_possible_sites_occupied_by_medium,fill = holc_grade, colour = holc_grade))+
  geom_density(alpha = 0.1) +
  xlim(0, 1)+ scale_colour_manual(values = holc_col)+scale_fill_manual(values = holc_col)+ theme(legend.position = "none")

#create density plot by proportion of possible site occupied by large trees
density_st_tree_hex_occupied_large<-ggplot(env,aes(proportion_of_possible_sites_occupied_by_large,fill = holc_grade, colour = holc_grade))+
  geom_density(alpha = 0.1) +
  xlim(0, 1)+ scale_colour_manual(values = holc_col)+scale_fill_manual(values = holc_col)+ theme(legend.position = "none")

#individual figures
#ggsave(file=sprintf("output/density_st_tree_hex_occupied_empty%s.pdf",date), plot=density_st_tree_hex_occupied_empty, width=5, height=5)
#ggsave(file=sprintf("output/density_st_tree_hex_occupied_small%s.pdf",date), plot=density_st_tree_hex_occupied_small, width=5, height=5)
#ggsave(file=sprintf("output/density_st_tree_hex_occupied_medium%s.pdf",date), plot=density_st_tree_hex_occupied_medium, width=5, height=5)
#ggsave(file=sprintf("output/density_st_tree_hex_occupied_large%s.pdf",date), plot=density_st_tree_hex_occupied_large, width=5, height=5)

#Create density Panel figure
Fig3densityplot<-plot_grid(density_st_tree_hex_occupied_empty, density_st_tree_hex_occupied_large, density_st_tree_hex_occupied_medium, density_st_tree_hex_occupied_small, 
                           labels = c('A.empty sites', 'B.Large trees','C.Medium trees','D. Small trees'), label_size = 10,
                           ncol = 4, nrow = 1)

ggsave(file=sprintf("output/Fig3density_proportion%s.pdf",date), plot=Fig3densityplot, width=7.5, height=3.5)



#######Create multipanel boxplot#######

#create plot of empty tree spots per hex by grade

tree_number_grade_hex_empty<-ggplot(env, aes(x=holc_grade, y=proportion_of_possible_sites_empty,fill=holc_grade))+ geom_jitter(width = 0.2, colour="black",alpha=.2)+
  geom_boxplot(notch=TRUE,outlier.shape = NA,weight=5)+scale_fill_manual(values = holc_col)+
  labs(x = "HOLC Grade", y = "proportion of potential sites in hexagon", fill = "HOLC grade", shape = "HOLC grade")+ theme(legend.position = "none")

#create plot of proportion small trees per hex by grade

tree_number_grade_hex_small<-ggplot(env, aes(x=holc_grade, y=proportion_of_possible_sites_occupied_by_small,fill=holc_grade))+ geom_jitter(width = 0.2, colour="black",alpha=.2)+
  geom_boxplot(notch=TRUE,outlier.shape = NA,weight=5)+scale_fill_manual(values = holc_col)+
  labs(x = "HOLC Grade", y = "proportion of potential sites in hexagon", colour = "HOLC grade", shape = "HOLC grade")+ theme(legend.position = "none")

#create plot of proportion medium trees per hex by grade

tree_number_grade_hex_medium<-ggplot(env, aes(x=holc_grade, y=proportion_of_possible_sites_occupied_by_medium,fill=holc_grade))+ geom_jitter(width = 0.2, colour="black",alpha=.2)+
  geom_boxplot(notch=TRUE,outlier.shape = NA,weight=5)+scale_fill_manual(values = holc_col)+
  labs(x = "HOLC Grade", y = "proportion of potenital sites in hexagon", colour = "HOLC grade", shape = "HOLC grade")+ theme(legend.position = "none")

#create plot of proportion large trees per hex by grade

tree_number_grade_hex_large<-ggplot(env, aes(x=holc_grade, y=proportion_of_possible_sites_occupied_by_large,fill=holc_grade))+ geom_jitter(width = 0.2, colour="black",alpha=.2)+
  geom_boxplot(notch=TRUE,outlier.shape = NA,weight=5)+scale_fill_manual(values = holc_col)+
  labs(x = "HOLC Grade", y = "proportion of potential sites in hexagon", colour = "HOLC grade", shape = "HOLC grade")+ theme(legend.position = "none")

#indivdual graphs

#ggsave(file=sprintf("output/tree_number_grade_hex_empty%s.pdf",date), plot=tree_number_grade_hex_empty, width=3, height=4)
#ggsave(file=sprintf("output/tree_number_grade_hex_small%s.pdf",date), plot=tree_number_grade_hex_small, width=3, height=4)
#ggsave(file=sprintf("output/tree_number_grade_hex_medium%s.pdf",date), plot=tree_number_grade_hex_medium, width=3, height=4)
#ggsave(file=sprintf("output/tree_number_grade_hex_large%s.pdf",date), plot=tree_number_grade_hex_large, width=3, height=4)


####create multipanel plot of proportion small trees per hex by grade

Fig3boxplot<-plot_grid(tree_number_grade_hex_empty, tree_number_grade_hex_large, tree_number_grade_hex_medium, tree_number_grade_hex_small, labels = c('I. No living tree', 'II. Large tree','III. Medium tree','IV. Small tree'), label_size = 12,ncol = 4, nrow = 1,hjust = 0.01, label_x = 0.24,vjust = -.2)+
  theme(plot.margin = unit(c(1,0,0,0), "lines"))

ggsave(file=sprintf("output/Fig3boxplot_proportion%s.pdf",date), plot=Fig3boxplot, width=7.5, height=3.5)


####### figures by totals not proportions#####


#create density plot by tree totals.living

density_st_tree_hex<-ggplot(env,aes(site.totals.live,fill = holc_grade, colour = holc_grade))+
  geom_density(alpha = 0.1) +
  xlim(0, 160)+ scale_colour_manual(values = holc_col)+scale_fill_manual(values = holc_col)

ggsave(file=sprintf("output/density_st_tree_hex%s.pdf",date), plot=density_st_tree_hex, width=5, height=5)



###### Occupancy models####

#check distribution of data- much more normal than I would think for count data!
ggplot(env,aes(site.totals.pot,fill = holc_grade, colour = holc_grade))+
  geom_density(alpha = 0.1) +
  xlim(0, 160)+ scale_colour_manual(values = holc_col)+scale_fill_manual(values = holc_col)

#linear mixed has better fit of residuals than glmer poisson

mod<-lmer(site.totals.pot~holc_grade+(1|holc_id),data=env)
summary(mod)
plot(mod)
drop1(mod, test="Chisq")
emms.pot<-emmeans(mod,~holc_grade,mode = "satterthwaite")
emm1df = as.data.frame(emms.pot)
plot(emms.pot,comparisons = TRUE) + theme_bw() + 
  labs(x = "Estimated marginal mean (potential tree locations)", y = "HOLC grade")

#B and D are only different pair

####

mod.large<-glmer(site.totals.large/site.totals.pot~holc_grade+(1|holc_id),weights=site.totals.pot,
           data=env,family="binomial")
summary(mod.large)
ranef(mod.large)
coef(mod.large)

#this is the intercept only model I am using drop1() Chisq to compare to:
#mod.large.nul<-glmer(site.totals.large/site.totals.pot~1+(1|holc_id),weights=site.totals.pot,data=env,family="binomial")

drop1(mod.large, test="Chisq")
#If holc_grade significant than further do paired comparisons with emmeans: it is
emms.large<-emmeans(mod.large,~holc_grade,mode = "satterthwaite",type = "response")
emms.large.df = as.data.frame(emms.large)
pairs(emms.large)
plot(emms.large,comparisons = TRUE) + theme_bw() + 
  labs(x = "Estimated marginal mean (Large tree in site:log odds)", y = "HOLC grade")

####(A&B;C;D) are the groups.

####

mod.medium<-glmer(site.totals.medium/site.totals.pot~holc_grade+(1|holc_id),weights=site.totals.pot,
                 data=env,family="binomial")
summary(mod.medium)
ranef(mod.medium)
coef(mod.medium)

#this is the intercept only model I am using drop1() Chisq to compare to:
#mod.medium.nul<-glmer(site.totals.medium/site.totals.pot~1+(1|holc_id),weights=site.totals.pot,data=env,family="binomial")

drop1(mod.medium, test="Chisq")

#HOLC_GRADE_NOT_SIG- do not run comparisions

#emms.medium<-emmeans(mod.medium,~holc_grade,mode = "satterthwaite",type = "response")
#emms.medium.df = as.data.frame(emms.medium)
#pairs(emms.medium)
#plot(emms.medium,comparisons = TRUE) + theme_bw() + 
 # labs(x = "Estimated marginal mean (medium tree in site:log odds)", y = "HOLC grade")

####

mod.empty<-glmer(site.totals.empty/site.totals.pot~holc_grade+(1|holc_id),weights=site.totals.pot,
                  data=env,family="binomial")
summary(mod.empty)
ranef(mod.empty)
coef(mod.empty)

#this is the intercept only model I am using drop1() Chisq to compare to:
#mod.empty.nul<-glmer(site.totals.empty/site.totals.pot~1+(1|holc_id),weights=site.totals.pot,data=env,family="binomial")

drop1(mod.empty, test="Chisq")
#If holc_grade significant than further do paired comparisons with emmeans: Not different.


####

mod.small<-glmer(site.totals.small/site.totals.pot~holc_grade+(1|holc_id),weights=site.totals.pot,
                 data=env,family="binomial")
summary(mod.small)
ranef(mod.small)
coef(mod.small)

#this is the intercept only model I am using drop1() Chisq to compare to:
#mod.small.nul<-glmer(site.totals.small/site.totals.pot~1+(1|holc_id),weights=site.totals.pot,data=env,family="binomial")

drop1(mod.small, test="Chisq")
#If holc_grade significant than further do paired comparisons with emmeans: It is!
emms.small<-emmeans(mod.small,~holc_grade,mode = "satterthwaite",type = "response")
emms.small.df = as.data.frame(emms.small)
pairs(emms.small)
plot(emms.small,comparisons = TRUE) + theme_bw() + 
  labs(x = "Estimated marginal mean (small tree in site:log odds)", y = "HOLC grade")


#B&D differ; no other pair-wise diffs



###### FIGURE 2: Species accumulation curves-extrapolation with I-Next package: 

st_tree_next<-st_tree2%>%
  filter(!SPP %in% undesired) %>% 
  filter(!is.na(SPP))%>% 
  filter(!SPP %in% undesiredSPP)

#create grade x species matrix for all trees
x<-pivot_wider(st_tree_next, id_cols=SPP, names_from = holc_grade, values_from = c(abundance),values_fn = list(abundance = sum),values_fill = 0)
x<-column_to_rownames(x, var = "SPP")

#ALLTREES FIGURE- rarefaction at all q levels

Holc_grade_est_S<-estimateD(x)

out.all<- iNEXT(x, q=c(0,1,2),datatype="abundance")
out.all$DataInfo # showing basic data information.
out.all$AsyEst # showing asymptotic diversity estimates.
out.all$iNextEst # showing diversity estimates with rarefied and extrapolated.

#summary rarefaction figure

qlabels <- c("0" = "richness (q=0)", "1" = "Shannon's EFN (q=1)","2" = "Simpson's EFN (q=2)")

accum_alltree_all_q<-ggiNEXT(out.all, type=1, facet.var="order") + theme_bw(base_size=10)+ theme_bw(base_size=10)+ xlim(c(0,20000))+scale_colour_manual(values = holc_col,name="HOLC grade")+scale_fill_manual(values = holc_col, name="HOLC grade")+
  labs(x = "Number of individual trees sampled", y = "Tree species diversity", colour = "HOLC grade", shape = "HOLC grade", fill="HOLC grade")+
  theme(legend.position="bottom",legend.title=element_blank())+facet_wrap(~order, scales="free",labeller=labeller(order = qlabels))

#save summary figure

ggsave(file=sprintf("output/accum_alltree_all_div_pooled_indv%s.pdf",date), plot=accum_alltree_all_q, width=7, height=4)


############FYI#######

#neighborhood accumulation curves- # not in current paper- need 20+ minutes to run

#create neighborhood x species matrix for all trees
y<-pivot_wider(st_tree_next, id_cols=SPP, names_from = holc_id, values_from = c(abundance),values_fn = list(abundance = sum),values_fill = 0)
y<-column_to_rownames(y, var = "SPP")

#ALLTREES FIGURE- rarefaction at all q levels

Holc_neigh_est_S<-estimateD(y)

out.all.neigh<- iNEXT(y, q=c(0,1,2),datatype="abundance")


out.all.neigh$DataInfo # showing basic data information.
out.all.neigh$AsyEst # showing asymptotic diversity estimates.
out.all.neigh$iNextEst # showing diversity estimates with rarefied and extrapolated.


#summary rarefaction figure

#36 neighborhoods- need vector to match

neigh_holc_grade<-c("A","A","A","A","B","B","B","B","B","B","B","B","B","B","B","B","C","C","C","C","C","C","C","C","C","C","C","C","C", "D","D","D","D","D","D","D")
neigh_holc_grade<-factor(neigh_holc_grade)
levels(neigh_holc_grade)<- c("#92BC6B" ,"#92C7C9", "#E7DC6B", "#E47D67")
neigh_col<-as.character(neigh_holc_grade)

accum_alltree_all_q_neigh<-ggiNEXT(out.all.neigh, type=1, facet.var="order",se="false") + theme_bw(base_size=10)+ theme_bw(base_size=10)+ xlim(c(0,6000))+scale_colour_manual(values = neigh_col,name="HOLC grade")+scale_fill_manual(values = neigh_col, name="HOLC grade")+
  labs(x = "Number of individual trees sampled", y = "Tree species diversity", colour = "HOLC grade", fill="HOLC grade")+
  theme(legend.position="right",legend.title=element_blank())+facet_wrap(~order, scales="free",labeller=labeller(order = qlabels))

ggiNEXT(out.all.neigh, facet.var="site", color.var="order")
