library(codyn)
library(tidyverse)
library(vegan)
library(lme4)
library(lmerTest)

library(gridExtra)
library(gtable)
library(grid)

theme_set(theme_bw(12))


setwd("C:\\Users\\mavolio2\\Dropbox\\BES Research\\redlining project")

HOLC<-read.csv("itree_w_HOLC_grades.csv")


plot_info<-HOLC%>%
  select(Plot, Year, Landuse, holc_id, holc_grade)%>%
  unique()

# How many plots? ---------------------------------------------------------

#Table 1 - number of plots through time


nplots_yr<-HOLC%>%
  select(Plot, Year, holc_grade)%>%
  unique()%>%
  group_by(Year, holc_grade)%>%
  summarize(n=length(Plot))

#i loose plot 149 in year 2014

#drop U plots
holc<-HOLC%>%
  filter(holc_grade!="U")

# describing the community structure --------------------------------------

#analysis of richness, evenness, and abundance
abund<-holc%>%
  group_by(Plot, Year, holc_grade)%>%
  summarise(abundance=sum(Abundance))

rich_even<-community_structure(df=holc, 
                               time.var = "Year", 
                               replicate.var = "Plot",
                               abundance.var="Abundance",
                               metric="SimpsonEvenness")

combine<-rich_even%>%
  left_join(abund)%>%
  mutate(evenness=ifelse(SimpsonEvenness==Inf, NA, SimpsonEvenness),
         holc_grade=as.factor(holc_grade),
         year=as.factor(Year))

#repeated measures model of holc and year
#Evenness
fit_Evenness<- lmer(evenness ~ holc_grade*year + (1|Plot), data=combine)
anova(fit_Evenness)
#no diff evenness

#Richness
fit_richness<- lmer(richness ~ holc_grade*year + (1|Plot), data=combine)
anova(fit_richness)
ls_means(fit_richness, which="year", pairwise = TRUE)
#diff by year only

#Abundance
fit_abund<- lmer(abundance ~ holc_grade*year + (1|Plot), data=combine)
anova(fit_abund)
#no diff abundance


plot<-combine%>%
  select(-SimpsonEvenness, -Year)%>%
  gather(key="measure", value="values", richness, evenness, abundance)%>%
  group_by(holc_grade, year, measure)%>%
  summarize(mean=mean(values, na.rm=T), sd=sd(values, na.rm=T), n=length(values))%>%
  mutate(se=sd/sqrt(n))

ggplot(data=plot, aes(x=year, y=mean, group=holc_grade, color=holc_grade))+
  geom_point(size=3)+
  geom_line(size=1)+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2)+
  scale_color_manual(name="HOLC Grade", values = c("olivedrab3", "lightcyan3", "khaki", "indianred3"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_blank())+
  xlab("Year")+
  ylab("Mean Value")+
  facet_wrap(~measure, ncol=1, scales = "free")

# Multivariate analyses ---------------------------------------------------

#doing this to compare landuse and time
#getting total inds for each species by landuse and year
data2<-holc%>%
  filter(Species!="#N/A")%>%
  select(-Native, -Commonname)%>%
  group_by(Year, holc_grade, Species)%>%
  summarize(abund=sum(Abundance))%>%
  ungroup()%>%
  spread(Species, abund, fill=0)

#do the NMDS  
plots<-data2[,1:2]
mds<-metaMDS(data2[,3:80], autotransform=FALSE, shrink=FALSE) 
mds #stress 0.035

# are there differences in communities by landuse
adonis(data2[,3:80]~as.factor(holc_grade), data2)
#sig diff communities by HOLC_Grade

#test whether Landuse have differences in dispersion
dist<-vegdist(data2[,3:80])
betadisp<-betadisper(dist,data2$holc_grade,type="centroid")
betadisp
permutest(betadisp)
#sig diff dispersion by holc_grade

scores <- data.frame(scores(mds, display="sites"))  # Extracts NMDS scores for each block
scores2<- cbind(plots, scores) # binds the NMDS scores landuse plot info

##plotting this
ggplot(scores2, aes(x=NMDS1, y=NMDS2, color=holc_grade, shape=as.factor(Year)))+
  geom_point(size=6.5)+
  scale_shape_manual(name="Year", values=c(15,17,10,19))+
  scale_color_manual(name="HOLC Grade", values = c("olivedrab3", "lightcyan3", "khaki", "indianred3"))+
  xlab("NMDS Axis 1")+
  ylab("NMDS Axis 2")+
  annotate("text", x=1.5, y=-1, label="stress = 0.04", size=4)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

#THis is not working and we drop most plots D b/c they have outlier species.
# #doing this to compare landuse beta diversity across all year NOT time.
# 
# 
# #getting total average inds for each species by landuse and plot
# data3<-holc%>%
#   filter(Species!="#N/A")%>%
#   select(-Native, -Commonname)%>%
#   group_by(Plot, holc_grade, Species)%>%
#   summarize(abund=mean(Abundance))
# 
# #getting data for wide format
# data4<-data3%>%
#   filter(Plot!=3&Plot!=16&Plot!=30&Plot!=100&Plot!=145&Plot!=147&Plot!=170&Plot!=173&Plot!=181&Plot!=183)%>%#dropping plots that are very dissimlar and throw off NMDS
#   ungroup()%>%
#   spread(Species, abund, fill=0)
# 
# #do the NMDS  - not good cannot find a convergent solution
# plots2<-data4[,1:2]
# mds2<-metaMDS(data4[,3:74], trymax=100) 
# 
# scores3 <- data.frame(scores(mds2, display="sites"))  # Extracts NMDS scores for each block
# scores4<- cbind(plots2, scores3)%>%# binds the NMDS scores landuse plot info
#   group_by(holc_grade)%>%
#   summarize(mNMDS1=mean(NMDS1), mNMDS2=mean(NMDS2), sd1=sd(NMDS1), sd2=sd(NMDS2), n=length(NMDS1))%>%
#   mutate(se1=sd1/sqrt(n), se2=sd2/sqrt(n))
# 
# ggplot(scores4, aes(x=NMDS1, y=NMDS2, color=holc_grade))+
#   geom_point(size=3)+
#   scale_color_manual(name="HOLC Grade", values = c("olivedrab3", "lightcyan3", "khaki", "indianred3"))+
#   #geom_errorbar(aes(ymin=mNMDS2-sd2, ymax=mNMDS2+sd2), width =0.01)+
#   #geom_errorbarh(aes(xmin=mNMDS1-sd1, xmax=mNMDS1+sd1))+
#   xlab("NMDS Axis 1")+
#   ylab("NMDS Axis 2")+
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
# 
# 
# # are there differences in communities by landuse
# adonis(data4[,3:110]~as.factor(Landuse), data4) # yes difference by landuse
# 
# #test whether Landuse have differences in dispersion
# dist<-vegdist(data4[,3:110])
# betadisp<-betadisper(dist,data4$Landuse,type="centroid")
# betadisp
# permutest(betadisp) #yes difference by landuse
# 
# #get differences among plots
# #this shows forest has least dispersion of all plots, lowest beta diversity
# 
# mult_diff<-multivariate_difference(data3,
#                                    species.var = "Species", 
#                                    abundance.var = "abund",
#                                    replicate.var = "Plot",
#                                    treatment.var = "Landuse")


# Rank Change analyses ----------------------------------------------------
holc2<-holc%>%
  filter(Species!="#N/A")

Rachange<-RAC_change(holc2, 
                     time.var = "Year", 
                     species.var = "Species", 
                     abundance.var = "Abundance", 
                     replicate.var = "Plot")

Turnover<-turnover(holc2, 
                   time.var = "Year", 
                   species.var = "Species", 
                   abundance.var = "Abundance", 
                   replicate.var = "Plot", metric="total")

rac2<-Rachange%>%
  left_join(plot_info)%>%
  mutate(holc_grade=as.factor(holc_grade),
         year=as.factor(Year))

#doing repeated measures anova's on changes over time
fit_rank<- lmer(rank_change ~ holc_grade*year + (1|Plot), data=rac2)
anova(fit_rank)

fit_gain<- lmer(gains ~ holc_grade*year + (1|Plot), data=rac2)
anova(fit_gain)

fit_loss<- lmer(losses ~ holc_grade*year + (1|Plot), data=rac2)#get an error message need to understand
anova(fit_loss)


plotrac<-rac2%>%
  select(-Year, -Year2, -richness_change, -evenness_change, -holc_id)%>%
  gather(key="measure", value="values", rank_change, gains, losses)%>%
  group_by(holc_grade, year, measure)%>%
  summarize(mean=mean(values, na.rm=T), sd=sd(values, na.rm=T), n=length(values))%>%
  mutate(se=sd/sqrt(n))

ggplot(data=plotrac, aes(x=as.factor(year), y=mean, group=holc_grade, color=holc_grade))+
  geom_point(size=3)+
  geom_line(size=1)+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2)+
  scale_color_manual(name="HOLC Grade", values = c("olivedrab3", "lightcyan3", "khaki", "indianred3"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  xlab("Year")+
  ylab("Mean Value")+
  scale_x_discrete(labels=c("1999-2004", "2004-2009", "2009-2014"))+
  facet_wrap(~measure, ncol=1, scales = "free")


