library(codyn)
library(tidyverse)
library(ggrepel)
library(vegan)

##do evenness in codyn

theme_set(theme_bw(12))

setwd("C:\\Users\\mavolio2\\Dropbox\\BES Research\\redlining project")

data<-read.csv("st_tree_inHEX2021_03_15.csv")

clean<-data%>%
  mutate(SPP2=ifelse(SPP=='Acer spp.'|SPP=='Carya spp.'|SPP=='Cornus spp.'|SPP=='Cornus x'|SPP=='Ficus spp.'|SPP=='Fraxinus spp.'|SPP=='Hydrangea spp.'|SPP=='Ilex spp.'|SPP=='Ilex x'|SPP=='Juniperus spp.'|SPP=='Magnolia spp.'|SPP=='Magnolia x'|SPP=='Photinia spp.'|SPP=='Picea spp.'|SPP=='Populus spp.'|SPP=='Quercus spp.'|SPP=='Quercus x'|SPP=='Rhus spp.'|SPP=='Salix spp.'|SPP=='Ulmus spp.', "unknown tree", SPP))%>%
  #mutate(SPP2=SPP)%>%
  filter(SPP2!='Vacant Site'&SPP2!='Vacant Potential'&SPP2!='Stump'&SPP2!='Vacant Site Not Suitable'&SPP2!='NA'&SPP2!='unknown shrub'&SPP2!="Z Add 01"&SPP2!=" "&SPP2!="Dead")%>%
  filter(LOC_TYPE!="Park")%>%
  filter(CONDITION!="Dead"&CONDITION!="Stump"&CONDITION!="Sprout")%>%
  select(UniqueID, Genus, SPP2, DBH, hex_id, holc_id, holc_grade)%>%
  mutate(present=1)%>%
  filter(SPP2!="unknown tree")

#what does condition absent mean?
#Making RACs for each holc grade

####no re-doing this but with same number of hexes per NB, now each NB has at least 5 hexes
numhex<-clean%>%
  select(holc_grade, holc_id, hex_id)%>%
  unique()%>%
  group_by(holc_grade, holc_id)%>%
  summarize(n=length(hex_id))%>%
  mutate(drop=ifelse(holc_id %in% c("B16", "B6", "B7"), 1, 0))%>%
  filter(drop==0)


clean2<-clean%>%
  right_join(numhex)%>%
  mutate(size=ifelse(DBH<=5, "S", ifelse(DBH>=25, "L", "drop")))%>%
  filter(size!="drop")

#total by each grade
abund<-clean2%>%
  group_by(holc_grade, size, SPP2)%>%
  summarize(abund=sum(present))%>%
  mutate(rank=rank(-abund, ties.method = "first"))%>%
  mutate(name=ifelse(rank<6, SPP2, ""))


ggplot(data=abund, aes(x=rank, y=abund, label=name))+
  geom_point()+
  facet_grid(size~holc_grade, scales="free")+
  theme(legend.position = "none")+
  geom_text_repel(max.overlaps = 100)

#average by each grade
mabund<-clean%>%
  group_by(holc_grade, holc_id, SPP2)%>%
  summarize(abund=sum(present))%>%
  group_by(holc_grade, SPP2)%>%
  summarize(mabund=mean(abund))%>%
  mutate(rank=rank(-mabund, ties.method = "first"))%>%
  mutate(name=ifelse(rank<6, SPP2, ""))%>%
  mutate(sp=ifelse(mabund>70, SPP2, ""))


ggplot(data=mabund, aes(x=rank, y=mabund, label=name, color=sp))+
  geom_point()+
  facet_wrap(~holc_grade)+
  geom_text_repel()+
  theme(legend.position = "none")

##doing beta diversity stuff
nbid<-clean%>%
  right_join(numhex)%>%
  group_by(holc_grade, holc_id, SPP2)%>%
  summarize(abund=sum(present))

##step 1 do an NMDS

nbid_wide<-nbid%>%
  spread(SPP2, abund, fill=0)

#do the NMDS  
plots<-nbid_wide[,1:2]
mds<-metaMDS(nbid_wide[,3:230], autotransform=FALSE, shrink=FALSE) 
mds #stress 0.09

# are there differences in communities by landuse
adonis(nbid_wide[,3:230]~as.factor(holc_grade), nbid_wide)
#not sig diff communities by HOLC_Grade

#test whether Landuse have differences in dispersion
dist<-vegdist(nbid_wide[,3:230])
betadisp<-betadisper(dist,nbid_wide$holc_grade,type="centroid")
betadisp
permutest(betadisp)
#not sig diff dispersion by holc_grade

scores <- data.frame(scores(mds, display="sites"))  # Extracts NMDS scores for each block
scores2<- cbind(plots, scores) # binds the NMDS scores landuse plot info

##plotting this
ggplot(scores2, aes(x=NMDS1, y=NMDS2, color=holc_grade))+
  geom_point(size=4)+
  scale_shape_manual(name="Year", values=c(15,17,10,19))+
  scale_color_manual(name="HOLC Grade", values = c("olivedrab3", "lightcyan3", "khaki", "indianred3"))+
  xlab("NMDS Axis 1")+
  ylab("NMDS Axis 2")+
  annotate("text", x=1.5, y=-1, label="stress = 0.08", size=4)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())  

#doing rac_difference
racdiff<-RAC_difference(df=nbid, species.var = "SPP2", abundance.var = "abund", replicate.var = "holc_id", treatment.var = "holc_grade")

racdiff_sub<-racdiff%>%
  filter(holc_grade==holc_grade2)

#do holc grades differ in species differences?
summary(aov(species_diff~holc_grade, data=racdiff_sub))
#yes, sig diff in species differences
TukeyHSD(aov(species_diff~holc_grade, data=racdiff_sub))

#do holc grades differ in rank differences?
summary(aov(rank_diff~holc_grade, data=racdiff_sub))
#yes, sig diff in rank differences
TukeyHSD(aov(rank_diff~holc_grade, data=racdiff_sub))

#do holc grades differ in richness differences?
summary(aov(abs(richness_diff)~holc_grade, data=racdiff_sub))
#no, sig diff in rich differences


#do holc grades differ in evenness differences?
summary(aov(abs(evenness_diff)~holc_grade, data=racdiff_sub))
#no,  sig diff in even differences


toplot<-racdiff_sub%>%
  gather(measure, value, richness_diff:species_diff)%>%
  group_by(holc_grade, measure)%>%
  summarize(mean=mean(abs(value)), sd=sd(value), n=length(value))%>%
  mutate(se=sd/sqrt(n))%>%
  mutate(text=ifelse(measure=="richness_diff"|measure=="evenness_diff", "",
              ifelse(measure=="rank_diff"&holc_grade=="A"|measure=="species_diff"&holc_grade=="A", "AB", 
       ifelse(measure=="rank_diff"&holc_grade=="B"|measure=="species_diff"&holc_grade=="C", "A", "B"))))

ggplot(data=toplot, aes(x=holc_grade, y=mean, fill=holc_grade, label=text))+
  geom_bar(stat="identity")+
  scale_fill_manual(name="HOLC Grade", values = c("olivedrab3", "lightcyan3", "khaki", "indianred3"))+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2)+
  facet_wrap(~measure)+
  geom_text(aes(y=(mean+se)+0.1))

####
##doing this to compare neighborhoods
racdiff_comp<-racdiff%>%
  filter(holc_grade!=holc_grade2)%>%
  mutate(comp=paste(holc_grade, holc_grade2, sep="-"))

#do holc grades differ in species differences?
summary(aov(species_diff~comp, data=racdiff_comp))
#no, no diff in how NB diff in species differences

#do holc grades differ in rank differences?
summary(aov(rank_diff~comp, data=racdiff_comp))
#yes, sig diff in rank differences
TukeyHSD(aov(rank_diff~comp, data=racdiff_comp))

#do holc grades differ in richness differences?
summary(aov(abs(richness_diff)~comp, data=racdiff_comp))
#no, sig diff in rich differences

#do holc grades differ in evenness differences?
summary(aov(abs(evenness_diff)~comp, data=racdiff_comp))
#yes,  sig diff in species differences
TukeyHSD(aov(abs(evenness_diff)~comp, data=racdiff_comp))

toplot2<-racdiff_comp%>%
  gather(measure, value, richness_diff:species_diff)%>%
  group_by(comp, measure)%>%
  summarize(mean=mean(abs(value)), sd=sd(value), n=length(value))%>%
  mutate(se=sd/sqrt(n))%>%
  mutate(text=ifelse(measure=="richness_diff"|measure=="species_diff", "", 
              ifelse(measure=="evenness_diff"&comp=="B-D"|measure=="rank_diff"&comp=="B-C"|measure=="evenness_diff"&comp=="C-D", "A",
              ifelse(measure=="evenness_diff"&comp=="A-B"|measure=="evenness_diff"&comp=="A-C"|measure=="evenness_diff"&comp=="B-C"|measure=="rank_diff"&comp=="C-D", "B", "AB"))))

ggplot(data=toplot2, aes(x=comp, y=mean, label=text))+
  geom_bar(stat="identity")+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2)+
  facet_wrap(~measure)+
  geom_text(aes(y=(mean+se)+0.05))

####no re-doing this but with same number of hexes per NB - need to select 5 hexes randomly from each nb

#subsetting so only hexes with 8 or more trees

subset<-clean%>%
  right_join(numhex)%>%
  group_by(holc_grade, holc_id, hex_id)%>%
  summarize(abund=sum(present))%>%
  filter(abund>7)

#doing rac_difference
racdiff<-RAC_difference(df=nbid, species.var = "SPP2", abundance.var = "abund", replicate.var = "holc_id", treatment.var = "holc_grade")

racdiff_sub<-racdiff%>%
  filter(holc_grade==holc_grade2)

#do holc grades differ in species differences?
summary(aov(species_diff~holc_grade, data=racdiff_sub))
#yes, sig diff in species differences
TukeyHSD(aov(species_diff~holc_grade, data=racdiff_sub))

#do holc grades differ in rank differences?
summary(aov(rank_diff~holc_grade, data=racdiff_sub))
#yes, sig diff in rank differences
TukeyHSD(aov(rank_diff~holc_grade, data=racdiff_sub))

#do holc grades differ in richness differences?
summary(aov(abs(richness_diff)~holc_grade, data=racdiff_sub))
#no, sig diff in rich differences


#do holc grades differ in evenness differences?
summary(aov(abs(evenness_diff)~holc_grade, data=racdiff_sub))
#no,  sig diff in even differences


toplot<-racdiff_sub%>%
  gather(measure, value, richness_diff:species_diff)%>%
  group_by(holc_grade, measure)%>%
  summarize(mean=mean(abs(value)), sd=sd(value), n=length(value))%>%
  mutate(se=sd/sqrt(n))%>%
  mutate(text=ifelse(measure=="richness_diff"|measure=="evenness_diff", "",
                     ifelse(measure=="rank_diff"&holc_grade=="A"|measure=="species_diff"&holc_grade=="A", "AB", 
                            ifelse(measure=="rank_diff"&holc_grade=="B"|measure=="species_diff"&holc_grade=="C", "A", "B"))))

ggplot(data=toplot, aes(x=holc_grade, y=mean, fill=holc_grade, label=text))+
  geom_bar(stat="identity")+
  scale_fill_manual(name="HOLC Grade", values = c("olivedrab3", "lightcyan3", "khaki", "indianred3"))+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2)+
  facet_wrap(~measure)+
  geom_text(aes(y=(mean+se)+0.1))

####
##doing this to compare neighborhoods
racdiff_comp<-racdiff%>%
  filter(holc_grade!=holc_grade2)%>%
  mutate(comp=paste(holc_grade, holc_grade2, sep="-"))

#do holc grades differ in species differences?
summary(aov(species_diff~comp, data=racdiff_comp))
#no, no diff in how NB diff in species differences

#do holc grades differ in rank differences?
summary(aov(rank_diff~comp, data=racdiff_comp))
#yes, sig diff in rank differences
TukeyHSD(aov(rank_diff~comp, data=racdiff_comp))

#do holc grades differ in richness differences?
summary(aov(abs(richness_diff)~comp, data=racdiff_comp))
#no, sig diff in rich differences

#do holc grades differ in evenness differences?
summary(aov(abs(evenness_diff)~comp, data=racdiff_comp))
#yes,  sig diff in species differences
TukeyHSD(aov(abs(evenness_diff)~comp, data=racdiff_comp))

toplot2<-racdiff_comp%>%
  gather(measure, value, richness_diff:species_diff)%>%
  group_by(comp, measure)%>%
  summarize(mean=mean(abs(value)), sd=sd(value), n=length(value))%>%
  mutate(se=sd/sqrt(n))%>%
  mutate(text=ifelse(measure=="richness_diff"|measure=="species_diff", "", 
                     ifelse(measure=="evenness_diff"&comp=="B-D"|measure=="rank_diff"&comp=="B-C"|measure=="evenness_diff"&comp=="C-D", "A",
                            ifelse(measure=="evenness_diff"&comp=="A-B"|measure=="evenness_diff"&comp=="A-C"|measure=="evenness_diff"&comp=="B-C"|measure=="rank_diff"&comp=="C-D", "B", "AB"))))

ggplot(data=toplot2, aes(x=comp, y=mean, label=text))+
  geom_bar(stat="identity")+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2)+
  facet_wrap(~measure)+
  geom_text(aes(y=(mean+se)+0.05))

####comparing hexes within a NB

#subsetting so only hexes with 8 or more trees. and for 5 hexes in each NB. not randomly picking hexes!

tosubset<-clean%>%
  right_join(numhex)%>%
  group_by(holc_grade, holc_id, hex_id)%>%
  summarize(abund=sum(present))%>%
  filter(abund>7)%>%
  select(-abund)

numhex_2<-tosubset%>%
  select(holc_grade, holc_id, hex_id)%>%
  unique()%>%
  group_by(holc_grade, holc_id)%>%
  mutate(num=1:n())%>%
  filter(num<6)%>%
  select(-num)


subset<-clean%>%
  right_join(numhex_2)%>%
  group_by(holc_grade, holc_id, hex_id, SPP2)%>%
  summarize(abund=sum(present))

racdiff<-RAC_difference(df=subset, species.var = "SPP2", abundance.var = "abund", replicate.var = "hex_id", treatment.var = "holc_id")

racdiff_sub<-racdiff%>%
  filter(holc_id==holc_id2)%>%
  mutate(holc_grade=substr(holc_id, 1, 1))%>%
  gather(measure, value, richness_diff:species_diff)%>%
  group_by(holc_grade, holc_id, measure)%>%
  summarize(mean=mean(abs(value)))%>%
  spread(measure, mean)
  

#do holc grades differ in species differences?
summary(aov(species_diff~holc_grade, data=racdiff_sub))
#no, sig diff in species differences

#do holc grades differ in rank differences?
summary(aov(rank_diff~holc_grade, data=racdiff_sub))
#no, sig diff in rank differences

#do holc grades differ in richness differences?
summary(aov(abs(richness_diff)~holc_grade, data=racdiff_sub))

#do holc grades differ in evenness differences?
summary(aov(abs(evenness_diff)~holc_grade, data=racdiff_sub))
#no,  sig diff in even differences


toplot<-racdiff_sub%>%
  gather(measure, value, evenness_diff:species_diff)%>%
  group_by(holc_grade, measure)%>%
  summarize(mean=mean(abs(value)), sd=sd(abs(value)), n=length(value))%>%
  mutate(se=sd/sqrt(n))

ggplot(data=toplot, aes(x=holc_grade, y=mean, fill=holc_grade))+
  geom_bar(stat="identity")+
  scale_fill_manual(name="HOLC Grade", values = c("olivedrab3", "lightcyan3", "khaki", "indianred3"))+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2)+
  facet_wrap(~measure)
