library(codyn)
library(tidyverse)
library(ggrepel)
library(vegan)
library(gridExtra)
library(grid)

#export all as TIFF

##do evenness in codyn
#test.
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
  mutate(size=ifelse(DBH<=5, "II. Small trees", ifelse(DBH>=20, "I. Large trees", "drop")))%>%
  filter(size!="drop")

#total by each grade
abund<-clean2%>%
  group_by(holc_grade, size, SPP2)%>%
  summarize(abund=sum(present))%>%
  mutate(rank=rank(-abund, ties.method = "first"))%>%
  mutate(name=ifelse(rank<6, SPP2, ""))%>%
  mutate(plotname=ifelse(name=="Acer rubrum", "A. rubrum", ifelse(name=="Liriodendron tulipifera", "L. tulipifera", ifelse(name=="Ulmus americana", "U. americana", ifelse(name=="Zelkova serrata", "Z. serrata", ifelse(name=="Cornus florida", "C. florida", ifelse(name=="Cupressocyparis leylandii", "C. lylandii", ifelse(name=="Lagerstroemia indica", "L. indica", ifelse(name=="Prunus spp.", "Prunus sp.", ifelse(name=="Acer saccharinum", "A. saccharinum", ifelse(name=="Tilia cordata", "T. cordata", ifelse(name=="Cercis canadensis", "C. canadensis", ifelse(name=="Gleditsia triacanthos", "G. triacanthos", ifelse(name=="Pyrus calleryana", "P. calleryana", ifelse(name=="Quercus palustris", "Q. palustris", ifelse(name=="Platanus x acerifolia", "P. x acerifolia", ""))))))))))))))))


rac<-
  ggplot(data=abund, aes(x=rank, y=abund, label=plotname))+
  geom_point(size=1)+
  geom_line()+
  facet_grid(holc_grade~size)+
  xlab("Rank")+
  ylab("Number of Trees")+
  theme(legend.position = "none")+
  geom_text_repel(max.overlaps = 50, size=2.5, nudge_x=15)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

g <- ggplot_gtable(ggplot_build(rac))
strip_both <- which(grepl('strip-', g$layout$name))
fills <- c("gray","gray","olivedrab3", "lightcyan3", "khaki", "indianred3")
k <- 1

for (i in strip_both) {
  j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}
grid.draw(g)


ggsave("RAC_figure.tiff", units="px", width=5200, height=5600, dpi=800)

tiff(filename="C:\\Users\\mavolio2\\Dropbox\\BES Research\\redlining project\\FigRAC.tif",height=5600,width=5200,units="px",res=800,compression="lzw")


# #average by each grade
# mabund<-clean%>%
#   group_by(holc_grade, holc_id, SPP2)%>%
#   summarize(abund=sum(present))%>%
#   group_by(holc_grade, SPP2)%>%
#   summarize(mabund=mean(abund))%>%
#   mutate(rank=rank(-mabund, ties.method = "first"))%>%
#   mutate(name=ifelse(rank<6, SPP2, ""))%>%
#   mutate(sp=ifelse(mabund>70, SPP2, ""))
# 
# 
# ggplot(data=mabund, aes(x=rank, y=mabund, label=name, color=sp))+
#   geom_point()+
#   facet_wrap(~holc_grade)+
#   geom_text_repel()+
#   theme(legend.position = "none")

##doing beta diversity stuff
nbid<-clean2%>%
  right_join(numhex)%>%
  group_by(holc_grade, holc_id, SPP2)%>%
  summarize(abund=sum(present))

nbid_n<-nbid%>%
  select(holc_grade, holc_id) %>% 
  unique()%>%
  group_by(holc_grade)%>%
  summarise(n=length(holc_id))

nbid_large<-clean2%>%
  right_join(numhex)%>%
  filter(size=="I. Large trees")%>%
  group_by(holc_grade, holc_id, SPP2)%>%
  summarize(abund=sum(present))

nbid_small<-clean2%>%
  right_join(numhex)%>%
  filter(size=="II. Small trees")%>%
  group_by(holc_grade, holc_id, SPP2)%>%
  summarize(abund=sum(present))

##step 1 do an NMDS
#change the numbers and dataset for each

nbid_wides<-nbid_small%>%
  spread(SPP2, abund, fill=0)

#do the NMDS  
plots<-nbid_wides[,1:2]
mdss<-metaMDS(nbid_wides[,3:200], distance="morisita", autotransform=FALSE, shrink=FALSE) 
mdss #stress all 0.09; stress large 0.11; stress small 0.14

# are there differences in communities by landuse
adonis(nbid_wides[,3:200]~as.factor(holc_grade), nbid_wides)
#not sig diff communities by HOLC_Grade; large trees, big sig diff holc grade# sig diff small trees

#test whether Landuse have differences in dispersion
dist<-vegdist(nbid_wides[,3:200])
betadisp<-betadisper(dist,nbid_wides$holc_grade,type="centroid")
betadisp
permutest(betadisp)
#not sig diff dispersion by holc_grade, no sig diff large trees; #no sig diff large trees

scores <- data.frame(scores(mdss, display="sites"))  # Extracts NMDS scores for each block
scores2<- cbind(plots, scores)%>%
  mutate(tree="II. Small trees")# binds the NMDS scores landuse plot info

##large trees
nbid_widel<-nbid_large%>%
  spread(SPP2, abund, fill=0)

#do the NMDS  
plotsl<-nbid_widel[,1:2]
mdsl<-metaMDS(nbid_widel[,3:125], autotransform=FALSE, shrink=FALSE) 
mdsl #stress all 0.09; stress large 0.11; stress small 0.14

# are there differences in communities by landuse
adonis(nbid_widel[,3:125]~as.factor(holc_grade), nbid_widel)

#not sig diff communities by HOLC_Grade; large trees, big sig diff holc grade# sig diff small trees

#test whether Landuse have differences in dispersion
dist<-vegdist(nbid_widel[,3:125])
betadisp<-betadisper(dist,nbid_widel$holc_grade,type="centroid")
betadisp
permutest(betadisp)
#not sig diff dispersion by holc_grade, no sig diff large trees; #no sig diff large trees

scores3 <- data.frame(scores(mdsl, display="sites"))  # Extracts NMDS scores for each block
scores4<- cbind(plotsl, scores3)%>%
  mutate(tree="I. Large trees")%>%
  bind_rows(scores2)

stress<-data.frame(tree=c("I. Large trees","II. Small trees"),
                      toplot=c("Stress = 0.11", "Stress = 0.14"))

##plotting this
#put figure legend on the bottom
nmds<-
  ggplot(scores4, aes(x=NMDS1, y=NMDS2))+
  geom_point(size=3, aes(color=holc_grade))+
  facet_wrap(~tree)+
  scale_shape_manual(name="Year", values=c(15,17,10,19))+
  scale_color_manual(name="HOLC\nGrade", values = c("olivedrab3", "lightcyan3", "khaki", "indianred3"))+
  xlab("NMDS Axis 1")+
  ylab("NMDS Axis 2")+
  stat_ellipse(size=1, aes(color=holc_grade))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_text(data=stress, aes(x=-1.5, y = -1.5, label = toplot), size=3)


#ggsave("NMDS_Size.pdf", nmds, units="mm", width=150, dpi=300)


###NMDS all NBs
nbid<-clean%>%
  right_join(numhex)%>%
  group_by(holc_grade, holc_id, SPP2)%>%
  summarize(abund=sum(present))

nbid_wide<-nbid%>%
  spread(SPP2, abund, fill=0)

#do the NMDS  - results are the same for bray, morisita and horn distance matrices.
plots<-nbid_wide[,1:2]
mds<-metaMDS(nbid_wide[,3:230], autotransform=FALSE, shrink=FALSE, trymax=100) 
mds #stress all 0.09; stress large 0.11; stress small 0.14

# are there differences in communities by landuse
adonis(nbid_wide[,3:230]~as.factor(holc_grade), nbid_wide)
#not sig diff communities by HOLC_Grade; large trees, big sig diff holc grade# sig diff small trees

#test whether Landuse have differences in dispersion
dist<-vegdist(nbid_wide[,3:230])
betadisp<-betadisper(dist,nbid_wide$holc_grade,type="centroid")
betadisp
permutest(betadisp)
#not sig diff dispersion by holc_grade, no sig diff large trees; #no sig diff large trees

scores5 <- data.frame(scores(mds, display="sites"))  # Extracts NMDS scores for each block
scores6<- cbind(plots, scores5)

#adjusting p-values for all
p<-c(0.912, 0.059, 0.015, 0.484, 0.001, 0.464)
p.adjust(p, "BH")


##plotting this
nmds_all<-
  ggplot(scores6, aes(x=NMDS1, y=NMDS2))+
  geom_point(size=3, aes(color=holc_grade))+
  scale_shape_manual(name="Year", values=c(15,17,10,19))+
  scale_color_manual(name="HOLC\nGrade", values = c("olivedrab3", "lightcyan3", "khaki", "indianred3"))+
  xlab("NMDS Axis 1")+
  ylab("NMDS Axis 2")+
  stat_ellipse(size=1, aes(color=holc_grade))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  annotate("text", label="Stress = 0.09", x=-1.8, y=-0.8, size=3)


ggsave("NMDS_All_AppendixFig.jpeg", nmds_all, units="mm", width=120, dpi=300)



#doing rac_difference for trees
#small

racdiffs<-RAC_difference(df=nbid_small, species.var = "SPP2", abundance.var = "abund", replicate.var = "holc_id", treatment.var = "holc_grade")

racdiff_subs<-racdiffs%>%
  filter(holc_grade==holc_grade2)%>%
  mutate(tree="II. Small trees")

#do holc grades differ in species differences?
summary(aov(species_diff~holc_grade, data=racdiff_subs))
#yes, sig diff in species differences
#no, no sp. diff for large trees
#sig sp diff small trees
TukeyHSD(aov(species_diff~holc_grade, data=racdiff_subs))

#do holc grades differ in rank differences?
summary(aov(rank_diff~holc_grade, data=racdiff_subs))
#yes, sig diff in rank differences, #yes sig RAC differences for large trees
#sig diff small trees
TukeyHSD(aov(rank_diff~holc_grade, data=racdiff_subs))

#do holc grades differ in richness differences?
summary(aov(abs(richness_diff)~holc_grade, data=racdiff_subs))
#no, sig diff in rich differences; no diff large trees, #no diff richness small trees

#do holc grades differ in evenness differences?
summary(aov(abs(evenness_diff)~holc_grade, data=racdiff_subs))
#no,  sig diff in even differences, #sig diff evenness large; no diff snall
TukeyHSD(aov(abs(evenness_diff)~holc_grade, data=racdiff_subs))

##large

racdiffl<-RAC_difference(df=nbid_large, species.var = "SPP2", abundance.var = "abund", replicate.var = "holc_id", treatment.var = "holc_grade")

racdiff_subl<-racdiffl%>%
  filter(holc_grade==holc_grade2)%>%
  mutate(tree="I. Large trees")

#do holc grades differ in species differences?
summary(aov(species_diff~holc_grade, data=racdiff_subl))
#yes, sig diff in species differences
#no, no sp. diff for large trees
#sig sp diff small trees
TukeyHSD(aov(species_diff~holc_grade, data=racdiff_subl))

#do holc grades differ in rank differences?
summary(aov(rank_diff~holc_grade, data=racdiff_subl))
#yes, sig diff in rank differences, #yes sig RAC differences for large trees
#sig diff small trees
TukeyHSD(aov(rank_diff~holc_grade, data=racdiff_subl))

#do holc grades differ in richness differences?
summary(aov(abs(richness_diff)~holc_grade, data=racdiff_subl))
#no, sig diff in rich differences; no diff large trees, #no diff richness small trees

#do holc grades differ in evenness differences?
summary(aov(abs(evenness_diff)~holc_grade, data=racdiff_subl))
#no,  sig diff in even differences, #sig diff evenness large; no diff snall
TukeyHSD(aov(abs(evenness_diff)~holc_grade, data=racdiff_subl))

###doing stats all at once.
racdiff_sub<-racdiff_subs%>%
  bind_rows(racdiff_subl)%>%
  gather(measure, value, richness_diff:species_diff)%>%
  group_by(measure, tree)%>%
  summarize(pval=summary(aov(value~holc_grade))[[1]][["Pr(>F)"]][1])%>%
  mutate(pad=p.adjust(pval, method = "BH"))

toplot<-racdiff_subs%>%
  bind_rows(racdiff_subl)%>%
  select(-richness_diff, -evenness_diff)%>%
  gather(measure, value, rank_diff:species_diff)%>%
  group_by(tree, holc_grade, measure)%>%
  summarize(mean=mean(abs(value)), sd=sd(value), n=length(value))%>%
  mutate(se=sd/sqrt(n))%>%
  mutate(text=ifelse(measure=="species_diff"&holc_grade=="A"&tree=="Small"|measure=="rank_diff"&holc_grade=="A"&tree=="I. Large trees", "xy", 
              ifelse(measure=="rank_diff"&holc_grade=="D"&tree=="II. Small trees"|measure=="species_diff"&holc_grade=="D"&tree=="II. Small trees"|measure=="rank_diff"&holc_grade=="B"&tree=="I. Large trees", "y",ifelse(tree=="I. Large trees"&measure=="species_diff", " ", "x"))))


###graphing this
parameter<-c(rank_diff="Species\nReordering", species_diff="Species\nTurnover")

betadiv<-
ggplot(data=toplot, aes(x=holc_grade, y=mean, fill=holc_grade, label=text))+
  facet_grid(measure~tree, labeller=labeller(measure = parameter))+
  geom_bar(stat="identity")+
  scale_fill_manual(name="HOLC Grade", values = c("olivedrab3", "lightcyan3", "khaki", "indianred3"))+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2, size=1)+
  geom_text(aes(y=(mean+se)+0.1))+
  xlab("HOLC Grade")+
  ylab("Mean")+
  theme(legend.position = "none")

fig4<-grid.arrange(nmds, betadiv, ncol=1)

ggsave("Fig4.tiff", fig4, units="px", width=5200, height=5600, dpi=800)

tiff(filename="C:\\Users\\mavolio2\\Dropbox\\BES Research\\redlining project\\FigBeta.tif",height=5600,width=5200,units="px",res=800,compression="lzw")

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
  mutate(se=sd/sqrt(n))
# %>%
#   mutate(text=ifelse(measure=="richness_diff"|measure=="species_diff", "", 
#               ifelse(measure=="evenness_diff"&comp=="B-D"|measure=="rank_diff"&comp=="B-C"|measure=="evenness_diff"&comp=="C-D", "A",
#               ifelse(measure=="evenness_diff"&comp=="A-B"|measure=="evenness_diff"&comp=="A-C"|measure=="evenness_diff"&comp=="B-C"|measure=="rank_diff"&comp=="C-D", "B", "AB"))))

ggplot(data=toplot2, aes(x=comp, y=mean))+#, label=text
  geom_bar(stat="identity")+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2)+
  facet_wrap(~measure)+
  ggtitle("Large trees (>20 DBH")
  #geom_text(aes(y=(mean+se)+0.05))

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




#####
#####NOT DOING THIS
#####
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
