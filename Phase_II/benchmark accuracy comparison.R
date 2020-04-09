algorithm=c("UFL_Fuzzy_ART","K-prototype","K-medoids","Ensemble","LCM","Improved_LCM","Proposed")
algorithm=factor(algorithm,levels=c("UFL_Fuzzy_ART","K-prototype","K-medoids","Ensemble","LCM","Improved_LCM","Proposed"))
accuracy=c(81.5,80.0,76.5,81.3,78.9,83.0,83.3)
benchmark=data.frame(algorithm,accuracy)
figure4=ggplot(data=benchmark,aes(x=algorithm,y=accuracy))+geom_bar(stat="identity",width=0.5)+
  scale_y_continuous()+coord_cartesian(ylim=c(70,86))+
  theme(axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  guides(shape=guide_legend(override.aes=list(color=override.color)))+
geom_text(aes(label=paste0(accuracy, '%')), vjust=-0.3, size=3.5)
figure4
tiff('/Users/yinlinfu/Dropbox/ESM Phase II/paper/figure/figure4.tiff', units="in", width=7, height=4, res=300)
figure4
dev.off()