#autocorrelation plot function script
autocorrelation.plot=function(Data,Time,Monkey){
  X=data.frame(Autocorrelation=NA,Months=NA,Monkey=NA)
  for(i in 1:length(Monkey)){
    time=Time[[i]]
    data=Data[[i]]
    #data=barcodetrackR::barcode_ggheatmap(data,printtable=TRUE,n_clones=500)
    time=time[-1]
    x=cor(data,method="spearman")
    x=x[-1,]
    x=diag(x)
    
    x=data.frame(Autocorrelation=x,Months=time,Monkey=Monkey[i])
    
    X=rbind(X,x)
  }
  print(ggplot(X[-1,],aes(x=Months,y=Autocorrelation,group=Monkey,color=Monkey))+
          geom_line(size=1)+
          geom_point(aes(shape=Monkey),size=3)+scale_y_continuous("Autocorrelation (Spearman)",limits=c(-0.1,1))+
          scale_color_manual(values=c(RColorBrewer::brewer.pal(3, "Set1"),rep("black",5)))+
          scale_x_continuous("Months post transplant", breaks=seq(5,25,by=5))+
          scale_shape_manual(values=c(19,19,19,19,17,15,18,0))+
          theme_classic() +
          theme(
            axis.text.x = element_text(size = 10, color = "black"),
            axis.text.y = element_text(size = 10, color = "black"),
            axis.title.x = element_text(size = 12, color = "black", margin = margin(t = 10, r = 0, b = 0, l = 0)),
            axis.title.y = element_text(size = 12, color = "black", margin = margin(t = 0, r = 10, b = 0, l = 0)),
            legend.text = element_text(size = 12, color = "black"),
            legend.title = element_text(size = 12, color = "black", face = "bold"),
          ))
}