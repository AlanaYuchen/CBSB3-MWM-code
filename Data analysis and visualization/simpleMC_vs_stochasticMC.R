data1<-read.csv('/Users/chenghui/Documents/CBSB3/ICA/Water_Maze_Modelling/data/best_simpleMC.csv',header = F)
data2<-read.csv('/Users/chenghui/Documents/CBSB3/ICA/Water_Maze_Modelling/data/best_S1T12.csv',header=F)
header=c('Group','Day','Trial','Wmult','sigma_pc','sigma_ac','Vdecay','ac_const','beta','etdecay','alpha','gamma','noise','wpun','gof','latency','dist','quadrant4','quadrant2','wall_zone',
         'speed_std','mean_angle','time_step')
colnames(data1)<-header
colnames(data2)<-header
data1$estimation<-'Simple MC'
data2$estimation<-'Iterative MC'
data=rbind(data1,data2)
data$day_trial=(data$Day-1)*4+data$Trial

library(ggplot2)
library(gridExtra)
library(emmeans)

plots=list()
count=0
for(j in 1:4){
for(i in c(1:6)){
  count=count+1
  datas=data[which(data$Group==j),]
  datass=datas[,c(i+15,24:25)]
  model=summary(lm(datass[,1]~datass$estimation+datass$day_trial))
  r=model$r.squared
  estimation.p=model$coefficients[2,4]
  trial.p=model$coefficients[3,4]
  #ph=emmeans(model,~estimation+day_trial)
  #pairs(ph)
  plots[[count]]=ggplot(datass,aes_string(x='day_trial',y=header[i+15],color='estimation'))+
           geom_line()+
    geom_point()+
    theme_classic()+
    labs(title = paste('Group',j,'R square',round(r,3),'est.p',round(estimation.p,3),'trial.p',round(trial.p,3)))
  
}}
library(ggpubr)
ggsave('/Users/chenghui/Documents/CBSB3/ICA/Water_Maze_Modelling/data/PM_plot_simple_vs_stochasticMC.png',
       plot=plotss)
plotss=ggarrange(plotlist = plots,ncol = 2,nrow = 2)


for(m in 0:3){
plotss1=grid.arrange(plots[[m*7+1]],plots[[m*7+2]],nrow=2,ncol=1)
ggsave(paste('/Users/chenghui/Documents/CBSB3/ICA/Water_Maze_Modelling/data/',m*4+1,'.png'),
       plot=plotss1)
plotss2=grid.arrange(plots[[m*7+3]],plots[[m*7+4]],nrow=2,ncol=1)
ggsave(paste('/Users/chenghui/Documents/CBSB3/ICA/Water_Maze_Modelling/data/',m*4+2,'.png'),
       plot=plotss2)
plotss3=grid.arrange(plots[[m*7+5]],plots[[m*7+6]],nrow=2,ncol=1)
ggsave(paste('/Users/chenghui/Documents/CBSB3/ICA/Water_Maze_Modelling/data/',m*4+3,'.png'),
       plot=plotss3)
plotss4=grid.arrange(plots[[m*7+7]],nrow=2,ncol=1)
ggsave(paste('/Users/chenghui/Documents/CBSB3/ICA/Water_Maze_Modelling/data/',m*4+4,'.png'),
       plot=plotss4)
}

for(m in 0:3){
  plotss1=grid.arrange(plots[[m*7+1]],plots[[m*7+2]],plots[[m*7+3]],plots[[m*7+4]],nrow=2,ncol=2)
  ggsave(paste('/Users/chenghui/Documents/CBSB3/ICA/Water_Maze_Modelling/data/',m*2+1,'.png'),
         plot=plotss1)
  plotss4=grid.arrange(plots[[m*7+5]],plots[[m*7+6]],plots[[m*7+7]])
  ggsave(paste('/Users/chenghui/Documents/CBSB3/ICA/Water_Maze_Modelling/data/',m*2+2,'.png'),
         plot=plotss4)
}

marrangeGrob(grobs = plots,ncol = 2,nrow = 3)%>%
ggsave(paste('/Users/chenghui/Documents/CBSB3/ICA/Water_Maze_Modelling/data/simplevsiter6.pdf'),
       plot=.,width = 12,height = 12,device = 'pdf',dpi = 300)

plots=list()
count=0
for(j in 1:4){
  for(i in 1:8){
    count=count+1
    datas=data[which(data$Group==j),]
    datass=datas[,c(i+15,24:25)]
    for(m in 1: length(datass$estimation)){
      datass[m,1]=Re(datass[m,1])
    }
    datass[,1]=as.numeric(datass[,1])
    model=summary(lm(datass[,1]~datass$estimation+datass$day_trial))
    r=model$r.squared
    estimation.p=model$coefficients[2,4]
    trial.p=model$coefficients[3,4]
    plots[[count]]=ggplot(datass,aes_string(x='day_trial',y=header[i+15],color='estimation'))+
      geom_line()+
      geom_point()+
    theme_classic()+
      labs(title = paste('Group',j,'R square',round(r,3),'est.p',round(estimation.p,3),'trial.p',round(trial.p,3)))
  }}
marrangeGrob(grobs = plots,ncol = 2,nrow = 4)%>%
  ggsave(paste('/Users/chenghui/Documents/CBSB3/ICA/Water_Maze_Modelling/data/simplevsiter8.pdf'),
         plot=.,width = 12,height = 16,device = 'pdf',dpi = 300)

plots=list()
count=0
for(j in 1:4){
  for(i in c(0)){
    count=count+1
    datas=data[which(data$Group==j),]
    datass=datas[,c(2,i+15,24:25)]
    for(m in 1: length(datass$estimation)){
      datass[m,2]=Re(datass[m,2])
    }
    datass[,2]=as.numeric(datass[,2])
    datass=aggregate(gof~Day+estimation,data=datass,FUN = 'mean')
    model=summary(lm(datass[,1]~datass$estimation+datass$Day))
    r=model$r.squared
    estimation.p=model$coefficients[2,4]
    trial.p=model$coefficients[3,4]
    plots[[count]]=ggplot(datass,aes_string(x='Day',y=header[i+15],color='estimation'))+
      geom_line()+
      geom_point()+
      theme_classic()+
      labs(title = paste('Group',j,'R square',round(r,3),'est.p',round(estimation.p,3),'trial.p',round(trial.p,3)))
  }}

# Plot parameters
plot_param=function(data){
  param=c('sigma_pc','Vdecay','ac_const','beta','alpha','gamma')
  plotlist=list()
  count=0
for(i in c(5,7:9,11,12)){
  count=count+1
  datas=data[,c(1,2,3,i)]
  datas$day_trial=(datas$Day-1)*4+datas$Trial
  datas$Group=as.factor(datas$Group)
  model=lm(datas[,4]~datas$Group+datas$day_trial)
  print(param[count])
  print(emmeans(model,pairwise~Group)$contrasts)
  plotlist[[count]]=ggplot(datas,aes_string(x='day_trial',y=param[count],color='Group'))+
    geom_point(alpha=0.7)+geom_smooth(alpha=0.5,method=loess,level=0)+theme_classic()+
    xlab('Trial')
}
  return(plotlist)
}
plot_param(data2)
marrangeGrob(grobs = plot_param(data1),ncol = 2,nrow = 3)%>%
  ggsave(paste('/Users/chenghui/Documents/CBSB3/ICA/Water_Maze_Modelling/data/param~trial_simpleMC.pdf'),
         plot=.,width = 12,height = 12,device = 'pdf',dpi = 300)

