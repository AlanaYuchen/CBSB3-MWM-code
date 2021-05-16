#======================================distance cell=========================================
setwd('/Users/chenghui/Documents/CBSB3/ICA/Water_Maze_Modelling/distance_based_model_random_simulation')
library(dplyr)
library(nlme)
PMs<-c( 'latency', 'distance', 'target quandrant', 'opposite quadrant','wall zone','speed_std','mean_angle','time step')

# fixed
fix0=read.csv('output_fixed_1000trial_alpha0.csv',header = F,col.names = c('PM','day','trial','rep','value'))
fix1=read.csv('output_fixed_1000trial.csv',header = F,col.names = c('PM','day','trial','rep','value'))


# variable
var0=read.csv('output_variable_1000trial_alpha0.csv',header = F,col.names = c('PM','day','trial','rep','value'))
var1=read.csv('output_variable_1000trial.csv',header = F,col.names = c('PM','day','trial','rep','value'))

mixed_effect=function(a,b){
  models<-list()
  pvalue=c()
  for(i in 1:8){
    a1=a[which(a$PM==i),]
    b1=b[which(b$PM==i),]
    c=rbind(a1,b1)
    d=cbind(c,group=c(rep('a',length(a1$PM)),rep('b',length(b1$PM))))
    model=lme(value~group+day*trial,data=d,random= ~1|rep)
    #print(emmeans(model,pairwise~group)$contrasts)
    #print(summary(model)$tTable[3,5])
    pvalue=c(pvalue,summary(model)$tTable[2,5])
    models[[i]]=model
  }
  return(pvalue)
}
mixed_effect(fix0,fix1)
DC_model_fix=mixed_effect(fix0,fix1)
mixed_effect(var0,var1)
DC_model_var=mixed_effect(var0,var1)

result=data.frame(PMs,DC_model_fix,DC_model_var)

#================================plotting======================================
# Regression on fixed platform with trial ~ behavioural variable
plotfun=function(data){
  mean_rep<-aggregate(value~PM+day+trial,data=data,FUN='mean')
  sd_rep<-aggregate(value~PM+day+trial,data=data,FUN='sd')
  mean_sd<-data.frame(mean_rep,sd_rep$value)
  names(mean_sd)[names(mean_sd) =="value"] <-"mean"
  names(mean_sd)[names(mean_sd) =="sd_rep.value"] <-"sd"
  
  bv<-c( 'latency', 'distance', ' quandrant4', 'quadrant2','wall zone','speed_std','mean_angle','time step')
  library(ggpubr) # for ggsave function
  library(ggplot2)
  
  # Calculate trial number by using 'day' and 'trial' in the output. E.g. day2-trial2 => the 6th trial
  day_trial<-(mean_sd$day-1)*4+mean_sd$trial
  mean_sd<-data.frame(mean_sd,day_trial)
  #  plotbox<-list()
  plotline<-list()
  # Boxplot with regression line for trial
  # for (i in 1:8){
  #   ss<-mean_sd[which(mean_sd$PM==i),]
  #   # linear regression
  #   l<-lm(ss$mean~ss$day_trial)
  #   plotbox[[i]]<-ggplot(ss,aes(x=day_trial,y=mean,group=day_trial))+
  #     geom_boxplot()+
  #     geom_abline(slope = l$coefficients[[2]],intercept = l$coefficients[[1]],color='red')+ # add regression line
  #     labs(x="trial",y=bv[i]) # add labels
  # }
  # 
  # calculate the means over replicates (the mean of data with same PM&day&trial)
  mm<-aggregate(mean~PM+day+trial,data = mean_sd,FUN='mean')
  # combine mean and sd into the same data frame
  mm<-data.frame(mm,mean_sd$sd)
  # change column name
  names(mm)[names(mm)=="mean_sd.sd"]<-'sd'
  
  # Calculate trial number by using 'day' and 'trial'
  day_trial2<-(mm$day-1)*4+mm$trial
  mm<-data.frame(mm,day_trial2)
  
  # plot line chart with regression line
  for (i in 1:8){
    ss2<-mm[which(mm$PM==i),]
    l1<-lm(ss2$mean~ss2$day_trial2)
    plotline[[i]]<-ggplot(ss2,aes(x=day_trial2,y=mean))+
      geom_line(key_glyph = draw_key_abline) + 
      geom_point(size=1, key_glyph = draw_key_point) + 
      geom_errorbar(aes(ymin = mean - sd/2, ymax = mean + sd/2), width = 0.1, key_glyph = draw_key_blank, linetype=2)+
      geom_abline(slope = l1$coefficients[[2]],intercept = l1$coefficients[[1]],color='red')+
      labs(x="trial",y=bv[i],subtitle =  paste("R2=",round(summary(l1)$r.squared,4),                                                "p.val=",round(summary(l1)$coefficients[8],4)))+
      theme_classic()
  }
  
  library(ggpubr)
  ggsave(filename=paste0(deparse(substitute(data)),"_trial~PM.pdf"),plot=ggarrange(plotlist = plotline),width = 12,height = 10,device = 'pdf',dpi = 300)
}

plotfun(fix0)
plotfun(fix1)
plotfun(var0)
plotfun(var1)
#====================================place cell===============================================
setwd('/Users/chenghui/Documents/CBSB3/ICA/Water_Maze_Modelling/place_cell_based_model_simulation')
library(dplyr)
library(nlme)
PMs<-c( 'latency', 'distance', 'target quandrant', 'opposite quadrant','wall zone','speed_std','mean_angle','time step')

# fixed
fix3=read.csv('output_fix_1000trial_alpha0.csv',header = F,col.names = c('PM','day','trial','rep','value'))
fix4=read.csv('output_fix_1000trial.csv',header = F,col.names = c('PM','day','trial','rep','value'))

# variable
var3=read.csv('output_variable_1000trial_alpha0.csv',header = F,col.names = c('PM','day','trial','rep','value'))
var4=read.csv('output_variable_1000trial.csv',header = F,col.names = c('PM','day','trial','rep','value'))

mixed_effect=function(a,b){
  models<-list()
  pvalue=c()
  for(i in c(1:8)){
    a1=a[which(a$PM==i),]
    b1=b[which(b$PM==i),]
    c=rbind(a1,b1)
    d=cbind(c,group=c(rep('a',length(a1$PM)),rep('b',length(b1$PM))))
    d=d[complete.cases(d),]
    model=lme(value~group+day*trial,data=d,random= ~1|rep)
    print(summary(model)$tTable[3,5])
    pvalue=c(pvalue,summary(model)$tTable[2,5])
    models[[i]]=model
  }
  return(pvalue)
}

mixed_effect(fix3,fix4)
PC_model_fix=mixed_effect(fix3,fix4)
mixed_effect(var0,var1)
PC_model_var=mixed_effect(var3,var4)

result=data.frame(PMs,PC_model_fix,PC_model_var,DC_model_fix,DC_model_var)

#================================plotting======================================
# Regression on fixed platform with trial ~ behavioural variable
plotfun=function(data){
  mean_rep<-aggregate(value~PM+day+trial,data=data,FUN='mean')
  sd_rep<-aggregate(value~PM+day+trial,data=data,FUN='sd')
  mean_sd<-data.frame(mean_rep,sd_rep$value)
  names(mean_sd)[names(mean_sd) =="value"] <-"mean"
  names(mean_sd)[names(mean_sd) =="sd_rep.value"] <-"sd"
  
  bv<-c( 'latency', 'distance', ' quandrant4', 'quadrant2','wall zone','speed_std','mean_angle','time step')
  library(ggpubr) # for ggsave function
  library(ggplot2)
  
  # Calculate trial number by using 'day' and 'trial' in the output. E.g. day2-trial2 => the 6th trial
  day_trial<-(mean_sd$day-1)*4+mean_sd$trial
  mean_sd<-data.frame(mean_sd,day_trial)
  #  plotbox<-list()
  plotline<-list()
  # Boxplot with regression line for trial
  # for (i in 1:8){
  #   ss<-mean_sd[which(mean_sd$PM==i),]
  #   # linear regression
  #   l<-lm(ss$mean~ss$day_trial)
  #   plotbox[[i]]<-ggplot(ss,aes(x=day_trial,y=mean,group=day_trial))+
  #     geom_boxplot()+
  #     geom_abline(slope = l$coefficients[[2]],intercept = l$coefficients[[1]],color='red')+ # add regression line
  #     labs(x="trial",y=bv[i]) # add labels
  # }
  # 
  # calculate the means over replicates (the mean of data with same PM&day&trial)
  mm<-aggregate(mean~PM+day+trial,data = mean_sd,FUN='mean')
  # combine mean and sd into the same data frame
  mm<-data.frame(mm,mean_sd$sd)
  # change column name
  names(mm)[names(mm)=="mean_sd.sd"]<-'sd'
  
  # Calculate trial number by using 'day' and 'trial'
  day_trial2<-(mm$day-1)*4+mm$trial
  mm<-data.frame(mm,day_trial2)
  
  # plot line chart with regression line
  for (i in 1:8){
    ss2<-mm[which(mm$PM==i),]
    l1<-lm(ss2$mean~ss2$day_trial2)
    plotline[[i]]<-ggplot(ss2,aes(x=day_trial2,y=mean))+
      geom_line(key_glyph = draw_key_abline) + 
      geom_point(size=1, key_glyph = draw_key_point) + 
      geom_errorbar(aes(ymin = mean - sd/2, ymax = mean + sd/2), width = 0.1, key_glyph = draw_key_blank, linetype=2)+
      geom_abline(slope = l1$coefficients[[2]],intercept = l1$coefficients[[1]],color='red')+
      labs(x="trial",y=bv[i],subtitle =  paste("R2=",round(summary(l1)$r.squared,4),                                                "p.val=",round(summary(l1)$coefficients[8],4)))+
      theme_classic()
  }
  
  library(ggpubr)
  library(gridExtra)
  ggsave(filename=paste0(deparse(substitute(data)),"_trial~PM.pdf"),plot=plotline,width = 12,height = 10,device = 'pdf',dpi = 300)
}


plotfun(fix3)
plotfun(fix4)
plotfun(var3)
plotfun(var4)
#===================================================================================
# adjust wall zone
wall08<-read.csv("output_fix_alpha001_w08.csv",header = F,col.names = c('PM','day','trial','rep','value'))
wall07<-read.csv("output_fix_alpha001_w07.csv",header = F,col.names = c('PM','day','trial','rep','value'))
plotfun(wall08)
plotfun(wall07)
