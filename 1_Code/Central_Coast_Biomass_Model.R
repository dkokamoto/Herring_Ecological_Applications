library(ggplot2);library(reshape2);
library(gridExtra);library(scales);library(grid)
library(zoo);library(dplyr);library(gdata);library(rstan)
library(lubridate)

biomass <- read.csv("2_Data/CC_BIOMASS_Data.csv")
data_mat <- acast(biomass, SECTION ~ YEAR ~ variable)

### create a numeric scaler so biomass is on a manageable scale 
### (mean sd in observed total biomass)
scaler <- mean(apply(data_mat[,,"TOT_BIOM"],2,sd,na.rm=T))

### create a numeric scaler so biomass is on a manageable scale
data_test <- t(data_mat[,,"SP_BIOMASS"]/scaler)
data_test2 <- t(data_mat[,,"TOTAL_CATCH"]/scaler)

### set observed NAs to zero (for catch) or the observed minimum
data_test[is.na(data_test)] <- min(data_test[data_test>0&!is.na(data_test)])
data_test2[is.na(data_test2)] <- 0

### inputs (p and q)
p=3 ### reproductive lag
q= 1 ### order of moving average in errors 
NS= length(unique(biomass$SECTION)) ### number of sites
T= length(unique(biomass$YEAR)) ### number of time periods

stan_data <- list(
  T=T,
  NS=NS,
  p=3,
  S_obs=data_test, 
  H= data_test2,
  q=1
)

fit_q <- stan(file= "1_Code/Varma_q.stan", 
                data=stan_data, iter=2000,warmup =1000,chains= 5,cores= 5)

 save(fit_q,data_test,data_test2,file="fit_q_CC.RData")
pars <- extract(fit_q)
p=3
dat <- melt(t(data_test),varnames= c("SECTION","TIME"),value.name = "obs")
dat <- join(dat,melt(t(data_test2),varnames= c("SECTION","TIME"),value.name = "H"))
dat$estimated <- melt(t(apply(exp(pars$S_hat),c(2,3),mean)),
                      varnames= c("SECTION","TIME"),value.name= "estimated")[,3]
dat$pred[dat$TIME>1988+p-1] <-melt(t(apply(exp(pars$S_mu),c(2,3),mean)),
                     varnames= c("SECTION","TIME"),value.name= "pred")[,3]
dat$CIL <- melt(t(apply(exp(pars$S_hat),c(2,3),quantile,probs= 0.05)),varnames= c("SECTION","TIME"),value.name= "CIL")[,3]
dat$CIU <- melt(t(apply(exp(pars$S_hat),c(2,3),quantile,probs= 0.95)),varnames= c("SECTION","TIME"),value.name= "CIU")[,3]
dat <- join(dat,ddply(dat,.(SECTION),summarise, mean= mean(estimated)))
dat <- join(dat,ddply(dat,.(SECTION),summarise, mean_tot= mean(estimated+H)))
dat <- join(dat,ddply(dat,.(SECTION),summarise, stdev2= sd(estimated+H)))
dat <- join(dat,ddply(dat,.(TIME),summarise, an_mean= mean(estimated)))
dat <- join(dat,ddply(dat,.(TIME),summarise, an_harv= mean(H)))
names(dat)
dat <- join(dat,unique(biomass[,c("SECTION","SECTION_NAME")]))
#levels(dat$SECTION_NAME)<- 
 # c("Area4","Upper Spiller","Englefield","Juan Perez","Kitkatla","Kwakshua/Kwakume","Kitasoo Bay","Laskeek","Louscoone Int","Higgins","Port Louis","Port Simpson","Powell/Lower Spiller","Rennel Sound","Skincuttle Inlt","Dundivan/Thompson")
#dat$SECTION_NAME <- factor(dat$SECTION_NAME,levels=  c("Port Louis","Rennel Sound","Englefield","Laskeek","Juan Perez","Skincuttle Inlt","Louscoone Int","Port Simpson","Area4","Kitkatla","Kitasoo Bay","Higgins","Upper Spiller","Powell/Lower Spiller","Dundivan/Thompson","Kwakshua/Kwakume"))    

plot.options <- theme(strip.background = element_rect(fill = NA, colour = NA), 
                      panel.grid.major.y =element_blank(), 
                      panel.grid.minor.y= element_blank(),
                      panel.grid.major.x =element_line(size=0.5), 
                      panel.grid.minor.x= element_line(size=0.25),
                      axis.text.x = element_text(size= 9),
                      strip.text.x = element_text(size= 9),
                      legend.key = element_blank(), 
                      legend.key.size = unit(0.4, "cm"),
                      legend.position = "right", 
                      legend.direction = "vertical", 
                      legend.key.width = unit(.5, "lines"), 
                      legend.title = element_blank(),
                      legend.background = element_rect(fill = NA, colour = NA),
                      panel.background = element_rect(fill = NA, colour = "black"),
                      plot.background = element_blank(), 
                      plot.margin = unit(c(0.05,0.05,0.05,0.05), "inches"))

dat$DATE <- parse_date_time(paste(dat$TIME,"3","15",sep= "-"),"%Y-%m-%d")
dat_breaks <- as.Date(parse_date_time(paste(seq(1990,2015,by= 5),"1","1",sep= "-"),"%Y-%m-%d"))
dat_limits <- as.Date(parse_date_time(c(paste(1988,"1","1",sep= "-"),paste(2015,"12","31",sep= "-")),"%Y-%m-%d"))
plot1 <- ggplot(aes(as.Date(DATE),(1/exp(mean(pars$q))*(obs)+H)*scaler/1000),data=dat)+
  facet_wrap(~SECTION_NAME)+
  scale_x_date(date_minor_breaks= "1 year",
               labels= c("'90","","'00","","'10",""),
               breaks= dat_breaks, 
               expand= c(0,0),
               limits= dat_limits)+
  geom_ribbon(aes(ymin= (CIL+H)*scaler/1000,ymax= (CIU+H)*scaler/1000),col="grey80",alpha=0.25)+
  ylab("Tonnes (X 1000)\n")+
  xlab("")+
  #geom_point(aes(y= pred*scaler/1000),col= "red",size= 0.5)+
  geom_line(aes(y= (estimated+H)*scaler/1000,colour= "Estimate"))+
  geom_hline(yintercept= 0)+
  geom_line(aes(y= (an_mean+an_harv)*scaler/1000,colour= "Stock Mean"))+
  geom_point(shape=21,fill= "white",size=0.8)+
  geom_segment(aes(y= H*scaler/1000,yend= 0.005*scaler/1000,
                   x= as.Date(DATE),xend=as.Date(DATE),colour= "Harvest"),size= 0.7)+
  scale_colour_manual(values= c("Estimate"="black","Harvest"="#d95f02","Stock Mean"="#7570b3"))+
  coord_cartesian(ylim=c(0,20))+
  theme_bw()+
  plot.options+
  theme(plot.margin = unit(c(0.05,0.05,0,0.05), "inches"),
        axis.title.x= element_blank(),
        axis.text.x= element_blank(),
        legend.position= c(0.925,0.8));plot1

plot2.q1 <- ggplot(aes(as.Date(DATE),H/(estimated+H)),data= dat)+
  geom_hline(yintercept= 0)+
  facet_wrap(~SECTION_NAME,ncol=6)+
  geom_linerange(aes(ymin= H/(CIL+H),ymax= H/(CIU+H)),size=0.4)+
  ylab("Exploitation \n Rate")+
  scale_x_date(date_minor_breaks= "1 year",
               labels= c("'90","","'00","","'10",""),
               breaks= dat_breaks,
               expand= c(0,0),
               limits= dat_limits)+
  xlab("")+
  geom_point(aes(y= H/(estimated+H),
                 shape= SECTION_NAME,fill= SECTION_NAME,colour= SECTION_NAME),size=1,alpha= 0.8)+
  scale_shape_manual(values= c(21,22,21:24))+
  scale_fill_manual(values= c("white","white",rep("black",2),rep("grey40",2)))+
  scale_colour_manual(values= c(1,1,"grey70","grey70","grey20","grey20"))+
  geom_hline(yintercept= 0.2,linetype= "dotted")+
  scale_y_continuous(breaks= seq(0,1,by= 0.2),labels= c("0","0.2","0.4","0.6","0.8","1"))+
  theme_bw()+
  plot.options+
  theme(axis.title.x= element_blank(),
        strip.text=element_blank(),
        legend.position= "none")

gAa <- ggplotGrob(plot1)
gBb <- ggplotGrob(plot2)
gCc <- ggplotGrob(plot2.q1)
gBb$widths <- gAa$widths
gCc$widths <- gAa$widths
pdf(width=8,height= 4.5,paste0("/Users/Dan/Dropbox/Post_Doc/MetapopulationModels/Metapopulation_Models/4_MS/Figures/combo_fig.pdf"),family= "Times")
a1 <- grid.arrange(gAa,gBb,gCc,heights= c(1.5,1.2,1.3))
dev.off()


plot3.q1 <-ggplot(aes(I(estimated+H)*scaler/1000,IDF),
                  data= subset(dat,TIME>1990&!(TIME%in%c(2008:2013))))+
  geom_ribbon(aes(ymin=0,ymax=IDF),fill= "grey70",alpha=0.5)+
 # geom_ribbon(aes(ymin=0,ymax=EH/(estimated+H)),fill= "blue",alpha=0.1)+
  #geom_line(aes(y=IDF_CIL),colour= "grey40",linetype= "dotted")+
  #geom_line(aes(y=IDF_CIU),colour= "grey40",linetype= "dotted")+
  geom_hline(yintercept= 0,size=0.25)+
  geom_vline(xintercept=0,size=0.25)+
  geom_vline(xintercept=20,size=0.25)+
  geom_hline(yintercept=1,size=0.25)+
  geom_line(aes(x= I(estimated+H)*scaler/1000,y=IDF,colour= "Ideal Free Distribution",size= "Ideal Free Distribution"))+
  geom_line(aes(y=EH/(estimated+H),colour= "Proportional Allocation",size="Proportional Allocation"))+
  geom_errorbar(aes(y= H/(estimated+H),ymin= H/(CIL+H),ymax= H/(CIU+H)),colour= "grey40",width=0.5,size= 0.25)+
  scale_colour_manual(values= c("Ideal Free Distribution"="grey45","Proportional Allocation"="blue",alpha= 0.5))+
  geom_rect(aes(ymax= H/(estimated+H),ymin=0,xmin= I(estimated+H)*scaler/1000-0.5,xmax= I(estimated+H)*scaler/1000+0.5,fill=SECTION_NAME),alpha=0.7,colour= NA)+
  geom_point(aes(y= H/(estimated+H),shape= SECTION_NAME,fill= SECTION_NAME),size=1.5)+
 
  geom_point(aes(y= H/(estimated+H),shape= SECTION_NAME,fill= SECTION_NAME),size=1.5,alpha= 0.8)+
  geom_text(aes(x=Inf, y = Inf, label = paste("'",substring(TIME, 3, 4),sep="")), hjust=1.05,vjust=1.5,
            data= subset(dat,SECTION==72&TIME>1990&!(TIME%in%c(2008:2013))),size=3.5)+
  facet_wrap(~TIME,ncol= 5)+
  scale_shape_manual(values= c(21,22,21,22,15,8))+
  scale_fill_manual(values= brewer.pal(11,"RdYlBu")[c(1,3,4,9,10,11)])+
  #scale_colour_manual(values= c("black","black",rep("grey60",2),rep("black",2)))+
  theme_bw()+
  scale_x_continuous(expand= c(0,0),limits= c(0,20),breaks = c(0,5,10,15))+
  scale_y_continuous(breaks= seq(0,1,by= 0.2),labels= c("0","20","40","60","80",""),
                     limits= c(0,1),expand= c(0,0))+
  ylab("Local Exploitation Rate (%)")+
  scale_size_manual(values= c("Ideal Free Distribution"=0.7,"Proportional Allocation"=0.4))+
  xlab("Estimated Pre-Spawning Biomass (x 1000)")+
  guides(size= FALSE)+
  plot.options+
  theme(axis.text.y= element_text(size=12),
        axis.text.x= element_text(size=12),
        axis.title.y= element_text(size=12),
        axis.title.x= element_text(size=12),
        legend.text= element_text(size=12),
        strip.text=element_blank(),
        legend.key.width= unit(1.5,"lines"),
        panel.grid.major.x= element_blank(),
        panel.grid.minor.x= element_blank(),
        panel.grid.major.y= element_blank(),
        panel.grid.minor.y=  element_blank(),
        legend.direction= "horizontal",
        legend.position= "top",legend.key.width=unit(2,"lines"),
        panel.background=element_rect(colour= NA,fil= NA),
        plot.background=element_rect(colour= NA,fil= NA),
        panel.border=element_rect(colour= NA,fil= NA))


gt <- ggplot_gtable(ggplot_build(plot3.q1))
gt$layout$clip[grep("panel",gt$layout$name)] <- "off"
pdf(width= 7,height= 5,"/Users/Dan/Dropbox/Post_Doc/MetapopulationModels/Metapopulation_Models/4_MS/Figures/Harvest_Rate_REF.pdf",family= "Times")
grid.draw(gt)
dev.off()

plot_H <- ggplot(aes(EH/(estimated+H),H/I(estimated+H)),data= subset(dat,TIME>1990&!(TIME%in%c(2008:2013))))+
  geom_point(aes(fill=((estimated+H)-mean_tot)/stdev2),shape=21,alpha=0.5)+
  scale_fill_gradientn(colours= rev(brewer.pal(10,"RdYlBu")),name= "")+
  geom_abline(aes(intercept= 0,slope=1))+
  scale_x_continuous(limits= c(0,0.91))+
  scale_y_continuous(limits= c(0,0.91))+
  scale_shape_manual(values= shapes)+
  scale_colour_manual(values= cols)+
  scale_fill_manual(values= fills)+
  coord_cartesian(ylim=c(0,1),
                  xlim=c(0,1))+
  theme_bw()+
  plot.options+
  theme(legend.position=c(0.6,0.8),
        legend.key.width= unit(1.5,"lines"),
        legend.direction= "vertical",
        axis.title.x= element_blank(),
        strip.text=element_blank(),
        panel.grid.major.y= element_line(size= 0.5),
        panel.grid.minor.y=  element_line(size= 0.5))+
  scale_size_continuous(limits= c(-1.7,3.7),range= c(0.1,6))


plot_H <- ggplot(aes(EH,H),data= subset(dat,TIME>1990&!(TIME%in%c(2008:2013))))+
  geom_point(aes(shape= SECTION_NAME,
                 fill= SECTION_NAME,colour= SECTION_NAME),size=1.5)+
  geom_point(aes(shape= SECTION_NAME,
                 fill= SECTION_NAME,colour= SECTION_NAME),
             size=1.5,alpha=0.8)+
  geom_abline(aes(intercept= 0,slope=1))+
  scale_shape_manual(values= shapes)+
  scale_colour_manual(values= cols)+
  scale_fill_manual(values= fills)+
  theme_bw()+
  plot.options+
  theme(legend.position="none",
        legend.key.width= unit(1.5,"lines"),
        legend.direction= "vertical",
        axis.title.x= element_blank(),
        strip.text=element_blank(),
        panel.grid.major.y= element_line(size= 0.5),
        panel.grid.minor.y=  element_line(size= 0.5))


plot_IDF <- ggplot(aes(IDF,H/I(estimated+H)),data= subset(dat,TIME>1990&!(TIME%in%c(2008:2013))))+
  geom_point()+
  geom_abline(aes(intercept= 0,slope=1))+
  scale_x_continuous(limits= c(0,0.91))+
  scale_y_continuous(limits= c(0,0.91))


pdf(width= 7,height= 6,"2_Figures/ARMA.pdf",family= "Times")
plot1
dev.off()


pdf(width= 7,height= 4,"2_Figures/Harvest_Rate.pdf",family= "Times")
plot2.q1
dev.off()

test <- ddply(dat,.(SECTION),summarize,std_ab = mean(estimated),SEC_CV = sd(estimated)/mean(estimated))
test2 <- join(test,dat)
test3 <- ddply(test2,.(TIME),summarize,spat_CV= sqrt(var(estimated/std_ab)*(6-1)/(6))/mean(estimated/std_ab),
               skew = mean((estimated-mean(estimated))^3)/(sum((estimated-mean(estimated))^2)/(6-1))^(3/2),
               spat_sd= sd(estimated),
               mean= mean(estimated))
test3$rolled_SD <- rollmean(test3$spat_sd, 3, na.pad=TRUE)
test3$rolled_CV <- rollmean(test3$spat_CV, 3, na.pad=TRUE)
test3$rolled_MU <- rollmean(test3$mean, 3, na.pad=TRUE)
test3$rolled_SKU <- rollmean(test3$skew, 3, na.pad=TRUE)

testb <- ddply(SHI2,.(SECTION),summarize,std_ab = mean(SHI_std),var)
test2b <- join(testb,SHI2)
test3b <- ddply(test2b,.(YEAR),summarize,spat_CV= sqrt(var(SHI_std,na.rm=T)),mean= mean(SHI_std,na.rm=T))
test3b$rolled_SD <- rollmean(test3b$spat_CV, 3, na.pad=TRUE)
test3b$rolled_CV <- rollmean(test3b$spat_CV/test3b$mean, 3, na.pad=TRUE)

fills <- c("white","white",brewer.pal(11,"RdYlBu")[c(4,9,10,11)])
cols <- brewer.pal(11,"RdYlBu")[c(1,3,4,9,10,11)]
shapes <- c(21,22,21,22,3,8)
  
plot4 <- ggplot(aes(TIME,(obs+H)*scaler/1000),data= subset(dat,SECTION%in%c(67,72,74,77,78,85),group= SECTION))+
  geom_rect(xmin= 2008,xmax= 2013,ymin= -Inf,ymax= Inf,fill= "grey90")+
  geom_path(aes(y= (an_mean+an_harv)*scaler/1000),col= "blue",size=2,alpha=0.5,lineend="round")+
  ylab("Tonnes (x 1000)")+
  xlab("")+
 # geom_segment(aes(y= mean(H*scaler/1000),yend= 0.005*scaler/1000,
#                   x= TIME,xend=TIME))+
  geom_line(aes(y= (estimated+H)*scaler/1000,colour= SECTION_NAME))+
  geom_point(aes(y= (estimated+H)*scaler/1000,shape= SECTION_NAME,fill= SECTION_NAME,colour= SECTION_NAME),size=1.5)+
  geom_line(aes(y= (estimated+H)*scaler/1000,colour= SECTION_NAME),alpha=0.8)+
  geom_point(aes(y= (estimated+H)*scaler/1000,shape= SECTION_NAME,fill= SECTION_NAME,colour= SECTION_NAME),size=1.5,alpha=0.8)+
  geom_hline(yintercept= 0.005*scaler/1000)+
  scale_shape_manual(values= shapes)+
  scale_fill_manual(values= fills)+
  scale_colour_manual(values= cols)+
  coord_cartesian(ylim=c(0,20),
                 xlim=c(1987.5,2015.5))+
  scale_x_continuous(expand= c(0,0))+
  scale_y_continuous(expand= c(0.01,0))+
  guides(fill= guide_legend(ncol= 2),colour= guide_legend(ncol=3),
         shape= guide_legend(ncol=3))+
  theme_bw()+
  plot.options+
  theme(legend.position="top",
        #axis.text.x= element_blank(),
        legend.key.width= unit(1.5,"lines"),
        legend.direction= "vertical",
        axis.title.x= element_blank(),
        strip.text=element_blank())

plot5 <- ggplot(aes(x=TIME,y=spat_CV),data=test3)+
  geom_rect(xmin= 2008,xmax= 2013,ymin= -Inf,ymax= Inf,fill= "grey90")+
  geom_line(aes(y=rolled_CV),size=1)+
  geom_point(size=1)+
  ylab("Spatial Variation (CV)")+
  xlab("")+
  theme_bw()+
  plot.options+
  coord_cartesian(ylim=c(0,1.4),
                  xlim=c(1987.5,2015.5))+
  scale_x_continuous(expand= c(0,0))+
  theme(strip.text=element_blank(),
        axis.text.x= element_blank(),
        axis.title.x= element_blank(),
        legend.position= c(0.4,0.9),
        legend.key.width = unit(1.25, "lines"))

plot6 <- ggplot(aes(TIME,(obs+H)*scaler/1000),data= subset(dat,SECTION%in%c(67,72,74,77,78,85),group= SECTION))+
  ylab("Exploitation Rate")+
  geom_rect(xmin= 2008,xmax= 2013,ymin= -Inf,ymax= Inf,fill= alpha("grey90",0.5))+
  geom_line(aes(y= an_harv/(an_mean+an_harv)),col= "blue",size=2,alpha=0.5)+
  xlab("")+
  # geom_segment(aes(y= mean(H*scaler/1000),yend= 0.005*scaler/1000,
  #                   x= TIME,xend=TIME))+
  geom_line(aes(y= H/(estimated+H),colour= SECTION_NAME))+
  geom_hline(yintercept= 0.2,linetype= "dotted")+
  geom_point(aes(y= H/(estimated+H),shape= SECTION_NAME,fill= SECTION_NAME,colour= SECTION_NAME),size=1.5)+
  geom_line(aes(y= H/(estimated+H),colour= SECTION_NAME),alpha=0.8)+
  geom_point(aes(y= H/(estimated+H),shape= SECTION_NAME,fill= SECTION_NAME,colour= SECTION_NAME),size=1.5,alpha=0.8)+
  #geom_point(aes(y= (estimated+H)*scaler/1000,shape= SECTION_NAME,fill= SECTION_NAME),size=0.25,alpha= 0.8)+
  geom_hline(yintercept= 0)+
  scale_shape_manual(values= shapes)+
  scale_colour_manual(values= cols)+
  scale_fill_manual(values= fills)+
  scale_y_continuous(breaks= c(0,0.25,0.5,0.75,1),labels= c(0,"",0.5,"",1))+
  coord_cartesian(ylim=c(0,1),
                  xlim=c(1987.5,2015.5))+
  scale_x_continuous(expand= c(0,0))+
  theme_bw()+
  plot.options+
  theme(legend.position="none",
        legend.key.width= unit(1.5,"lines"),
        legend.direction= "vertical",
        axis.title.x= element_blank(),
        strip.text=element_blank(),
        panel.grid.major.y= element_line(size= 0.5),
        panel.grid.minor.y=  element_line(size= 0.5))
gSa <- ggplotGrob(plot4)
gSb <- ggplotGrob(plot5)
gSc <- ggplotGrob(plot6)
gSb$widths <- gSa$widths
gSc$widths <- gSa$widths

pdf(width= 5,height= 6,"/Users/Dan/Dropbox/Post_Doc/MetapopulationModels/Metapopulation_Models/4_MS/Figures/Figure_4.pdf",family= "sans",pointsize= 9)
grid.arrange(gSa,gSb,gSc,ncol= 1,heights= c(1.5,1,1))
dev.off()


a <- array(NA,dim=c(4500,25,6))               
for(i in 1:4500){a[i,,] <- pars$S_hat_pre[i,,]-pars$S_mu[i,,]}
Cor <- matrix(apply(apply(pars$L_Sigma,1,function(x) diag(1/sqrt(diag(x%*%t(x))))%*%(x%*%t(x))%*%diag(1/sqrt(diag(x%*%t(x))))),1,mean),ncol=6)


dat2 <- dat
dat2 <- drop.levels(dat2)
a <- melt(Cor)
a$Source  <- "Demographic Volatility"
a2 <- melt(cor(apply(exp(pars$S_hat),c(2,3),mean)))

test <- array(NA,dim= c(4500,6))


test <- dat%>%
  group_by(SECTION_NAME,SECTION)%>%
  do(data.frame(spec.pgram(log(.$estimated),plot= FALSE,spans= 10)[c("freq","spec")]))%>%
  mutate(mu= sum(spec))
test2 <- data.frame(spec.pgram(log(unique(dat$an_mean)),plot=FALSE,spans= 10)[c("freq","spec")])%>%
  mutate(mu= sum(spec))

spec1 <- ggplot(aes(freq,spec),data= test)+
  geom_line(aes(colour= SECTION_NAME))+
  geom_ribbon(aes(ymax=spec,ymin= 0.01,freq),data= test2,fill= "grey70",alpha=0.8)+
  geom_path(aes(freq,spec),data= test2,colour= "grey20",size= 1.5,lineend= "round",alpha=0.8)+
  theme_bw()+
  geom_point(aes(shape= SECTION_NAME,fill= SECTION_NAME),size=2)+
  scale_shape_manual(values= c(21,22,21:24))+
  scale_fill_manual(values= fills)+
  scale_colour_manual(values= cols)+
  ylab("power")+
  plot.options+
  guides(col=guide_legend(ncol= 2),
         linetype= guide_legend(ncol= 2),
         shape=guide_legend(ncol= 2))+
  theme(legend.position="top",
        legend.key.width= unit(2,"lines"),
        legend.direction= "horizontal",
        axis.title.x= element_blank(),
        axis.text.x= element_blank(),
        plot.margin = unit(c(0.05,0.1,0.05,0.05), "inches"))+
  annotation_logticks(base= 10)+
  scale_x_continuous(limits= c(1/20,0.51),expand= c(0,0),breaks= c(1/8,1/4,1/2),labels= c("1/8","1/4","1/2"))+
  coord_cartesian(xlim= c(1/10,1/2),ylim= c(0.01,5),expand =c(0,0))+
  scale_y_continuous(trans= 'log10')

spec2 <- ggplot(aes(freq,spec/mu),data= test)+
  geom_ribbon(aes(ymax=spec/mu,ymin= 0.004,freq),data= test2,fill= "grey70",alpha=0.8,col= NA)+
  geom_path(aes(freq,spec/mu),data= test2,colour= "grey30",size= 1.5,alpha=0.8,lineend="round")+
  geom_line(aes(colour= SECTION_NAME))+
  theme_bw()+
  geom_point(aes(shape= SECTION_NAME,fill= SECTION_NAME),size=2)+
  scale_shape_manual(values= c(21,22,21:24))+
  scale_fill_manual(values= fills)+
  scale_colour_manual(values= cols)+
  ylab("normalized power")+
  plot.options+
  theme(legend.position= "none",
        plot.margin = unit(c(0.05,0.1,0.05,0.05), "inches"))+
  annotation_logticks(base= 10)+
  xlab(expression(paste("frequency ", (Years)^-1)))+
  scale_x_continuous(limits= c(1/20,0.51),expand= c(0,0),breaks= c(1/8,1/4,1/2),labels= c("1/8","1/4","1/2"))+
  scale_y_continuous(trans= 'log10', breaks= c(0.01,0.1),labels= c("0.01","0.1"))+coord_cartesian(xlim= c(1/10,1/2),ylim= c(0.3,0.006),expand =c(0,0))+
  xlab(expression(paste("frequency ", (Years)^-1)))

gSa <- ggplotGrob(spec1)
gSb <- ggplotGrob(spec2)
gSb$widths <- gSa$widths
pdf(width= 4,height= 6,"/Users/Dan/Dropbox/Post_Doc/MetapopulationModels/Metapopulation_Models/4_MS/Figures/spec.pdf",family= "Times")
grid.arrange(gSa,gSb,heights= c(1.1,1))
dev.off()



a2$Source <- "Estimated Abundance"
a3 <- rbind(a,a2)
a3$SEC <- factor(7-a3$Var1,labels=rev(c("Kitasoo Bay", "Higgins","Upper Spiller","Powell/Lower Spiller","Dundivan/Thompson","Kwakshua/Kwakume")))
a3$SEC2 <- factor(7-a3$Var2,labels=rev(c("Kitasoo Bay", "Higgins","Upper Spiller","Powell/Lower Spiller","Dundivan/Thompson","Kwakshua/Kwakume")))
a3$SEC2 <- factor(a3$SEC2,levels= rev(levels(a3$SEC2)))

pdf(width= 4,height= 6,"2_Figures/cors.pdf",family= "Times")
ggplot(aes(SEC2,SEC),data= a3)+
  facet_grid(Source~.)+
  geom_tile(aes(fill= value))+
  scale_fill_gradientn(colours= rev(brewer.pal(9,"RdYlBu")),name= expression(paste(rho," (correlation)")),hight= 1,low=-1)+
  scale_x_discrete(expand= c(0,0))+
  scale_y_discrete(expand= c(0,0))+
  theme(axis.text.x= element_text(angle= 30,hjust=1),
        axis.title=element_blank(),
        legend.direction= "horizontal",
        legend.position= "top",
        legend.key.height= unit(0.5,"lines"),
        legend.key.width= unit(2.5,"lines"),
        strip.background= element_blank())+
  guides(fill= guide_colourbar(title.position= "top",title.hjust= 0.5))
dev.off()



dat2 <- dat
dat2$TIME <- dat2$TIME-1
dat2$estimated2 <- dat2$estimated
dat2$H2 <- dat2$H

dat3 <- join(dat,dat2[,c("SECTION","TIME","estimated2","H2")])

ggplot(aes(x=log(estimated),y=estimated2<0.25*(mean_tot)),data= subset(dat3,I(obs+H)*scaler/1000>1))+
  geom_point(aes(fill=SECTION_NAME,shape= SECTION_NAME))+
  stat_smooth()
  theme_bw()+geom_hline(yintercept= 0)+geom_vline(xintercept= 0)+
  geom_point(aes(x=an_harv/(an_harv+an_mean),y=an_mean/mean(an_mean)),shape= "+",size= 3,col= "red")+
  scale_shape_manual(values= c(21,22,21:24))+ylab("standardized biomass")+
  scale_fill_manual(values= c("white","white",rep("black",2),rep("grey40",2)))+plot.options+
  scale_colour_manual(values= c("black","black",rep("black",2),rep("grey40",2)))+xlab("Exploitation Rate")

fit1 <- lm(log(estimated/estimated2)~)

ggplot(aes(log(H/(estimated+H)),log(e/mean)),data= dat3)+geom_point(aes(colour= log(estimated/estimated2)))+scale_colour_gradient(low= "red",high= "blue")

ggplot(aes((x=H/(estimated+H)),y=estimated/mean),data= subset(dat,H>0&estimated>1))+
  geom_point(aes(colour= SECTION_NAME))+
  theme_bw()+scale_y_continuous(trans= "log10")+
  geom_text(aes(x=an_harv/(an_harv+an_mean),y=an_mean/mean(an_mean),label= TIME))

head(dat)
