
devtools::load_all(".")

#Get data from Google Drive if more recent than on Desktop
require(pacman)
pacman::p_load(googledrive,googlesheets4,this.path,tidyverse)

main_dir<-here(..=2) #This is "C:/Users/John.Syslo/Documents"

#read in data files from desktop
#catch
d7_catch<-read.csv("C:\\Users\\John.Syslo\\Documents\\Deep-7-Benchmark-2024\\Data\\Total_catch.csv") 
d7_catch[,2]<-d7_catch[,2]/1000000 #Catch in millions of lbs
#remove 1948 
d7_catch<-d7_catch[-1,]

#cpue
d7_cpue<-read.csv("C:\\Users\\John.Syslo\\Documents\\Deep-7-Benchmark-2024\\Data\\CPUE.csv") 
#remove 1948 
d7_cpue<-d7_cpue[-1,]

#cpue se - commercial is given as SE, BFISH is given as CV
#need to be converted to log scale? unclear whether it is log(se) or CV
d7_se<-read.csv("C:\\Users\\John.Syslo\\Documents\\Deep-7-Benchmark-2024\\Data\\SE.csv") 
#remove 1948 
d7_se<-d7_se[-1,]

d7_se_log<-cbind(d7_se[,1],d7_se[,2]/d7_cpue[,2],d7_se[,3]) #here we have converted the se to a cv
#the FRS CV is too low, because of the sample size.
#Do we need to make relative and multiply by obs error as in BUGS model?
#BFISH SE already entered as CV - may need to adjust these 


#make FRS CV relative
#d7_se_log[,2]<-d7_se_log[,2]/min(d7_se_log[,2],na.rm=T)

#multiply FRS by 10 to scale with BFISH
#d7_se_log[,2]<-d7_se_log[,2]*10

#Reduce survey SE by 10 to fit tighter
#d7_se_log[,3]<-d7_se_log[,3]/1000

#The fixed.obsE is added to both CPUE series in JABBA
#add here for FRS ONLY and then set to 0.001 
d7_se_log[,2]<-(sqrt(d7_se_log[,2]^2+0.2^2)) #this will be nearly the same as multiplying it by 10 above


??JABBA
?fit_jabba                                                  
?jabba2jags
?mp_jabba
?build_jabba #this is where everything happens
?jabba_plots
?jbplot_ensemble


setwd(file.path(main_dir,"Results")) 


jbinput = build_jabba(catch=d7_catch,
                      cpue=d7_cpue,se=d7_se_log,
                      assessment="Deep7",scenario = "TestRun",
                      model.type ="Pella_m",
                      catch.cv = 0.15, 
                      r.prior = c(0.1,0.25),  
                      K.prior=c(29,0.5), 
                      psi.prior = c(0.5,0.5), 
                      psi.dist = "beta",
                      sigma.proc=T,  
                      igamma=c(0.001, 0.001),           
                      proc.dev.all = T,               
                      BmsyK =  0.5,  
                      shape.CV = 1.0,
                      sigma.est=TRUE,
                      fixed.obsE = 0.001,             
                      catch.metric = "Million lb",
                      projection = FALSE,
                      # TACs = seq(0.01,1.0,0.1),
                      # imp.yr=2024,
                      #  pyrs=5,
                      verbose=TRUE)

#plotting Indices
jbplot_indices(jbinput,legend.loc = "topleft")

#fit                      
fit_test = fit_jabba(jbinput,ni=1000,nt=10,nb=200,nc=3,verbose=TRUE,save.csvs = TRUE,do.ppc=TRUE,save.all = TRUE) #,output.dir=file.path(main_dir,"Results",scenario)) 

#save.jabba==TRUE #I believe this saves the posteriors

summary(fit_test)
fit_test$pars
fit_test$estimates
fit_test$stats