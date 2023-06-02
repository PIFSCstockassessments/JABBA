
#' Prior and posterior distributions of derived quantities: MSY, BMSY, HMSY, and BMSY/K
#'
#' @param jabba output list from fit_jabba
#' @param quants vector of derived quantities to plot ("MSY", "BMsy", "Fmsy", "BmsyK")
#' @param output.dir directory to save plots
#' @param as.png save as png file of TRUE
#' @param add if true don't call par() to allow construction of multiplots
#' @param width plot width
#' @param height plot hight
#' @param mfrow set up plot frame
#' @param addPP show PPMR and PPVR
#' @param verbose silent option  
#' @export
jbplot_ppderived <- function(jabba, quants = c("MSY", "Bmsy", "Fmsy", "BmsyK"), output.dir=getwd(),as.png = FALSE,mfrow=c(round((ncol(jabba$refpts_posterior))/3+0.33,0),3),width  = 8, height = 2.5*round(ncol(jabba$refpts_posterior)/3,0),cex=0.8,verbose=TRUE,addPP=TRUE){

if(verbose) cat(paste0("\n","><> jbplot_ppderived() - prior and posterior distributions of derived quantities <><","\n"))

Prs = as.matrix(cbind(jabba$settings$K.pr,jabba$settings$r.pr,c(jabba$settings$mu.m, jabba$settings$m.CV)))
out = jabba$refpts_posterior
node_id = names(out) 
node_id = node_id[node_id %in% quants]

Par = list(mfrow=mfrow,mai=c(0.4,0.1,0,.1),omi = c(0.3,0.5,0.1,0) + 0.1,mgp=c(1,0.1,0), tck = -0.02,cex=0.8)
  if(as.png){png(file = paste0(output.dir,"/Posteriors_",jabba$assessment,"_",jabba$scenario,".png"),width  = width, height = height,
      res = 200, units = "in")}
  par(Par)

 rpr =  rlnorm(10000,log(Prs[1,2]),Prs[2,2]) #column 2 is r prior
 mpr = rlnorm(10000,log(Prs[1,3]),Prs[2,3]) #column 3 is m prior
 kpr = rlnorm(10000,log(Prs[1,1]),Prs[2,1]) #column 1 is K prior
# Derived quantities
prior_bmsyK = mpr^(-(1/(mpr-1)))
prior_bmsy = kpr*mpr^(-1/(mpr-1))
prior_fmsy = (rpr/(mpr-1))*(1-1/mpr)
prior_msy = prior_fmsy*prior_bmsy


  for(i in 1:length(quants)){
    quant = quants[i]
    post.par = as.numeric(unlist(out[paste(quant)]))

    if(quant == "MSY"){
      pdf = stats::density(post.par,adjust=2)
      prior = stats::density(prior_msy)
      plot(pdf,type="l",ylim=range(0,pdf$y),
      xlim = range(c(pdf$x,quantile(prior_msy,c(0.0001,0.95)))),
      yaxt="n",xlab="MSY",ylab="",xaxs="i",yaxs="i",main="")
      polygon(prior, col = grey(0.4,1))
      polygon(pdf, col=gray(0.7,0.7))
      legend('right',c("Prior","Posterior"),pch=22,cex=cex+0.1,pt.cex=cex+0.1,pt.bg = c(grey(0.4,1),grey(0.8,0.6)),bty="n")
      PPVR = round((sd(post.par)/mean(post.par))^2/(sd(prior_msy)/mean(prior_msy))^2,3)
      PPVM = round(mean(post.par)/mean(prior_msy),3)
      legend("topright",c(paste("PPMR =",PPVM),paste("PPVR =",PPVR)),cex=cex,bty="n", x.intersp=0.5,y.intersp=0.85)

    }

    if(quant == "Bmsy"){
      
      pdf = stats::density(post.par,adjust=2)
      prior = stats::density(prior_bmsy)
      plot(pdf, type="l",ylim=range(0,pdf$y),
      xlim=range(c(pdf$x,quantile(prior_bmsy,c(0.0001,0.95)))),
      yaxt="n",xlab=expression(B[MSY]),ylab="",xaxs="i",yaxs="i",main="")
      polygon(prior, col = grey(0.4,1))
      polygon(pdf,col=gray(0.7,0.7))
      PPVR = round((sd(post.par)/mean(post.par))^2/(sd(prior_bmsy)/mean(prior_bmsy))^2,3)
      PPVM = round(mean(post.par)/mean(prior_bmsy),3)
     if(addPP) legend("topright",c(paste("PPMR =",PPVM),paste("PPVR =",PPVR)),cex=cex,bty="n", x.intersp=0.5,y.intersp=0.85)


    }

    if(quant == "Fmsy"){

        pdf = stats::density(post.par,adjust=2)
        prior = stats::density(prior_fmsy)
        plot(pdf$x,pdf$y,type="l",ylim=range(0,pdf$y),
        xlim=range(c(post.par,quantile(prior_fmsy,c(0.0001,0.95)))),
        yaxt="n",xlab=expression(F[MSY]),ylab="",xaxs="i",yaxs="i",main="")

        polygon(prior,col=gray(0.4,1))
        polygon(pdf,col=gray(0.7,0.7))
        PPVR = round((sd(post.par)/mean(post.par))^2/(sd(prior_fmsy)/mean(prior_fmsy))^2,3)
        PPVM = round(mean(post.par)/mean(prior_fmsy),3)
       if(addPP) legend("topright",c(paste("PPMR =",PPVM),paste("PPVR =",PPVR)),cex=cex,bty="n", x.intersp=0.5,y.intersp=0.85)

    }

    if(quant == "BmsyK"){
      
      pdf = stats::density(post.par,adjust=2)
      prior = stats::density(prior_bmsyK)

      plot(pdf,type="l",ylim=range(0,pdf$y),
      xlim=range(c(post.par,quantile(prior_bmsyK,c(0.001,0.95)))),
      yaxt="n",xlab=expression(frac(B[MSY],K)),ylab="",xaxs="i",yaxs="i",main="")
      polygon(prior, col=gray(0.4,1))
      polygon(pdf, col=gray(0.7,0.7))
      PPVR = round((sd(post.par)/mean(post.par))^2/(sd(prior_bmsyK)/mean(prior_bmsyK))^2,3)
      PPVM = round(mean(post.par)/mean(prior_bmsyK),3)
      if(addPP) legend("topright",c(paste("PPMR =",PPVM),paste("PPVR =",PPVR)),cex=cex,bty="n", x.intersp=0.5,y.intersp=0.85)

    }

  }
  mtext(paste("Density"), side=2, outer=TRUE, at=0.5,line=1,cex=0.9)
  if(as.png){dev.off()}

}