# JABBA: Just Another Bayesian Biomass Assessment
The materials in this repository present the stock assessment tool ‘Just Another Bayesian Biomass Assessment’ JABBA. The motivation for developing JABBA was to provide a user-friendly R to JAGS (Plummer) interface for fitting generalized Bayesian State-Space SPMs with the aim to generate reproducible stock status estimates and diagnostics. Building on recent advances in optimizing the fitting procedures through the development of Bayesian state-space modelling approaches, JABBA originates from a continuous development process of a Bayesian State-Space SPM tool that has been applied and tested in many assessments across oceans. JABBA was conceived in the Hawaiian Summer of 2015 as a collaboration between young researchers from South Africa and the Pacific Islands Fisheries Science Center (NOAA) in Honolulu, HI USA. The goal was to provide a bridge between age-structured and biomass dynamic models, which are still widely used. JABBA runs quickly and by default generates many useful plots and diagnosic tools for stock assessments.

Inbuilt JABBA features include:

+ Integrated state-space tool for averaging multiple CPUE series (+SE) for optional use in assessments
+ Automatic fitting of multiple CPUE time series and associated standard errors
+ Fox, Schaefer or Pella Tomlinson production function (optional as input Bmsy/K)
+ Kobe-type biplot plotting functions 
+ Forecasting for alternative TACs 
+ Residual and MCMC diagnostics 
+ Estimating or fixing the process variance
+ Optional estimation additional observation variance for individual or grouped CPUE time series
+ Easy implementation of time-block changes in selectivity
+ New: Estimating Catch with Error 
+ New: Estimating shape with prior
+ New: Catch-Only option with additional relative biomass priors
+ New: Inbuilt retrospective and hindcasting run and plotting options 

## Installing JABBA as R package

`library(remotes)` <br>
`install_github("jabbamodel/JABBA")`

To install the PIFSC-dev branch use:
`remotes::install_github("PIFSCstockassessments/JABBA", ref = "PIFSC-dev")`

## Test-drive JABBA

`library(JABBA)`<br>
`data(iccat)`

`# Compile JABBA JAGS model and input object for bigeye tuna (bet)`<br>
`jbinput = build_jabba(catch=iccat$bet$catch,cpue=iccat$bet$cpue,se=iccat$bet$se,assessment="BET",scenario = "TestRun",model.type = "Fox",sigma.est = FALSE,fixed.obsE = 0.01)`

`# Fit JABBA (here mostly default value - careful)`<br>
`bet1 = fit_jabba(jbinput,quickmcmc=TRUE)`

`# Make individual plots` <br>
`jbplot_catcherror(bet1)` <br>
`jbplot_ppdist(bet1)`     <br>
`jbplot_cpuefits(bet1)`<br>
`jbplot_logfits(bet1)`<br>

`# Plot Status Summary` <br>
`par(mfrow=c(3,2),mar = c(3.5, 3.5, 0.5, 0.1))` <br>
`jbplot_trj(bet1,type="B",add=T)` <br>
`jbplot_trj(bet1,type="F",add=T)` <br>
`jbplot_trj(bet1,type="BBmsy",add=T)`<br>
`jbplot_trj(bet1,type="FFmsy",add=T)`<br>
`jbplot_spphase(bet1,add=T)`<br>
`jbplot_kobe(bet1,add=T)`<br>
`# Test run end`

More detailed examples of JABBA applications to ICCAT stocks can be found in [Examples_iccat.R](https://github.com/jabbamodel/JABBA/blob/master/Example/Examples_iccat.R).

An updated vignette is currently under development...

## PIFSC-dev additions 

In `build_jabba()`:  
*   catch.error - added the option to use a uniform distribution for catch error for the first 55 years of the model, then lognormal distribution is used for the remaining years. To use this options set argument to "deep7".  
*   catch.adj - argument needed if catch.error = "deep7". Input a vector of the adjustment (up and down) to catch for years of uniform distribution eg c(0.8,1.2)
*   cpue_lamda - added for data weighting, it takes a vector of weights for the cpue series. If you want to weight all cpue equally, set all to 1
*   index_type - added to indicate if an absoulte abundance is used for bridging between 2019 BFISH assessment and 2023. If an absolute abundance index is used, set the argument equal to a vector the length of n.cpue.series and indicate either "relative" or "absolute" for each cpue eg c("relative", "absolute"). If you are not using an absolute abundance, can ignore (default is NULL). 
    *   When you indicate an absolute abundance, inside `build_jabba` function, it will calculate the number of estimated q's (nran.q) as `length(unqiue(sets.q[-abs.ind]))` where abs.ind is the index of cpue that are absolute abundances. This variable is used inside jabba2jags code to set the q priors. If there are fixed q's then it will set the priors differently than if they are all estimated.
*   rad.prior - if including an absolute abundance index and using a radius prior (for camera survey) include a vector of radius mean and CV
*   n.grid - needed if using rad.prior, indicate the number of sampling grids in a domain
*   a.grid - needed if using rad.prior, indicate the area within a sampling grid
*   s_lambda - needed if using rad.prior, used to indicate uncertainty with BFISH survey 
*   nsig.off.ind - added to allow estimating observation error for some indices but not others. The argument takes either 0 or 1 (defaults to 0), and this indicates how many indices you don't want to estimate observation error for. Note that it is assumed the index which obsE is not estimated for is the last column in cpue.
*   n.var - adjust this argument when using nsig.off.ind to indicate how many indices you will be estimating obsE for

`jabba_plots()`:  
*   jbplot_ppdist - added radius prior/posterior distributions to the plot if they are specified
*   jbplot_kobe_bfrac - kobe plot with state-dependent reference point (MSST), in `jabba_plots()` can change `statusplot = "bfrac"` to generate this version of kobe plot.
*   jbplot_TOE - stacked total observation error plot

`jbplot_ppderived()`: 
*   function to plot prior and posterior distributions of derived quantities (MSY, BMsy, Fmsy, and BmsyK)

**Reference**  
[Winker, H., Carvalho, F., Kapur, M. (2018) <U>JABBA: Just Another Bayesian Biomass Assessment.</U> *Fisheries 
Research* **204**: 275-288.](https://www.sciencedirect.com/science/article/pii/S0165783618300845)   


--------------------------------------------------------------------------------

#### [JABBAbeta](https://github.com/Henning-Winker/JABBAbeta) GitHub repository
JABBA development version for testing new JABBA features and stock assessment examples 
