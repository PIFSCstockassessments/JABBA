
test_that("build_jabba works", {
  #load in test data, should include dfs for catch, cpue, and se
  load(test_path("fixtures", "test_data.RData"))
  #check that dfs were loaded in properly
  expect_equal(min(test.catch$Year), 1950)
  expect_equal(max(test.catch$Year), 2015)
  
  jabba_obj <- build_jabba(catch = test.catch, 
                           cpue = test.cpue,
                           se = test.se,
                           assessment = "test",
                           model.type = "Pella_m",
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
                           verbose = FALSE)
  
  #check some settings
  expect_equal(jabba_obj$settings$psi.dist, "beta")
  expect_equal(jabba_obj$settings$BmsyK, 0.5)
  expect_equal(jabba_obj$settings$add.catch.CV, TRUE)
  #should be 66 years of data
  expect_equal(jabba_obj$jagsdata$N, 66)
  expect_equal(length(jabba_obj$data$yr), 66)
  #shouldn't be any 0s in se
  expect_equal(length(which(is.na(jabba_obj$jagsdata$SE2) == TRUE)), 0) 
  
  })

