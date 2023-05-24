test_that("fit_jabba works", {
  #load in test data, should include dfs for catch, cpue, and se
  load(test_path("fixtures", "test_data.RData"))
  #create simple jabba input object
  jabba_obj <- build_jabba(catch = test.catch, 
                           cpue = test.cpue,
                           se = test.se, verbose = FALSE)
  #fit jabba model
  fit_test <- fit_jabba(jabba_obj, ni = 100000, nt = 5, nb = 7000, nc = 2, save.csvs = TRUE, output.dir = test_path("Fixtures"), jagsdir = test_path("Fixtures"))
  expect_equal(fit_test$pars$Median[2], 0.3838237) #r
  expect_equal(fit_test$pars$Median[5], 0.88276642020) #psi
  expect_equal(fit_test$pars$Median[9], 2.00000000000) #m
  expect_equal(fit_test$est.catch[66,2], 10256.308461) #last year of catch (mu)
  expect_equal(fit_test$stats$Value[5], 35.90) #RMSE
  expect_equal(fit_test$refpts$bmsy[1], 147956.9034898) #fit_test$refpts[1,]
  
  #test there are 4 files saved
  nfiles <- length(list.files(test_path("fixtures"), pattern = ".csv"))
  expect_equal(nfiles, 5)
  #delete all csv files created
  unlink(test_path("fixtures", "*.csv"))
  
})
