test_that("fit_jabba works", {
  #load in test data, should include dfs for catch, cpue, and se
  load(test_path("fixtures", "test_data.RData"))
  #create simple jabba input object
  jabba_obj_abs <- build_jabba(catch = test.catch, 
                           cpue = test.cpue,
                           se = test.se, verbose = FALSE)
  #fit jabba model
  fit_test_abs <- fit_jabba(jabba_obj_abs, ni = 30000, nt = 5, nb = 5000, nc = 2, save.csvs = TRUE, output.dir = test_path("Fixtures"), jagsdir = test_path("Fixtures"))
  expect_equal(fit_test$pars$Median[2], 0.29887245227) #r
  expect_equal(fit_test$pars$Median[5], 0.88068436754) #psi
  expect_equal(fit_test$pars$Median[9], 2.00000000000) #m
  expect_equal(fit_test$est.catch[66,2], 10258.2488575) #last year of catch (mu)
  expect_equal(fit_test$stats$Value[5], 36.80) #RMSE
  expect_equal(fit_test$refpts$bmsy[1], 179934.6423743) #fit_test$refpts[1,]
  
  #test there are 4 files saved
  nfiles <- length(list.files(test_path("fixtures"), pattern = ".csv"))
  expect_equal(nfiles, 5)
  #delete all csv files created
  unlink(test_path("fixtures", "*.csv"))
  
})

