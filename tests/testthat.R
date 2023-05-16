library(testthat)
library(JABBA)

Sys.setenv(JAGS_HOME="C:/JAGS-4.3.1")
test_check("JABBA")