## Code to set up test fixtures
test.catch <- read.csv(file = file.path(getwd(), "Version1.1_files", "SWO_SA", "CatchSWO_SA.csv"))
test.cpue <- read.csv(file = file.path(getwd(), "Version1.1_files", "SWO_SA", "cpueSWO_SA.csv"))
test.se <- read.csv(file = file.path(getwd(), "Version1.1_files", "SWO_SA", "seSWO_SA.csv"))

test.cpue <- test.cpue[,c(1,2,5)] #keep Year, BRA_LL1, and SPA_LL2
test.se <- test.se[,c(1,2,5)] #keep Year, BRA_LL1, and SPA_LL2

dir.create(file.path(getwd(), "tests", "testthat", "fixtures"))
save(list=c("test.catch", "test.cpue", "test.se"), file = file.path(getwd(), "tests", "testthat", "fixtures", "test_data.RData"))


