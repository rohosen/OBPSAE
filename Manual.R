pack <- "OBPSAE"
path <- find.package(pack)
system(paste(shQuote(file.path(R.home("bin"), "R")),
             "CMD", "Rd2pdf", shQuote(path)))


###
library(tools)
Rd2latex('man/bpeFH.Rd')
Rd2latex('man/fitSpline.cv.Rd')
Rd2latex('man/Hospital.Rd')
Rd2latex('man/MSPE_boot.Rd')
Rd2latex('man/MSPE_boot_adjusted.Rd')
Rd2latex('man/MSPE_boot_augmented.Rd')
Rd2latex('man/MSPE_JNR.Rd')
Rd2latex('man/MSPE_McJack.Rd')
Rd2latex('man/MSPE_McJack_Benchmark.Rd')
Rd2latex('man/MSPE_McSpline.Rd')
Rd2latex('man/MSPE_McSpline_Benchmark.Rd')
Rd2latex('man/MSPE_naive.Rd')
Rd2latex('man/obpFH.Rd')
Rd2latex('man/obpFHbenchmark.Rd')
Rd2latex('man/obpFH_adjusted.Rd')
Rd2latex('man/obpFH_augmented.Rd')
