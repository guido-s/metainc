#
# (1) Make R packages available
#

library("devtools")
library("roxygen2")


#
# (2) Create documentation file(s)
#

document()


#
# (3) Build R package and PDF file with help pages
#

build()
build_manual()


#
# (4) Install R package
#

install()


#
# (5) Check R package
#

check()


#
# (6) Check examples
#

setwd("..")
run_examples("metainc", run_dontrun = TRUE, run_donttest = TRUE)
warnings()

