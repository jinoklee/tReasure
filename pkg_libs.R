withr::with <- makevars(c(PKG <- LIBS="-liconv"), install.packages("tidyverse"), assignment="+=")
