
df <- c()
for (x in list.files(path =Sys.getenv("USERPROFILE") , pattern="tReasure_0.3.1.tar", recursive=TRUE, full.names =TRUE )) {
  print(paste("install", x, sep = " "))
}


install.packages(x, repos = NULL, type = "source")
