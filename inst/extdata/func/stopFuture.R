stopFuture <- function(x){
  tools::pskill(x,signal = tools::SIGTERM)
  tools::pskill(x,signal = tools::SIGKILL)
}
