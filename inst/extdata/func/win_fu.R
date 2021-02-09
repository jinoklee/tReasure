win_fu <- function(h,...) {
  stopFuture <- function(x){
    tools::pskill(x,signal = tools::SIGTERM)
    tools::pskill(x,signal = tools::SIGKILL)
  }
  val <- gconfirm("Really close window", parent=h$obj)
  if(as.logical(val)){
    if(dir.exists("apid")){
      load("apid")
      stopFuture(apid)}
    if(dir.exists("bpid")){
      load("bpid")
      stopFuture(bpid)}
    dispose(window)
    gtkMainQuit()
  }else{
    return(TRUE)
  }
}
