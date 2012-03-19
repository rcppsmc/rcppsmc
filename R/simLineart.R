simLineart <- function(len = 250)
{
  sim <- c();
  
  statex <- cumsum(rnorm(len));
  statey <- cumsum(rnorm(len));

  datax  <- statex + rt(len,df=4);
  datay  <- statey + rt(len,df=4);

  sim$state <- c(statex,statey);
  dim(sim$state) <- c(2,2);
  sim$data  <- c(datax,datay);
  dim(sim$data) <- c(2,2);

  invisible(sim);
}
