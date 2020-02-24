## Data generation functions

generate_data_binary <- function(Nabc,spec){
  N     = Nabc[1] #Number of samples
  Z <- rnorm(N)
  Z2 <- Z^2
  X = rbinom(N,1,1/(exp(-Z -spec[1]*Z2)+1))
  M = Nabc[2]*X + Z+spec[2]*Z2 +rnorm(N)
  Y = Nabc[3]*M + Nabc[4]*X + Z +spec[3]*Z2 +rnorm(N)
  return(data.frame(Z=Z,X=X,M=M,Y=Y))
}

gd_bin    <- function(Nabc){generate_data_binary(Nabc,c(0,0,0))}
gd_bin_XY <- function(Nabc){generate_data_binary(Nabc,c(0,1,0))}
gd_bin_MX <- function(Nabc){generate_data_binary(Nabc,c(0,0,1))}
gd_bin_X  <- function(Nabc){generate_data_binary(Nabc,c(0,1,1))}
gd_bin_M  <- function(Nabc){generate_data_binary(Nabc,c(1,0,1))}
gd_bin_Y  <- function(Nabc){generate_data_binary(Nabc,c(1,1,0))}










