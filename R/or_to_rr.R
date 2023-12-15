or_to_rr <- function(sims,br){
  sims_rr <- (sims/(1-br*(1-sims)))
  sims_rr
}