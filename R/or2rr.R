# Zhang J and Yu KF (1998):
# What’s the Relative Risk?: A Method of Correcting the Odds Ratio in
# Cohort Studies of Common Outcomes.
# JAMA 280, 1690-91.

or2rr <- function(or, br)
  or / (1 - br + br * or)
