//local drive Z:/
local drive /Users/Michael/Documents/merlin
cd "`drive'"
adopath ++ "`drive'"
adopath ++ "`drive'/merlin"
adopath ++ "`drive'/stmerlin"
clear all

tr:do ./build/buildmlib.do
mata mata clear

//requires merlin v1.13 
// net install merlin, from(https://www.mjcrowther.co.uk/code/merlin)

//simulate and fit a piecewise-exponential model with a piecewise time-dependent effect
pr drop _all
clear
set seed 249857
set obs 1000
gen trt = runiform()>0.5

//simulate from a piecewise baseline, and a piecewise time-dependent treatment effects
//_> will simulate with a hazard ratio of exp(-0.5) before 3 years, and exp(0.1) after 3 years
survsim stime died, maxtime(5) 									///
					hazard(	(0.1:*({t}:>0 :& {t}:<1) 			///
							:+ 0.2:*({t}:<2 :& {t}:>=1) 		///
							:+ 0.15:*({t}:>=2)) 				///
							:* exp(								///
							-0.5:*trt:*({t}:<3)					///
							:+ 0.1:*trt:*({t}:>=3))				///
							)

stset stime , f(died)

//fit the true model
merlin (stime trt /*#pc(stime, knots(3))*/ 	///
				, family(pwexp, knots(1 2) failure(died)))

// stmerlin , dist(pwe) knots(1 2)

//note - no intercept is contained, as pc() creates binary indicator variables for each of the 2 windows in 
//this case, so each of them will represent the (log) hazard ratio of the treatment effect in each 

/*
Fitting full model:

Iteration 0:   log likelihood = -37667.174  
Iteration 1:   log likelihood = -15115.742  
Iteration 2:   log likelihood = -14658.303  
Iteration 3:   log likelihood = -14547.777  
Iteration 4:   log likelihood = -14547.736  
Iteration 5:   log likelihood = -14547.736  

Mixed effects regression model                  Number of obs     =     10,000
Log likelihood = -14547.736
------------------------------------------------------------------------------
             |      Coef.   Std. Err.      z    P>|z|     [95% Conf. Interval]
-------------+----------------------------------------------------------------
stime:       |            
  trt#pc():1 |  -.5261864   .0357518   -14.72   0.000    -.5962586   -.4561143
  trt#pc():2 |   .1122879   .0395306     2.84   0.005     .0348093    .1897665
       dap:1 |  -2.351055   .0390934   -60.14   0.000    -2.427676   -2.274433
       dap:2 |  -1.570912   .0306238   -51.30   0.000    -1.630934   -1.510891
       dap:3 |  -1.898588   .0254346   -74.65   0.000    -1.948439   -1.848738
------------------------------------------------------------------------------
*/
