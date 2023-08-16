
//local drive Z:/
local drive /Users/Michael/Documents
cd "`drive'/merlin"
adopath ++ "./merlin"
clear all

do ./build/buildmlib.do
mata mata clear

set seed 72549

pr drop _all
clear

// In this example, I will first show how to simulate interval censored survival times, and then show how to use `merlin` 
// to fit an interval censored flexible parametric survival model. There are two ways to simulate interval censored survival 
// times - either under a discrete time setting, or a continuous time setting. I will focus on the second, as I wish to 
// fit a continuous time survival model, whilst accounting for interval censoring.

// First we'll simulate a dataset of 1000 observations, and include a binary variable which we will assume is 
// a treatment group indicator coded 0 for control, and 1 for active treatment. I'll assume allocation to each 
// group is 50% at random.

set obs 100000
gen trt = runiform()>0.5

// To keep things nice and simple I will use [`survsim`]() to simulate survival times from a continuous time, Weibull 
// distribution with shape and scale parameters, 1.2 and 0.1. I'll assume a log hazard ratio due to treatment of -0.5. 
// I'll apply administrative censoring at 5 years.

survsim stime event , dist(weib) lambda(0.1) gamma(1.2) covariates(trt -0.5) maxt(5)

// I've simulated continuous, exacty observed survival times, and will now apply interval censoring. I'm going to assume 
// observations were only made at annual intervals, i.e. every year.
local maxt 5
//So any events will be rounded down to the nearest year and the left interval taken as the previoys year
gen st1 = floor(stime) if stime < `maxt'
replace stime = st1 + 1 if st1 < `maxt'

// merlin identifies interval censored oberservations through the `failure()` indicator being coded as a 2. 
replace event = 2 if event == 1

// In this example I'm assuming all observations are either interval censored events, or right censored. `merlin` 
// allows any combination of exactly observed events, interval censoring, right censoring and left truncation.

// Now we are all setup and ready to fit an interval censored survival model with `merlin`. My approach to methods 
// development these days focusses on general code, in that by adding interval censoring to merlin, it means that it will 
// work with any of the inbuilt survival distributions, including the user-defined ones. However, here I will keep it 
// relatively simple and show how to fit a Royston-Parmar flexible parametric survival model, allowing for interval 
// censoring.

// The new option added to `merlin` is the `linterval()` option, which lets you pass the variable which contains the 
// left interval time for those observations which are interval censored, which must have their associated event indicator 
// specified as a 2.

merlin 	(stime trt , family(rp, df(3) failure(event) linterval(st1))) , evaltype(gf1)

// replace st1 = stime if event==0
// replace stime = . if event==0
// stintreg trt, dist(weib) interval(st1 stime) nohr
// Our model fits smoothly and we recover our log hazard ratio. For comparison, we can also fit the same model using 
// Patrick Royston's `stpm` command (this has been superceded by `stpm2`, but `stpm2` does not support interval 
// censoring).

// stpm requires that for non interval censored observations, the left interval contains the event time, so we 
// make that adjustment first. Our data must also be stset for `stpm`.
gen died = event>0
stset stime, f(died)
replace st1 = _t if st1==.
stpm trt, scale(h) df(3) left(st1)

// We obtain identical estimates and log-likelihood, which is reassuring!

