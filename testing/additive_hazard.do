//local drive Z:/
local drive /Users/Michael/Documents/merlin
cd "`drive'"
adopath ++ "`drive'"
adopath ++ "`drive'/merlin"
adopath ++ "`drive'/stmerlin"
clear all

tr:do ./build/buildmlib.do
mata mata clear

pr drop _all
set seed 7254
clear
set obs 10000
gen id1 = _n
expand 5
bys id1: gen id2 = _n
gen trt = runiform()>0.5

survsim stime dead , hazard(0.001 :+ 0.01 :* {t} :+ 0.05 :* trt) maxt(5)
stset stime, f(dead)
mat startvals = (0.04,0.01,0.1)

timer clear
timer on 1
merlin (stime trt rcs(stime, df(3) log orthog event) , family(hazard, failure(dead))), 

stmerlin trt , dist(addrcs) df(3)

timer off 1
timer list

//fit null exp and use as intercept starting value
