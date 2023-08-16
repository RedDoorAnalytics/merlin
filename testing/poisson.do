//local drive Z:/
local drive /Users/Michael/Documents
cd "`drive'/merlin"
adopath ++ "./merlin"
clear all

tr:do ./build/buildmlib.do
mata mata clear

set seed 72549

pr drop _all
clear
set obs 100000
gen id1 = _n
gen trt = runiform()>0.5
gen age = rnormal(55,5)
gen y = rpoisson(0 + 0.1 * age -0.5 * trt)

timer clear
timer on 1
merlin (y trt age , family(poisson))
timer off 1
timer on 2
merlin (y trt age , family(poisson)), evaltype(gf1)
timer off 2
timer on 3
merlin (y trt age , family(poisson)), evaltype(gf2)
timer off 3
timer list
