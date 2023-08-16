//local drive Z:/
local drive /Users/Michael/Documents
cd "`drive'/merlin"
adopath ++ "./merlin"
clear all

tr:do ./build/buildmlib.do
mata mata clear

set seed 725488
pr drop _all
clear
set obs 1000
gen id1 = _n
expand 5
bys id1: gen id2 = _n
gen trt = runiform()>0.5
mat cor1 = (1,0.25\0.25,1)
drawnorm u1 u2, means(0 0) sds(1 0.5) corr(cor1)
bys id1 (id2) : replace u1 = u1[1]
bys id1 (id2) : replace u2 = u2[1]
gen trtui = (-0.5+u2) * trt

survsim stime dead , dist(weib) lambda(0.1) gamma(1.2) cov(trtui 1 u1 1) maxt(5) 
stset stime, f(dead)

timer clear
timer on 1
mestreg trt || id1: trt, dist(weib) cov(unstr)
timer off 1
predict refs1*, reffects

timer on 2
merlin (stime trt trt#M1[id1]@1 M2[id1]@1, family(weib, failure(dead))), cov(unstr)
timer off 2
timer list
predict refs21 refs22, reffects
