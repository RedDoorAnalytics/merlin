//local drive Z:/
local drive /Users/Michael/Documents
cd "`drive'/merlin"
adopath ++ "./merlin"
clear all

do ./build/buildmlib.do
mata mata clear

use "/Users/Michael/Documents/merlin/data/stjm_pbc_example_data",clear
stset stop, enter(start) f(event=1) id(id)

cap bys id: gen stime2=stime if _n==1
cap bys id: gen died2=died if _n==1

merlin (stime2 trt EV[logb] trt#fp(stime2,pow(0)), family(weibull, failure(died2)) timevar(stime2)) ///
       (logb fp(time,pow(1)) fp(time,pow(1))#M2[id]@1 M1[id]@1, family(gaussian) timevar(time)), restartvalues(M2 0.1)
