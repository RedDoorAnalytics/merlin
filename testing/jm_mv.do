//local drive Z:\
local drive /Users/Michael/Documents
cd "`drive'/stdev/gml"
adopath ++ "./gml"
clear all

do ./build/buildmlib.do
mata mata clear

use "`drive'/stdev/stjm/data/stjm_pbc_example_data",clear
stset stop, enter(start) f(event=1) id(id)


cap bys id: gen stime2=stime if _n==1
cap bys id: gen died2=died if _n==1


megenreg 	(logb fp(1)@lam1 fp(1)#M1[id] M2[id], family(gaussian) timevar(time)) 	///
			(pro  fp(1)@lam2 fp(1)#M3[id] M4[id], family(gaussian) timevar(time))	///
			(stime2 trt EV[logb]@alpha1 EV[pro]@alpha2 , family(weib, failure(died2))) , ///
			intmethod(mc)
