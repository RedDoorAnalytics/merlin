//local drive Z:/
local drive /Users/Michael/Documents
cd "`drive'/forge"
adopath ++ "./forge"
clear all

do ./build/buildmlib.do
mata mata clear

use "./data/orth",clear
gen agec = age - 11
gen sex = Sex==1
//drop if sex==1

//sort distance
//gsem (distance <- agec M1[id] , family(gaussian)), evaltype(gf0) intmethod(gh)

/*
megenreg (distance agec sex M1[id] , family(gaussian))
mat svs = e(b)
mat list svs
megenreg 	(distance agec M1[id] , family(user, llfunction(logl) nap(1))) ///
			, intmethod(gh) from(svs) intpoints(7)
mat svs = e(b)
*/
local q = 0.5

egen id1 = group(id)

forge 	(distance agec M1[id1]@1 , family(lquant, q(`q'))) ///
			, intmethod(gh) intpoints(25) //diff 
est store m1

sort id1
qui forge 	(distance agec M1[id1]@1 , family(lquant, q(`q'))) ///
			, intmethod(gh) intpoints(25) //diff 
est store m2

sort id1 distance
qui forge 	(distance agec M1[id1]@1 , family(lquant, q(`q'))) ///
			, intmethod(gh) intpoints(25) //diff 
est store m3

est tab m*
