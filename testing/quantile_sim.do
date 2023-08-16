//local drive Z:\
local drive /Users/Michael/Documents
cd "`drive'/megenreg/Stata"
adopath ++ "./megenreg"
clear all

do ./build/buildmlib.do
mata mata clear

set seed 2487

local M = 300
local n = 5

set obs `M'
gen id = _n
gen b0 = rnormal(0,1)
expand `n'

gen x = runiform()

gen y = b0 + 10 * x + rchi2(3)

megenreg 	(y x M1[id] , family(gaussian))
mat svs = e(b)
megenreg 	(y x M1[id] , family(lquant, q(0.5))) 	///
			, from(svs) //intmethod(gh) intpoints(11) //diff 
