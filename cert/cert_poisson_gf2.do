
set seed 725498

clear
set obs 1000
gen trt = runiform()>0.5

gen time = runiform()

gen y = rpoisson(0.4*time + 5*trt)

		
merlin (y time trt, family(poisson)) , evaltype(gf0)
mat b1 = e(b)
mat v1 = e(V)
local b1 : display %9.5f b1[1,1]
local v1 : display %9.5f v1[1,1]

merlin (y time trt, family(poisson)) , evaltype(gf1)
mat b2 = e(b)
mat v2 = e(V)
local b2 : display %9.5f b2[1,1]
local v2 : display %9.5f v2[1,1]

merlin (y time trt, family(poisson)) , evaltype(gf2)
mat b3 = e(b)
mat v3 = e(V)
local b3 : display %9.5f b3[1,1]
local v3 : display %9.5f v3[1,1]

assert `b3'==`b2'
assert `v3'==`v2'

merlin (y rcs(time, df(3) orthog) trt, family(poisson)) , evaltype(gf0)
mat b1 = e(b)
mat v1 = e(V)
local b1 : display %9.5f b1[1,1]
local v1 : display %9.5f v1[1,1]

merlin (y rcs(time, df(3) orthog) trt, family(poisson)) , evaltype(gf1)
mat b2 = e(b)
mat v2 = e(V)
local b2 : display %9.5f b2[1,1]
local v2 : display %9.5f v2[1,1]

merlin (y rcs(time, df(3) orthog) trt, family(poisson)) , evaltype(gf2)
mat b3 = e(b)
mat v3 = e(V)
local b3 : display %9.5f b3[1,1]
local v3 : display %9.5f v3[1,1]

assert `b3'==`b2'
assert `v3'==`v2'
