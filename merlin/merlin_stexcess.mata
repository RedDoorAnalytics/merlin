*! version 1.0.0 ?????2016

local gml 	struct merlin_struct scalar
local pgml	pointer(struct merlin_struct scalar) scalar
local TR 	transmorphic
local RS 	real scalar
local RC 	real colvector
local SS 	string scalar
local PS 	pointer scalar
local RR 	real rowvector
local RM 	real matrix
local PC 	pointer colvector
local PM 	pointer matrix
local SC 	string colvector

version 15.1

mata:

real matrix merlin_stexcess_logh(`gml' gml, real matrix t)
{
	haz_expect = exp(merlin_util_xzb(gml,t))
        haz_excess = exp(merlin_util_xzb_mod(gml,2,t))
	return(log(haz_expect :+ 
                        gml.indicator[merlin_get_index(gml)] :* 
                        haz_excess))
}

end
