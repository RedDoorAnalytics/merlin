*! version 1.0.0 ?????2016

local gml 		struct merlin_struct scalar
local pgml		pointer(struct merlin_struct scalar) scalar
local Egml		struct merlin_ereturn_struct scalar
local RS 		real scalar
local SS 		string scalar
local PS 		pointer scalar
local RM 		real matrix
local SM		string matrix
local PC 		pointer colvector
local SC 		string colvector
local SR		string rowvector
local TR		transmorphic
local RC		real colvector
local RR		real rowvector
local PM		pointer matrix

version 14.2

mata:

void merlin_error_checks(`gml' gml)
{
	for (i=1; i<=gml.Nmodels; i++) {
		
		if (st_local("timevar"+strofreal(i))!="" & (gml.familys[i]=="gamma")) {
			merlin_error("timevar() not currently supported with family(gamma)")
		}
	
	}
}
end
