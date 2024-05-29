*! version 1.0.0 ?????2016

local gml 		struct merlin_struct scalar
local pgml		pointer(struct merlin_struct scalar) scalar
local TR 		transmorphic
local RS 		real scalar
local RC 		real colvector
local SS 		string scalar
local PS 		pointer scalar
local RR 		real rowvector
local RM 		real matrix
local PC 		pointer colvector
local PM 		pointer matrix
local SC 		string colvector

version 14.1

mata:

`RM' merlin_logl_ologit(`gml' gml,| 	///
                         `RM' G, 	///
                         `RM' H)
{
	mod 	= gml.model
	z 	= merlin_util_xzb(gml)
	cuts 	= asarray(gml.distancb,(mod,1))

	//handle 1 and ncuts separately	
	oindex		= asarray(gml.OrdIndexes,(mod,4))
	oneindex 	= asarray(gml.OrdIndexes,(mod,1))
	ncutsindex 	= asarray(gml.OrdIndexes,(mod,3))
	rest 		= asarray(gml.OrdIndexes,(mod,2))

	z[rest,] 	= 1:/ (1 :+ exp(-cuts[oindex[rest,1]] :+ z[rest,])) :- 1:/ (1 :+ exp(-cuts[oindex[rest,2]] :+ z[rest,]))
	z[oneindex,] 	= 1:/ (1 :+ exp(-cuts[oindex[oneindex,1]] :+ z[oneindex,]))
	z[ncutsindex,] 	= 1 :- 1:/ (1 :+ exp(-cuts[oindex[ncutsindex,2]] :+ z[ncutsindex,]))
	
	if (gml.todo==0) return(log(z))
}

`RM' merlin_logl_oprobit(`gml' gml,| 	///
                         `RM' G, 	///
                         `RM' H)
{
	mod 	= gml.model
	z 	= merlin_util_xzb(gml)
	cuts 	= asarray(gml.distancb,(mod,1))
	
	//handle 1 and ncuts separately	
	oindex		= asarray(gml.OrdIndexes,(mod,4))
	oneindex 	= asarray(gml.OrdIndexes,(mod,1))
	ncutsindex 	= asarray(gml.OrdIndexes,(mod,3))
	rest 		= asarray(gml.OrdIndexes,(mod,2))

	z[rest,] 	= normal(cuts[oindex[rest,1]] :- z[rest,]) :- normal(cuts[oindex[rest,2]] :- z[rest,])
	z[oneindex,] 	= normal(cuts[oindex[oneindex,1]] :- z[oneindex,])
	z[ncutsindex,] 	= 1 :- normal(cuts[oindex[ncutsindex,2]] :- z[ncutsindex,])
	
	if (gml.todo==0) return(log(z))
}

end

