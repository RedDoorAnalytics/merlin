
clear 
set seed 24575240
set obs 300
gen b21 = rnormal(0,0.5) //+  0.1*trt
gen sdre = sqrt(exp(b21))
gen b11 = rnormal()
gen id = _n
expand 5
bys id: gen time = _n-1
gen y = 0.5 + b11 + 0.1*time + rnormal(0,sdre)

mata:
real matrix lev1_logl(gml,| G,H)
{
     y        	= merlin_util_depvar(gml)               //response
     xb 		= merlin_util_xzb(gml)                  //main lin. pred.
     residxb 	= exp(merlin_util_xzb_mod(gml,2))       //lev1 lin. pred
	 logl		= lnnormalden(y,xb,sqrt(residxb))		//logl
     return(logl)
}
end

merlin 	(y time M1[id]@1, family(user, llfunction(lev1_logl)))  ///
		(M2[id]@1, family(null))
do ./cert/predictions.do
		
