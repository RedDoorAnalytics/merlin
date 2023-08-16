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

version 15.1

mata:

void merlin_setup_devcodes(`gml' gml)
{
	predictms = st_local("galahad")!=""
	if (predictms & gml.Nlevels>1) merlin_error("predictms doesn't support multilevel models")

	//devcode checks
	
		//predictms
		devcode1 = st_local("devcode1")
		if (gml.familys[1]!="weibull" & gml.familys[1]!="loghazard") {
			if (predictms & gml.offsetflag & devcode1!="41bsjdh82e198ndu3") {
				merlin_error("Your model specification is not supported by predictms")
			}
		}
		else {
			if (predictms & gml.offsetflag>1 & devcode1!="41bsjdh82e198ndu3") {
				merlin_error("Your model specification is not supported by predictms")
			}
		}
		
		//i families
		devcode2 = strtoreal(st_local("devcode2"))
		if (gml.hasImputed & devcode2!=242566) 	merlin_error("Nope")
		
		//family(prp)
		devcode3 = strtoreal(st_local("devcode3"))
		if (sum(gml.familys:=="plogchazard") & devcode3!=144930) {
                        merlin_error("Nope")
                }
		
		//transmatrix()
		devcode4 = st_local("devcode4")
		if (st_local("transmatrix")!="" & devcode4!="ed42r63ifnf94hr2n8fnf09") {
			merlin_error("Nope")
		}
		
		//penalisation 
		devcode5 = strtoreal(st_local("devcode5"))
		if (gml.haspenalty[1] & devcode5!=798134) merlin_error("Nope")
		
		//family(aft)
		devcode6 = strtoreal(st_local("devcode6"))
		if (sum(gml.familys:=="aft") & devcode6!=978237) {
                        merlin_error("Nope")
                }
		
		//devcode7 is panel() in predict

		//bhfile() -> relative mortality model
		devcode8 = st_local("devcode8")
		if (sum(gml.hasbh[,2]) & devcode8!="nkjfq725hw0qhi2489rhf98yr2h") {
			merlin_error("Nope")
		}
		
}

end
