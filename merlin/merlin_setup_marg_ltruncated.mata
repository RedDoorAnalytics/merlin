*! version 1.0.0 ?????2016

local gml 	struct merlin_struct scalar
local pgml	pointer(struct merlin_struct scalar) scalar
local Egml	struct merlin_ereturn_struct scalar
local RS 	real scalar
local SS 	string scalar
local PS 	pointer scalar
local RM 	real matrix
local SM	string matrix
local PC 	pointer colvector
local SC 	string colvector
local SR	string rowvector
local TR	transmorphic
local RC	real colvector
local RR	real rowvector
local PM	pointer matrix

version 15.1

mata:

void merlin_setup_marg_ltruncated(`gml' gml)
{	
        gml.hasmargltrunc = 0
        if (!gml.Nrelevels) return

        gml.survind	= 0
	familys 	= gml.familys
	
	mltmodcheck = J(gml.Nmodels,1,0)
	for (mod=1; mod<=gml.Nmodels; mod++) {
		mltmodcheck[mod] = strtoreal(st_local("margltruncated"+strofreal(mod)))
	}
	
	if (sum(mltmodcheck)==0) return
	
	if (sum(mltmodcheck)>1) {
		errprintf("Only one model can have ltruncated(, marginal)\n")
		exit(198)
	}
	
	if (gml.Nrelevels>1) {
		errprintf("marginal left truncation not supported in models with more than 2 levels\n")
                exit(198)
	}
	
	mod = gml.mltmodel = selectindex(mltmodcheck:==1)

	// if no obs with t0>0 then bail out
	if (!gml.Nsurv[mod,4]) return
	
	gml.model = gml.modtoind = gml.mltmodel
	
	// no get the info
	
	y = merlin_util_depvar(gml)
	
	// if multilevel model, see if all obs within a cluster have _t0>0 and if so
	// index the earliest _t0 within each cluster to do the marginal survival calc.

	// level i , model j
	// asarray(gml.panelindexes,(i,j))
	index = J(0,1,.)
	iindex = asarray(gml.panelindexes,(1,mod))
	mltllindex = J(gml.Nobs[1,mod],1,.)
	for (i=1;i<=gml.Nobs[1,mod];i++) {
		
		t0 = panelsubmatrix(y[,3],i,iindex)
		mint0 = min(t0)
		if (mint0>0) {
			mltllindex[i] = i
			//get index for this ob
			t0index = selectindex(t0:==mint0)
			if (rows(t0index)>1) {
				errprintf("multiple obs with same ltruncated() time\n")
				exit(1986)
			}
			index = index\(iindex[i,1]:+t0index:-1)
		}
	}
	gml.mltllindex = select(mltllindex,mltllindex:!=.)

	if (rows(index) & cols(index)) {
		gml.hasmargltrunc = 1
		gml.hasltrunc[mod] = 0			//override to marginal
		gml.hasanyltrunc = sum(gml.hasltrunc)
		asarray(gml.surv_index,(mod,7),index)
		gml.Nsurv[mod,7] = rows(index)
	}
	else return
	
        //warning
        if (gml.hasmargltrunc) {
                printf("note: marginal left truncation with random effects is only\n")
                printf("      valid if no gaps within a cluster or one \n")
                printf("      observation per cluster\n")
        }
        
}

end

