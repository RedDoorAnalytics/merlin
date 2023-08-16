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

void merlin_setup_survival(`gml' gml)
{	
	gml.survind	= 0
	gml.Nsurv 	= J(gml.Nmodels,7,0) 
	gml.surv_index 	= asarray_create("real",2)
	familys 	= gml.familys
	gml.NI 		= J(gml.Nmodels,1,0)
	
	for (mod=1; mod<=gml.Nmodels; mod++) {
		fam = familys[mod]
		if (fam=="pwexponential" | 
                        fam=="exponential" | fam=="weibull" | 
                        fam=="gompertz" | (fam=="user" & 
                        st_local("hazfunction"+strofreal(mod))!="") | 
                        (fam=="user" & 
                        st_local("loghazfunction"+strofreal(mod))!="") | 
                        (fam=="user" & 
                        st_local("failure"+strofreal(mod))!="")) {
			if (gml.istimedep[mod,1]) gml.NI[mod] = 1
		}
		if (fam=="user" & (st_local("hazfunction"+strofreal(mod))!="" | 
                        st_local("loghazfunction"+strofreal(mod))!="") & 
                        st_local("chazfunction"+strofreal(mod))=="") {
			gml.NI[mod] = 1
		}
		if (fam=="loghazard" | fam=="addhazard") {
			gml.NI[mod] = 1
		}
		if (gml.hasbh[mod,2]) gml.NI[mod] = 1
	}
	
	if (sum(gml.hasbh[,1])) gml.bhazards = asarray_create("real",1)
	if (sum(gml.hasbh[,2])) {
		printf("Setting up expected rates...")
		displayflush()
		gml.bhazards = asarray_create("real",2)
	}
	
	if (st_local("pchintpoints")!="") {
                gml.chip = strtoreal(st_local("pchintpoints"))
        }
	else {
		chipstr = st_local("chintpoints")
		if (chipstr=="") gml.chip = 30
		else 		 gml.chip = strtoreal(chipstr)
	}
	
	for (mod=1; mod<=gml.Nmodels; mod++) {
		
		fam = familys[mod]
		if (gml.failures[1,mod]!="") {
                        
                gml.model = gml.modtoind = mod
                y = merlin_util_depvar(gml)

                if (fam=="cox") {
                        
                        //store index for sorted unique event times
                        index = selectindex(y[,2]:==1)
                        index = sort((y[index,1],index),1)[,2]
                        asarray(gml.surv_index,(mod,1),index)
                        Ns = rows_cols(index)
                        if (Ns[1] & Ns[2]) gml.Nsurv[mod,1] = Ns[1]
                        
                }
                else {
                        
                        if (gml.hasbh[mod,2]) {
                                merlin_setup_check_expratevars(gml)
                                merlin_setup_get_expratefile(gml)
                        }
                        
                        //exactly observed events
                        //-> hazard function	
                        index 	= selectindex(y[,2]:==1)
                        Ns 	= rows_cols(index)
                        if (Ns[1] & Ns[2]) {
                                asarray(gml.surv_index,(mod,1),index)
                                gml.Nsurv[mod,1] = Ns[1]
                                
                                //expected hazard at event times
                                if (gml.hasbh[mod,1] | gml.hasbh[mod,2]) {
                                        merlin_setup_exph1(gml)
                                }
                        }
                        
                        //exactly observed events and/or right censoring
                        //-> survival function
                        index 	= selectindex(y[,2]:<2)
                        Ns 	= rows_cols(index)
                        if (Ns[1] & Ns[2]) {
                                asarray(gml.surv_index,(mod,2),index)
                                gml.Nsurv[mod,2] = Ns[1]
                                if (gml.hasbh[mod,2]) merlin_setup_exph2(gml)
                        }
                        
                        //interval censoring
                        //-> cdf function
                        //-> left interval handled separately in case of 0s
                        if (gml.haslint[mod]) {
                                index 	= selectindex(y[,2]:==2)
                                Ns 	= rows_cols(index)
                                if (Ns[1] & Ns[2]) {
                                        asarray(gml.surv_index,(mod,3),index)
                                        gml.Nsurv[mod,3] = Ns[1]
                                }
                                ind2 = 3
                                if (gml.hasltrunc[mod]) ind2 = 4
                                index 	= selectindex((y[,2]:==2) :* (y[,ind2]:>0))
                                Ns 	= rows_cols(index)
                                if (Ns[1] & Ns[2]) {
                                        asarray(gml.surv_index,(mod,5),index)
                                        gml.Nsurv[mod,5] = Ns[1]
                                }
                        }

                        //left truncation
                        //-> survival function
                        if (gml.hasltrunc[mod]) {
                                index 	= selectindex(y[,3]:>0)
                                Ns 	= rows_cols(index)
                                if (Ns[1] & Ns[2]) {
                                        asarray(gml.surv_index,(mod,4),index)
                                        gml.Nsurv[mod,4] = Ns[1]
                                        
                                }
                                if (!gml.Nsurv[mod,4]) {
                                        gml.hasltrunc[mod] = 0
                                }
                        }
                        
                        //right censored - needed for generalised gamma/log normal/ log logistic
                        if (gml.familys[mod]=="ggamma" | 
                                gml.familys[mod]=="lognormal" | 
                                gml.familys[mod]=="loglogistic") {
                                
                                index = selectindex(y[,2]:==0)
                                Ns = rows_cols(index)
                                if (Ns[1] & Ns[2]) {
                                        asarray(gml.surv_index,(mod,6),index)	
                                        gml.Nsurv[mod,6] = Ns[1]
                                }
                        }
                        
                }
		
		} //endif
	}
	
	if (sum(gml.hasbh[,2])) {
		printf("done\n")
		displayflush()
	}
	
	gml.hasanyltrunc = sum(gml.hasltrunc)   //update in case ltrunc 
                                                //spec'd but all 0
		
	if (gml.hasImputed) {
	
		gml.imputing 		= 1
		gml.NsurvImp 		= J(gml.Nmodels,6,0) 
		gml.surv_index_imp 	= asarray_create("real",2)
		
		for (i=1;i<=gml.Nmodels;i++) {	
			if (gml.failures[1,i]!="") {
				gml.model = gml.modtoind = i
				y = merlin_util_depvar(gml)

				if (gml.familys[i]=="cox") {
					
					//store index for sorted unique event times
					index = selectindex(y[,2]:==1)
					index = sort((y[index,1],index),1)[,2]
					asarray(gml.surv_index_imp,(i,1),index)
					gml.NsurvImp[i,1] = rows(index)

				}
				else {

					//exactly observed events
					//-> hazard function	
					index = selectindex(y[,2]:==1)
					asarray(gml.surv_index_imp,(i,1),index)
					gml.NsurvImp[i,1] = rows(index)

					//exactly observed events and/or right censoring
					//-> survival function
					index = selectindex(y[,2]:<2)
					asarray(gml.surv_index_imp,(i,2),index)
					gml.NsurvImp[i,2] = rows(index)

					//interval censoring
					//-> cdf function
					//-> left interval handled separately in case of 0s
					if (gml.haslint[i]) {
						index = selectindex(y[,2]:==2)
						asarray(gml.surv_index_imp,(i,3),index)
						gml.NsurvImp[i,3] = rows(index)
						ind2 = 3
						if (gml.hasltrunc[i]) ind2 = 4
						index = selectindex((y[,2]:==2) :* (y[,ind2]:>0))
						asarray(gml.surv_index_imp,(i,5),index)
						gml.NsurvImp[i,5] = rows(index)
					}

					//left truncation
					//-> survival function
					if (gml.hasltrunc[i]) {
						index = selectindex(y[,3]:>0)
						asarray(gml.surv_index_imp,(i,4),index)
						gml.NsurvImp[i,4] = rows(index)
					}
					
					//right censored - needed for generalised gamma etc
					if (gml.familys[i]=="ggamma" | gml.familys[i]=="lognormal" | gml.familys[i]=="loglogistic") {
						index = selectindex(y[,2]:==0)
						if (rows(index) & cols(index)) {
							asarray(gml.surv_index_imp,(i,6),index)
							gml.NsurvImp[i,6] = rows(index)
						}
						
					}
					
				}
			
			}
		}

		gml.imputing = 0
		
	}
	
}

end

