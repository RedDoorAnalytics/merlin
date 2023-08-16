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

void merlin_setup_latents(`gml' gml)
{
	merlin_parse_latents(gml)
	merlin_init_vcvs(gml)
	merlin_parse_covstructures(gml)
	merlin_parse_vcv_eqns(gml)
	if (gml.Nrelevels) {
		merlin_init_integration(gml)
		if (gml.hasmargltrunc) merlin_init_ip_ltrunc(gml)
		st_local("mlprolog","derivprolog(merlin_prolog())")
	}
}

void merlin_parse_latents(`gml' gml)
{
	//notes; not effected by @'s
	//	     not affected by EV[]

	//unique latents at each level
	gml.latlevs 	= asarray_create("real",1)
	gml.Nres 	= J(gml.Nlevels,1,0)
	
	for (i=1;i<=gml.Nlevels;i++) {
		lats = J(0,1,"")
		//go through everything for each level
		for (j=1;j<=gml.Nmodels;j++) {		
			depvars = merlin_get_indepvars(j)
			Ndv 	= cols(depvars)
			for (k=1;k<=Ndv;k++) {
				dv 	= strtrim(depvars[1,k])
				pos = 1	
				while (pos) {
					pos = strpos(dv,"#")
					if (pos) {
						dv2 = substr(dv,1,pos-1)
						dv 	= substr(dv,pos+1,.)
					}
					else dv2 = dv
					posid = strpos(dv2,gml.levelvars[2,i])
					if (posid) lats = lats\substr(dv2,1,posid-1)
				}	
			}
		}
		
		if (rows(lats)) {
			//they could appear in multiple models and multiple times in same model
			lats 	= uniqrows(lats)			
			//now check that if a, b, c is specified, then its core random effect is also specified
			Nlats 	= rows(lats)
			len 	= strlen(lats)
			last 	= substr(lats,len,.)
			for (p=1;p<=Nlats;p++) {
				if (strtoreal(last[p])==.) {
					//strip off a, b to match
					core = substr(lats[p],1,len[p]-1)	
					if (!sum(core:==lats)) merlin_error("random effect "+lats[p]+" requires "+core)
				}
			}
			
			//now go through elements, expanding any for design matrices
			latsup = J(0,1,"")
			for (j=1;j<=gml.Nmodels;j++) {		
				depvars 	= merlin_get_indepvars(j)
				Ncmp 		= cols(depvars)
				cmpxindex 	= asarray(gml.CmpXBIndex,(j,1))
				for (c=1;c<=Ncmp;c++) {
					dv 	= strtrim(depvars[1,c])
					el  = 1
					pos = 1	
					while (pos) {
						pos = strpos(dv,"#")
						if (pos) {
							dv2 = substr(dv,1,pos-1)
							dv 	= substr(dv,pos+1,.)
						}
						else dv2 = dv
						posid = strpos(dv2,gml.levelvars[2,i])
						if (posid) {							
							Nvars = cmpxindex[c,2] :- cmpxindex[c,1] :+ 1
							if (Nvars>1) {
								for (let=1;let<=Nvars;let++) {
									latsup = latsup\(substr(dv2,1,posid-1)+merlin_get_letter(let))
								}
							}
							else latsup = latsup\substr(dv2,1,posid-1)
							//store level
							//store re index(es)
							asarray(gml.elinfo,(j,c,el),i)
						}
						el++
					}	
				}
			}
			uniqlats 	= uniqrows(latsup)
			gml.Nres[i]	= rows(uniqlats)
			
			//now go through elements, storing level and random effect index 
			for (j=1;j<=gml.Nmodels;j++) {		
				depvars 	= merlin_get_indepvars(j)
				Ncmp 		= cols(depvars)
				cmpxindex 	= asarray(gml.CmpXBIndex,(j,1))
				for (c=1;c<=Ncmp;c++) {
					dv 	= strtrim(depvars[1,c])
					el  = 1
					pos = 1	
					while (pos) {
						pos = strpos(dv,"#")
						if (pos) {
							dv2 = substr(dv,1,pos-1)
							dv 	= substr(dv,pos+1,.)
						}
						else dv2 = dv
						posid = strpos(dv2,gml.levelvars[2,i])
						if (posid) {							
							Nvars 		= cmpxindex[c,2] :- cmpxindex[c,1] :+ 1
							ellatinfo 	= asarray_create("real",1)
							asarray(ellatinfo,1,i)								//store level
							reindex = J(0,1,.)
							if (Nvars>1) {
								for (let=1;let<=Nvars;let++) {
									lat = substr(dv2,1,posid-1)+merlin_get_letter(let)
									reindex = reindex\selectindex(lat:==uniqlats)
								}
							}
							else {
								lat = substr(dv2,1,posid-1)
								reindex = reindex\selectindex(lat:==uniqlats)
							}
							asarray(ellatinfo,2,reindex)	//store re index(es)
							asarray(gml.elinfo,(j,c,el),ellatinfo)
						}
						el++
					}	
				}
			}
			asarray(gml.latlevs,i,uniqlats)
		}
	}
}

end
