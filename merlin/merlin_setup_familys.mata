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

void merlin_setup_rp(`gml' gml)
{
	for (i=1;i<=gml.Nmodels;i++) {
		if (gml.familys[i]=="rp") {
			gml.model = gml.modtoind = i
			rcsopts = st_local("rcsopts"+strofreal(i))
			if (rcsopts!="") merlin_rcs_build(gml,rcsopts)
		}
	}
}

void merlin_rcs_build(`gml' gml, `SS' rcsopts)
{	
	gml.survind = 0
	y = merlin_util_depvar(gml)

	if (gml.predict) { 													//called from predict
		strk 		= strofreal(gml.model)
		knots 		= strtoreal(tokens(st_global("e(knots"+strk+")")))
		hasorthog 	= st_global("e(orthog"+strk+")")!=""
		if (hasorthog) {
			rmat 	= st_matrix("e(rcsrmat_"+strk+")")
		}
	}
	else {

		//parse options
		stata("local 0 , "+rcsopts)
		stata("syntax , [DF(string) Knots(string) NOORTHog]")
		
		knots 		= st_local("knots")
		df 		= st_local("df")
		hasorthog 	= st_local("noorthog")==""
		
		if (knots!="" & df!="") merlin_error("Can't specify both df() and knots()")
		
		if (df!="") {
			df 		= strtoreal(df)
			tv 		= select(y[,1],y[,2])
			tv 		= log(tv)
			tv	 	= sort(tv,1)
			nrows 	= rows(tv)
			if (df==1) 			index = 1\nrows
			else {
				if 		(df==2) index = 50
				else if (df==3) index = 33.3333333333\66.66666666
				else if (df==4) index = 25\50\75
				else if (df==5) index = 20\40\60\80
				else if (df==6) index = 17\33\50\67\83
				else if (df==7) index = 14\29\43\57\71\86
				else if (df==8) index = 12.5\25\37.5\50\62.5\75\87.5
				else if (df==9) index = 11.1\22.2\33.3\44.4\55.6\66.7\77.8\88.9
				else if (df==10) index = 10\20\30\40\50\60\70\80\90 
				index = 1\round(index :/100 :* nrows)\nrows
			}
			knots = tv[index]'
		}
		else knots = strtoreal(tokens(knots))

		if (hasorthog)	rmat 	= merlin_orthog(merlin_rcs(log(y[,1]),knots)) 
	}

	if (hasorthog) 	rcsvars = merlin_rcs(log(y[,1]),knots,0,rmat)
	else 		rcsvars = merlin_rcs(log(y[,1]),knots)

	//store stuff
	asarray(gml.distancb,(gml.model,2),rcsvars)
	asarray(gml.distancb,(gml.model,3),knots)
	asarray(gml.distancb,(gml.model,4),hasorthog)
	if (hasorthog) asarray(gml.distancb,(gml.model,5),rmat)

	if (gml.familys[gml.model]=="rp" | gml.familys[gml.model]=="prp") {
		if (hasorthog) 	drcsvars = merlin_rcs(log(y[,1]),knots,1,rmat)
		else 		drcsvars = merlin_rcs(log(y[,1]),knots,1)
		asarray(gml.distancb,(gml.model,6),drcsvars) 

		//if left trunacted, need baseline splines at ltruncated()
		if (gml.hasltrunc[gml.model]) {
			if (hasorthog) 	t0rcsvars = merlin_rcs(log(y[asarray(gml.surv_index,(gml.model,4)),3]),knots,0,rmat)
			else 		t0rcsvars = merlin_rcs(log(y[asarray(gml.surv_index,(gml.model,4)),3]),knots,0)
			asarray(gml.distancb,(gml.model,8),t0rcsvars) 
		}

		//if interval censoring, need baseline splines at left interval as well
		if (gml.haslint[gml.model]) {
			if (gml.hasltrunc[gml.model]) 	tcol = 4
			else 							tcol = 3
			if (hasorthog) 	l0rcsvars = merlin_rcs(log(y[asarray(gml.surv_index,(gml.model,5)),tcol]),knots,0,rmat)
			else 		l0rcsvars = merlin_rcs(log(y[asarray(gml.surv_index,(gml.model,5)),tcol]),knots,0)
			asarray(gml.distancb,(gml.model,7),l0rcsvars)
		}
	}

	//post to Stata and ml equation local
	stub = "_rcs"+strofreal(gml.model)+"_"
	stata("cap drop "+stub+"*")
	Nvars = cols(rcsvars)
	names = J(1,0,"")
	eqnames = J(1,0,"")
	for (r=1;r<=Nvars;r++) {
		names = names,(stub+strofreal(r))
	}

	if (!gml.predict) {
		id = st_addvar("double",names)
		st_store(.,id,gml.modeltouses[gml.model],rcsvars)
		if (!gml.nogen) printf("variables created: "+stub+"1 to "+stub+strofreal(Nvars)+"\n")
	}
	
	st_local("rcsvars"+strofreal(gml.model),invtokens(names))
	
}

void merlin_setup_ordinal(`gml' gml)
{
	if (sum(gml.familys:=="ordinal")) {
		gml.OrdIndexes = asarray_create("real",2)
		for (i=1;i<=gml.Nmodels;i++)  {
			if (gml.familys[i]=="ordinal") {
			
				//get Ndf from response var i.e. # of cut points
				y 		= asarray(gml.y,i)
				ks 		= uniqrows(y)
				nks 	= rows(ks)
				kindex 	= 1::nks
				oindex 	= J(gml.Nobs[gml.Nlevels,i],2,.)
				for (j=1;j<=gml.Nobs[gml.Nlevels,i];j++) {
					oindex[j,1] = selectindex(ks:==y[j])
				}
				oindex[,2] 		= oindex[,1] :- 1
				Ndf 			= rows(ks) - 1
				gml.Ndistancp[i] 	= Ndf
				
				asarray(gml.OrdIndexes,(i,4),oindex)
				asarray(gml.OrdIndexes,(i,1),selectindex(oindex[,1]:==1))
				asarray(gml.OrdIndexes,(i,2),selectindex((oindex[,1]:!=1) :* (oindex[,2]:!=Ndf)))
				asarray(gml.OrdIndexes,(i,3),selectindex(oindex[,2]:==Ndf))

			}
		}
	}
}

void merlin_setup_cox(`gml' gml)
{	
	gml.survind		= 0
	gml.nocox		= 1	
	gml.firth		= st_local("firth")!=""

	for (mod=1;mod<=gml.Nmodels;mod++) {	
	
		if (gml.failures[1,mod]!="") {
			gml.model = gml.modtoind = mod
			y = merlin_util_depvar(gml)

			if (gml.familys[mod]=="cox") {
				index = asarray(gml.surv_index,(mod,1))
				if (!gml.istimedep[mod,1] & !gml.hasltrunc[mod]) {
					ord = order(y[,1],-1)
					asarray(gml.surv_index,(mod,2),ord)		//index to arrange all times descending sorted order
					asarray(gml.surv_index,(mod,3),invorder(ord))	//to put back to original order
				}
				else {
					if (gml.nocox) gml.cox_index = asarray_create("real",3)
					if (gml.hasltrunc[mod]) {
						for (j=1;j<=gml.Nsurv[mod,1];j++) {
							atrindex = selectindex((y[index[j],1] :<= y[,1]) :& (y[index[j],1] :> y[,3]))
							asarray(gml.cox_index,(mod,1,j+6),rows(atrindex))
							asarray(gml.cox_index,(mod,2,j+6),atrindex)
						}
					}
					else {	
						for (j=1;j<=gml.Nsurv[mod,1];j++) {
							atrindex = selectindex(y[index[j],1] :<= y[,1])
							asarray(gml.cox_index,(mod,1,j+6),rows(atrindex))
							asarray(gml.cox_index,(mod,2,j+6),atrindex)
						}
					}
				}
				gml.nocox = 0
			}
		
		}
	}
}

void merlin_setup_pwexp(`gml' gml)
{
	for (i=1;i<=gml.Nmodels;i++) {
		if (gml.familys[i]=="pwexponential") {
			gml.model = gml.modtoind = i
			opts = st_local("rcsopts"+strofreal(i))
			merlin_pwexp_build(gml,opts)
		}
	}
}

void merlin_pwexp_build(`gml' gml, `SS' opts)
{	
	gml.survind = 0
	y = merlin_util_depvar(gml)

	if (gml.predict) { 													//called from predict
		strk 		= strofreal(gml.model)
		knots 		= strtoreal(tokens(st_global("e(knots"+strk+")")))
	}
	else {
		//parse options
		stata("local 0 , "+opts)
		stata("syntax , Knots(string)")
		knots = st_local("knots")
		knots = strtoreal(tokens(knots))
	}
	
	//store stuff
	Nknots = cols(knots)
	gml.Ndistancp[gml.model] = Nknots + 1
	asarray(gml.distancb,(gml.model,2),knots)
	asarray(gml.distancb,(gml.model,3),Nknots)
	
	gml.survind = 0
	baseindex 	= J(merlin_get_nobs(gml),1,.)
	index		= merlin_get_index(gml)
	
	//first interval
	ix 	= selectindex(y[index,1]:<knots[1])
	nrc = rows(ix)\cols(ix)
	if (nrc[1] & nrc[2]) baseindex[ix] = J(nrc[1],1,1)
	//mid
	if (Nknots>1) {
		for (i=1;i<Nknots;i++) {
			ix 	= selectindex((knots[i]:<=y[index,1]) :& (y[index,1]:<knots[i+1]))
			nrc = rows(ix)\cols(ix)
			if (nrc[1] & nrc[2]) baseindex[ix] = J(nrc[1],1,i+1)
		} 
	}
	//last interval
	ix 	= selectindex(knots[Nknots]:<=y[index,1])
	nrc = rows(ix)\cols(ix)
	if (nrc[1] & nrc[2]) baseindex[ix] = J(nrc[1],1,Nknots+1)
	asarray(gml.distancb,(gml.model,4),baseindex)
	
}

end
