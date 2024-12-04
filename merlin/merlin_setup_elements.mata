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

// 	prep for varname element
`RC' merlin_setup_var(`gml' gml, `RS' mod, `RS' i, `RS' Nels, `SS' el)
{
	asarray(gml.elvarname,(mod,i,Nels),el)
	varinfo = asarray_create("real",1)
	istvar	= strtrim(el)==gml.tvarnames[mod]
	ist0var	= strtrim(el)==gml.ltruncated[mod]
	asarray(varinfo,1,istvar\ist0var)
	asarray(varinfo,2,st_data(.,el,gml.touse))
	asarray(gml.elinfo,(mod,i,Nels),varinfo)
	return(asarray(varinfo,2))
}

// random effect
`RC' merlin_setup_re(`gml' gml, `RS' mod, `RS' i, `RS' Nels, `SS' el)
{
// 	name = substr(el,1,(strpos(el,"[")-1))										//get name of random effect
// 	info = J(2,1,0)																//get level and dimension index
// 	for (l=1;l<=gml.Nlevels;l++) {
// 		if (sum(name:==asarray(gml.latlevs,l))) {
// 			info[1] = l															//level
// 		}
// 	}
// 	info[2] = (1..gml.Nres[info[1]]) * (name:==asarray(gml.latlevs,info[1]))	//dimensionindex
// //!!
// 	asarray(gml.elinfo,(mod,i,1),info)										//store info
// //!!
	return(J(gml.N /*obs[gml.Nlevels,gml.model]*/,1,1))
}

// 	prep for rcs() element
// 	- stores appropriate info 
`RM' merlin_setup_rcs(`gml' gml, `RS' mod, `RS' i, `RS' Nels, `SS' el)
{
	pos1 = strpos(el,"(")
	pos2 = strlen(el)
	synt = substr(el,pos1+1,pos2-pos1-1)

	//parse rcs() syntax

	if (strpos(synt,",")) 	stata("local 0 "+synt)
	else 					stata("local 0 , "+synt)
	stata("syntax anything, [DF(string) KNOTS(numlist asc) LOG ORTHog EVent OFFset(varname) MOFFset(varname)]")

	anything	= st_local("anything")
	rcsvar 		= st_data(.,anything,gml.touse)
	
	islog 		= st_local("log")!=""
	orth 		= st_local("orthog")!=""

	//offset
	hasoffset 	= st_local("offset")!=""
	hasmoffset	= st_local("moffset")!=""
	offsetist0var 	= .
	moffsetist0var 	= .

	if (hasoffset) {
		offsetvar	= strtrim(st_local("offset"))
		offsetist0var 	= offsetvar==gml.ltruncated[mod]
		if (!offsetist0var) {
                        offset = st_data(.,st_local("offset"),gml.touse)
                }
		else offset = asarray(gml.y,mod)[,3]
	}
	if (hasmoffset) {
		moffsetvar	= strtrim(st_local("moffset"))
		moffsetist0var 	= moffsetvar==gml.ltruncated[mod]
		if (!moffsetist0var) {
                        moffset = -st_data(.,st_local("moffset"),gml.touse)
                }
		else 	moffset = -asarray(gml.y,mod)[,3]
	}
	
	istvar	= strtrim(anything)==gml.tvarnames[mod]
	ist0var	= strtrim(anything)==gml.ltruncated[mod]

	//get element info
	
	if (gml.predict) { //called from predict
		index 	= strofreal(mod)+"_"+strofreal(i)+"_"+strofreal(Nels)
		knots 	= strtoreal(tokens(st_global("e(knots_"+index+")")))
		if (orth) rmat 	= st_matrix("e(rmat_"+index+")")
	}
	else {
		
		modindex = asarray(gml.xbindex,mod)
		
		if (st_local("df")!="") {
			df = strtoreal(st_local("df"))	
			
			tv = rcsvar[modindex]
			if (hasoffset) 		tv = tv :+ offset[modindex]
			if (hasmoffset) 	tv = tv :+ moffset[modindex]
			if (st_local("event")!="")  {
				gml.survind = 0
                                if (st_local("reffailure"+strofreal(mod))!="") {
                                        refmod = strtoreal(st_local("reffailure"+strofreal(mod)))
                                        tv = select(tv,merlin_util_depvar_mod(gml,refmod)[,2])
                                }
                                else {                            
                                        tv = select(tv,merlin_util_depvar_mod(gml,mod)[,2])
                                }
			}
			if (islog) tv = log(tv)

			tv	= sort(tv,1)
			nrows 	= rows(tv)
			if (df==1) 			index = 1\nrows
			else {
				if (df==2) 		index = 50
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
			if (orth) {
				tv = rcsvar[modindex]
				if (hasoffset)  tv = tv :+ offset[modindex]
				if (hasmoffset) tv = tv :+ moffset[modindex]
				if (islog) 	tv = log(tv)
				rmat = merlin_orthog(merlin_rcs(tv,knots)) 
			}
		}
		else {
			knots = strtoreal(tokens(st_local("knots")))
			if (orth) {
				tv = rcsvar[modindex]
				if (hasoffset)  tv = tv :+ offset[modindex]
				if (hasmoffset) tv = tv :+ moffset[modindex]
				if (islog)	rmat = merlin_orthog(merlin_rcs(log(tv),knots)) 
				else 		rmat = merlin_orthog(merlin_rcs(tv,knots)) 
			}
		}
	
	}

        //check knots
        if (cols(knots)!=rows(uniqrows(knots'))) {
                errprintf("Knot locations not unique\n")
                exit(198)
        }
        
	//store element info
	info = asarray_create("real",1)
	asarray(info,2,rcsvar)				//store input var
	asarray(info,3,knots)				//knots
	asarray(info,4,islog)				//log time or time
	asarray(info,5,orth)				//orthog or not
	if (orth) asarray(info,6,rmat)			//orthog matrix
	asarray(info,7,(hasoffset\offsetist0var))	//has offset
	if (hasoffset & !offsetist0var) {
		asarray(info,8,offset)			//offset
	}
	asarray(info,9,(hasmoffset\moffsetist0var))	//has moffset
	if (hasmoffset & !moffsetist0var) {
		asarray(info,10,moffset)		//moffset
	}
		
	//store info
	asarray(gml.elinfo,(mod,i,Nels),info)
	
	if (hasoffset) 	rcsvar = rcsvar :+ offset
	if (hasmoffset) rcsvar = rcsvar :+ moffset
	if (orth) {
		if (islog)      return(merlin_rcs(log(rcsvar),knots,0,rmat))
		else 		return(merlin_rcs(rcsvar,knots,0,rmat))
	}
	else {
		if (islog) 	return(merlin_rcs(log(rcsvar),knots,0))
		else 		return(merlin_rcs(rcsvar,knots,0))
	}
}

`RM' merlin_setup_fp(`gml' gml, `RS' mod, `RS' i, `RS' Nels, `SS' el)
{
	pos1 = strpos(el,"(")
	pos2 = strlen(el)
	synt = substr(el,pos1+1,pos2-pos1-1)
	
	//parse fp() syntax
	
	stata("local 0 "+synt)
	stata("syntax varlist(max=1) , POWers(numlist max=2) [OFFset(varname) MOFFset(varname)]")
	
	varname			= st_local("varlist")
	fpvar			= st_data(.,varname,gml.touse)
	powers 			= strtoreal(tokens(st_local("powers")))
	Nfp 			= cols(powers)
	
	//offset
	hasoffset 		= st_local("offset")!=""
	hasmoffset		= st_local("moffset")!=""
	offsetist0var 	= .
	moffsetist0var 	= .

	if (hasoffset) {
		offsetvar	= strtrim(st_local("offset"))
		offsetist0var 	= offsetvar==gml.ltruncated[mod]
		if (!offsetist0var) {
                        offset = st_data(.,st_local("offset"),gml.touse)
                }
		else    offset = asarray(gml.y,mod)[,3]
	}
	if (hasmoffset) {
		moffsetvar	= strtrim(st_local("moffset"))
		moffsetist0var 	= moffsetvar==gml.ltruncated[mod]
		if (!moffsetist0var) 	moffset = -st_data(.,st_local("moffset"),gml.touse)					
		else 			moffset = -asarray(gml.y,mod)[,3]
	}
	
	istvar			= varname==gml.tvarnames[mod]
	ist0var			= varname==gml.ltruncated[mod]
	
	//store element info
	fpinfo = asarray_create("real",1)
	asarray(fpinfo,1,istvar\ist0var)
	if (!istvar & !ist0var) asarray(fpinfo,2,fpvar)	//only need to store if not timevar()
	asarray(fpinfo,3,powers)		
	asarray(fpinfo,4,Nfp)
	
	asarray(fpinfo,5,(hasoffset\offsetist0var))	//has offset
	if (hasoffset & !offsetist0var) {
		asarray(fpinfo,7,offset)		//offset
	}
	asarray(fpinfo,6,(hasmoffset\moffsetist0var))	//has moffset
	if (hasmoffset & !moffsetist0var) {
		asarray(fpinfo,8,moffset)		//moffset
	}
	
	asarray(gml.elinfo,(mod,i,Nels),fpinfo)		
	if (hasoffset) 	fpvar = fpvar :+ offset
	if (hasmoffset) fpvar = fpvar :+ moffset
	return(merlin_fp(fpvar,powers))
}

/*
	prep for bs() element
	- stores appropriate info 
*/

`RM' merlin_setup_bs(`gml' gml, `RS' mod, `RS' i, `RS' Nels, `SS' el)
{

	pos1 = strpos(el,"(")
	pos2 = strlen(el)
	synt = substr(el,pos1+1,pos2-pos1-1)
	
	//parse bs() syntax
	
	if (strpos(synt,",")) 	stata("local 0 "+synt)
	else 			stata("local 0 , "+synt)
	stata("syntax varlist(max=1) , [Order(string) Degree(string) DF(string) Knots(numlist asc) BKnots(numlist asc max=2 min=2) LOG EVent INTercept OFFset(varname) MOFFset(varname)]")
	
	varname			= st_local("varlist")
	rcsvar			= st_data(.,varname,gml.touse)
	islog 			= st_local("log")!=""
	hasint 			= st_local("intercept")!=""
	
	//offset
	hasoffset 		= st_local("offset")!=""
	hasmoffset		= st_local("moffset")!=""
	offsetist0var 	= .
	moffsetist0var 	= .

	if (hasoffset) {
		offsetvar		= strtrim(st_local("offset"))
		offsetist0var 	= offsetvar==gml.ltruncated[mod]
		if (!offsetist0var) offset = st_data(.,st_local("offset"),gml.touse)
		else 				offset = asarray(gml.y,mod)[,3]
	}
	if (hasmoffset) {
		moffsetvar		= strtrim(st_local("moffset"))
		moffsetist0var 	= moffsetvar==gml.ltruncated[mod]
		if (!moffsetist0var) 	moffset = -st_data(.,st_local("moffset"),gml.touse)					
		else 					moffset = -asarray(gml.y,mod)[,3]
	}
	
	istvar			= varname==gml.tvarnames[mod]
	ist0var			= varname==gml.ltruncated[mod]
	
	//get element info
	
	if (gml.predict) { 													//called from predict
		index 			= strofreal(mod)+"_"+strofreal(i)+"_"+strofreal(Nels)
		knots 			= strtoreal(tokens(st_global("e(knots_"+index+")")))
		degree			= strtoreal(tokens(st_global("e(degree_"+index+")")))
		Nbasis			= strtoreal(tokens(st_global("e(Nbasis_"+index+")")))
	}
	else {
		
		modindex = asarray(gml.xbindex,mod)

		if (st_local("order")!="") 	degree = strtoreal(st_local("order"))-1
		else 						degree = strtoreal(st_local("degree"))
		ord = degree + 1
		
		if (st_local("bknots")=="") {
			tv = rcsvar[modindex]
			if (hasoffset)  tv = tv :+ offset[modindex]
			if (hasmoffset) tv = tv :+ moffset[modindex]
			if (islog) tv = log(tv)
			bknots = minmax(tv)
		}
		else bknots = strtoreal(tokens(st_local("bknots")))
		
		if (st_local("knots")!="") {
			knots = strtoreal(tokens(st_local("knots")))
			Niknots = cols(knots)
			knots = J(1,ord,bknots[1,1]),knots,J(1,ord,bknots[1,2])		
		}
		else if (st_local("df")!="") {
			
			df = strtoreal(st_local("df"))	
			if (df==1) 			{
				knots = J(1,ord,bknots[1,1]),J(1,ord,bknots[1,2])		
				Niknots = 0
			}
			else {
				tv = rcsvar[modindex]
				if (hasoffset)  tv = tv :+ offset[modindex]
				if (hasmoffset) tv = tv :+ moffset[modindex]
				if (st_local("event")!="") {
					tv = select(tv,asarray(gml.y,mod)[,2])
				} 
				if (islog) tv = log(tv)
				
				tv	 	= sort(tv,1)
				nrows 	= rows(tv)
			
				if (df==2) 		index = 50
				else if (df==3) index = 33.3333333333\66.66666666
				else if (df==4) index = 25\50\75
				else if (df==5) index = 20\40\60\80
				else if (df==6) index = 17\33\50\67\83
				else if (df==7) index = 14\29\43\57\71\86
				else if (df==8) index = 12.5\25\37.5\50\62.5\75\87.5
				else if (df==9) index = 11.1\22.2\33.3\44.4\55.6\66.7\77.8\88.9
				else if (df==10) index = 10\20\30\40\50\60\70\80\90 
				index = round(index :/100 :* nrows)
				knots = tv[index]'
				Niknots = cols(knots)
				knots = J(1,ord,bknots[1,1]),knots,J(1,ord,bknots[1,2])
			}			
		}
		else {
			knots = J(1,ord,bknots[1,1]),J(1,ord,bknots[1,2])		
			Niknots = 0
		}
		
		Nbasis = Niknots + degree + 1 
		
	}

	//store element info
	info = asarray_create("real",1)
	asarray(info,1,istvar\ist0var)
	if (!istvar & !ist0var) asarray(info,2,rcsvar)		//only need to store if not timevar()
	asarray(info,3,knots)								//knots
	asarray(info,4,degree)								//degree
	asarray(info,5,islog)								//log time or time
	
	asarray(info,6,(hasoffset\offsetist0var))				//has offset
	if (hasoffset & !offsetist0var) {
		asarray(info,7,offset)								//offset
	}
	asarray(info,11,(hasmoffset\moffsetist0var))				//has moffset
	if (hasmoffset & !moffsetist0var) {
		asarray(info,12,moffset)								//moffset
	}
	
	asarray(info,8,Nbasis)								//N basis functions
	asarray(info,9,st_varformat(st_local("varlist")))	//to post knots so equality check matches for upper boundary 
	asarray(info,10,hasint)								//has intercept
	
	asarray(gml.elinfo,(mod,i,Nels),info)
	
	if (hasoffset) 	rcsvar = rcsvar :+ offset
	if (hasmoffset) rcsvar = rcsvar :+ moffset
	if (islog) 	return(merlin_bs(log(rcsvar),knots,degree,Nbasis,hasint))
	else 		return(merlin_bs(rcsvar,knots,degree,Nbasis,hasint))
}

// 	prep for pc() element
// 	- stores appropriate info 
`RM' merlin_setup_pc(`gml' gml, `RS' mod, `RS' i, `RS' Nels, `SS' el)
{
	pos1 = strpos(el,"(")
	pos2 = strlen(el)
	synt = substr(el,pos1+1,pos2-pos1-1)

	//parse pc() syntax
	if (strpos(synt,",")) 	stata("local 0 "+synt)
	else 					stata("local 0 , "+synt)
	stata("syntax anything, KNOTS(numlist asc min=1) [OFFset(varname) MOFFset(varname) NOREFerence]")

	varname			= st_local("anything")
	pcvar 			= st_data(.,varname,gml.touse)
	//offset
	hasoffset 		= st_local("offset")!=""
	hasmoffset		= st_local("moffset")!=""
	offsetist0var 	= .
	moffsetist0var 	= .
	
	if (hasoffset) {
		offsetvar		= strtrim(st_local("offset"))
		offsetist0var 	= offsetvar==gml.ltruncated[mod]
		if (!offsetist0var) offset = st_data(.,st_local("offset"),gml.touse)
		else 				offset = asarray(gml.y,mod)[,3]
	}
	if (hasmoffset) {
		moffsetvar		= strtrim(st_local("moffset"))
		moffsetist0var 	= moffsetvar==gml.ltruncated[mod]
		if (!moffsetist0var) 	moffset = -st_data(.,st_local("moffset"),gml.touse)					
		else 					moffset = -asarray(gml.y,mod)[,3]
	}

	istvar			= varname==gml.tvarnames[mod]
	ist0var			= varname==gml.ltruncated[mod]

	//get element info
	knots = strtoreal(tokens(st_local("knots")))

	//store element info
	info = asarray_create("real",1)
	asarray(info,1,istvar\ist0var)
	if (!istvar & !ist0var) asarray(info,2,pcvar)		//only need to store if not timevar()
	asarray(info,3,knots)								//knots
	asarray(info,4,(hasoffset\offsetist0var))				//has offset
	if (hasoffset & !offsetist0var) {
		asarray(info,5,offset)								//offset
	}
	asarray(info,6,(hasmoffset\moffsetist0var))				//has moffset
	if (hasmoffset & !moffsetist0var) {
		asarray(info,7,moffset)								//moffset
	}
	asarray(info,8,st_local("noreference")!="")
	
	//store info
	asarray(gml.elinfo,(mod,i,Nels),info)
	
	if (hasoffset) 	pcvar = pcvar :+ offset
	if (hasmoffset) pcvar = pcvar :+ moffset
	return(merlin_pc(pcvar,knots,asarray(info,8)))
}

/*
	prep for mf() element
	- stores appropriate info 
*/

`RC' merlin_setup_mf(`gml' gml, `RS' mod, `RS' i, `RS' Nels, `SS' el)
{
	start 	= strpos(el,"(") + 1												//get function name
	len 	= strpos(el,")") - start
	fname 	= substr(el,start,len)

	stata("mata: pf = &"+fname+"()")											//get pointer
	external pf
	asarray(gml.elinfo,(mod,i,Nels),pf)											//store info
	return(merlin_xz_mf(gml,i,Nels))
}

/*
	prep for ?EV element
	- stores appropriate info 
*/

`RC' merlin_setup_expval(`gml' gml, `RS' mod, `RS' i, `RS' Nels, `SS' el)
{	
	info 	= asarray_create("real",1)
	kpos 	= strpos(el,"[")	//get response varname
	kpos2 	= strpos(el,"]")
	y 	= substr(el,kpos+1,kpos2-kpos-1)
	if (strpos(y,",")) {
		opts = substr(y,strpos(y,",")+1,kpos2)
		y = substr(y,1,strpos(y,",")-1) //strip off , options if there
		//get time()
		stata("local 0 , "+opts)
		stata("syntax , [Time(string) IMPUTE]")
		if (st_local("time")!="") {
			time = strtoreal(st_local("time"))
			asarray(info,2,1) //flag
			asarray(info,3,time)
		}
	}
	else asarray(info,2,0) //time() flag

	for (k=1;k<=gml.Nmodels;k++) {	//get response varname index and store
		if (st_local("response"+strofreal(k))==y | strofreal(k)==y) {
			asarray(info,1,k) //response index
		}
	}
	asarray(gml.elinfo,(mod,i,Nels),info) //store info
	return(J(gml.N /*obs[gml.Nlevels,gml.model]*/,1,1))
}

end
