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

/*
	prep for varname element
	- stores appropriate info 
*/

`RS' merlin_setup_var(`gml' gml, `RS' mod, `RS' i, `RS' Nels, `SS' el)
{
	asarray(gml.elinfo,(mod,i,Nels),st_data(.,el,gml.touse))
	return(1)
}

/*
	prep for random effect element
	- stores appropriate info 
*/

`RS' merlin_setup_re(`gml' gml, `RS' mod, `RS' i, `RS' Nels, `SS' el)
{
	name = substr(el,1,(strpos(el,"[")-1))										//get name of random effect
	
	info = J(2,1,0)																//get level and dimension index
	for (l=1;l<=gml.Nlevels;l++) {
		if (sum(name:==asarray(gml.latlevs,l))) {
			info[1] = l															//level
		}
	}
	info[2] = (1..gml.Nres[info[1]]) * (name:==asarray(gml.latlevs,info[1]))	//dimensionindex
	asarray(gml.elinfo,(mod,i,Nels),info)										//store info
	return(2)
}

/*
	prep for ?EV element
	- stores appropriate info 
*/

`RS' merlin_setup_expval(`gml' gml, `RS' mod, `RS' i, `RS' Nels, `SS' el)
{	
	kpos 	= strpos(el,"[")													//get response varname
	kpos2 	= strpos(el,"]")
	y 		= substr(el,kpos+1,kpos2-kpos-1)
	if (strpos(y,",")) y = substr(y,1,strpos(y,",")-1)							//strip off , options if there

	for (k=1;k<=gml.Nmodels;k++) {												//get response varname index and store
		if (st_local("response"+strofreal(k))==y | strofreal(k)==y) {
			asarray(gml.elinfo,(mod,i,Nels),k)									//store info
		}
	}

	if (substr(el,1,2)=="EV")		return(4)
	else if (strpos(el,"iEV")) 		return(5)
	else if (strpos(el,"dEV")) 		return(6)
	else if (strpos(el,"d2EV")) 	return(7)
}

/*
	prep for mf() element
	- stores appropriate info 
*/

`RS' merlin_setup_mf(`gml' gml, `RS' mod, `RS' i, `RS' Nels, `SS' el)
{
	start 	= strpos(el,"(") + 1												//get function name
	len 	= strpos(el,")") - start
	fname 	= substr(el,start,len)
	fname
	stata("mata: pf = &"+fname+"()")											//get pointer
	external pf
	asarray(gml.elinfo,(mod,i,Nels),pf)											//store info
	return(3)	
}

/*
	prep for fp() element
	- stores appropriate info 
*/

`RS' merlin_setup_fp(`gml' gml, `RS' mod, `RS' i, `RS' Nels, `SS' el)
{

	pos1 = strpos(el,"(")
	pos2 = strlen(el)
	synt = substr(el,pos1+1,pos2-pos1-1)
	
	//parse fp() syntax
	
	stata("local 0 "+synt)
	stata("syntax varlist(max=1) , POWers(numlist max=2) [OFFset(varname) MOFFset(varname)]")
	
	fpvar			= st_data(.,st_local("varlist"),gml.touse)
	powers 			= strtoreal(tokens(st_local("powers")))
	Nfp 			= cols(powers)
	hasoffset 		= st_local("offset")!=""
	hasmoffset		= st_local("moffset")!=""
	hasanyoffset 	= hasoffset | hasmoffset
	if (hasoffset) 	offset = st_data(.,st_local("offset"),gml.touse)
	if (hasmoffset) offset = -st_data(.,st_local("moffset"),gml.touse)
	istvar			= st_local("varlist")==st_local("timevar"+strofreal(mod))
	
	//store element info
	
	fpinfo = asarray_create("real",1)
	asarray(fpinfo,1,istvar)
	if (!istvar) asarray(fpinfo,2,fpvar)							//only need to store if not timevar()
	asarray(fpinfo,3,powers)		
	asarray(fpinfo,4,Nfp)
	asarray(fpinfo,5,hasanyoffset)
	if (hasanyoffset) asarray(fpinfo,6,offset)
	
	asarray(gml.elinfo,(mod,i,Nels),fpinfo)		
	
	return(9)
}

/*
	prep for rcs() element
	- stores appropriate info 
*/

`RS' merlin_setup_rcs(`gml' gml, `RS' mod, `RS' i, `RS' Nels, `SS' el)
{

	pos1 = strpos(el,"(")
	pos2 = strlen(el)
	synt = substr(el,pos1+1,pos2-pos1-1)
	
	//parse rcs() syntax
	
	if (strpos(synt,",")) 	stata("local 0 "+synt)
	else 					stata("local 0 , "+synt)
	stata("syntax varlist(max=1) , [DF(string) KNOTS(numlist asc) LOG ORTHog EVent OFFset(varname) MOFFset(varname)]")
	
	rcsvar			= st_data(.,st_local("varlist"),gml.touse)
	islog 			= st_local("log")!=""
	orth 			= st_local("orthog")!=""
	hasoffset 		= st_local("offset")!=""
	hasmoffset		= st_local("moffset")!=""
	hasanyoffset 	= hasoffset | hasmoffset
	if (hasoffset) 	offset = st_data(.,st_local("offset"),gml.touse)
	if (hasmoffset) offset = -st_data(.,st_local("moffset"),gml.touse)
	istvar			= st_local("varlist")==st_local("timevar"+strofreal(mod))
	
	//get element info
	
	if (gml.predict) { 													//called from predict
		index 			= strofreal(mod)+"_"+strofreal(i)+"_"+strofreal(Nels)
		knots 			= strtoreal(tokens(st_global("e(knots_"+index+")")))
		if (orth) rmat 	= st_matrix("e(rmat_"+index+")")
	}
	else {
		
		modindex = asarray(gml.xbindex,mod)
		
		if (st_local("df")!="") {
			df = strtoreal(st_local("df"))	
			tv = rcsvar[modindex]
			if (hasanyoffset) tv = tv :+ offset[modindex]
			if (st_local("event")!="") {
				tv = select(tv,asarray(gml.y,mod)[,2])
			} 
			if (islog) tv = log(tv)
			
			tv	 	= sort(tv,1)
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
				if (hasanyoffset) tv = tv :+ offset[modindex]
				if (islog) 	tv = log(tv)
				rmat = merlin_orthog(merlin_rcs(tv,knots)) 
			}
		}
		else {
			knots = strtoreal(tokens(st_local("knots")))
			if (orth) {
				tv = rcsvar[modindex]
				if (hasanyoffset) tv = tv :+ offset[modindex]
				if (islog)	rmat = merlin_orthog(merlin_rcs(log(tv),knots)) 
				else 		rmat = merlin_orthog(merlin_rcs(tv,knots)) 
			}
		}
	
	}
	
	//store element info
	
	rcsinfo = asarray_create("real",1)
	asarray(rcsinfo,1,istvar)
	if (!istvar) asarray(rcsinfo,2,rcsvar)			//only need to store if not timevar()
	asarray(rcsinfo,3,knots)						//knots
	asarray(rcsinfo,4,islog)						//log time or time
	asarray(rcsinfo,5,orth)							//orthog or not
	if (orth) asarray(rcsinfo,6,rmat)				//orthog matrix
	asarray(rcsinfo,7,hasanyoffset)					//has offset
	if (hasanyoffset) asarray(rcsinfo,8,offset)		//offset
	
	asarray(gml.elinfo,(mod,i,Nels),rcsinfo)

	return(8)										//rcs()

}

end
