*! version 1.0.0 MJC

local TS 	transmorphic scalar
local RS 	real scalar
local RC 	real colvector
local RM 	real matrix
local Pgml 	pointer(struct merlin_struct scalar) scalar
local gml 	struct merlin_struct scalar

version 14.2
mata:

/*
- one function for complex predictor
- one function for each distributional ancillary parameter
*/

`RM' merlin_hessian_ob(`gml' gml)
{
	index = 1
	for (j=1;j<=gml.Nmodels;j++) {
		
		gml.model = gml.modtoind = j
		nb1 = asarray(gml.NHbs,j)
		
		// different model VCV - block 0
		k=1
		while (k<j) {
			
			nj = nk = 0
			nb2 = asarray(gml.NHbs,k)
			for (p=1;p<=gml.NHeqns[j];p++) nj = nj + nb1[p]
			for (p=1;p<=gml.NHeqns[k];p++) nk = nk + nb2[p]
			null = J(nk,nj,0)
			k++
		}
	
		//model specific VCV
		for (p1=1;p1<=gml.NHeqns[j];p1++) {
			
			p2 		= 1
			Href	= J(nb1[p1],0,.)

			while (p2<p1) {

				d2 		= quadcolsum((*gml.hessP[index++])(gml),1)
				Hsub 	= J(nb1[p1],nb1[p2],0)

				el = 1				
				for (e1=1;e1<=nb1[p1];e1++) {
					for (e2=1;e2<=nb1[p2];e2++) {
						Hsub[e1,e2] = d2[el++]
					}					
				} 
				Href = Href,Hsub
				p2++
			}

			//p1==p2
			d2 		= quadcolsum((*gml.hessP[index++])(gml),1)
			Hsub 	= J(nb1[p1],nb1[p2],0)	
			el 		= 1
			for (e1=1;e1<=nb1[p1];e1++) {
				e2 = 1
				while (e2<=e1) {
					if (e1==e2) Hsub[e1,e1] = d2[el++]
					else 		Hsub[e2,e1] = Hsub[e1,e2] = d2[el++]
					e2++
				}
			
			}

			if (p1==1 & p2==1) H = Hsub
			else H = H,Href'\Href,Hsub

		}

		if (k==1 & j==1) mainH = H
		else mainH = mainH,null\null',H
		
	}
	return(mainH)
		
}

end
