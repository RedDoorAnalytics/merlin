
local gmlStruct struct GML_struct scalar
local TR transmorphic
local RS real scalar
local RC real colvector
local SS string scalar
local PS pointer scalar
local RR real rowvector
local RM real matrix
local PC pointer colvector
local PM pointer matrix
local SC string colvector

mata:



transmorphic scalar gml_init()
{
	`gmlStruct' S
	
	gml_Nlevels(S)
	gml_panelinfo(S)
	gml_Nres_at_levels(S)
	
	return(S)
}

function gml_Nmodels(`gmlStruct' S, | `RS' Nmods)
{
	if (args()==1) return
	S.Nmodels = Nmods
}

function gml_Nobs(`gmlStruct' S, | `RS' kobs)
{
	if (args()==1) return
	S.Nobs = kobs
}

function gml_Nlevels(`gmlStruct' S, | `RS' k)
{
	if (args()==1) return
	S.Nlevels = k
}

function gml_panelinfo(`gmlStruct' S, | `SS' levelvarstub, `SS' tousevarstub)
{
	if (args()==1) {
		S.panelindexes = asarray_create("real",1)
		return
	}
	if (S.Nlevels==0) {
		errprintf("Nlevels = 0\n")
		error(1986)
	}
	S.Npanels = J(S.Nlevels,1,.)
	for (i=1;i<=S.Nlevels;i++) {
		pansetup = panelsetup(st_data(.,st_local(levelvarstub+strofreal(i)),st_local(tousevarstub+strofreal(i+1))),1)
		asarray(S.panelindexes,i,pansetup)
		S.Npanels[i] = panelstats(pansetup)[1]
	}
}

function gml_Nres_at_levels(`gmlStruct' S, | `RC' Nreslevs)
{
	if (args()==1) return
	if (rows(Nreslevs)!=S.Nlevels) {
		errprintf("Check dimensions\n")
		exit(1986)
	}
	S.Nres = Nreslevs
}

function gml_p_logl(`gmlStruct' S, | `PC' Plnlfunc)
{
	if (args()==1) return
	S.Plnl = Plnlfunc
}

function gml_init_lnli(`gmlStruct' S)
{
	S.lnfi1 = J(S.Npanels[1],1,.)
	S.lhd = asarray_create("real",S.Nlevels)
}

function gml_vcv_mats(`gmlStruct' S)
{
	S.vcvs = asarray_create("real",1)
	for (i=1;i<=S.Nlevels;i++) {
		asarray(S.vcvs,i,I(S.Nres[i]))
	}
}

function gml_vcv_structures(`gmlStruct' S, | `SS' covs)
{
	S.covariances = tokens(covs):=="independent"\tokens(covs):=="exchangeable"\tokens(covs):=="unstructured"
}



function gml_todo(`gmlStruct' S, | `RS' tod)
{
	S.todo = tod
}

function gml_init_Z(`gmlStruct' S, | `SS' Zstub, `SS' touse)
{
	S.Z = S.Zb = asarray_create("real",1)
	if (sum(S.covariances[1,])>0 & S.todo>0) S.Zb_re = asarray_create("real",2)
	if (sum(S.covariances[3,])>0 & S.todo>0) S.Zb_dre = asarray_create("real",2)
	if (sum(S.covariances[3,])>0 & S.todo>1) S.Zb_d2re = asarray_create("real",3)
	for (i=1;i<=S.Nlevels;i++) {
		asarray(S.Z,i,st_data(.,tokens(st_local(Zstub+strofreal(i))),touse))
		asarray(S.Zb,i,J(S.Nobs,S.ndim[i],.))
	}	
}

function gml_adaptiveGH(`gmlStruct' S, | `RS' ad)
{
	S.adapt = ad
}

function gml_NmleqnsFEs(`gmlStruct' S, | `RS' Nmleqns)
{
	S.NmleqnsFEs = Nmleqns
}

function gml_init_X(`gmlStruct' S, | `SS' Xstub, `SS' touse)
{
	S.Ncoefeq = J(S.NmleqnsFEs,1,.)
	for (i=1;i<=S.NmleqnsFEs;i++) {
		S.Ncoefeq[i] = cols(tokens(st_local(Xstub+strofreal(i))))
	}
	S.Ncoefs = sum(S.Ncoefeq)
	//if (S.todo==0) return

	S.X = asarray_create("real",1)
	for (i=1;i<=S.NmleqnsFEs;i++) {
		asarray(S.X,i,st_data(.,tokens(st_local(Xstub+strofreal(i))),touse))
	}
}

function gml_P_scoreFE(`gmlStruct' S, | `PC' Pscorefuncs)
{
	if (args()==1 | S.todo==0) {
		return
	}
	S.PscoreFEs = Pscorefuncs
}

function gml_init_G(`gmlStruct' S)
{
	if (S.todo==0) return
	S.G = asarray_create("real",S.Nlevels+2)		//store score across nodes at lower level, each equation, each coefficient
	S.Gres = asarray_create("real",S.Nlevels+2)		//+ level + re at that level
}

function gml_P_hessian(`gmlStruct' S, | `PM' Phesfuncs)
{
	if (args()==1 | S.todo<2) {
		return
	}
	S.Phessian = Phesfuncs
}

function gml_init_xb(`gmlStruct' S)
{
	S.xb = J(S.Nobs,S.NmleqnsFEs,.)
}

function gml_init_y(`gmlStruct' S, | `RM' y)
{
	if (args()==1) return
	S.y = y
}

function gml_bhazard(`gmlStruct' S, | `RC' bhaz)
{
	S.bhazard = bhaz
}

function gml_P_scoreSigma(`gmlStruct' S, | `PC' Pscorefunc)
{
	if (args()==1 | S.todo==0) {
		return
	}
	S.PscoreSigma = Pscorefunc
}

function gml_chainrule(`gmlStruct' S, | `RC' usecr)
{
	if (args()==1) S.chainrule = J(S.NmleqnsFEs,1,1)
	else S.chainrule = usecr
}
end
