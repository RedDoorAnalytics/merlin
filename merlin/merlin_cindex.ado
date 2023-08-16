*! version 1.0.0 14may2018 MJC

//based on PR's stcstat2

program define merlin_cindex, rclass sortpreserve
	version 14.2
	
	syntax	, [outcome(string)]
	
	if "`outcome'"=="" {
		local outcome 1
	}
	
	if "`e(failure`outcome')'"=="" {
		di as error "outcome(`outcome') not a survival model"
		exit 198
	}
	
	local t : word 1 of `e(response`outcome')'
	local d : word 2 of `e(response`outcome')'
	
	//need linear predictor for outcome
	tempvar h Dv
	predict `h' if e(sample), eta outcome(`outcome')
	
	marksample touse
	markout `touse' `h'
	
	sort `touse' `h'
	qui count if `touse'

	local obs = r(N)

	quietly {
		local D 0
		local N 0   /* N = # as expected; on output we will call
		N E */
		local T 0
		local i = _N - `obs' + 1
		while `i' < _N {
			local j = `i' + 1
				
			gen byte `Dv' = `d'[`i'] & `d' in `j'/l

			replace `Dv' = 2 /*
			*/ if (!`Dv') & `d'[`i'] & `t'[`i']<=`t' in `j'/l
			replace `Dv' = 3 /* 
			*/ if (!`Dv') & `d' & `t'[`i']>=`t' in `j'/l
			replace `Dv'=0 /*
			*/ if (abs(`t'-`t'[`i'])<1e-12)& (`Dv'==1) in `j'/l

			count if `Dv' in `j'/l
			local D = `D' + r(N)

			count if `Dv' & `h'[`i']==`h' in `j'/l
			local T = `T' + r(N)

			count if `Dv'==1 & `h'[`i']!=`h' & `t'[`i']>`t' /*
			*/ in `j'/l
			local N = `N' + r(N)

			count if `Dv'==3 & `h'[`i']!=`h' & `t'[`i']>=`t' /*
			*/ in `j'/l
			local N = `N' + r(N)

			drop `Dv'
			local i = `i' + 1
		}
	}
	
	ret scalar C = (`N'+`T'/2)/`D'
	
end

version 14.2

mata:
void merlin_cindex()
{
	//linear predictor
	st_view(eta=.,.,st_local("eta"),st_local("touse"))
	
	//survival times
	y = tokens(st_global("e(response"+st_local("outcome")+")"))
	st_view(stime=.,.,y[1],st_local("touse"))
	st_view(event=.,.,y[2],st_local("touse"))
	
	N = rows(eta)
	
	num = 0
	denom = 0
	
	for (i=1;i<=N;i++) {
		for (j=1;j<=N;j++) {
			if (stime[i]>stime[j] & eta[i]>eta[j] & event[j]) num++
			if (stime[i]>stime[j] & event[j]) denom++
		}
	}

	st_numscalar(st_local("cindex"),num/denom)
	
}
	
end

