
local hasres = "`e(levelvars)'"!=""

local ind = 1

forval mod = 1/`e(Nmodels)' {
	n di "model `mod'"
	cap drop _pv*
	cap drop _refs*
	cap drop _serefs*
	n di "`ind++'"
	predict _pv1, eta fixedonly outcome(`mod') ci

        n di "`ind++'"
	if `hasres' {
		cap predict _pv2, eta marginal  outcome(`mod')
		if _rc {
			di as txt "captured eta marginal error"
		}
	}
	else {
		predict _pv2, eta standardise outcome(`mod')
	}

	
	if "`e(failure`mod')'"!="" & "`e(llfunction`mod')'"=="" {
		n di "`ind++'"
		predict _pv3, survival fixedonly   outcome(`mod')
		n di "`ind++'"
		predict _pv4, hazard fixedonly  outcome(`mod')
		n di "`ind++'"
		predict _pv5, rmst fixedonly outcome(`mod')
		n di "`ind++'"
		predict _pv6, chazard fixedonly  outcome(`mod')
		n di "`ind++'"
		predict _pv7, cif fixedonly  outcome(`mod')
		n di "`ind++'"
		predict _pv9, timelost fixedonly  outcome(`mod')
		n di "`ind++'"
                
                predict _pv99, density fixedonly  outcome(`mod')
		n di "`ind++'"
		
		//marginal or standardise predictions
		if `hasres' {
			local predtype marginal
			
			predict _pv10, hazard outcome(`mod') `predtype' 
			n di "`ind++'"
			predict _pv11, rmst outcome(`mod') `predtype' 
			n di "`ind++'"
			predict _pv12, chazard outcome(`mod') `predtype' 
			n di "`ind++'"
			predict _pv13, cif outcome(`mod') `predtype' 
			n di "`ind++'"
			predict _pv14, rmst outcome(`mod') `predtype' 
			n di "`ind++'"
			predict _pv15, timelost outcome(`mod') `predtype' 
			n di "`ind++'"
			predict _pv16, survival outcome(`mod') `predtype' 
			n di "`ind++'"
                        predict _pv96, density outcome(`mod') `predtype' 
			n di "`ind++'"
			if `e(Nlevels)'==2 & "`e(re_dist1)'"=="normal" {
				predict _refs`mod'*, reffects
				predict _serefs`mod'*, reses
			}
		}
		else {
			local predtype standardise
			
			predict _pv11, rmst outcome(`mod') `predtype' 
			n di "`ind++'"
			predict _pv13, cif outcome(`mod') `predtype' 
			n di "`ind++'"
			predict _pv14, rmst outcome(`mod') `predtype' 
			n di "`ind++'"
			predict _pv15, timelost outcome(`mod') `predtype' 
			n di "`ind++'"
			predict _pv16, survival outcome(`mod') `predtype' 
			n di "`ind++'"
                        predict _pv996, density outcome(`mod') `predtype' 
			n di "`ind++'"
		}
		
		
		//check timevar()
		
			cap range tvar 0 10 30
			cap gen t0 = tvar*0.01
			predict _pv17, survival fixedonly  timevar(tvar) outcome(`mod')
			n di "`ind++'"
			predict _pv33, survival fixedonly  timevar(tvar) outcome(`mod') ltruncated(t0)
			n di "`ind++'"
			predict _pv18, hazard fixedonly timevar(tvar) outcome(`mod')
			n di "`ind++'"
			predict _pv19, rmst fixedonly timevar(tvar) outcome(`mod')
			n di "`ind++'"
			predict _pv20, chazard fixedonly timevar(tvar) outcome(`mod')
			n di "`ind++'"
			predict _pv21, eta fixedonly timevar(tvar) outcome(`mod')
			n di "`ind++'"
			predict _pv22, cif fixedonly timevar(tvar) outcome(`mod')
			n di "`ind++'"
			predict _pv23, rmst fixedonly timevar(tvar) outcome(`mod')
			n di "`ind++'"
			predict _pv24, timelost fixedonly timevar(tvar) outcome(`mod')
			n di "`ind++'"
                        predict _pv9924, density fixedonly timevar(tvar) outcome(`mod')
			n di "`ind++'"
		
			//marginal or standardise
			if `hasres' {
				local predtype marginal
				
				predict _pv25, hazard timevar(tvar) outcome(`mod') `predtype' 
				n di "`ind++'"
				predict _pv26, rmst timevar(tvar) outcome(`mod') `predtype' 
				n di "`ind++'"
				predict _pv27, chazard timevar(tvar) outcome(`mod') `predtype' 
				n di "`ind++'"
				cap predict _pv28, eta timevar(tvar) outcome(`mod') `predtype' 
				n di "`ind++'"
				predict _pv29, cif timevar(tvar) outcome(`mod') `predtype' 
				n di "`ind++'"
				predict _pv30, rmst timevar(tvar) outcome(`mod') `predtype' 
				n di "`ind++'"
				predict _pv31, timelost timevar(tvar) outcome(`mod') `predtype' 
				n di "`ind++'"
				predict _pv32, survival timevar(tvar) outcome(`mod') `predtype' 
				n di "`ind++'"
                                predict _pv32a, density timevar(tvar) outcome(`mod') `predtype' 
				n di "`ind++'"
			}
			else {
				local predtype standardise
				
				predict _pv26, rmst timevar(tvar) outcome(`mod') `predtype' 
				n di "`ind++'"
				cap predict _pv28, eta timevar(tvar) outcome(`mod') `predtype' 
				n di "`ind++'"
				predict _pv29, cif timevar(tvar) outcome(`mod') `predtype' 
				n di "`ind++'"
				predict _pv30, rmst timevar(tvar) outcome(`mod') `predtype' 
				n di "`ind++'"
				predict _pv31, timelost timevar(tvar) outcome(`mod') `predtype' 
				n di "`ind++'"
				predict _pv32, survival timevar(tvar) outcome(`mod') `predtype' 
				n di "`ind++'"
                                predict _pv32a, density timevar(tvar) outcome(`mod') `predtype' 
				n di "`ind++'"
				predict _pv34, survival timevar(tvar) outcome(`mod') `predtype' ltruncated(t0)
				n di "`ind++'"
			}
		
	}

}
