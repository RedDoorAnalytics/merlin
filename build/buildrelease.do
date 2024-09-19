//build new version of merlin
//
// --> run whole do file

//!!
// BUILD in 17
//!!

local drive /Users/Michael
cd "`drive'/My Drive/software/merlin/"

local includemata       = 0

//=======================================================================================================================//

//build new release -> current version up is 2.3.0
local newversion 2_4_5
if `includemata' {
        local newversion `newversion'_mata
}
cap mkdir ./release/version_`newversion'
local fdir ./release/version_`newversion'/

//=======================================================================================================================//

//pkg files
copy ./build/merlin_details.txt `fdir', replace
	
//=======================================================================================================================//

//merlin

	copy ./merlin/merlin.ado `fdir', replace
	copy ./merlin/merlin_parse.ado `fdir', replace
	copy ./merlin/merlin_p.ado `fdir', replace

	//help files
	copy ./merlin/merlin.sthlp `fdir', replace
	copy ./merlin/merlin_estimation.sthlp `fdir', replace
	copy ./merlin/merlin_model_options.sthlp `fdir', replace
	copy ./merlin/merlin_models.sthlp `fdir', replace
	copy ./merlin/merlin_user.sthlp `fdir', replace
	copy ./merlin/merlin_reporting.sthlp `fdir', replace
	copy ./merlin/merlin_postestimation.sthlp `fdir', replace
	
//exptorcs

	copy ./merlin/exptorcs.ado `fdir', replace
	copy ./merlin/exptorcs.sthlp `fdir', replace


//mlib
cap erase `fdir'lmerlin.mlib
copy ./lmerlin.mlib `fdir', replace
