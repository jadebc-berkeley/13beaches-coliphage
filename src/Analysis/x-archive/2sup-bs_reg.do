*---------------------------------------------
* Program      : bs_reg.do
* Programmer   : Ben Arnold
* Date         : 12 Mar 2010
* Description  :
/*
This program is a "lean" version of the original
backward deletion algorithm that Tim Wade wrote
for the beaches analyses: "backselect_reg.ado". 


It select variables for adjusted analyses using 
a minimum change in estimate criteria.

This code assumes a single exposure variable
It returns two scalar macros:
r(covarlist) : string list of selected covariates
r(sparse)    : equal to 1 if there were problems
               with the selection regressions 
               (collinearity or non-convergence)
               
               
Updates
30 Dec 2011 : commented out printing to reduce log output

*/
*---------------------------------------------

*--------------------------------------------
* short routine to identify collinear results
* as an argument, requires the "e(rules)" 
* matrix that is stored after running a logit model
* returns two scalar macros:
* r(Nvar_dropped): num of variables dropped
* r(Nobs_dropped): num of observations dropped
*--------------------------------------------
capture mata: mata drop iscollinear()
version 9
mata:
void iscollinear(string scalar r)
{
	real matrix R, Rsum, nV
	R = st_matrix(r)
	Rsum = colsum(R)
	if (Rsum[,1]>0) nV = rows(R)
	else nV = 0
	st_numscalar("r(Nvar_dropped)",nV)
	st_numscalar("r(Nobs_dropped)",Rsum[,4])
}
end

*---------------------------------------------
* bs_reg
*---------------------------------------------
capture program drop bs_reg
program bs_reg, rclass
syntax  [if] [in], model(string) outcome(string) exposure(string) [keepvar(string)] [covar(string)] [criteria(real 0.10)] [unit(real 1)] [options(string)]
	* model    : model call (e.g., -logit- -regress-)
	* outcome  : outcome of interest (assumed binary)
	* swimvar  : swimming exposure level
	* keepvar  : a string of variables besides swim exposure to force into models (not subject to model selection)
	* covar    : a string of variables that are candidates for removal
	* unit     : units of change on the log10 scale in level of indicator exposure
	* criteria : change in estimate criteria used by the algorithm to remove variables from the specification
	* options  : options to pass to the model call
	

	display in white "`outcome' `exposure' "
	
	set more off
	marksample touse
	
	local bswarn "Backward selection impossible for `outcome' : `exposure'. Data are too sparse!"
	local dropwarn "Warning: at least 1 candidate covariate had collinearity problems"
	local ndropped = 0
	
	local covar2 "`covar'"
	
	if index("`covar2'", "i.")!=0 {
		local covar2=subinstr("`covar2'", "i.", "", .) 
		}
	
	capture xi: `model' `outcome' `exposure' `keepvar' `covar' if `touse', `options'
	*if (_rc==0) mata: iscollinear("e(rules)")
	if (_rc!=0) {
		noisily di in red _n "`bswarn'"
		noisily di in red _n "Error code:" _rc
		return local sparse    = 1
		return local covarlist = ""
		exit
	}
	if (_rc==0 & r(Nvar_dropped)>0 ) local ndropped= r(Nvar_dropped)
	capture lincom `unit'*`exposure' 
	if (_rc!=0) {
		noisily di in red _n "`bswarn'"
		noisily di in red _n "Error code:" _rc
		return local sparse    = 1
		return local covarlist = ""
		exit
	}
	local full=r(estimate)
	global sparse = 0
	
	local minest=0
	while `minest'<`criteria' {

		if ltrim("`covar'")=="" {
			noisily di in yel "No variables meet criteria"
			if (`ndropped'>0) di in red "`dropwarn', N vars dropped: `ndropped'"
			return local sparse    = 0
			return local covarlist = ""
			exit 		
		}
	
		local i = 1
		foreach var in `covar'{
			local varlist: subinstr local covar "`var'" ""
			*display in yellow "-`var'"
		
			capture xi: `model' `outcome' `exposure' `keepvar' `varlist' if `touse', `options'
			*if (_rc==0) mata: iscollinear("e(rules)")
			if (_rc!=0) {
				noisily di in red _n "`bswarn'"
				noisily di in red _n "Error code:" _rc
				return local sparse    = 1
				return local covarlist = ""
				exit
			}
			if (_rc==0 & r(Nvar_dropped)>0 ) local ndropped= r(Nvar_dropped) 
			capture lincom `unit'*`exposure' 
			if (_rc!=0) {
				noisily di in red _n "`bswarn'"
				noisily di in red _n "Error code:" _rc
				return local sparse    = 1
				return local covarlist = ""
				exit
			}
			local est = r(estimate)
			local chg = abs((`full'-`est')/(`full'))
			if (`i'==1) {
				mata: n = "`var'"
				mata: x = `chg'
			}
			else {
				mata: n = n \ "`var'"
				mata: x = x \ `chg'			
			}
			local i = `i'+1
	
		}
	
		mata: o = order(x,1)
		mata: st_global("r(minname)",n[o][1])
		mata: st_numscalar("r(minest)",x[o][1])
	
		local minname = r(minname)
		local minest = r(minest)
		mata: (n[o], strofreal(x[o]))
		
		*di in white "REMOVE `minname'"
		local oldlist "`covar'"
		local covar: subinstr local covar "`minname'" ""
		*di in white "new varlist is `covar'"
		
	}
	
	noisily di "Variables that meet criteria include: `oldlist'"
	*if (`ndropped'>0) noisily di in red "`dropwarn', N vars dropped: `ndropped'"
	return local sparse    = 0
	return local covarlist = "`oldlist'"

noi di as res _n "--------------------------------------------" _n "Program bs_reg complete" _n  "--------------------------------------------" 

end bs_reg




