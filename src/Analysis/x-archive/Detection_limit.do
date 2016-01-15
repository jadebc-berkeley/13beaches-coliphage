*Get detection limits for coliphage assays

* ---------- Doheny -----------

use "/Users/jadebc/Dropbox/Coliphage/Data/Untouched/doheny_inddata.dta", clear

tabstat result if groupindex=="14FMC051DO071" & result!=-88, by(qualifier) stats(n min max mean)
tabstat result if groupindex=="16FMC081DO071" & result!=-88, by(qualifier) stats(n min max mean)
tabstat result if groupindex=="14FPC051DO071" & result!=-88, by(qualifier) stats(n min max mean)
tabstat result if groupindex=="16FPC081DO071" & result!=-88, by(qualifier) stats(n min max mean)


log10	result
3.0791812	1200
3.146128	1400
5.1760913	150000

* ---------- Malibu -----------

use "/Users/jadebc/Dropbox/Coliphage/Data/Untouched/malibu_inddata.dta", clear

tabstat result if groupindex=="16FMC081MA091" & result!=-88, by(qualifier) stats(n min max mean)
tabstat result if groupindex=="14FPC051MA091" & result!=-88, by(qualifier) stats(n min max mean)
tabstat result if groupindex=="16FPC081MA091" & result!=-88, by(qualifier) stats(n min max mean)
