** ADO FILE FOR EFFECTIVENESS INDICATORS STAND ALONE COMMANDS

** VERSION AND NOTES (changes between versions described under CHANGES)
** v1.0 29jun2019 
*! (beta version please report any bugs)

** CHANGES



cap program drop ceqtaxstar
program define ceqtaxstar, rclass 
	#delimit;
	syntax [if] [in] [pw aw iw fw/] [,
			/*Incomes*/
			startinc(varname)
			endinc(varname)
			tax(varname)
			]
		;	
		#delimit cr	
		
		if "`exp'" == "" {
			local exp = 1
			}
			
		cap drop  ___startinc
		gen double ___startinc = `startinc'
		cap drop ___tax
		gen double ___tax = `tax'
		gsort -___startinc

		replace ___tax = abs(___tax)
		qui sum ___tax [aw=`exp'] 
		local tot=r(sum) //total amount to redistribute
		qui sum ___startinc [aw=`exp'] 
		
		if (`tot' > r(sum) | `tot' == 0 ) {
			if `tot' > r(sum) return scalar t_gr = 1
			if `tot' == 0 return scalar t_0 = 1
			exit
		}
		else {
			// Taking the difference of income b/w one person and then previous
			gen double ___diff_y = ___startinc - ___startinc[_n-1]
			recode ___diff_y (.=0)
			gen double ___cum_w = sum(`exp')
			gen double ___diff_y_i = ___diff_y*___cum_w[_n-1]
			recode ___diff_y_i (.=0)
			gen double ___cum_diff_y_i = sum(___diff_y_i)
			// Determining who is taxed.
			gen ___taxed = (abs(___cum_diff_y_i) < `tot')
			gen ___last_taxed = ___taxed==1 & ___taxed[_n+1]==0

			gen ___id = _n
			sum ___id if ___last_taxed== 1
			assert r(min)==r(max) 
			local which = r(mean) // Giving observation of which person is taxed.
			// Generating optimal tax
			gen double ____taxstar = 0
			replace ____taxstar = ___startinc - ___startinc[`which'] if ___taxed==1

			local still_tax =  `tot' - abs(___cum_diff_y_i[`which'])
			
			local still_tax_per_person = `still_tax' / ___cum_w[`which']

			replace ____taxstar = ____taxstar + `still_tax_per_person' if ___taxed==1

			sum ____taxstar [aw=`exp']
			// Ensuring that we how allocating the exact amount of tax available
			*assert round(`tot',.1) == round(r(sum),.1)
			nois assert abs((`tot' - r(sum))/`tot') < 0.0001
			// Generating optimal income 
			cap drop ____id_taxstar
			cap drop ____ytaxstar
			gen double ____ytaxstar = ___startinc
			replace ____ytaxstar = ____ytaxstar - ____taxstar if ___taxed == 1
			gen ____id_taxstar = ___taxed

			return scalar twarn = 0

			drop ___taxed ___last_taxed ____taxstar ___diff_y ___diff_y_i ///
				 ___cum_diff_y_i ___cum_w ___id 
		}
		cap drop ___startinc ___taxes
		end
		
		
		
* Program to compute ideal transfer for impact effectiveness
*generates income variable (ybenstar) whith ideal transfers 
*rest of observations have no transfer income
*var benstar has the ideal transfers for those obs. rest is missing
cap program drop ceqbenstar
program define ceqbenstar, rclass 
	#delimit;
	syntax [if] [in] [pw aw iw fw/] [,
			/*Incomes*/
			startinc(varname)
			endinc(varname)
			ben(varname)
			]
		;	
			
		#delimit cr
		if "`exp'" == "" {
			local exp = 1
			}
			
		cap drop  ___startinc
		gen double ___startinc = `startinc'
		cap drop  ___ben
		gen double ___ben = `ben'
		sort  ___startinc
		
		qui sum ___ben [aw=`exp'] 
		local tot=r(sum) //total amount to redistribute
		qui sum ___startinc [aw=`exp'] 
		
		if (`tot' > r(sum) | `tot' == 0 ) {
			if `tot' > r(sum) return scalar b_gr = 1
			if `tot' == 0 return scalar b_0 = 1
			exit
		}
		else {
			// Taking the difference of income b/w one person and then previous
			gen double ___diff_y = ___startinc - ___startinc[_n-1]
			recode ___diff_y (.=0)
			gen double ___cum_w = sum(`exp')
			gen double ___diff_y_i = ___diff_y*___cum_w[_n-1]
			recode ___diff_y_i (.=0)
			gen double ___cum_diff_y_i = sum(___diff_y_i)
			// Determining who is taxed.
			gen ___gets_ben = (abs(___cum_diff_y_i) < `tot')
			gen ___last_rec = ___gets_ben==1 & ___gets_ben[_n+1]==0

			gen ___id = _n
			sum ___id if ___last_rec== 1
			assert r(min)==r(max) 
			local which = r(mean) // Giving observation of which person is taxed.

			gen double ____benstar = 0
			replace ____benstar = `startinc'[`which'] - `startinc'   if ___gets_ben==1

			local still_ben =  `tot' - abs(___cum_diff_y_i[`which'])
			local still_ben_per_person = `still_ben' / ___cum_w[`which']
			
			replace ____benstar = ____benstar + `still_ben_per_person' if ___gets_ben==1

			sum ____benstar [aw=`exp']
			// Ensuring that we how allocating the exact amount of benefits available
			*nois assert round(`tot',.1) ==  round(r(sum),.1)
			nois assert abs((`tot' - r(sum))/`tot') < 0.0001
			
			
			cap drop ____id_benstar
			cap drop ____ybenstar
			gen double ____ybenstar = `startinc'
			replace ____ybenstar = ____ybenstar + ____benstar if ___gets_ben == 1
			gen ____id_benstar = ___gets_ben

			return scalar bwarn = 0

			drop ___gets_ben ___last_rec ____benstar ___diff_y ___diff_y_i ///
				 ___cum_diff_y_i ___cum_w ___id 
			}

		drop  ___startinc ___ben

		end
		
* Program to compute harm tax formula for poverty impact effectiveness
*generates income variable (ytaxharm) whith ideal taxes 
*rest of observations have no tax income
*var taxharm has the harm tax for those obs. rest is missing

cap program drop ceqtaxharm
program define ceqtaxharm, rclass 
	#delimit;
	syntax [if] [in] [pw aw iw fw/] [,
			/*Incomes*/
			endinc(varname)
			tax(varname)
			]
		;	
			#delimit cr
			
		if "`exp'"!="" {
			local aw = "[aw = `exp']" //weights
			local pw = "[pw = `exp']"
		}
		else {
			tempvar ww
			gen `ww' = 1
			}
		*no taxes income
		cap drop ___ww
		gen double ___ww = `exp'
		cap drop ___notax
		gen double ___notax =`endinc'+abs(`tax')
		cap drop ___taxes
		gen double ___taxes =  `tax' 
		sort ___notax

		qui sum ___taxes `aw'
		local tot=r(sum) //total amount to redistribute

		qui sum ___notax `aw'

		if `tot' > r(sum) { 
			return scalar thwarn = 1
			exit
		}
		else {

			gen double ___inc_wght= ___notax*___ww
			gen double ___cum_inc_wght = sum(___inc_wght)
			gen ___taxed = ___cum_inc_wght < `tot'
			gen ___new_inc = ___notax
			replace ___new_inc = 0 if ___taxed == 1
			gen ___last_taxed = 1 if ___taxed == 1 & ___taxed[_n+1] == 0
			gen ___id = _n
			sum ___id if ___last_taxed == 1
			assert r(N) ==1
			local which = r(mean)

			scalar remainder = `tot' - ___cum_inc_wght[`which']
			assert remainder > 0 & remainder < ___cum_inc_wght[`which' + 1] 
			replace ___new_inc = ___notax - (remainder/___ww) in `=`which' + 1'
			cap drop ___taxstar

			gen double ___taxstar = ___notax - ___new_inc

			sum ___taxstar `aw'
			local rsum = r(sum)
			// Ensuring harm tax is equal to amount avaible to tax
			local tax1 = round(`tot',.1)
			local tax2 = round(`rsum',.1) 
			assert (abs((`tax1'-`tax2')/`tax1') < 0.00001)

			cap drop ____id_taxharm
			cap drop ____ytaxharm
			gen double ____id_taxharm=___taxed
			gen double ____ytaxharm=___new_inc
			return scalar thwarn = 0

			cap drop ___id ___notax ___inc_wght ___cum_inc_wght ___taxed ___last_taxed ___new_inc

		}
		cap drop  ___taxes ___notax ___ww
		end
		
				***Marginal contribution ID
	
	cap program drop _ceqmcid
program define _ceqmcid, rclass 
	#delimit;
	syntax [if] [in] [pw aw iw fw/] [,
			/*Incomes*/
			inc(varname)
			sptax(varname)
			spben(varname) 
			pline(string)  
			]
			
		;	
			#delimit cr
		if "`exp'" !="" {
			local aw = "[aw = `exp']" //weights
			local pw = "[pw = `exp']"
		}
		local id_tax=0
		local id_ben=0
		
		tempvar inter
		*See if we are dealing with taxes or transfers
		if wordcount("`sptax'")>0{
		local id_tax=1
		gen double `inter'=abs(`sptax')
		
		}
		if wordcount("`spben'")>0{
		local id_ben=1
		gen double `inter'=abs(`spben')
		}
		
		
		if `id_tax'==1{
		tempvar o_inc
		gen double `o_inc'=`inc'+`inter'
		}
		if `id_ben'==1{
		tempvar o_inc
		gen double `o_inc'=`inc'-`inter'
		}
		
		
			*gini final income
			covconc `inc' `pw'
			local g_f=r(gini)
			covconc `o_inc' `pw'
			local g_o=r(gini)
		
			local mc=`g_o'-`g_f'
			return scalar mc_ineq =  `mc' //Marginal Contribution
		
		if wordcount("`pline'")>0{
			tempvar pov0_o pov1_o pov2_o
			tempvar pov0_f pov1_f pov2_f
			gen `pov0_o'=(`o_inc'<`pline')
			gen `pov0_f'=(`inc'<`pline')
			
			qui gen `pov1_o' = max((`pline'-`o_inc')/`pline',0) // normalized povety gap of each individual
			qui gen `pov2_o' = `pov1_o'^2 
			qui gen `pov1_f' = max((`pline'-`inc')/`pline',0) // normalized povety gap of each individual
			qui gen `pov2_f' = `pov1_f'^2 
			forvalues f=0/2{
			qui sum `pov`f'_o' `aw'
			local p`f'_o=r(mean)
			qui sum `pov`f'_f' `aw'
			local p`f'_f=r(mean)
			local mc`f'=`p`f'_o'-`p`f'_f'
			return scalar mc_p`f'=`mc`f'' //Marginal Contribution of poverty fgt:`f'
			}
		}
		end
					
		
****Spending Effectiveness program
cap program drop _ceqspend
program define _ceqspend, rclass 
	#delimit;
	syntax [if] [in] [pw aw iw fw/] [,
			/*Incomes*/
			inc(varname)
			sptax(varname)
			spben(varname)
			]
		;	
			
		#delimit cr
		
		// note variable names represent column names in Effectives_SE Cal.xlsx in Methods folder of CEQStataPackage shared dropbox folder
		local id_tax=0
		local id_ben=0

		set type double 
		*See if we are dealing with taxes or transfers


		tempvar G
		qui gen double `G' = `inc'

		if wordcount("`sptax'")>0{

		local id_tax=1
		tempvar F E
		qui gen double `F'=abs(`sptax')
		qui gen double `E' = `G'+ `F'

		} 


		if wordcount("`spben'")>0{
		local id_ben=1
		tempvar F E
		qui gen double `F' =abs(`spben')
		qui gen double `E' = `G'-`F'
		
		}

		
	

		if `id_ben' == 1 {
			sort `E'  // Because we are dealing with benefits
			tempvar C D
			qui sum `exp' 
			qui gen double `C' = `exp'/r(sum) // normalizing the weights
			qui gen double `D' = sum(`C') // culmative sum of weights

			tempvar B
			qui gen int `B' = _n // generating observation number, ith person
			 
			tempvar S T 
			qui gen double `S' = `E'[_n+1] - `E' // diff in starting income between i and ith person
			qui gen double `T' = `S'*`D' in 1
			qui replace `T' = `S'*`D' + `T'[_n-1] if _n>1 // Sum of differnce in OI.

			tempvar F_hat w_F
			qui gen double `F_hat' = `D'[_n-1] + `C'/2 // empirical distribution / rankings
			qui replace `F_hat' = `C'/2 in 1
			qui gen double `w_F' = `C'*`F_hat' // weighted rankings
			qui summ  `w_F'
			scalar Fbar = r(sum) // average weighted rankings

			// Based on theory this should always be .5
			assert abs(Fbar-.5) <= 0.0001

		    tempvar W X 
		    qui gen double `W' = (`F_hat' - Fbar)*`C' // difference with avg. weighted rankings
			qui gen double `X' = `W' in 1
			qui replace  `X' = `W' + `X'[_n-1] if _n>1  // culmative sum of diff with avg weighted rankings

			tempvar Y Z 
			qui gen double `Y' = `X'*`S' // Multiplying X with difference in original income
			qui gen double `Z' = `Y' in 1
			qui replace `Z' = `Y' + `Z'[_n-1] if _n>1 // Culmative sum of difference times culmative sum of weight rankings difference

			tempvar AD 
			qui qui gen double `AD' = 2*`Z'

			qui covconc `E' [pw=`C']
			scalar Gini_OI = r(gini)

			qui summarize `E' [aw=`C'], meanonly
			scalar mu_OI = r(mean)

			tempvar Gini_star

			qui gen double `Gini_star' = (mu_OI*Gini_OI + `AD')/(mu_OI + `T')
			

			// suppose the observed ending income Gini is saved in scalar Gini_EI
			qui covconc `G' [pw=`C']
			scalar Gini_EI = r(gini)

			tempvar higher_gini
			qui gen byte `higher_gini' = (`Gini_star' < Gini_EI) & ///
								   (`Gini_star'[_n-1] > Gini_EI ) 
			


				// note higher_gini is column AG
				
			qui summ `B' if `higher_gini'==1
			assert r(min)==r(max) // just one obs
			local which = r(mean) 

			// Composition of s
			// signs change depending on tax or benefit
			scalar  AT2 = (Gini_OI*mu_OI + `AD') in `=`which'-1'
			
			scalar AV2 = (mu_OI + `T') in `=`which'-1'
			
			scalar AR2 = `D' in `which'

			scalar AX2 = 2*`X' in `which'

			*scalar s = (AT2 - Gini_EI*AV2)/(Gini_EI*AR2 - AX2) // scalar wasn't storing value
			local s = (AT2 - Gini_EI*AV2)/(Gini_EI*AR2 - AX2)
			
			tempvar AJ AK AZ BA 
			qui gen double `AJ' = (`E'[`which'] - `E')*(`B' <= `which')

			qui gen double `AK' = `E' + `AJ'

			qui gen double `AZ' = (`AJ' +`s')*(`B'<= `which')

			qui gen double `BA' = `AZ' + `E'

			covconc `BA' [pw=`C']
			 
			assert abs(r(gini) - Gini_EI)/Gini_EI < 0.0001


			qui summ `AZ' [aw=`C']

			scalar tot_EHB = r(sum)

			qui summ `F' [aw=`C']

			scalar tot_TB = r(sum)

			return scalar sp_ef = tot_EHB / tot_TB
			
		}

		if `id_tax' == 1 {
			*set trace on
			sort `E' // even though we are dealing with taxes the orderings are still intended to be from least to greatest, the reverse
			tempvar B // comes with a second wieght variable we create
			gen `B' = _n
			

			sum `exp'
			tempvar C D
			gen double `C' = `exp'/r(sum)
			gen double `D' = sum(`C')

			
			// Storing value for starting income Gini
			covconc `E' [pw=`C']
			scalar Gini_OI = r(gini)
			
			// suppose the observed ending income Gini is saved in scalar Gini_EI
			covconc `G' [pw=`C']
			scalar Gini_EI = r(gini)

			tempvar F_hat
			gen double `F_hat' = `D'[_n-1] + `C'/2
			replace `F_hat' = `C'/2 in 1
			tempvar w_F
			gen double `w_F' = `C' * `F_hat'
			summ `w_F'
			scalar Fbar = r(sum)

			assert abs(r(sum)-.5)/.5 < 0.00001
			tempvar W
			gen double `W' = (`F_hat' - Fbar)*`C'

			tempvar  S 
			gen double `S' =  `E' - `E'[_n-1] 

			gsort - `E'
			tempvar Drev T
			gen double `Drev' = sum(`C')
			gen double `T' = `S'*`Drev' in 1
			replace `T' = `S'*`Drev' + `T'[_n-1] if _n>1

			tempvar X 
			gen double `X' = `W' in 1
			replace `X' = `W' + `X'[_n-1] if _n>1

			tempvar Y Z
			gen double `Y' = `X'*`S'
			gen double `Z' = `Y' in 1
			replace `Z' = `Y' + `Z'[_n-1] if _n>1

			tempvar AD
			gen double `AD' = 2*`Z'

			// Resort based on starting income to calculate Ginis 

			sort `E'



			summarize `E' [aw=`C'], meanonly
			scalar mu_OI = r(mean)
			
			tempvar Gini_star
			gen double `Gini_star' = (mu_OI*Gini_OI - `AD')/(mu_OI - `T')
				// remember to change to negative signs in both numerator and denom
				//  when doing it for taxes


			tempvar higher_gini
		       gen byte `higher_gini' = (`Gini_star' < Gini_EI) & (`Gini_star'[_n+1] > Gini_EI)
				// note higher_gini is column AG
			*qui gen byte `higher_gini' = (`Gini_star' < Gini_EI) & ///
							(`Gini_star'[_n-1] > Gini_EI ) 

			summ `B' if `higher_gini'==1
			assert r(min)==r(max) // just one obs
			local which = r(mean)

			scalar  AT2 = (Gini_OI*mu_OI - `AD') in `=`which'+1'
			scalar AV2 = (mu_OI - `T') in `=`which'+1'
			scalar AR2 = `Drev' in `which'
			scalar AX2 = 2*`X' in `which'

			*scalar s = (AT2 - Gini_EI*AV2)/(AX2-Gini_EI*AR2) // scalar not returning value
			local s = (AT2 - Gini_EI*AV2)/(AX2-Gini_EI*AR2)			

			tempvar AJ AK AZ BA
			gen double `AJ' = (`E'[`which'] - `E')*(`B' >=`which')



			gen double `AK' = `E' + `AJ'

			gen double `AZ' = (`AJ' - `s')*(`B'>=`which')

			gen double `BA' = `AZ' + `E'

 			covconc `BA' [pw=`C']

			assert abs(r(gini) - Gini_EI)/Gini_EI < 0.0001


			summ `AZ' [aw=`C']

			scalar tot_EHB = -1*r(sum)

			summ `F' [aw=`C']

			scalar tot_TB = r(sum)

			return scalar se_ind = tot_EHB / tot_TB
			*set trace off
			
		}
		

		scalar drop tot_EHB tot_TB Gini_EI Gini_OI  ///
		  AR2 AX2 AT2 AV2 mu_OI Fbar
end

	
	****Poverty Spending effectiveness
* Program to computepoverty spending effectiveness for FGT1 and FGT2
*generates scalars  sp_ef_pov_1 and sp_ef_pov_2
cap program drop ceqbensp
program define ceqbensp, rclass 
	#delimit;
	syntax [if] [in] [pw aw iw fw/] [,
			/*Incomes*/
			startinc(varname)
			endinc(varname)
			ben(varname)
			zz(string)
			obj1(string)
			obj2(string)

			*
			]
		;
		#delimit cr
	
		
		scalar A = `zz'
		scalar B = _N
		tempvar D C
		gen double `D' = `exp'
		gen int `C' = _n

		tempvar G H I
		gen double `G' = `startinc'
		gen double `H' = `ben'
		gen double `I' = `G' + `H'

		if 	"`obj1'" != "" {
			// Poverty gap calculation
			tempvar opt_ben
			gen double `opt_ben' = min(`H',A-`G')*(`G'<A) // optimal benefits according to alis agorithm in Spending Effectiveness for the poverty indicators_Nov 24, 2017.

			sum `opt_ben' [aw=`D']
			scalar prime1 = r(sum)

			sum  `H' [aw=`D']
			scalar tot_ben1 = r(sum)

			return scalar sp_ef_pov_1 = prime1/tot_ben1
			scalar drop prime1 tot_ben1
		}

		if "`obj2'" != "" {
			sum `D'
			scalar E = r(sum)
			tempvar F
			gen double `F' = sum(`D') // Cumulative weight 

			tempvar J K L M
			gen double `J' = ((A - `G')^2/A^2)*(`G'<A) // SGR for OI

			gen double `L' =  ((A - `I')^2/A^2)*(`I'<A) // SGR EI

			sum `J' [aw=`D']
			scalar N = r(mean) // Squared poverty gap ratio for OI

			sum `L' [aw=`D'] // Squared poverty gap ratio for EI
			scalar O = r(mean)

			scalar P = N-O // Marg contribution between OI and EI
			local P = N-O  // Also using local because there have been issues with using scalar.
			

			scalar Q = A^2*E*P // Target Value 1

			tempvar R S T U V W X Y
			gen double  `R' = min(`G'[_n+1] - `G', A -`G')*(`G'<A) // Diff in OI bw I and I+1th person for poor

			gen double `S' = (A-`G')*(`G'<A) // z - y_i

			gen double `T' = `R'^2 // DIfference squared.

			gen double `U' = (2*`R'*`S') - `T' // 2*(diff)*(z-y_i) - (diff)^2

			gen double `V' = `U'*`F' // U times cumulative sum of weights

			gen double `W' = `V'/(A^2*E) // V scaled by pov line and sum of weights

			gen double `X' = `W' in 1
			replace `X' = `W' + `X'[_n-1] if `C' > 1

			gen byte `Y' = (`X' <= P & `X'[_n+1] > P)

			summ `C' if `Y' == 1
			assert r(min) == r(max)
			local which = r(mean)
			scalar AA = (`C'+1)  in `which' // person after J

			scalar AC = P - `X' in `which' // Remained MC

			scalar AD = A^2*E*AC
			tempvar AE
			gen `AE' = `X' > P & `X'[_n-1] < P // Adjusted Y column
			local p = P

		    summ `C' if `AE' == 1
			assert r(min) == r(max)
			local which2 = r(mean)
			scalar AG = `S' in `which2' // z-y_J
			scalar AH = AG^2 // (z-y_J)^2
			scalar AJ = `F' in `which2'

			scalar AK = (2*AG+sqrt(4*AH - 4*(AD/AJ)))/2 // First root of equation

			scalar AL = (2*AG - sqrt(4*AH - 4*(AD/AJ)))/2 // second root of equation

			scalar AM = min(max(AL,0),max(AK,0)) // taking smallest root greter than 0.

			scalar AO = `G' in `which2'
			tempvar AP AQ AR
			gen double `AP' = ((AO-`G') + AM)*(`G'<=AO) // Optimum benefit

			gen double `AQ' = `G' + `AP' // New EI

			gen double `AR' = ((A - `AQ')/A)^2*(`AQ'<A)

			sum `AR' [aw=`D']
			scalar AT= r(mean)

			*scalar AU = N - AT
			local AU = N - AT

		        *assert (`AU' == `P')
			assert abs(((`AU' - `P')/`P')) < 0.0001

			summ `AP' [aw=`D']
			scalar prime2 = r(sum)

			summ `H' [aw=`D']
			scalar tot_ben2 = r(sum)

			return scalar sp_ef_pov_2 = prime2/tot_ben2

			scalar drop E N O P Q AA AC AD AG AH AJ AK AL AM AO AT prime2 tot_ben2

		}

		scalar drop A B 
end





// BEGIN ceqeffect 
//  Calculates spending and impact effectiveness indicators of a specific fiscal intervention
//   measures for a specific poverty line and inequality 
* general programming locals

	local dit display as text in smcl
	local die display as error in smcl
	
capture program drop ceqeffect
program define ceqeffect, rclass
	
	#delimit; 
	syntax [if] [in] [pweight/] , 
		startinc(string)
		endinc(string)
		
		[
		tax(string)
		ben(string)
		z(string)
		noineq 
		nopov 
		]
	
	;
	
	#delimit cr
	
	
	preserve
	if wordcount("`if' `in'")!=0 quietly keep `if' `in'
	if "`exp'" !="" {
			local aw = "[aw = `exp']" //weights
			local pw = "[pw = `exp']"
		}
	
	
	if "`tax'"=="" local id_tax=0
	if "`tax'"!="" local id_tax=1
	if "`ben'"=="" local id_ben=0
	if "`ben'"!="" local id_ben=1
	
	local taxben=`id_tax'+`id_ben'
	
	//TAXES
	qui {
	if `id_tax'==1 { 
	
		// Inequality
		*TAX STAR
		ceqtaxstar `pw' , startinc(`startinc') tax(`tax') 
		local twarn = 0 
		if  r(t_gr) ==1 { 
			nois `dit' "Sum of `tname_`rw'_`cc'' exceed ``rw'', so impact effectiveness indicator not produced from ``rw'' to ``cc''" 
			local warning `warning' "Sum of `tname_`rw'_`cc'' exceed ``rw'', so impact effectiveness indicator not produced from ``rw'' to ``cc''" 						
			local twarn = r(t_gr) 
		} 
		else if  r(t_0) ==1 { 
			nois `dit' "Sum of `tname_`rw'_`cc'' is 0, so impact effectiveness indicator not produced from ``rw'' to ``cc''" 
			local warning `warning' "Sum of `tname_`rw'_`cc'' is 0, so impact effectiveness indicator not produced from ``rw'' to ``cc''" 					
			local twarn = r(t_0) 
		} 
		else  {
			tempvar ystar
			gen double `ystar'=____ytaxstar
			cap drop ____ytaxstar ____ybenstar ____id_benstar ____id_taxstar
			if "`noineq'"=="" & `id_tax'==1 {	
		
		// Impact Effectiveness
		
		covconc `endinc' `pw' //gini of column income
		local g1_end=r(gini)
		covconc `startinc' `pw' //gini of row income
		local g2_start=r(gini)
		covconc `ystar' `pw' //gini of star income
		local g_star=r(gini)
		local ie_gini=(`g2_start'-`g1_end')/(`g2_start'-`g_star')
		return scalar ie_gini = `ie_gini'

		// Spending Effectiveness
		_ceqmcid `pw', inc(`endinc') sptax(`tax') 
		*If Marg. Cont. is negative, SE is missing
		if r(mc_ineq)<0{
			return scalar se_gini = .
		}
		else{
			cap drop ___t 
			gen double ___t = `endinc' + abs(`tax') 
			covconc ___t `pw' 
			local gini1 = r(gini) 
			covconc `endinc' `pw'  
			local gini2 = r(gini) 
			if abs((`gini1' - `gini2')/`gini1') < 0.009 { 
				nois `dit' "Difference beween starting and ending Ginis is too small. Spending effectiveness indicator for ``y'' considering ``y'_`ext'' is not produced" 
				local warning `warning' `dit' "Difference beween starting and ending Ginis is too small. Spending effectiveness indicator for ``y'' considering ``y'_`ext'' is not produced" 
			} 
			else { 
				_ceqspend `pw',inc(`endinc') sptax(`tax') 
				return scalar se_gini = r(se_ind)	
				local se_gini=r(se_ind)
			} 
		}
		}
		
		// Poverty
		if ("`nopov'"=="") & `id_tax'==1 {
		// Impact Effectiveness
		
			tempvar normalized1  // create temporary variable end income
			tempvar normalized2  // create temporary variable start income...
			tempvar normalized3  // create temporary variable ystar income...
			qui gen `normalized1' = `endinc'/`z' // normalized by pov line 
			qui gen `normalized2' = `startinc'/`z' // normalized by pov line
			qui gen `normalized3' = `ystar'/`z' // normalized by pov line
			local _pline = 1                       // and normalized pov line is 1
			local vtouse1 `normalized1' // use normalized income in the calculations
			local vtouse2 `normalized2' // use normalized income in the calculations
			local vtouse3 `normalized3' // use normalized income in the calculations
			
			tempvar fgt1_1 fgt2_1 fgt1_2 fgt2_2 fgt1_3 fgt2_3    
			qui gen `fgt1_1' = max((`_pline'-`vtouse1')/`_pline',0) // normalized povety gap of each individual
			qui gen `fgt2_1' = `fgt1_1'^2   // square of normalized poverty gap
			qui gen `fgt1_2' = max((`_pline'-`vtouse2')/`_pline',0) // normalized povety gap of each individual
			qui gen `fgt2_2' = `fgt1_2'^2                            // square of normalized poverty gap
			qui gen `fgt1_3' = max((`_pline'-`vtouse3')/`_pline',0) // normalized povety gap of each individual
			qui gen `fgt2_3' = `fgt1_3'^2                            // square of normalized poverty gap
			forval i=1/2 {
				qui summ `fgt`i'_1' `aw', meanonly // `if' `in' restrictions already taken care of by `touse' above	
				local p`i'_orig=r(mean)
				qui summ `zyzfgt`i'_2' `aw', meanonly // `if' `in' restrictions already taken care of by `touse' above
				local p`i'_wo=r(mean)
				qui summ `zyzfgt`i'_3' `aw', meanonly // `if' `in' restrictions already taken care of by `touse' above
				local p`i'_3=r(mean)
				drop  `fgt`i'_3' `fgt`i'_2' `fgt`i'_1'
				****Poverty Impact effectiveness
				//Marginal contributions for fgt 1,2
				local mp_`i'=`p`i'_wo'-`p`i'_orig' //Observed MC
				****For Impact effectiveness there can only be a negative effect, we use the harm formula Ch. 5 CEQ Handbook
				// Harm formula (taxes)
				if `mp_`i''<0{
					ceqtaxharm `pw', endinc(`endinc') tax(`tax')
					if r(thwarn) == 1 { 
						nois `dit' "Sum of `tname_`rw'_`cc'' exceed ``rw'', so impact effectiveness indicator not produced from ``rw'' to ``cc''" 
						local warning `warning' "Sum of `tname_`rw'_`cc'' exceed ``rw'', so impact effectiveness indicator not produced from ``rw'' to ``cc''" 
					}
							
					tempvar yharm
						gen double `yharm'=____ytaxharm
						cap drop ____ytaxharm   ____id_taxharm
						tempvar vtousehnorm
						gen `vtousehnorm'=(`yharm'/`z')
						local _pline = 1                       // and normalized pov line is 1
						local vtouseh `vtousehnorm'			
						tempvar fgt1_h fgt2_h
						qui gen `fgt1_h' = max((`_pline'-`vtouseh')/`_pline',0) // normalized povety gap of each individual
						qui gen `fgt2_h' = `fgt1_h'^2                            // square of normalized poverty gap							
						qui summ `fgt`i'_h' `aw', meanonly // `if' `in' restrictions already taken care of by `touse' above	
						local p`i'_h=r(mean)
						local mst_p`i'_h=`p`i'_3' - `p`i'_h' //Ideal MC with tax formula
						local ie_p`i'=(`mp_`i''/`mst_p`i'_h')*(-1)
						
						return scalar ie_p`i'  = `ie_p`i''
						local se_p`i'=.
						return scalar se_p`i'  = .	
				}
				else{
					local ie_p`i'=.
					return scalar ie_p`i'  = .	
					local se_p`i'=.
					return scalar se_p`i'  = .	
				}
			}
		}
				
			// Spending Effectiveness
		*No spending effectiveness for taxes on poverty
		//from Enami (2008): The spending effectiveness indicator can only be calculated for the taxes and transfers
		// with a positive MC (and as a result, the spending effectiveness of taxes on poverty
		//reduction is undefined).
	}
}
									
		
		
		

		
		
	
	//Transfers
	
	if `id_ben'==1 { 
		*Ben STAR
		ceqbenstar `pw', startinc(`startinc') ben(`ben')
		local bwarn = 0 
		if r(b_gr) == 1 { 
			nois `dit' "Sum of `bname_`rw'_`cc'' is 0, so impact effectiveness indicator not produced from ``rw'' to ``cc''" 
			local warning `warning' "Sum of `bname_`rw'_`cc'' exceed ``rw'', so impact effectiveness indicator not produced from ``rw'' to ``cc''" 
			local bwarn = r(b_gr) 
		} 
		else if r(b_0) == 1 { 
			nois `dit' "Sum of `bname_`rw'_`cc'' is 0, so impact effectiveness indicator not produced from ``rw'' to ``cc''" 
			local warning `warning' "Sum of `bname_`rw'_`cc'' exceed ``rw'', so impact effectiveness indicator not produced from ``rw'' to ``cc''" 
			local bwarn = r(b_0) 
		} 
		else  {		
			tempvar ystar
			gen double `ystar'=____ybenstar
			cap drop ____ytaxstar ____ybenstar ____id_benstar ____id_taxstar
		}
								
		
		// Inequality
		if ("`noineq'"=="") & `id_ben'==1{	
		// Impact Effectiveness
			if !( "`bwarn'" == "1" & "`twarn'" == "1" ) { 
				covconc `endinc' `pw' //gini of end income
				local g1_`endinc'=r(gini)
				covconc `startinc' `pw' //gini of orig income
				local g2_`startinc'=r(gini)
				covconc `ystar' `pw' //gini of star income
				local g_star=r(gini)
				local ie_gini=(`g2_`startinc''-`g1_`endinc'')/(`g2_`startinc''-`g_star')
				
				 return scalar ie_gini=`ie_gini'
			}
		// Spending Effectiveness
		
		_ceqmcid `pw', inc(`endinc') spben(`ben') 
			*If Marg. Cont. is negative, SE is missing
			if r(mc_ineq)<0{
				 return scalar se_gini=.
				local se_gini=.
			}
			else{
				cap drop ___t  
				gen double ___t = `endinc' - abs(`ben')  
				covconc ___t `pw'  
				local gini1 = r(gini)  
				covconc `endinc' `pw'  
				local gini2 = r(gini) 
				if abs((`gini1' - `gini2')/`gini1') < 0.009 { 
					nois `dit' "Difference beween starting and ending Ginis is too small. Spending effectiveness indicator for `ben' considering `endinc' is not produced" 
					local warning `warning' `dit' "Difference beween starting and ending Ginis is too small. Spending effectiveness indicator for `ben' considering `endinc' is not produced" 
								
				} 

				else { 
					_ceqspend `pw',inc(`endinc') spben(`ben') // ! Remove capture once debugged.
						*local spef=r(sp_ef)
					local se_gini=r(sp_ef)
					return scalar se_gini=r(sp_ef)

				}  
			}
		
		}
		
		// Poverty
		if ("`nopov'"=="") & `id_ben'==1{
			// Impact Effectiveness
			
			tempvar normalized4 normalized5 normalized6
			qui gen `normalized4' = `endinc'/`z' // normalized by pov line  
			qui gen `normalized5' = `startinc'/`z' // normalized by pov line
			if `bwarn' == 0 qui gen `normalized6' = `ystar'/`z' // normalized by pov line
			local _pline = 1                       // and normalized pov line is 1
			local vtouse1 `normalized4' // use normalized income in the calculations 
			local vtouse2 `normalized5' // use normalized income in the calculations 
			local vtouse3 `normalized6' // use normalized income in the calculations 
			
			qui gen double ___zyzfgt1_1 = max((`_pline'-`vtouse1')/`_pline',0) // normalized povety gap of each individual
			qui gen double ___zyzfgt2_1 = ___zyzfgt1_1^2                            // square of normalized poverty gap
			qui gen double ___zyzfgt1_2 = max((`_pline'-`vtouse2')/`_pline',0) // normalized povety gap of each individual
			qui gen double ___zyzfgt2_2 = ___zyzfgt1_2^2                            // square of normalized poverty gap
			if `bwarn' == 0 qui gen double ___zyzfgt1_3 = max((`_pline'-`vtouse3')/`_pline',0) // normalized povety gap of each individual
			if `bwarn' == 0 qui gen double ___zyzfgt2_3 = ___zyzfgt1_3^2                            // square of normalized poverty gap
						
			qui summ ___zyzfgt1_1 `aw', meanonly // `if' `in' restrictions already taken care of by `touse' above	
			local p1_orig=r(mean)
			qui summ ___zyzfgt1_2 `aw', meanonly // `if' `in' restrictions already taken care of by `touse' above	
			local p1_2=r(mean)
			if `bwarn' == 0 qui summ ___zyzfgt1_3 `aw', meanonly // `if' `in' restrictions already taken care of by `touse' above	
			local p1_3_st=r(mean)
			qui summ ___zyzfgt2_1 `aw', meanonly // `if' `in' restrictions already taken care of by `touse' above	
			local p2_orig=r(mean)
			qui summ ___zyzfgt2_2 `aw', meanonly // `if' `in' restrictions already taken care of by `touse' above	
			local p2_2=r(mean)
			if `bwarn' == 0 qui summ ___zyzfgt2_3 `aw', meanonly // `if' `in' restrictions already taken care of by `touse' above	
			local p2_3_st=r(mean)
						
			//Marginal contributions for fgt 1,2
			local mp_1_`p'=`p1_2'-`p1_orig' //Observed MC
			if `bwarn' == 0 local mp_1_`p'_s=`p1_2'-`p1_3_st' //Star MC
			local mp_2_`p'=`p2_2'-`p2_orig' //Observed MC
			if `bwarn' == 0 local mp_2_`p'_s=`p2_2'-`p2_3_st' //Star MC						
			// If warning was produced we should skip this section.
			if `bwarn' != 0 local mp_1_`p'_s= 1 
			if `bwarn' != 0 local mp_2_`p'_s= 1 
			forval i=1/2 {
				****Poverty Impact effectiveness
				****For Impact effectiveness with Transfers there can only be a positive effect
				if `mp_`i'_`p''>0{
					*Impact effectiveness
					return scalar ie_p`i' =  (`mp_`i'_`p''/`mp_`i'_`p'_s') //MC/MCstar
					local ie_p`i' = 	(`mp_`i'_`p''/`mp_`i'_`p'_s') //MC/MCstar
					// Spending Effectiveness
					cap drop ___bentouse
					gen double ___bentouse=abs(`vtouse1'-`vtouse2')
					covconc `vtouse1' `pw'  
					local gini1 = r(gini)  
					covconc `vtouse2' `pw'  
					local gini2 = r(gini) 
					if abs((`gini1' - `gini2')/`gini1') < 0.009 { 
						nois `dit' "Difference beween starting and ending Ginis is too small. Poverty spending effectiveness indicator for `ben' considering `startinc' is not produced" 
					    local warning `warning' `dit' "Difference beween starting and ending Ginis is too small. Poverty spending effectiveness indicator for `ben' considering `startinc' is not produced" 
					} 
					else { 
						ceqbensp  `pw', startinc(`vtouse2') ben(___bentouse) zz(`_pline') obj1(`p1_orig') obj2(`p2_orig')	
						local se_p`i'= r(sp_ef_pov_`i')    							
						return scalar se_p`i'=r(sp_ef_pov_`i')    							
					} 

				}
				else{
					local ie_p`i'= .   							
					return scalar ie_p`i'=.    		
					local se_p`i'= .   							
					return scalar se_p`i'=.   		
				}
							
			}
		

		}
	}
	
	}

	tempname table
			.`table'  = ._tab.new, col(3)  separator(1) lmargin(0)
			.`table'.width  20 20 20  
			.`table'.strcolor yellow yellow yellow  
			.`table'.numcolor yellow yellow yellow    
			.`table'.numfmt %16s  %16.5f %16.5f

	       
	      	.`table'.sep, top
	      	.`table'.titles "Indicator" "Impact Eff." "Spending Eff." 
			
			scalar r1 = "Inequality - Gini"
			scalar r2 = "Poverty Gap"
			scalar r3 = "Squared Poverty Gap"

			.`table'.sep, mid
			.`table'.row r1 `ie_gini' `se_gini'
			.`table'.row r2 `ie_p1' `se_p1'
			.`table'.row r3 `ie_p2' `se_p2'
			
		   .`table'.sep,bot
restore			
end // END ceqeffect
