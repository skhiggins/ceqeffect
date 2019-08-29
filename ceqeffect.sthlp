{smcl}
{* 30jun2019}{...}
{cmd:help ceqeffect} (beta version; please report bugs) {right: Rodrigo Aranda}
{hline}

{title:Title}

{p 4 11 2}
{hi: ceqeffect} {hline 2} Calculates effectiveness indicators of fiscal interventions on poverty and inequality based on Enami (2016).


{pstd}

{title:Syntax}

{p 8 11 2}
    {cmd:ceqeffect} {ifin} {weight}  [{cmd:,} {startinc(varname)} {endinc(varname)} {tax(varname)} {ben(varname)} {z(string)} {it:options}]{break}

{synoptset 29 tabbed}{...}
{synopthdr}
{synoptline}


{syntab:Options}
{synopt :{opt z(string)}}Poverty Line (can be a number or a variable, has to be in the same unit as income and intervention){p_end}
{synopt :{opt noineq}}Do not run results with inequality{p_end}		
{synopt :{opt nopov}}Do not run results with poverty{p_end}		

   

{synoptline}		
{p 4 6 2}
{cmd:pweight} allowed; see {help weights}. 


{title:Description}

{pstd} 
The command {cmd:ceqeffect} estimates impact effectiveness and spending effectiveness indicators for Gini, poverty gap and square poverty gap.
To be able to identify the effectiveness of a fiscal intervention, {cmd:ceqeffect} needs {opth startinc(varname)} (starting income),
 {opth endinc(varname)} and the intervention ({opth tax(varname)} or {opth ben(varname)}). Only one intervention can be measured at a time. 
 The poverty line can be a number or a variable but has to be in the same unit as the incomes and the intervention.
 

{title:Examples}

{pstd}Effectiveness indicators for a transfer{p_end}
{phang} {cmd:. ceqeffect [pw=weight], startinc(y0) endinc(y1) ben(transf1) z(5000)}{p_end}

{pstd}Effectiveness indicators for a tax{p_end}
{phang} {cmd:. ceqeffect [pw=weight], startinc(y_beforetax) endinc(y_aftertax) tax(tax1) z(5000)}{p_end}

{title:Saved results}

 {phang} {cmd:.r(ie_gini) - Impact Effectiveness Indicator for Gini}{p_end}
 {phang} {cmd:.r(se_gini) - Spending Effectiveness Indicator for Gini}{p_end}
 {phang} {cmd:.r(ie_p1) - Impact Effectiveness Indicator for Poverty Gap}{p_end}
 {phang} {cmd:.r(ie_p2) - Impact Effectiveness Indicator for Squared Poverty Gap}{p_end}
 {phang} {cmd:.r(se_p1) - Spending Effectiveness Indicator for Poverty Gap}{p_end}
 {phang} {cmd:.r(se_p2) - Spending Effectiveness Indicator for Squared Poverty Gap}{p_end}


{title:Author}

{p 4 4 2}Rodrigo Aranda,  rarandabal@gmail.edu


{title:References}

{phang}
Cite
Enami, Ali. "Measuring the effectiveness of taxes and transfers in fighting inequality and poverty." Commitment to equity handbook: Estimating the impact of fiscal policy on inequality and poverty (2018): 207-216.


