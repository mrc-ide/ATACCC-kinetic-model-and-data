Code and data files for VL kinetic models in Singanayagam et al, Lancet ID, 2021
================================================================================

Code
====

All code (c) Neil Ferguson, Imperial College London, licensed under Apache v2 (http://www.apache.org/licenses/LICENSE-2.0), as specified in each source file. neil.ferguson@imperial.ac.uk

1. To use, unpack all source files and the two data files into a single source directory. Load R 4.1 or later and setwd to that directory.

2. Linux recommended as Stan can be less stable under Windows

3. 64GB RAM needed, 12+ core system recommended

4. Check R packages listed in libs.r are installed

5. Run "run_fits.r" to fit the 6 model variants listed in the paper for each of the ORF1 and E gene Ct value datasets

6. Each fit will take 4-12h, including running loo_cv, so better to run each fit on a cluster

7. Run "process.r" to process the results and generate Fig 1 from the paper and the results listed in the tables

8. Sorry, no automated generation of tables included - need to look at/save the interactive R session console output

9. Code has bare bones comments - enough to follow, with effort. Email me if help needed (no promises!)


Data file structure
===================

- Long format used. Summarise by person_id to look at subject specific variables
- 20210920_ATACCC_L.csv - ORF1ab Ct value measurements, 20210920_ATACCC_E_L.csv - E gene CT value measurements
- Fields:
	PersonID - numerical identifier for subject
	Subject - original ATACCC alphanumeric identifier of subject
	Age - subject's age in years
	AgeGT35 - 1 if age>=35, 0 otherwise
	SexFemale - dummy (always 1 and not used)
	ContactType - Household or Non-Household
	Vaccinated - 0 =  unvaccinated at time of enrolment, 1 =  vaccinated
	VaccInclude - 1 =  included in analysis (unvaccinated or fully vaccinated at enrolment), 0 =  not included (partially vaccinated)
	VacVOC = Group ID: 0 = Pre-Alpha (unvacc), 1 = Alpha (unvacc), 2 = Delta (unvacc), 3 = Delta (vacc)
	DateEnrolment - Date of first swab
	DaysTrueExp2FirstSwab - Delay from reported date of exposure (of subject) to enrolment (non-household contacts only). This is used to generate pseudo-absence data points (a "measurement" of undetectable VL, Ct=40, on the day of exposure)
	FirstSwabNotMax - 1 if first measurement is also the max VL recorded, 0 otherwise
	FirstSwabHalfMaxVL - 1 if first VL measurement is < 1/2 of max VL, 0 otherwise
	FullFollowup - 1 if at least 12 days of follow-up, 0 otherwise
	fullProfile - 1 if first VL measurement is undetectable (Ct=40), meaning the peak VL has been resolved, 0 otherwise 
	B117Status - 0 for Pre-Alpha, 1 for Alpha, 2 for Delta
	Symptomatic - 1 if there were symptoms meeting WHO definition on day of enrolment, 0 otherwise
	MinDay - 10 + Day after enrolment min Ct (max VL) was recorded. Hence 11 implies day of enrolment and first swab 
	MinCT - Minimum Ct value recorded
	TestDateIndex - day of swab, on a scale where 0 is always the day of max VL (min CT)
	CtT1 -  recorded Ct value, where 40 represents undetectable
	

