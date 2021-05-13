# head_motion
pediatric head motion characterization and analysis

## dataset
- Data from the Healthy Brain Network biobank (HBN) were used for all analyses (N=1388).
- http://fcon_1000.projects.nitrc.org/indi/cmi_healthy_brain_network/
	> Alexander LM, Escalera J, Ai L, Andreotti C, Febre K, Mangone A, et al.  
	An open resource for transdiagnostic research in pediatric mental health and learning disorders. Sci Data. 2017;4: 1–26.  
	doi:10.1038/sdata.2017.181

- Note: In this release, MCFLIRT .par files and HBN demographic and behavioural data are not included as per our data usage agreement.
- Subject lists: 
- - HBN-1388: `~\data\subjectList_hbn1388.txt`
- - HBN-865:  `~\data\subjectList_hbn865.txt` 

## preprocessing
Run hm_import.m then hm_preprocess.m to generate data structures

### 1: hm_import.m
- Import all MCFLIRT .par files into MATLAB structure "hm_data"
- Import behavioural and demographic data from COINS into "hm_data"

### 2: hm_prepreocess.m
- calculate FD, mean FD, and composition of FD
- generate motion cohort indices
- FD calculated as in: 
	> Power JD, Barnes KA, Snyder AZ, Schlaggar BL, Petersen SE.  
	Spurious but systematic correlations in functional connectivity MRI networks arise from subject motion. NeuroImage. 2012;59: 2142–2154.  
	doi:10.1016/j.neuroimage.2011.10.018

### hm_fft.m
- calculates the spectral power of raw displacement data
- used in `fig_frequency.m`
- closely follows:
	> Fair DA, Miranda-Dominguez O, Snyder AZ, Perrone A, Earl EA, Van AN, et al.   
	Correction of respiratory artifacts in MRI head motion estimates. NeuroImage. 2020;208: 116400.  
	doi:10.1016/j.neuroimage.2019.116400

## analyses
The remaining scripts visualize or export data for further analysis in PRISM
 		
***
Naturalistic Neuroimaging Lab  
BC Children's Research Institute  
University of British Columbia   
https://www.headspacestudios.org/

