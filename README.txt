	 "Spatial Source Phase: A New Feature for Identifying Spatial Differences Based on Complex-Valued Resting-State fMRI Data," submitted to Human Brain Mapping

	Requirements:   
        Matlab R2014a   GroupICATv4.0b

	Datatsets:
        Our results data can be downloaded from http://pan.dlut.edu.cn/share?id=cag3ass9t7ec
        The results data are obtained by using entropy bound minimization (EBM) algorithm of ICA (number of components N = 120, 10 runs, 82 subjects), and provided in data\res_SM_TC\res_SM_TC_k_*.mat. 
        We also provide the index for in-brain voxels in data\mask_ind.mat, and the spatial reference for DMN [1] and Auditory [2] in data\ref_DMN.mat and data\ref_AUD.mat.


	Experimental steps:     
        (1) Run data_prepare.m to perform component selection,phase processing including de-ambiguity and phase de-noising, and best run selection. 
        Note: Our results of spatial maps estimates for DMN and auditory cortex are stored in data\sm_EBM_DMN.mat and data\sm_EBM_AUD.mat.

        (2) Run data_validation.m to generate bootstrap samples and resampling samples for intra- and inter- difference.
        Note: Our results of resampling samples for DMN and auditory cortex are stored in data\bs_DMN.mat ,data\bs_AUD.mat, data\intra_inter_DMN.mat and data\intra_inter_AUD.mat.

        (3) Run main.m to perform Variance difference analysis and validation. 
        Note: Results for DMN and auditory cortex in the manuscript "Spatial Source Phase: A New Feature for Identifying Spatial Differences Based on Complex-Valued Resting-State fMRI" (submitted to Human Brain Mapping) can be obtained from demo_for_DMN.m and demo_for_Auditory.m.


    Functions: 
        (1) barweb.m: Plot histograms with standard deviation.
        (2) Bca_bootrap.m: Computes the bias-corrected and accelerated bootstrap estimate [3].
        (3) bestRunSelection.m: The best run selection by combining cross-run averaging and a one-sample t-test proposed in [4].
        (4) componentSelection.m: Component selection based on reference.
        (5) diffInNetwork.m: Calculate q and m in the network.
        (6) groupMask.m: Mask the spatial maps for shared voxels which are activated in more than half of all the subjects.
        (7) loo_pre.m: Contruct vectors of length N of the leave-one-out values of the statistic of interest to prepare for the bias-corrected and accelerated bootstrap estimate [3].
        (8) phase_de_ambiguity.m: Phase de-ambiguity based on maximizing the energy of the TC real part, which was proposed in [5].
        (9) phase_denoising.m: Phase de-denoising based on thresholds for spatial source phase and magnitude, which was proposed in [5].
        (10) SMdisplay.m: Show the saptial maps. The function is adapted from the GIFT toolbox (http://mialab.mrn.org/software/gift/index.html) and some functions and data involved are provided in the folder named SMshow.
        (11) tRef.m: Subfunction of bestRunSelection.m, which generates reference by combining cross-run averaging and a one-sample t-test.
        (12) varAnalysis.m: Calculate variance difference using F-test and FDR correction in terms of difference variance maps of SZ-HC and m.

    References:
 
        [1] Smith, S. M., Fox, P. T., Miller, K. L., Glahn, D. C., Fox, P. M., & Mackay, C. E., et al. (2009). Correspondence of the brain's functional architecture during activation and rest. Proceedings of the National Academy of Sciences of the United States of America, 106(31), 13040-13045.   
        [2] Allen, E. A., Erhardt, E. B., Damaraju, E., Gruner, W., Segall, J. M., & Silva, R. F., et al. (2011). A baseline for the multivariate comparison of resting-state networks. Frontiers in Systems Neuroscience, 5(2), 1¨C23. 
        [3] Efron, B., Tibshirani R. J. (1993). An Introduction to the Bootstrap. New York, Chapman & Hall/CRC.
        [4] Kuang,L.D., Lin,Q.H., Gong,X.F., Cong,F., Sui,J., & Calhoun,V.D.(2018). Model order effects on ICA of resting-state complex-valued fMRI data: application to schizophrenia. Journal of Neuroscience Methods, 304, 24-38. 
        [5] Yu, M. C., Lin, Q. H., Kuang, L. D., Gong, X. F., Cong, F., & Calhoun, V. D. (2015). ICA of full complex-valued fMRI data using phase information of spatial maps. Journal of Neuroscience Methods, 249, 75¨C91. 


