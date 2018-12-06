%*****************************************************************
% Description: calculate variance difference using F-test and FDR correction
%             in terms of difference variance maps of SZ-HC and m (number of voxels survived the voxel-wise F-test)
% Usage:
%   [diff, m, diff_FDR, m_FDR] = varAnalysis(data_N, data_A, alpha)
% Input: 
    % data_N: data for healthy controls with dimension k1 x V, k1 is the number of healty controls
	% data_A: data for patients with dimension k2 x V, k2 is the number of patients
    % alpha: a value ALPHA between 0 and 1 specifying the significance level as (100*ALPHA)%.
% Output: 
    % diff: the difference variance maps of SZ-HC (F-test)
    % m: number of voxels survived the voxel-wise F-test
    % diff_FDR: the difference variance maps of SZ-HC (F-test with FDR correction)
    % m_FDR: number of voxels survived the voxel-wise F-test with FDR correction
% Reference: "Spatial Source Phase: A New Feature for Identifying Spatial Differences 
%            Based on Complex-Valued Resting-State fMRI Data," submitted to Human Brain Mapping
% Date: October 2017
% Author: Yue Qiu
%*****************************************************************
function [diff, m, diff_FDR, m_FDR] = varAnalysis(data_N, data_A, alpha)
varpS = var(data_N);csz_varpS = var(data_A);diff=csz_varpS-varpS;
% F-test
[h,p,~,~]=vartest2(data_N, data_A, alpha); 
m = length(find(p < alpha));
diff = diff.*h;
% FDR correction
p_FDR = mafdr(p); m_FDR = length(find(p_FDR < alpha));
h(find(p_FDR > alpha)) = 0;
diff_FDR = diff .* h;