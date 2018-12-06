%*****************************************************************
% Description: contruct vectors of length N of the leave-one-out values of the statistic
% of interest to prepare for the bias-corrected and accelerated bootstrap estimate of Efron,
% B., & Tibshirani, R. J. (1993). An Introduction to the Bootstrap, Chapman & Hall/CRC: New York.
% Usage:
%   [loo,loo_FDR] = loo_pre(data_N, data_A, half_ind, area, alpha) 
% Input: 
    % data_N: data for healthy controls with dimension k1 x V, k1 is the number of healty controls
	% data_A: data for patients with dimension k2 x V, k2 is the number of patients
    % half_ind: the index for voxels survived magnitude and spatial sourse phase thresholds in more than half of subjects
    % area: the sub-network, e.g.,'PCC', 'ACC'
    % alpha: a value specifying the significance level as (100*ALPHA)%
% Output: 
    % loo: a vector of length N of the leave-one-out values of the statistic of interest.
    % loo_FDR: a vector of length N of the leave-one-out values of the statistic of interest with FDR correction
% Date: December 2018
% Author: Yue Qiu
%*****************************************************************
function [loo,loo_FDR] = loo_pre(data_N, data_A, half_ind, area, alpha) 
    for i = 1:size(data_N, 1) + size(data_A, 1)
        a = data_N; b = data_A;
        if(i <= size(data_N, 1)) 
            a(i,:) = []; 
        else
            b(i - size(data_N, 1), :) = [];
        end   
        [diff(i, :), m(i), diff_FDR(i, :), m_FDR(i)] = varAnalysis(a, b, alpha);
        [~, loo(i)] = diffInNetwork(diff(i, :), half_ind, area, 0);
        [~, loo_FDR(i)] = diffInNetwork(diff_FDR(i, :), half_ind, area, 0);
    end

