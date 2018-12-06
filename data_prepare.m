%*****************************************************************
% Description: component selection and phase processing
% Input: (1) mask_ind.mat: index for in-brain voxels, including mask_ind
%        (2) ref_*.mat: saptial reference for network *, including Sr
%        (3) res_SM_TC_k_*.mat: multiple spatial maps and associated time courses
%        of subject * obtained from EBM, including At and S
% Output: sm_EBM_*.mat: index (half_ind) of shared voxels in all the subjects and complex-valued spatial maps (Sm) constructed by these voxels
% Reference:"Spatial Source Phase: A New Feature for Identifying Spatial Differences 
%            Based on Complex-Valued Resting-State fMRI Data," submitted to Human Brain Mapping
% Note£ºThe following is an example for analyzing the DMN
% Date: October 2017
% Author: Yue Qiu
%*****************************************************************

clc;clear all;

%% ---------------------Initialize parameters------------------ %%
K = 82; % number of subjects
C = 10; % number of candidate components for component of interest
Zr = 0.5;  % magnitude threshold
alpha = 0.05;  % a value specifying the significance level as (100*ALPHA)%

%% ---------------------------Input data----------------------- %%    
load data\mask_ind.mat
load data\ref_DMN.mat  % "load data\ref_AUD.mat" for auditory

for k = 1 : K
    disp(['k=',num2str(k)]);
    eval(['load D:\Mydata\rest_SM_TC\res_SM_TC_k_' num2str(k) '.mat At S'])  
    R = size(S,1);
    for r = 1 : R
        % Select DMN from N components
        inn(r, k) = componentSelection(S(r,:,:), Sr, C); 
        % Perform phase processing
        [cS(k,:,r),cA(:,k,r)] = phase_de_ambiguity(At(r,:,inn(r,k)).',S(r,inn(r,k),:),Sr);
        pS(k,:,r) = phase_denoising(cS(k,:,r));
    end
    % Best run selection among R runs
    [cSm(k,:), pSm(k,:), rr(k), in(k)] = bestRunSelection(cS(k,:,:),pS(k,:,:), mask_ind, alpha);
end

% Obtain shared voxels(Sm) in all the subjects and their associated index(half_ind)
[Sm, half_ind] = groupMask(cSm, mask_ind, Zr);

% Save index and spatial source phase and magnitude vectors
save data\sm_EBM_DMN.mat Sm half_ind   % "save data\sm_EBM_AUD.mat Sm half_ind" for auditory

