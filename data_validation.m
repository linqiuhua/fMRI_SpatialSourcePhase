%*****************************************************************
% Description: generate data for resampling validation
%              including bootstrap samples and resampling for intra- and inter- differences
% Input: sm_EBM_*.mat: saptial maps by EBM for network *  with component selection, 
%        phase processing and best-run seletction,including half_ind and Sm.
% Output: (1) bs_*.mat: bootstrap samples for network *, including bs_diff, bs_diff_FDR, bs_mag_diff and bs_mag_diff_FDR.
%         (2) intra_inter_*.mat: resampling samples for intra- and inter- difference in network *, 
%            including HC_diff, HC_diff_FDR, HS_diff, HS_diff_FDR, SZ_diff, SZ_diff_FDR,
%            HC_mag_diff, HC_mag_diff_FDR, HS_mag_diff, HS_mag_diff_FDR, SZ_mag_diff and SZ_mag_diff_FDR.
% Reference: "Spatial Source Phase: A New Feature for Identifying Spatial Differences 
%            Based on Complex-Valued Resting-State fMRI Data," submitted to Human Brain Mapping
% Note£ºThe following is an example for analyzing the DMN
% Date: August 2018
% Author: Yue Qiu
%*****************************************************************

clc;clear all;

%% -------------------Initialize parameters-------------------- %%
L = 1000;  % the number of replications
alpha = 0.05  % a value specifying the significance level as (100*ALPHA)%
k1 = 40;  % number of healthy controls
k2 = 42;  % number of patients
kk = 20;  % number of subjects of G1 or G2

%% -------------------------Input data------------------------- %%    
load data\sm_EBM_DMN.mat half_ind Sm
% "load data\sm_EBM_AUD.mat" for auditory

phase = angle(Sm); magnitude = abs(Sm);

% I: bootstrap resampling:
[~,bootsam] = bootstrp(L,@mean,[1:k1 + k2]);
for l = 1 : L
    disp(['l=',num2str(l)]);
    x = bootsam(bootsam(:,l) <= k1, l);
    y = bootsam(bootsam(:,l) > k1, l);
    ph_N = phase(x,:);ph_A = phase(y,:);
    [bs_diff(l, :), ~, bs_diff_FDR(l, :),~] = varAnalysis(ph_N, ph_A, alpha);
    m_N=magnitude(x,:);m_A=magnitude(y,:);
    [bs_mag_diff(l, :), ~, bs_mag_diff_FDR(l, :), ~] = varAnalysis(m_N, m_A,alpha);
end
save data\bs_DMN.mat bs_diff bs_diff_FDR bs_mag_diff bs_mag_diff_FDR
%save data\bs_AUD.mat bs_diff bs_diff_FDR bs_mag_diff bs_mag_diff_FDR



% II: intra- and inter- resampling:
phase_N = phase(1:k1,:);phase_A = phase(k1+1:k1+k2,:); 
mag_N = magnitude(1:k1,:); mag_A = magnitude(k1+1:k1+k2,:);

for l = 1 : L
    disp(['l=',num2str(l)]);
    x = randperm(k1, kk); flag = 1;
    y = zeros(1, kk);
    for i = 1 : k1
        if(ismember(x, i) == 0)
            y(flag) = i;
            flag = flag+1;
        end
    end
    ph_N1 = phase_N(x,:); ph_N2 = phase_N(y,:); ph_A1 = phase_A(x,:); ph_A2 = phase_A(y,:);
    [HC_diff(l, :), ~, HC_diff_FDR(l, :), ~] = varAnalysis(ph_N1, ph_N2, alpha);
    [HS_diff(l, :), ~, HS_diff_FDR(l, :), ~] = varAnalysis(ph_N1, ph_A1, alpha);
    [SZ_diff(l, :), ~, SZ_diff_FDR(l, :), ~] = varAnalysis(ph_A1, ph_A2, alpha);
    m_N1 = mag_N(x,:);m_N2 = mag_N(y,:); m_A1 = mag_A(x,:);m_A2 = mag_A(y,:);
    [HC_mag_diff(l, :), ~, HC_mag_diff_FDR(l, :), ~] = varAnalysis(m_N1, m_N2, alpha);
    [HS_mag_diff(l, :), ~, HS_mag_diff_FDR(l, :), ~] = varAnalysis(m_N1, m_A1, alpha);
    [SZ_mag_diff(l, :), ~, SZ_mag_diff_FDR(l, :), ~] = varAnalysis(m_A1, m_A2, alpha);
end
save data\intra_inter_DMN.mat HC_diff HC_diff_FDR HS_diff HS_diff_FDR SZ_diff SZ_diff_FDR
save data\intra_inter_DMN.mat HC_mag_diff HC_mag_diff_FDR HS_mag_diff HS_mag_diff_FDR SZ_mag_diff SZ_mag_diff_FDR -append
% "save data\intra_inter_AUD.mat HC_diff HC_diff_FDR HS_diff HS_diff_FDR SZ_diff SZ_diff_FDR" for Auditory
% "save data\intra_inter_AUD.mat HC_mag_diff HC_mag_diff_FDR HS_mag_diff HS_mag_diff_FDR SZ_mag_diff SZ_mag_diff_FDR -append" for Auditory
