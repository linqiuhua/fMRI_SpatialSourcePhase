%**************************************************************************
% Description: demo for getting results of Auditory in the submission
% Input: (1) mask_ind.mat: index for in-brain voxels, including mask_ind.
%        (2) sm_EBM_AUD.mat: saptial maps estimate for Auditory by EBM with component 
%            selection, phase processing and best-run seletction,including half_ind and Sm.
%        (3) bs_AUD.mat: 1000 bootstrap samples for Auditory, including bs_diff, bs_diff_FDR, bs_mag_diff and bs_mag_diff_FDR.
%        (4) intra_inter_AUD.mat: 1000 resampling results for intra- and inter- difference in Auditory, 
%            including HC_diff, HC_diff_FDR, HS_diff, HS_diff_FDR, SZ_diff, SZ_diff_FDR,
%            HC_mag_diff, HC_mag_diff_FDR, HS_mag_diff, HS_mag_diff_FDR, SZ_mag_diff and SZ_mag_diff_FDR.
% Output: variance difference of SZ-HC on spatial source phase and magnitude in Auditory. 
% Reference: "Spatial Source Phase: A New Feature for Identifying Spatial Differences 
%            Based on Complex-Valued Resting-State fMRI Data," submitted to Human Brain Mapping
% Note: details of component selection and phase processing, including phase de-amiguity and phase de-noising, are introduced in data_prepare.m
% Note: details of resample techniques are shown in data_validation.m
% Date: December 2018
% Author: Yue Qiu
%**************************************************************************

clc;clear all;

%% ----------------------Initialize parameters----------------- %%
sub_networks = {'AUDL','AUDR'};  % sub-networks of DMN
location = [145,110,80];  % slice coordinates for SMdisplay.m
k1 = 40;   % number of healthy controls
k2 = 42;  % number of schizophrenia patients
alpha = 0.05;  % a value specifying the significance level as (100*ALPHA)%
null = 0.5;  % the value of the statistic of interest under the null hypothesis
color = 5;  % colormap for SMdisplay.m
width = 0.5;  % width of the bars for Histogram
error_sides = 2; % parameter for plots +/- std
color_map = autumn;  % colormap for histogram
color_line = {'r', 'g', 'b', 'y'}; % line color
close all;

%% ---------------------------Input data----------------------- %%      
load data\mask_ind.mat
load data\sm_EBM_AUD.mat
load data\bs_AUD.mat
load data\intra_inter_AUD.mat
load data\intra_inter_AUD.mat

%% ----Construct spatial source phase and magnitude vectors---- %%
phase = angle(Sm); phase_N = phase(1:k1, :); phase_A = phase(k1+1:k1+k2, :);
magnitude = abs(Sm); mag_N = magnitude(1:k1, :); mag_A = magnitude(k1+1:k1+k2, :);

%% Calculate voxel-wise differences in variance 
 
        % ----------------(1) Spatial source phase---------------%
% Variance difference obtained by F-test (with and without FDR correction)
[diff, m_total, diff_FDR, m_total_FDR] = varAnalysis(phase_N, phase_A, alpha);

% Show the difference variance maps of SZ-HC
SMdisplay(diff, half_ind,location, color);
SMdisplay(diff_FDR, half_ind,location, color);

% Calculate m and signed q values for sub-networks 
for ind = 1 : length(sub_networks)
    [m(ind), q(ind)] = diffInNetwork(diff, half_ind, sub_networks(ind), 1);
    [m_FDR(ind), q_FDR(ind)] = diffInNetwork(diff_FDR, half_ind, sub_networks(ind), 1);
end

% Show histogram for m and q
mMin = 0; mMax = 1600; qMin = -1.2; qMax = 1.2;
fig_m = [m'  m_FDR']; fig_q = [q'  q_FDR'];
figure,bar(fig_m, width), ylim([mMin ,mMax]),colormap(color_map),
set(gca,'XTickLabel',sub_networks)
title('Phase: m'),legend('Without FDR','With FDR')
figure,bar(fig_q, width), ylim([qMin ,qMax]),colormap(color_map),
set(gca,'XTickLabel',sub_networks)
title('Phase: signed q'),legend('Without FDR','With FDR')
%-------------------------------------------------------------------------%

          % -------------------(2) Magnitude------------------%
% Variance difference obtained by F-test (with and without FDR correction)         
[mag_diff, mag_m_total, mag_diff_FDR,  mag_m_total_FDR] = varAnalysis(mag_N, mag_A, alpha);

% Show the difference variance maps of SZ-HC
SMdisplay(mag_diff,half_ind,location,color);
SMdisplay(mag_diff_FDR,half_ind,location,color);

% Calculate m and signed q values for sub-networks 
for ind= 1 : length(sub_networks)
    [mag_m(ind), mag_q(ind)] = diffInNetwork(mag_diff, half_ind, sub_networks(ind), 1);
    [mag_m_FDR(ind), mag_q_FDR(ind)] = diffInNetwork(mag_diff_FDR, half_ind, sub_networks(ind), 1);
end

% Show histogram for m and q
figm = [mag_m' mag_m_FDR']; fig_q = [mag_q' mag_q_FDR'];
figure,bar(fig_m, width), ylim([mMin,mMax]),colormap(color_map)
set(gca,'XTickLabel',sub_networks)
title('Magnitude: m'),legend('Without FDR','With FDR')
figure,bar(fig_q, width), ylim([qMin,qMax]),colormap(color_map)
set(gca,'XTickLabel',sub_networks)
title('Magnitude: signed q'),legend('Without FDR','With FDR')
%-------------------------------------------------------------------------%

%% Resampling validation

%---------------------------I: Bootstrap sampling-------------------------%

        % ----------------(1) Spatial source phase---------------%
for ind = 1 : length(sub_networks)
    % Calculate unsigned q values for actual data
    [m(ind), q(ind)] = diffInNetwork(diff, half_ind, sub_networks(ind), 0);
    [m_FDR(ind), q_FDR(ind)] = diffInNetwork(diff_FDR, half_ind, sub_networks(ind), 0); 
    
    % Calculate m and unsigned q values for bootstrap samples
    [bs_m_sub(ind,:), bs_q_sub(ind,:), bs_mG_sub(ind), bs_qG_sub(ind), bs_mstd_sub(ind), bs_qstd_sub(ind)] = diffInNetwork(bs_diff, half_ind, sub_networks(ind), 0);
    [bs_mFDR_sub(ind,:), bs_qFDR_sub(ind,:), bs_mFDRG_sub(ind), bs_qFDRG_sub(ind), bs_mFDRstd_sub(ind), bs_qFDRstd_sub(ind)] = diffInNetwork(bs_diff_FDR, half_ind, sub_networks(ind), 0);
    
    % Calculate the p-values
    [loo(ind, :), loo_FDR(ind, :)] = loo_pre(phase_N, phase_A, half_ind, sub_networks(ind), alpha);
    [p_sub(ind),~] = BCa_bootstrap(q(ind),loo(ind, :),bs_q_sub(ind,:), null);
    [pFDR_sub(ind),~] = BCa_bootstrap(q_FDR(ind),loo_FDR(ind, :),bs_qFDR_sub(ind,:), null);
end

% Show probability density estimate for q
figure,
for ind = 1 : length(sub_networks)
    hold on,
    [f,xi] = ksdensity(bs_q_sub(ind,:));
    plot(xi,f, cell2mat(strcat('-', color_line(ind)))); 
end
xlabel('q'),ylabel('probability density estimate');
title('Without FDR'), legend(sub_networks,'Location','best');
figure,
for ind = 1 : length(sub_networks)
    hold on,
    [f,xi] = ksdensity(bs_qFDR_sub(ind,:));
    plot(xi,f, cell2mat(strcat('-', color_line(ind)))); 
end
xlabel('q'),ylabel('probability density estimate');
title('With FDR'), legend(sub_networks,'Location','best');
%-------------------------------------------------------------------------%

          % -------------------(2) Magnitude------------------%
% for ind = 1 : length(sub_networks)
%     % Calculate unsigned q values for actual data
%     [mag_m(ind), mag_q(ind)] = diffInNetwork(mag_diff, half_ind, sub_networks(ind), 0);
%     [mag_m_FDR(ind), mag_q_FDR(ind)] = diffInNetwork(mag_diff_FDR, half_ind, sub_networks(ind), 0);
%     
%     % Calculate m and unsigned q values for bootstrap samples
%     [bs_mag_m_sub(ind,:), bs_mag_q_sub(ind,:), bs_mag_mG_sub(ind), bs_mag_qG_sub(ind), bs_mag_mstd_sub(ind), bs_mag_qstd_sub(ind)] = diffInNetwork(bs_mag_diff, half_ind, sub_networks(ind), 0);
%     [bs_mag_mFDR_sub(ind,:), bs_mag_qFDR_sub(ind,:), bs_mag_mFDRG_sub(ind), bs_mag_qFDRG_sub(ind), bs_mag_mFDRstd_sub(ind), bs_mag_qFDRstd_sub(ind)] = diffInNetwork(bs_mag_diff_FDR, half_ind, sub_networks(ind), 0);
%     
%     % Calculate the p-values
%     % [loo(ind, :), loo_FDR(ind, :)] = loo_pre(mag_N, mag_A, half_ind, sub_networks(ind), alpha);
%     % [mag_p_sub(ind),~] = BCa_bootstrap(mag_q(ind),loo(ind, :),bs_mag_q_sub(ind,:), null);
%     % [mag_pFDR_sub(ind),~] = BCa_bootstrap(mag_q_FDR(ind),loo_FDR(ind, :),bs_mag_qFDR_sub(ind,:), null);
% end
% 
% % Show probability density estimate for q
% figure,
% for ind = 1 : length(sub_networks)
%     hold on,
%     [f,xi] = ksdensity(bs_mag_q_sub(ind,:));
%     plot(xi,f, cell2mat(strcat('-', color_line(ind)))); 
% end
% xlabel('q'),ylabel('probability density estimate');
% title('Without FDR'), legend(sub_networks,'Location','best');
% figure,
% for ind = 1 : length(sub_networks)
%     hold on,
%     [f,xi] = ksdensity(bs_mag_qFDR_sub(ind,:));
%     plot(xi,f, cell2mat(strcat('-', color_line(ind)))); 
% end
% xlabel('q'),ylabel('probability density estimate');
% title('With FDR'), legend(sub_networks,'Location','best');
%-------------------------------------------------------------------------%

%-----------------------II: intra- and inter- sampling--------------------%

% Initialize parameters
group = {'SZ','HC','HS'};  % HS_ for inter-group HC_ for intra-HC and SZ_ for intra-SZ
mTick = [0:100:1500]; qTick = [-1.4:0.2:1.4];

        % ----------------(1) Spatial source phase---------------%
for gind = 1:length(group)
    for ind = 1 : length(sub_networks)
        [~, ~, mG_sub(ind), qG_sub(ind), mstd_sub(ind), qstd_sub(ind)] = diffInNetwork(eval([cell2mat(group(gind)), '_diff']), half_ind, sub_networks(ind), 1);
        [~, ~, mFDRG_sub(ind), qFDRG_sub(ind), mFDRstd_sub(ind), qFDRstd_sub(ind)] = diffInNetwork(eval([cell2mat(group(gind)), '_diff_FDR']), half_ind, sub_networks(ind), 1);
    end
    % Show histogram of m and q
    fig_m = [mG_sub' mFDRG_sub']; err_m = [mstd_sub' mFDRstd_sub'];
    fig_q = [qG_sub' qFDRG_sub']; err_q = [qstd_sub' qFDRstd_sub'];
    leg={'Without FDR','With FDR'};  XTickLabel = sub_networks;
    figure, barweb(fig_m,err_m,width,XTickLabel,strcat('Phase ',group(gind),': m'), ' ', 'm',color_map,'none',leg,error_sides,'plot');
    set(gca,'yTick', mTick); ylim([0 1500]);
    figure, barweb(fig_q,err_q,width,XTickLabel,strcat('Phase ',group(gind),': q'), ' ', 'q',color_map,'none',leg,error_sides,'plot');
    set(gca,'yTick', qTick); ylim([-1.4 1.4]);
end
%-------------------------------------------------------------------------%

          % -------------------(2) Magnitude------------------%
for gind = 1:length(group)
    for ind = 1 : length(sub_networks)
        [~, ~, mG_sub(ind), qG_sub(ind), mstd_sub(ind), qstd_sub(ind)] = diffInNetwork(eval([cell2mat(group(gind)), '_mag_diff']), half_ind, sub_networks(ind), 1);
        [~, ~, mFDRG_sub(ind), qFDRG_sub(ind), mFDRstd_sub(ind), qFDRstd_sub(ind)] = diffInNetwork(eval([cell2mat(group(gind)), '_mag_diff_FDR']), half_ind, sub_networks(ind), 1);
    end
    % Show histogram of m and q
    fig_m = [mG_sub' mFDRG_sub']; err_m = [mstd_sub' mFDRstd_sub'];
    fig_q = [qG_sub' qFDRG_sub']; err_q = [qstd_sub' qFDRstd_sub'];
    leg={'Without FDR','With FDR'};  XTickLabel = sub_networks;
    figure, barweb(fig_m,err_m,width,XTickLabel,strcat('Magnitude ',group(gind),': m'), ' ', 'm',color_map,'none',leg,error_sides,'plot');
    set(gca,'yTick', mTick); ylim([0 1500]);
    figure, barweb(fig_q,err_q,width,XTickLabel,strcat('Magnitude ',group(gind),': q'), ' ', 'q',color_map,'none',leg,error_sides,'plot');
    set(gca,'yTick', qTick); ylim([-1.4 1.4]);
end

