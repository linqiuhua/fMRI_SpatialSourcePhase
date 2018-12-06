%*****************************************************************
% Description: The best run selection by combining cross-run averaging and a one-sample t-test
% Usage:
%   [cSm, pSm, r, in] = bestRunSelection(cS,pS,mask_ind,alpha)
% Input: 
    % cS: the corrected spatial map estimate with dimension 1 x V.
    % pS: the corrected and denoised spatial map estimate with dimension 1 x V.
    % mask_ind: mask for in-brain voxels with dimension 1 x V.
    % alpha: A value specifying the significance level as (100*ALPHA)%.
% Output: 
    % cSm: the corrected spatial map estimate at the best run with dimension 1 x V.
    % pSm: the corrected and denoised spatial map estimate at the best run with dimension 1 x V.
    % r: the correlation coefficient between the component at the best run and mean reference.
    % in: the best run index.
% Reference: Kuang,L.D., Lin,Q.H., Gong,X.F., Cong,F., Sui,J., & Calhoun,V.D.(2018). 
%            Model order effects on ICA of resting-state complex-valued fMRI data: application to schizophrenia. 
%            Journal of Neuroscience Methods, 304, 24¨C38. 
% Date: May 2017
% Author: Li-Dan Kuang   
%*****************************************************************
function [cSm, pSm, r, in] = bestRunSelection(cS,pS,mask_ind,alpha)
cS = squeeze(cS); pS = squeeze(pS);
R = size(cS, 2);
% Calculate the t-value
pm = mean(pS , 2);   
varpS = zeros(length(pm), 1);
for run = 1 : R
    varpS = varpS + (pS(:,run) - pm).^2;  
end
stdpS = sqrt(varpS / R - 1);
tmap_ind = find(stdpS > eps);
divisor = squeeze(stdpS(tmap_ind))./sqrt(R - 1); 
tpS = zeros(1,length(mask_ind));
tpS(tmap_ind) = squeeze(pm(tmap_ind))./ divisor;  
% Generate reference by combining cross-run averaging and a one-sample t-test
St = mean(pS ,2);
Stt = tRef(St,tpS,mask_ind,alpha,R-1);
for run=1:R
    p = corrcoef(abs(Stt),abs(pS(:,run))); Rm(run) = (p(1,2));
end
[r,in] = max(Rm); 
cSm = cS(:,in); pSm = pS(:,in);
