%*****************************************************************
% Description: mask the spatial maps for shared voxels activated in more than half of all the subjects
% Usage:
%   [Sm, half_ind] = groupMask(cSm, mask_ind, Zr)
% Input: 
    % cSm: the corrected spatial map estimate at the best run with dimension 1 x V
	% mask_ind: mask for in-brain voxels with dimension 1 x V 
    % Zr: magnitude threshold
% Output: 
    % Sm: voxels survived magnitude and spatial sourse phase thresholds in more than half of subjects
    % half_ind: the index for voxels survived magnitude and spatial sourse phase thresholds in more than half of subjects
% Reference: "Spatial Source Phase: A New Feature for Identifying Spatial Differences 
%            Based on Complex-Valued Resting-State fMRI Data," submitted to Human Brain Mapping
% Date: October 2017
% Author: Yue Qiu
%*****************************************************************
function [Sm, half_ind] = groupMask(cSm, mask_ind, Zr)
K = size(cSm, 1);
S = ones(62336,82);
for k = 1 : K
    S(find(abs(angle(cSm(k,:))) > pi/4) ,k) = 0;
    S(find(abs(cSm(k,:)) < Zr), k) = 0;
end   
l = mean(S,2);
in = find(l >= 0.5); half_ind = mask_ind(in);
Sm = cSm(:,in);