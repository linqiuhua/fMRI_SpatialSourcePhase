%*****************************************************************
% Description: calculate q and m in the network
% Usage:
%   [m, q, mG, qG, m_std, q_std] = diffInNetwork(diff, half_ind, area, signed) 
% Input: 
    % diff£ºthe difference variance maps of SZ-HC with dimension N x length(half_ind)
	% half_ind: the index for voxels survived magnitude and spatial sourse phase thresholds in more than half of subjects
    % area: the sub-network, e.g.,'PCC', 'ACC'
    % signed: 1 means signed q and 0 means unsigned q
% Output: 
    % m: number of voxels with significant difference in the area
    % q: variance difference index in the area
    % mG: the mean of m
    % qG: the mean of unsigned q
    % m_std: the standard deviation of m 
    % q_std: the standard deviation of unsigned q 
% Reference: "Spatial Source Phase: A New Feature for Identifying Spatial Differences 
%            Based on Complex-Valued Resting-State fMRI Data," submitted to Human Brain Mapping
% Date: October 2017
% Author: Yue Qiu
%*****************************************************************
function [m, q, mG, qG, m_std, q_std] = diffInNetwork(diff, half_ind, area, signed) 
if strcmp(area, 'PCC')
   indx = 19 : 34; indy = 1 : 40; indz = 1 : 46;
elseif strcmp(area, 'IPL')
   indx = [1 : 18 35 : 53]; indy = 1 : 40; indz = 1 : 46;
elseif strcmp(area, 'ACC')
   indx = 1 : 53; indy = 41 : 63; indz = 1 : 46;
elseif strcmp(area, 'AUDL')
   indx = 1 : 26; indy = 1 : 63; indz = 1 : 46;
elseif strcmp(area, 'AUDR')
   indx = 27 : 53; indy = 1 : 63; indz = 1 : 46;
end
    
Smm=zeros(size(diff, 1),53*63*46);Smm(:,half_ind) = diff;
Smm_sub=zeros(1,153594);
for i = 1 : size(diff, 1)
    out=zeros(53,63,46);out1=zeros(53,63,46);
    temp = Smm(i, :); out(:)=temp(:); 
    out1(indx, indy, indz)=out(indx, indy, indz); Smm_sub(:)=out1(:);
    m(i) = length(find(Smm_sub~=0));
    q(i) = length(find(Smm_sub > 0)) / length(find(Smm_sub~=0));
end
mG = mean(m); qG = mean(q); m_std = std(m); q_std = std(q);

if(signed == 1)
    q(q < 0.5) = q(q < 0.5) - 1;
    qG(qG < 0.5) = qG - 1;
end
