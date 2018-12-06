%*****************************************************************
%Description: component selection based on the reference network
% Usage:
%   inn = componentSelection(S, Sr, C)
%Input: 
    % S: the spatial map estimate with dimension N x V£¬N is the number of components and V is the number of in-brain voxels
	% Sr: the spatial reference network with dimension 1 x V
    % C: number of candidate components with top C correlation coefficients with Sr
%Output: 
    % inn: the index of the selected component, and inn is less than N
% Note: The component selection is complicated, and there are a few mis-selections, so we manually checked the results and made a few corrections.
% Reference: "Spatial Source Phase: A New Feature for Identifying Spatial Differences 
%            Based on Complex-Valued Resting-State fMRI Data," submitted to Human Brain Mapping
% Date: August 2018
% Author: Yue Qiu
%*****************************************************************
function inn = componentSelection(S, Sr, C)
S = squeeze(S);
mask = Sr; mask(mask > 0) = 1;
N = size(S, 1);
for n = 1:N
    p = corrcoef(abs(S(n,:)),Sr); 
    bb(n) = abs(p(1,2));
end
[~,index] = sort(bb, 'descend');
ind = 1;
for n = index(1 : C)
    temp = abs(S(n,:)); temp_mask = mask .* temp;
    cc(ind) = length(find(temp_mask > 0.5)) / length(find(mask > 0)) * length(find(temp_mask > 0.5)) / length(find(temp > 0.5));
    ind = ind + 1;
end
[~,in] = max(cc);
inn = index(in);
