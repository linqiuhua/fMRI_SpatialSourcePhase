%*****************************************************************
% Description: phase de-noising based on thresholds for spatial source phase
% Usage:
%   pS = phase_denoising(cS)
% Input: 
    % cS: the corrected spatial map estimate with dimension 1 x V
% Output: 
    % pS: the corrected and denoised spatial map estimate with dimension 1 x V
% Reference: Yu,M.C., Lin,Q.H., Kuang,L.D., Gong,X.F., Cong,F., & Calhoun,V.D. (2015). 
%            ICA of full complex-valued fMRI data using phase information of spatial maps. 
%            Journal of Neuroscience Methods, 249, 75¨C91. 
% Date: December 2014
% Author: Mou-Chuan Yu
%*****************************************************************
function pS = phase_denoising(cS)
pS=zeros(size(cS));
out_p=angle(cS);
ind=find(out_p<(pi/4)&out_p>(-pi/4));
pS(ind)=cS(ind);

