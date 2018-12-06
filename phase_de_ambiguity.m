%*****************************************************************
% Description: phase de-ambiguity based on maximizing the energy of the TC real part
% Usage:
%   [cS,cA] = phase_de_ambiguity(At,S,Sr)
% Input: 
    % At: the time course estimate with dimension T x 1£¬T is the number of time points
    % S: the spatial map estimate with dimension 1 x V£¬V is the number of in-brain voxels
	% Sr: the spatial reference network with dimension 1 x V
% Output: 
    % cA: the corrected time course estimate with dimension T x 1
    % cS: the corrected spatial map estimate with dimension 1 x V
% Reference: Yu,M.C., Lin,Q.H., Kuang,L.D., Gong,X.F., Cong,F., & Calhoun,V.D. (2015). 
%            ICA of full complex-valued fMRI data using phase information of spatial maps. 
%            Journal of Neuroscience Methods, 249, 75¨C91. 
% Date: December 2014
% Author: Mou-Chuan Yu
%*****************************************************************
function [cS,cA] = phase_de_ambiguity(At,S,Sr)
i=sqrt(-1); [p,q]=size(S);
Re=real(At)';
Im=imag(At)';
A0=[Re;Im];
[V, Lambda] = icatb_v_pca(A0',1,2,0);
V=circshift(V,1);
AA=V*A0;
Ar=AA(1,:)';
Ai=AA(2,:)';
cA=Ar+i*Ai;
b=acos(V(1,1));
A_new1=At.*exp(i*b);
if ~eq((round(A_new1(1,1).*10000))./10000,(round(cA(1,1).*10000))./10000)
    b = -b;
    A_new1=At.*exp(i*b);
end
cS=S.*exp(i*(-b));

temp=real(cS);
temp1=corrcoef(temp,Sr);
corA=temp1(1,2);
if (corA>0)
    cA=cA;
    cS=cS;
else
    cA=-cA;
    cS=-cS;
end
