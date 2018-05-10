function L = DLT11Calibration(P2,P3)
%% function for calculating the DLT parameters in step 1
%
% INPUT:
% P2: a Nx2 array of N 2D image points
% P3: a Nx3 array of N 3D world points
% The points must be corresponded (the indeces of the points have to match)
%
% OUTPUT:
% L: an 11x1 array representing the 11 DLT parameters
% These parameters are to be used for stereo calibration of corresponded image points

%%
if size(P2,1)~=size(P3,1)
    error('Number of points in both matrices must match');
end
if size(P2,2)~=2
    error('Size of first input must be Npx2');
end
if size(P3,2)~=3
    error('Size of first input must be Npx3');
end

N=size(P3,1);

P2array(1:2:2*N-1,1)=P2(:,1); % a vector containing both coordiantes of the image points [u1, v1, u2, v2, u3, v3,...]
P2array(2:2:2*N,1)=P2(:,2);

% Matrix for solving DLT calibration parameters (2N x 11)
M(1:2:2*N-1,:)=[P3 ones(N,1) zeros(N,4) -P2(:,1).*P3];
M(2:2:2*N,:)=[zeros(N,4) P3 ones(N,1)  -P2(:,2).*P3];

L=M\P2array; % The solution for the 11 DLT parameters


end

 
%% 
% MultiDIC: a MATLAB Toolbox for Multi-View 3D Digital Image Correlation
% 
% License: <https://github.com/MultiDIC/MultiDIC/blob/master/LICENSE.txt>
% 
% Copyright (C) 2018  Dana Solav
% 
% If you use the toolbox/function for your research, please cite our paper:
% <https://engrxiv.org/fv47e>