function P3D = DLT11Reconstruction(P1,P2,L1,L2)
%% function for reconstructing 2d points into 3d using the DLT parameters in step 1p and step 3
%
% INPUT:
% P1: a Nx2 array of N 2D image points from camera 1
% P2: a Nx2 array of N 2D image points from camera 2. The points must be corresponded (the indeces of the points have to match)
% L1: a 11x1 array representing the 11 DLT parameters of camera 1
% L2: a 11x1 array representing the 11 DLT parameters of camera 2
%
% OUTPUT:
% P3D:  a Nx3 array of N reconstructed 3D world points
%
%%

P3D=zeros(size(P1,1),3);

for ii=1:size(P1,1) % loop over number of mutual points
    if isnan(P1(ii,1)) || isnan(P2(ii,1))
        P3D(ii,:)= [NaN NaN NaN];
    else
    M =[P1(ii,1)*L1(9)-L1(1)  P1(ii,1)*L1(10)-L1(2)   P1(ii,1)*L1(11)-L1(3);
        P1(ii,2)*L1(9)-L1(5)  P1(ii,2)*L1(10)-L1(6)   P1(ii,2)*L1(11)-L1(7);
        P2(ii,1)*L2(9)-L2(1)  P2(ii,1)*L2(10)-L2(2)   P2(ii,1)*L2(11)-L2(3);
        P2(ii,2)*L2(9)-L2(5)  P2(ii,2)*L2(10)-L2(6)   P2(ii,2)*L2(11)-L2(7)];
    
    V= [L1(4)-P1(ii,1)
        L1(8)-P1(ii,2)
        L2(4)-P2(ii,1)
        L2(8)-P2(ii,2)];
       
    P3D(ii,:)=M\V;
    end
end

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