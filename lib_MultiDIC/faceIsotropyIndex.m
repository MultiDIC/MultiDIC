function [FisoInd]=faceIsotropyIndex(F,V)
%% function for computing  triangular faces isotropy index. 
% This is an index of the  "regularity" of the triangle. Equilateral triangle (perfectly regular) 
%
% INPUTS:
% * F: nFaces-by-3 array representing the list of vertices for triangular faces
% * V: nVertices-by-3 array representing the 3d positions of the vertices of F in the reference positions
%
% OUTPUTS: 
% * FisoInd: isotropy index for each face. value of 1 means perfectly
% equilateral. value of 0 means vertices are aligned (on a line)
%
%From the paper: Surface-Marker Cluster Design Criteria for 3-D Bone Movement Reconstruction (1997)
% Aurelio Cappozzo, Angelo Cappello, Ugo Della Croce, and Francesco Pensalfini
%%

FisoInd=zeros(size(F,1),1);
for iface=1:size(F,1)
    
    x1=V(F(iface,1),:)';
    x2=V(F(iface,2),:)';
    x3=V(F(iface,3),:)';
    
    if any(any(isnan([x1 x2 x3]))) % if any nan, iso index=nan
        FisoInd(iface)=NaN;
    else
        xa=(x1+x2+x3)/3;
        X=[x1-xa,x2-xa,x3-xa]; %cluster position model
        K=X*X'/3;
        Keig=real(eig(K)); %eigenvalues of X
        FisoInd(iface)=2*Keig(2)/(Keig(2)+Keig(3)); % isotropy index
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