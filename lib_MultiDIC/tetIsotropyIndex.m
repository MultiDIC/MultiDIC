function [tetIsoIndex]=tetIsotropyIndex(E,V)
%% function for computing  triangular faces isotropy index. 
% This is an index of the  "regularity" of the triangle. Equilateral triangle (perfectly regular) 
%
% INPUTS:
% * E: nElements-by-4 array representing the list of vertices for tetrahedral elements
% * V: nVertices-by-3 array representing the 3D positions of the vertices of E
%
% OUTPUTS: 
% * tetIsoIndex: isotropy index for each element. value of 1 means
% perfectly isotropic (equilateral tetrahedron). value of 0 means vertices are coplanar.
%
%From the paper: Surface-Marker Cluster Design Criteria for 3-D Bone Movement Reconstruction (1997)
% Aurelio Cappozzo, Angelo Cappello, Ugo Della Croce, and Francesco Pensalfini
%%

tetIsoIndex=zeros(size(E,1),1);
for itet=1:size(E,1)
    
    x1=V(E(itet,1),:)';
    x2=V(E(itet,2),:)';
    x3=V(E(itet,3),:)';
    x4=V(E(itet,3),:)';
    
    if any(any(isnan([x1 x2 x3 x4]))) % if any nan, iso index=nan
        tetIsoIndex(itet)=NaN;
    else
        xa=(x1+x2+x3+x4)/3;
        X=[x1-xa,x2-xa,x3-xa,x4-xa]; %cluster position model
        K=X*X'/4;
        Keig=real(eig(K)); %eigenvalues of X
        tetIsoIndex(itet)=3*Keig(1)/(Keig(1)+Keig(2)+Keig(3)); % isotropy index
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