function [D1]=triSurfSetDist(F1,V1,F2,V2,distMetric)

switch distMetric
    case 'dist'
        D1=minDist(V1,V2);
    case 'ray'
        [~,~,N1]=patchNormal(F1,V1);
        
        %Ray trace 1 onto 2
        optStruct.eps      = 1e-6;
        optStruct.triangle = 'two sided';
        optStruct.ray      = 'ray';
        optStruct.border   = 'normal';
        
        numSteps=size(V1,1);
        D1=nan(numSteps,1);
        c=1;
        hw=waitbar(c/numSteps,['Ray tracing...',num2str(round(100.*c/numSteps)),'%']);
        for q=1:1:numSteps
            v1=V1(q,:);
            n1=N1(q,:);
            [V_intersect,L_intersect,~] = triangleRayIntersection(v1(ones(size(F2,1),1),:),n1(ones(size(F2,1),1),:),V2,F2,optStruct);
            V_intersect=V_intersect(L_intersect,:);
            if nnz(L_intersect)>0
                d=dist(v1,V_intersect');
                d1=min(d,[],2);
                D1(q)=d1(:);
            end
            waitbar(c/numSteps,hw,['Ray tracing...',num2str(round(100.*c/numSteps)),'%']);
            c=c+1;
        end
        close(hw);
    case 'dist-ray'
        [D1d]=triSurfSetDist(F1,V1,F2,V2,'dist');
        [D1r]=triSurfSetDist(F1,V1,F2,V2,'ray');
        D1=nanmin([D1d D1r],[],2);     
    case 'near-norm'
        [~,indMin]=minDist(V1,V2);
        [~,~,Nv]=patchNormal(F2,V2);
        N=Nv(indMin,:);
        W=V1-V2(indMin,:);
        D1=abs(dot(N,W,2));
end
 
%% 
% _*GIBBON footer text*_ 
% 
% License: <https://github.com/gibbonCode/GIBBON/blob/master/LICENSE>
% 
% GIBBON: The Geometry and Image-based Bioengineering add-On. A toolbox for
% image segmentation, image-based modeling, meshing, and finite element
% analysis.
% 
% Copyright (C) 2018  Kevin Mattheus Moerman
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
