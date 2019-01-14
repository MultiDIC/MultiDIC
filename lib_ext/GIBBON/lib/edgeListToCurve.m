function [indList]=edgeListToCurve(E)

% function [indList]=edgeListToCurve(E)
% ------------------------------------------------------------------------
%
% ------------------------------------------------------------------------

%% Check for empty case

if isempty(E)
    indList=[];
    return;
end

%% Check for errors
if any(E(:,1)==E(:,2))
    ind=find(E(:,1)==E(:,2));
    error(['Invalid edges encountered, edge(s) "',sprintf('%d, ',ind),'" have the same start and end point']);
end

%%
if size(E,1)==1
    indList=E;
else
    [~,indV,~]=tesIND(E,[],0);
    logicEndPoints=sum(indV>0,2)==1;
    
    if nnz(logicEndPoints)==0 %closed loop
        indStartPoint=E(1,1);
    else
        indStartPoint=find(logicEndPoints,1);
    end
    indStartEdge=find(any(E==indStartPoint,2),1);
    E_now=E(indStartEdge,:);
    
    indUni=unique(E(:));
    
    ind1=E(indStartEdge,E_now==indStartPoint);
    ind2=E(indStartEdge,E_now~=indStartPoint);
    indList=ones(1,numel(indUni));
    indList(1)=ind1;
    indList(2)=ind2;
    q=3;
    Es=E;
    Es(indStartEdge,:)=NaN;
    while 1
        [indE_next,~]=find(Es==ind2,1);
        E_now=Es(indE_next,:);
        Es(indE_next,:)=NaN;
        ind3=E_now(E_now~=ind2);
        indList(q)=ind3;
        ind2=ind3;
        if nnz(isnan(Es))==numel(Es)
            break
        end
        q=q+1;
    end
    
    %% Fix order
    logicFlip=(E(:,1)==indList(1)) & (E(:,2)==indList(2));
    if ~any(logicFlip)
        indList=flipud(indList); %invert curve to conform to edge directions
    end
end

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
