function varargout=tesgroup(varargin)

% function [G,G_iter]=tesgroup(F,optionStruct)
% ------------------------------------------------------------------------
%
% This function finds groups in the tesselation defined by F. F may
% represent patch type faces or for instances node indices for
% tetrehedrons, hexahedrons. Row entries in F (e.g. tetrahedron vertex
% indices) which are "connected" (sharing vertex indices with other row
% entries in F) are grouped together. The output G is a logic matrix of
% size(F,1) rows and "number of groups" columns. Each column represents a
% group and ones appear in the column for each face belonging to the group.
%
%
% Kevin Mattheus Moerman
% gibbon.toolbox@gmail.com
%
% 2010/07/15 Created
% 2018/11/05 Cleaned up some comments
% 2018/11/05 Added exclude points 
% 2018/11/21 Added option structure input
% 2018/11/21 Added label ouput type which is more efficient when there are
% many groups as it avoids the creation of large arrays. 
%------------------------------------------------------------------------

%% Parse input

switch nargin
    case 1
        F=varargin{1};
        optionStruct=[];
    case 2
        F=varargin{1};
        optionStruct=varargin{2};
end

%Check optionStruct against default
defaultOptionStruct.indExclude=[];
defaultOptionStruct.outputType='array';
[optionStruct]=structComplete(optionStruct,defaultOptionStruct,1); %Complement provided with default if missing or empty

%Get variables from structure
indExclude=optionStruct.indExclude;
outputType=optionStruct.outputType;

%%

if isempty(indExclude)
    excludeOption=0; 
else
    excludeOption=1; 
end

IND_F=(1:1:size(F,1))';
IND_F_search=IND_F;
G=zeros(size(F,1),1);
L=false(size(F,1),1);
L_previous=false(size(F,1),1);

G_ind=nan(size(F,1),1);
v_search=[ ];
done=0;
num_v_search=0;
group_found=1;
group_n=0;
q=1; %Counter
while done==0 
    L=false(size(F,1),1);
    if group_found==1
        indNext=find(IND_F_search>0,1); %next un-grouped element
        L(indNext)=1;
        IND_F_search(L)=0; %Setting found to zero
        group_found=0;
    else
        L = any(ismember(F,v_search), 2);
        IND_F_search(L)=0; %Setting found to zero
    end
    v_new=F(L,:);
    
    %Remove exclude points
    if excludeOption        
        v_new=v_new(~ismember(v_new,indExclude)); 
    end
    
    v_search=unique([v_search; v_new(:)]); %Growing number of search vertices    
    
    G_ind(L&~L_previous)=q;
    L_previous=L; 
    
    if numel(v_search)==num_v_search %If the group has not grown
        group_found=1;        
        group_n=group_n+1;
        switch outputType
            case 'array'
                G(:,group_n)=L;        
            case 'label'
                G(L)=group_n;        
        end
        v_search=[ ];
    end
    
    if all(IND_F_search==0)
        done=1;
        group_found=1;        
        group_n=group_n+1;
        if ~all(any(G,2))%if not all points are grouped keep remainder
            G(:,group_n)=L;
        end        
        v_search=[ ];
    end    
    num_v_search=numel(v_search);
    
    q=q+1; %Increment counter
end

switch outputType
    case 'array'
        G=G>0;
    case 'label'
        
end
%%

G_iter=G_ind(:,ones(size(G,2),1));
G_iter(~G)=NaN;
G_iter_min=nanmin(G_iter,[],1);
G_iter=G_iter-G_iter_min(ones(size(G_iter,1),1),:)+1;

%% Prepare output
switch nargout
    case 1
        varargout{1}=G;
    case 2
        varargout{1}=G;
        varargout{2}=G_iter;
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
