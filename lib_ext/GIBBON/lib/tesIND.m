function [varargout]=tesIND(varargin)

% function [IND_F,IND_V,IND_FF]=tesIND(F,V,sparseOpt);
% ------------------------------------------------------------------------
%
%
% Kevin Mattheus Moerman
% gibbon.toolbox@gmail.com
% 
% 2015/02/23
%------------------------------------------------------------------------

%% PARSE INPUT

switch nargin
    case 1 
        TES=varargin{1};
        sparseOpt=0;
        V=[];
    case 2 %Get size from vertices
        TES=varargin{1};
        V=varargin{2};        
        sparseOpt=0;               
    case 3
        TES=varargin{1};
        V=varargin{2};
        sparseOpt=varargin{3};                
    otherwise
        error('Wrong number of input arguments');
end

%Check numVertices
if isempty(V)
    numV=max(TES(:)); %Assuming tesselation is "clean"
else
    numV=size(V,1);
end

sizTES2=size(TES,2);

%% COMPUTE VERTEX-FACE CONNECTIVITY

I=TES(:); %The vertex indices in the face matrix as column
J=(1:1:size(TES,1))';
J=J(:,ones(1,sizTES2)); %Face indices copied for each vertex entry
J=J(:); %Face indices as column

%The above generates a column of face numbers indices which can be compared
%to a column of corresponding vertex indices.  

%Creating a sparse array where on the rows I and columns J we place the
%number J. Rows indicate vertex index and the colum the face index
IND_F=sparse(I,J,J,numV,size(TES,1));

%% COMPUTE VERTEX-VERTEX CONNECTIVITY

IND_V=[];
if nargout>=2
    %Get edges matrix
    E=patchEdges(TES,1);
    E=[E;fliplr(E)];
    
    %Create IND_V
    IND_V=sparse(E(:,1),E(:,2),E(:,2),numV,numV);
end

%% COMPUTE FACE-FACE CONNECTIVITY

IND_FF=[];
if nargout==3
    
    % IND_FF=[];
    % for q=1:1:sizTES2;
    %     IND_FF=[IND_FF IND_F(TES(:,q),:)];
    % end
    
    nnzMax=size(TES,2).*full(nnz(IND_F));
    
    IND_FF=sparse([],[],[],size(TES,1),size(IND_F,2)*size(TES,2),nnzMax);
    
    for q=1:1:sizTES2
        ind=1+(q-1)*size(IND_F,2);
        IND_FF(:,ind:(ind+size(IND_F,2)-1))=IND_F(TES(:,q),:);
    end
    
    [I,~,S] = find(IND_FF);
    P=[I(:) S(:)];
    [~,indUni,~]=unique(P,'rows');
    I=I(indUni);
    S=S(indUni);
    L=I~=S;
    I=I(L);
    S=S(L);
    IND_FF = sparse(I,S,S,size(TES,1),size(TES,1));
end

%% COLLECT OUTPUT

if ~isempty(IND_F)
    if sparseOpt==0 %Crop and convert to full array
        IND_F=sort(IND_F,2);
        [~,J,~] = find(IND_F);
        IND_F=IND_F(:,min(J):end);
        IND_F=full(IND_F);
        
        %Move possible zeros to the final columns
        L=IND_F==0; 
        if any(L(:))            
            maxLevel=max(IND_F(:))+1;
            IND_F(L)=maxLevel;
            IND_F=sort(IND_F,2);
            IND_F(IND_F==maxLevel)=0;
        end
    end
end

if ~isempty(IND_V)
    if sparseOpt==0 %Crop and convert to full array
        IND_V=sort(IND_V,2);
        [~,J,~] = find(IND_V);
        IND_V=IND_V(:,min(J):end);
        IND_V=full(IND_V);
    end
end

if ~isempty(IND_FF)
    if sparseOpt==0 %Crop and convert to full array
        IND_FF=sort(IND_FF,2);
        [~,J,~] = find(IND_FF);
        IND_FF=IND_FF(:,min(J):end);
        IND_FF=full(IND_FF);
    end
end

varargout{1}=IND_F;
varargout{2}=IND_V;
varargout{3}=IND_FF;


 
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
