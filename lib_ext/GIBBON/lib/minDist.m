function varargout=minDist(varargin)

% function [D1,minIND]=minDist(V1,V2,maxVarSize,selfAvoid,numFreeBytes)

%% Parse input
if nargin<2
    error('Insufficient input arguments');
end

V1=varargin{1};
V2=varargin{2};
switch nargin
    case 2        
        maxVarSize=[]; %Empty will force calcucation below
        selfAvoid=0; 
        numFreeBytes=[];
    case 3
        maxVarSize=varargin{3};
        selfAvoid=0; 
        numFreeBytes=[];
    case 4
        maxVarSize=varargin{3};
        selfAvoid=varargin{4}; 
        numFreeBytes=[];
    case 5
        maxVarSize=varargin{3};
        selfAvoid=varargin{4};
        numFreeBytes=varargin{5};        
end

%Get free memory
if isempty(numFreeBytes)
    [numFreeBytes]=freeMemory;
end

%Get max variable size available        
if isempty(maxVarSize)    
    maxVarSize=numFreeBytes/2;
end

if isnan(maxVarSize)
    numSteps=1;
else
    %Derive class dependent variable size
    [~,b1]=maxnumel(V1(1),numFreeBytes);
    [~,b2]=maxnumel(V2(1),numFreeBytes);
    b=max([b1 b2]);
    numelVar=numel(V1)*numel(V2);
    varSize=numelVar*b;
    
    numSteps=ceil(varSize/maxVarSize);
    indSteps=round(linspace(0,size(V1,1),numSteps));
    indSteps=sort(unique(indSteps));
    numSteps=numel(indSteps);
end

if numSteps>1 %In steps
    D1=zeros(size(V1,1),1);
    minIND=zeros(size(V1,1),1);
    for q=1:1:numSteps-1
        v1=V1(indSteps(q)+1:indSteps(q+1),:);
        try 
            d=dist(v1,V2'); %dist from Neural network toolbox
        catch
            d=distND(v1,V2); %GIBBON's dist function
        end
        if selfAvoid
            %Set "diagonal" to something too large so self is avoided in
            %minimum (could use NaN and nanmin but the latter is a toolbox
            %function)
            I=1:size(v1,1);
            J=indSteps(q)+1:indSteps(q+1);
            ind=sub2ind(size(d),I,J); %Indices of selfies            
            d(ind)=1+max(d(:)); %Overide selfies
        end
        
        [min_d,min_ind]=min(d,[],2);
        D1(indSteps(q)+1:indSteps(q+1))=min_d;
        minIND(indSteps(q)+1:indSteps(q+1))=min_ind;        
    end
else %In one go
    try
        D=dist(V1,V2'); %dist from Neural network toolbox
    catch
        D=distND(V1,V2); %GIBBON's dist function
    end    
    if selfAvoid
        %Set "diagonal" to something too large so self is avoided in
        %minimum (could use NaN and nanmin but the latter is a toolbox
        %function) 
        L=eye(size(D))>0;
        D(L)=1+max(D(:));
    end
    [D1,minIND]=min(D,[],2);         
    D1=D1(:);
    minIND=minIND(:);
end

switch nargout
    case 1
        varargout{1}=D1;
    case 2
        varargout{1}=D1;
        varargout{2}=minIND;
    otherwise
        error('wrong number of output arguments');
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
