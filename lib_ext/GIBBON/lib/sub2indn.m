function ind = sub2indn(varargin)

% function ind = sub2indn(siz,A,nanOutBound)
% ------------------------------------------------------------------------
% This function is similar to MATLAB's sub2ind function. However the input
% subscrip indices are provided as an nxm array, such that numel(siz)=m and
% numel(ind)=n. The output are the equivalent linear indices
%
%   See also: ind2subn ind2sub sub2ind.
%
% Change log
% 2018/05/17 Optionally output NaN values for out of bound indices
% ------------------------------------------------------------------------

%% Parse input

switch nargin
    case 2
        siz=varargin{1};
        A=varargin{2};
        nanOutBound=0;
    case 3
        siz=varargin{1};
        A=varargin{2};
        nanOutBound=varargin{3};
end

%%
siz = double(siz);
numDim = numel(siz);
k = cumprod(siz);

if numDim < 2
    error('Invalid size specified. The number of dimensions should be equal or larger than 2');
end

if any(size(A)==0)
    A=[];
end

if ~isempty(A)    
    
    if size(A,2) ~= numDim
        error('The specified array size and number of subscript index columns do not match');
    end

    if ~nanOutBound
        if any(A(:)<1) || any(max(A,[],1)>siz)
            %Verify subscripts are within range
            error('Index out of range');
        end
    else
        A(A<1)=NaN;
        A(A>siz(ones(size(A,1),1),:))=NaN;
    end
    
    ind=A(:,1);
    for q=2:1:numDim
        ind = ind + (double(A(:,q))-1)*k(q-1);
    end    
else
    ind=[];
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
% Copyright (C) 2019  Kevin Mattheus Moerman
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
