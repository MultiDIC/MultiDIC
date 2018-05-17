function [C]=structComplete(A,B,emptyFixOpt)

% function [C]=structComplete(A,B,emptyFixOpt)
% ------------------------------------------------------------------------
% This function fills in the missing data in the structure A with the
% content from B. The structure B can be seen as a default input structure
% for instance and A might be an incomplete set of inputs. The optional
% parameter emptyFixOpt (0=no, 1=yes) determines whether empty entries in A
% are overwritten by the corresponding default values in B. 
%
% Kevin Mattheus Moerman
% gibbon.toolbox@gmail.com
% 
% 2017/11/20
%------------------------------------------------------------------------

%%

fieldNameSet=fieldnames(B);

C=A; %Initialize C as the same as A
for q=1:1:numel(fieldNameSet) %Loop over field names
    fieldNameNow=fieldNameSet{q}; %Current field name
    if isfield(A,fieldNameNow) %If A contains the field in B                
        
        if isempty(A.(fieldNameNow)) %if the field in A is empty
            if emptyFixOpt==1 %if empty field fixing is on
                C.(fieldNameNow)=B.(fieldNameNow); %Replace empty in C by default
            end
        end    
        
        if isstruct(A.(fieldNameNow)) && isstruct(B.(fieldNameNow)) %If the field in A is a structure, check structure recursively
            [C.(fieldNameNow)]=structComplete(A.(fieldNameNow),B.(fieldNameNow),emptyFixOpt);
        end
        
    else %If the field is missing, add it
        C.(fieldNameNow)=B.(fieldNameNow);            
    end
end

%%
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
