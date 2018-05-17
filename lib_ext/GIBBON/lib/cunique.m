function [varargout]=cunique(A)

% function [A_uni,ind1,ind2,Ac]=cunique(A)
%-------------------------------------------------------------------------
%This function is similar to MATLAB's unique function. There are three
%differences: 1) An additional 4th optional output is available providing
%the count, or number of occurances, for each element in the input array.
%2) The 2nd output mathces the size of the first input, 3) The 3rd output
%is reshaped to be the size of the input variable. 
%
% See also: unique
%
% Kevin Mattheus Moerman
% gibbon.toolbox@gmail.com
%
% 2018/03/21: Created
%-------------------------------------------------------------------------

%%

[A_uni,ind1,ind2]=unique(A);

varargout{1}=A_uni;
varargout{2}=reshape(ind1,size(A_uni));
varargout{3}=reshape(ind2,size(A));

if nargout==4
    [subInd] = ind2subn(size(A_uni),ind2);
    Ac=accumarray(subInd,ones(numel(ind2),1),size(A_uni));
    Ac=reshape(Ac(ind2),size(A));
    varargout{4}=Ac;
end
