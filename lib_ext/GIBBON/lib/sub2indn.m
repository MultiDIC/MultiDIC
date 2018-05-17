function ind = sub2indn(siz,A)

% function ind = sub2indn(siz,A)
% ------------------------------------------------------------------------
% This function is similar to MATLAB's sub2ind function. However the input
% subscrip indices are provided as an nxm array, such that numel(siz)=m and
% numel(ind)=n. The output are the equivalent linear indices
%
%   See also: ind2subn ind2sub sub2ind.
% ------------------------------------------------------------------------

%%

siz = double(siz);
numDim = numel(siz);
k = cumprod(siz);

if numDim < 2
    error('Invalid size specified. The number of dimensions should be equal or larger than 2');
end

if size(A,2) ~= numDim
    error('The specified array size and number of subscript index columns do not match');
end

if any(A(:)<1) || any(max(A,[],1)>siz)
    %Verify subscripts are within range
    error('Index out of range');
end

ind=A(:,1);
for q=2:1:numDim    
    ind = ind + (double(A(:,q))-1)*k(q-1);
end
