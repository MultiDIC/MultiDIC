function [ind]=inddiag(A)

i=1:size(A,2);
if size(A,2)>size(A,1)
    i=i(1:size(A,2));
end
ind=sub2ind(size(A),i,i);