function P3D = DLT11Reconstruction(P1,P2,L1,L2)
% this function returns the reconstructed 3D coordinates of points as
% measured by 2 cameras P1 and P2, with the respective DLT
% 11 parameters vectors L1 and L2

P3D=zeros(size(P1,1),3);

for ii=1:size(P1,1) % loop over number of mutual points
    if isnan(P1(ii,1)) || isnan(P2(ii,1))
        P3D(ii,:)= [NaN NaN NaN];
    else
    M =[P1(ii,1)*L1(9)-L1(1)  P1(ii,1)*L1(10)-L1(2)   P1(ii,1)*L1(11)-L1(3);
        P1(ii,2)*L1(9)-L1(5)  P1(ii,2)*L1(10)-L1(6)   P1(ii,2)*L1(11)-L1(7);
        P2(ii,1)*L2(9)-L2(1)  P2(ii,1)*L2(10)-L2(2)   P2(ii,1)*L2(11)-L2(3);
        P2(ii,2)*L2(9)-L2(5)  P2(ii,2)*L2(10)-L2(6)   P2(ii,2)*L2(11)-L2(7)];
    
    V= [L1(4)-P1(ii,1)
        L1(8)-P1(ii,2)
        L2(4)-P2(ii,1)
        L2(8)-P2(ii,2)];
       
    P3D(ii,:)=M\V;
    end
end

end