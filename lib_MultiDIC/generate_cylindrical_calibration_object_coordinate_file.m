% generate a 3d coordinate file for a cylindrical calibration object

prompt = {'Cylinder diameter [mm]','Number of rows:','Number of clumns:','Distance between rows [mm]:','Distance between columns [mm]:'};
title = 'Input calibration object parameters';
dims = [1 35];
% definput = {'20','hsv'};
answer = inputdlg(prompt,title);

D=str2num(answer{1});
Nr=str2num(answer{2});
Nc=str2num(answer{3});
dr=str2num(answer{4});
dc=str2num(answer{5});

% total number of point
Np=Nr*Nc;
% angle between columns
angInc=2*dc/D;

P3d=zeros(Nr,Nc,3);
for ir=1:Nr
    z=dr*(ir-1);
    for ic=1:Nc
        ang=angInc*(ic-1);
        x=.5*D*cos(ang);
        y=.5*D*sin(ang);
        P3d(ir,ic,:)=[x y z];
    end
end

P3dVec=reshape(P3d,Np,3);
%plot
cFigure; hold all; axisGeom;
plotV(P3dVec,'sb','MarkerFaceColor','b');

% save file
uisave('P3d','myCylinderCoordinates');


