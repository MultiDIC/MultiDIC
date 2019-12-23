function [R,t,PfromTransformed]=rigidTransformation(Pfrom,Pto,varargin)
%% Find the optimal rigid body transformation between two point clouds %%%%%%%%%
% The optimal solution is the least square minimization solution between
% the two point clouds. See (Rubin & Solav 2016)
% The two point clouds must be corresponded (indeces refer to the same point)

%% INPUT:
% Pfrom=  point cloud (to be transformed into TruP)
% Pto= true 3D point cloud
% * optional input: 
% estimation method: 'rigid' (default) or 'affine'

%% OUTPUT:
% R = estimated rotation matrix from Pfrom to Pto (3-by-3)
% t = estimated translation vector (3-by-1)
% PfromTransformed = Pfrom transformed using R and t

%%
nArg=numel(varargin);

switch nArg
    case 0
        estimationMethod='rigid';
    case 1
        estimationMethod=varargin{1};
end



%%
a=Pfrom;
b=Pto; 

% find NaNs
aNanInd=find(isnan(a(:,1)));
bNanInd=find(isnan(b(:,1)));
abNanInd=[aNanInd; bNanInd];
% delete Nans
aNoNan=a;
aNoNan(abNanInd,:)=[];
bNoNan=b;
bNoNan(abNanInd,:)=[];

ac=(1/length(aNoNan))*sum(aNoNan); % centroid of  RecP.
bc=(1/length(bNoNan))*sum(bNoNan); % centroid of  TruP.

da=aNoNan-ac;
db=bNoNan-bc;

switch estimationMethod
    
    case 'rigid'
        M=db'*da;
        [U,~,V] = svd(M);
        S=[1 0 0; 0 1 0; 0 0 det(U*V')];
        R=U*S*V'; %rotation matrix from reconstructed coordinates to true
        t=bc'-R*ac'; %translation vector from centroid of reconstructed coordinates to true- this is the translation of the origin of the coordinate system!
    
    case 'affine'
        M=db'*da;
        A=da'*da;
        F=M/A;
%         R=F*(F'*F)^-.5; % this gives the same solution as the svd
        [U,~,V] = svd(F);
        S=[1 0 0; 0 1 0; 0 0 det(U*V')];
        R=U*S*V'; %rotation matrix from reconstructed coordinates to true
        t=bc'-R*ac'; %translation vector from centroid of reconstructed coordinates to true- this is the translation of the origin of the coordinate system!
        
end

PfromTransformed=(R*a'+t)'; %transformed constructed coordinates


% residuals=PfromTransformed-Pto;
% residualsMgn=sqrt(sum(residuals.^2,2));
% residulasMgnRMS=sqrt(nansum(residualsMgn.^2)/sum(~isnan(residualsMgn)));


end
 
%% 
% MultiDIC: a MATLAB Toolbox for Multi-View 3D Digital Image Correlation
% 
% License: <https://github.com/MultiDIC/MultiDIC/blob/master/LICENSE.txt>
% 
% Copyright (C) 2018  Dana Solav
% 
% If you use the toolbox/function for your research, please cite our paper:
% <https://engrxiv.org/fv47e>