function [C,J] = findCentroids(IMbw,Np)
%% function to compute the centroids of black dots on white background from an image in step 1
% used inside plot_centroids_threshold.
% The function receives a black&white image and a number of dots and
% returns the dots' centroid positions.
% 
% INPUT:
% * IMbw: black&white (binary) image
% * Np: number of dots (black areas) on the image to be identified
% 
% OUTPUT:
% * C: the dots' centroid positions
% * J: The negative image with the dots regions

%%
%remove salt and pepper noise
J = medfilt2(IMbw,[3 3]);
%remove black pixels at the image corners
J(1,1)=1; J(1,end)=1; J(end,end)=1; J(end,1)=1;
% turn 0 into 1 and viceversa
J=~J;

%label the regions in the image (each "dot" is a region labeled with a different number in Regions)
[Regions,RegionsNum] = bwlabel(J);

%find the geometrical properties (in our case the centroids and areas) of regions in the image
Centroids = regionprops(Regions,'Centroid'); % a struct of num fields, each one is a 2D centroid position
Areas = regionprops(Regions,'Area');
C=zeros(numel(Centroids),2);
for k=1:numel(Centroids)
    C(k,:)=Centroids(k).Centroid;
end

% sort by area of regions
AreasArray=struct2array(Areas);
[~,indA]=sort(AreasArray);

% delete the centroids with the smallets area
if RegionsNum~=Np
%     warning('wrong number of points');
    C(indA(1:RegionsNum-Np),:)=[];
end


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