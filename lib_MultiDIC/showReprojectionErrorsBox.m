function []=showReprojectionErrorsBox(camera_parameters)
%% function for plotting the reprojection error statistics as bloxplot, in STEP0
% The function is used inside the function plot_camera_parameters_2tabs 
%
% INPUTS:
% * camera_parameters:  a camera_parameters class saved inside cameraCBparameters structure
%%
errors=camera_parameters.ReprojectionErrors;

errorsAll=reshape(errors,size(errors,1)*size(errors,3),2);
errorsAll=[errorsAll sqrt(sum(errorsAll.^2,2))];
boxplot(errorsAll,'Labels',{'X','Y','Mgn'});

title('Reprojection errors statistics');
ylabel('Error [pixel]','fontsize',12);

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