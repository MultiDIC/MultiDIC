function []=showReprojectionPointsErrors(camera_parameters)
%% function for plotting the reprojection error of all corner points as scatter, in STEP0
% The function is used inside the function plot_camera_parameters_2tabs 
%
% INPUTS:
% * camera_parameters:  a camera_parameters class saved inside cameraCBparameters structure
%%
errors=camera_parameters.ReprojectionErrors;
for ii=1:size(camera_parameters.ReprojectionErrors,3)   
    scatter(errors(:,1,ii),errors(:,2,ii),'+'); hold on
end

axis equal
eMax= max([max(max(abs(errors(:,1,:)))) max(max(abs(errors(:,2,:))))]);
xlim([-eMax eMax]);
ylim([-eMax eMax]);
refline(0,0);
line([0,0],ylim)
xlabel('x [pixel]','fontsize',14);
ylabel('y [pixel]','fontsize',14);

title('Reprojection errors distribution');

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