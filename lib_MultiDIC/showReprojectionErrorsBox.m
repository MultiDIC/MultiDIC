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