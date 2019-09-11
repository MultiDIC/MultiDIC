function [xl,yl,zl]=axesLimits(V)
%% calculate axes limits for plotting
% V can be a nX3 matrix of 3D vertices or a cell array which contains nX3
% matrices representing different time frames.
%%
if ~iscell(V)
    Vtemp=V;
    V=cell(1);
    V{1}=Vtemp;
end


% xl=[0 0]; yl=[0 0]; zl=[0 0];
xl=[min(V{1}(:,1)) max(V{1}(:,1))]; yl=[min(V{1}(:,2)) max(V{1}(:,2))]; zl=[min(V{1}(:,3)) max(V{1}(:,3))];
for it=1:numel(V)
    xl(1)=min([min(V{it}(:,1)) xl(1)]);
    xl(2)=max([max(V{it}(:,1)) xl(2)]);
    yl(1)=min([min(V{it}(:,2)) yl(1)]);
    yl(2)=max([max(V{it}(:,2)) yl(2)]);
    zl(1)=min([min(V{it}(:,3)) zl(1)]);
    zl(2)=max([max(V{it}(:,3)) zl(2)]);
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