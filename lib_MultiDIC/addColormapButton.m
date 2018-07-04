function []=addColormapButton(varargin)
%% Viewing figures: changing colormap
% addColormapButton
% addColormapButton(hf)

% This script creates a new pushtool in a figure's toolbar to manipulate
% the color scheme displayed in the figure images
%%
switch nargin
    case 0
        hf=gcf;
    case 1
        hf=varargin{1};
end

% Create icon
S=zeros(16,16,3);
S(1,:,:)=parula(16);
S(2,:,:)=parula(16);
S(3,:,:)=parula(16);
S(4,:,:)=parula(16);
S(5,:,:)=jet(16);
S(6,:,:)=jet(16);
S(7,:,:)=jet(16);
S(8,:,:)=jet(16);
S(9,:,:)=gray(16);
S(10,:,:)=gray(16);
S(11,:,:)=gray(16);
S(12,:,:)=gray(16);
S(13,:,:)=coldwarm(16);
S(14,:,:)=coldwarm(16);
S(15,:,:)=coldwarm(16);
S(16,:,:)=coldwarm(16);

hb = findall(hf,'Type','uitoolbar');

% Create a uipushtool in the toolbar
uipushtool(hb(1),'TooltipString','colormap','CData',S,'Tag','colormap_button','ClickedCallback',{@colormapFunc,{hf}});

end

%% Colormap function  colormapFunc

function colormapFunc(~,~,inputCell)

hf = inputCell{1};
prompt = {'Colormap'};
title = '';
format = struct('type','list');
format.style = 'popupmenu';
format.items = {' -- Select a color scheme -- ','coldwarm','parula','jet','hsv','hot','cool','gray','bone'};
default = cell(size(prompt,1),1);
default{1,1} = 1; % Default selection will always be instruction item
default = cell2struct(default,prompt,1);
prompt = repmat(prompt,1,2);
options.AlignControls = 'on';
choice = inputsdlg(prompt,title,format,default,options);

% Change color scheme for all existing images accordingly
if choice.Colormap == 2
    colormap(0.8*coldwarm);
elseif choice.Colormap == 3
    colormap(parula);
elseif choice.Colormap == 4
    colormap(jet);
elseif choice.Colormap == 5
    colormap(hsv);
elseif choice.Colormap == 6
    colormap(hot);
elseif choice.Colormap == 7
    colormap(cool);
elseif choice.Colormap == 8
    colormap(gray);
elseif choice.Colormap == 9
    colormap(bone);
end
end

%% 
% MultiDIC: a MATLAB Toolbox for Multi-View 3D Digital Image Correlation
% 
% License: <https://github.com/MultiDIC/MultiDIC/blob/master/LICENSE.txt>
% 
% Copyright (C) 2018  Dana Solav
% 
% Modified by Rana Odabas 2018
% 
% If you use the toolbox/function for your research, please cite our paper:
% <https://engrxiv.org/fv47e>