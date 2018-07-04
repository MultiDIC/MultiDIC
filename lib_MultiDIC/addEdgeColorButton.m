function []=addEdgeColorButton(varargin)
%% Viewing figures: setting edge color
% addEdgeColorButton
% addEdgeColorButton(hf)

% This script creates a new pushtool in a figure's toolbar to manipulate
% the edge color of all patch objects on a displayed figure
%%
switch nargin
    case 0
        hf=gcf;
    case 1
        hf=varargin{1};
end

% Get icon
filePath=mfilename('fullpath');
toolboxPath=fileparts(fileparts(filePath));
iconPath=fullfile(toolboxPath,'lib_ext','GIBBON','icons');

hb = findall(hf,'Type','uitoolbar');
D=imread(fullfile(iconPath,'spline.png'));
S=double(D);
S=S-min(S(:));
S=S./max(S(:));
if size(S,3)==1
    S=repmat(S,[1 1 3]);
end

% Create a uipushtool in the toolbar
uipushtool(hb(1),'TooltipString','edge color','CData',S,'Tag','edgeColor_button','ClickedCallback',{@edgeColorFunc,{hf}});
end

%% Edge color function  edgeColorFunc

function edgeColorFunc(~,~,inputCell)

hf = inputCell{1};

% Change edge color for all existing images accordingly
hpatches = findobj(hf,'type','patch');
if ~isempty(hpatches)
    
    prompt = {'EdgeColor'};
    title = '';
    format = struct('type','list');
    format.style = 'popupmenu';
    format.items = {' -- Select an edge color -- ','none','black','white'};
    default = cell(size(prompt,1),1);
    default{1,1} = 1; % Default selection will always be instruction item
    default = cell2struct(default,prompt,1);
    prompt = repmat(prompt,1,2);
    options.AlignControls = 'on';
    choice = inputsdlg(prompt,title,format,default,options);
    
    numPatches = size(hpatches,1);
        if choice.EdgeColor == 2
            for i = 1:numPatches
                hpatches(i).EdgeColor = 'none';
            end
        elseif choice.EdgeColor == 3
            for i = 1:numPatches
                hpatches(i).EdgeColor = 'k';
            end
        elseif choice.EdgeColor == 4
            for i = 1:numPatches
                hpatches(i).EdgeColor = 'w';
            end
        end
else
    msgbox('There are no patch objects in the figure');
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