function []=addFaceLightingButton(varargin)
%% Viewing figures: setting Specular strength
% addSpecularStrengthButton
% addSpecularStrengthButton(hf)

% This script creates a new pushtool in a figure's toolbar to manipulate
% the Specular strength of a given patch of light on a displayed axis
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
D=imread(fullfile(iconPath,'smooth.png'));
S=double(D);
S=S-min(S(:));
S=S./max(S(:));
if size(S,3)==1
    S=repmat(S,[1 1 3]);
end

% Create a uipushtool in the toolbar
uipushtool(hb(1),'TooltipString','Face Lighting','CData',S,'Tag','FaceLighting_button','ClickedCallback',{@FaceLightingFunc,{hf}});
end

%% Face Lighting function  FaceLightingFunc


function FaceLightingFunc(~,~,inputCell)

hf = inputCell{1};

hpatches = findobj(hf,'type','patch');
if ~isempty(hpatches)
    numPatches = size(hpatches,1);
    for ii = 1:numPatches
        default = hpatches(ii).FaceLighting;
    end
    
    prompt = {'FaceLight'};
    title = '';
    format = struct('type','list');
    format.style = 'popupmenu';
    format.items = {' -- Select a lighting option -- ','flat','gouraud','none'};
    default = cell(size(prompt,1),1);
    default{1,1} = 1; % Default selection will always be instruction item
    default = cell2struct(default,prompt,1);
    prompt = repmat(prompt,1,2);
    options.AlignControls = 'on';
    choice = inputsdlg(prompt,title,format,default,options);
    
    % Unless the user closes the pop up window, clear all existing light
    if choice.FaceLight ~= 1
        % Apply selection
        if choice.FaceLight == 2
            for ii = 1:numPatches
                hpatches(ii).FaceLighting = 'flat';
            end
        elseif choice.FaceLight == 3
            for ii = 1:numPatches
                hpatches(ii).FaceLighting = 'gouraud';
            end
        elseif choice.FaceLight == 4
            for ii = 1:numPatches
                hpatches(ii).FaceLighting = 'none';
            end
        end
    end
else
    msgbox('There are no light patches in the figure');
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