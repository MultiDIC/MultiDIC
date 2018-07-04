function []=addQuiverFactorButton(varargin)
%% Viewing figures: changing quiver factor
% addQuiverFactorButton
% addQuiverFactorButton(hf)

% This script creates a new pushtool in a figure's toolbar to manipulate
% the quiver factor of the figure images
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

D=importdata(fullfile(iconPath,'quiver.png'));
S=double(D);
S=S-min(S(:));
S=S./max(S(:));
S(S==1)=NaN;
if size(S,3)==1
    S=repmat(S,[1 1 3]);
end

% Create a uipushtool in the toolbar
uipushtool(hb(1),'TooltipString','quiver','CData',S,'Tag','quiverFactor_button','ClickedCallback',{@quiverFactorFunc,{hf}});

end

%% Quiver factor function  quiverFactorFunc

function quiverFactorFunc(~,~,inputCell)

hf = inputCell{1};

prompt = {'Set quiver scale factor'};
dlg_title = 'Set quiver scale factor';

allQuiver = findall(hf,'Type','quiver');
if ~isempty(allQuiver)
    for ic = 1:size(allQuiver,1)
        allQuiver(ic).AutoScale = 'on';
    end
    currentFactor = allQuiver(1).AutoScaleFactor;
    defaultOptions = {num2str(currentFactor)};
    s = 25+max([cellfun(@numel,prompt) cellfun(@numel,defaultOptions)]);
    
    Q = inputdlg(prompt,dlg_title,[1 s],defaultOptions);
    if ~isempty(Q)
        newFactor = str2double(Q{1});
        
        for ic = 1:size(allQuiver,1)
            allQuiver(ic).AutoScaleFactor = newFactor;
        end
    end
else
    msgbox('There is no quiver in the figure');
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