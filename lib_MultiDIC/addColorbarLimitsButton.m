function []=addColorbarLimitsButton(varargin)
%% Viewing figures: changing colorbar limits
% addColorbarLimitsButton
% addColorbarLimitsButton(hf)

% This script creates a new pushtool in a figure's toolbar to manipulate
% the limits of all colorbar axes displayed in the figure
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

D=importdata(fullfile(iconPath,'colorbar.jpg'));
S=double(D);
S=S-min(S(:));
S=S./max(S(:));
S(S==1)=NaN;
if size(S,3)==1
    S=repmat(S,[1 1 3]);
end

% Create a uipushtool in the toolbar
uipushtool(hb(1),'TooltipString','colorbar','CData',S,'Tag','colorbar_button','ClickedCallback',{@colorbarFunc,{hf}});

end

%% Colorbar function  colorbarFunc

function colorbarFunc(~,~,inputCell)

hf = inputCell{1};

prompt = {'Colorbar minimum value:','Colorbar maximum value:'};
dlg_title = 'Set colorbar limits';

allAxes = findall(hf,'Type','axes');
if ~isempty(allAxes)
    currentLimits = caxis(allAxes(1));
    defaultOptions = {num2str(currentLimits(1)),num2str(currentLimits(2))};
    s = 25+max([cellfun(@numel,prompt) cellfun(@numel,defaultOptions)]);

    Q = inputdlg(prompt,dlg_title,[1 s],defaultOptions);
    if ~isempty(Q)
        minC = str2double(Q{1});
        maxC = str2double(Q{2});

        for ic = 1:size(allAxes,1)
            caxis(allAxes(ic),[minC maxC]);
        end
    end
else
    msgbox('There are no axes in the figure');
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