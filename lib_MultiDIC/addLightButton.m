function []=addLightButton(varargin)
%% Viewing figures: projecting light
% addLightButton
% addLightButton(hf)

% This script creates a new pushtool in a figure's toolbar to manipulate
% the light projected on a displayed axis
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

D=imread(fullfile(iconPath,'lightbulb.jpg'));
S=double(D);
S=S-min(S(:));
S=S./max(S(:));
if size(S,3)==1
    S=repmat(S,[1 1 3]);
end

% Create a uipushtool in the toolbar
uipushtool(hb(1),'TooltipString','light','CData',S,'Tag','light_button','ClickedCallback',{@lightFunc,{hf}});

end

%% Light function  lightFunc

function lightFunc(~,~,inputCell)

hf = inputCell{1};
ha = findobj(hf,'type','axes');

if isempty(ha)
    msgbox('There are no axes in the figure');
else
    prompt = {'Light'};
    title = '';
    format = struct('type','list');
    format.style = 'popupmenu';
    format.items = {' -- Select a lighting option -- ','none','left','right','headlight'};
    default = cell(size(prompt,1),1);
    default{1,1} = 1; % Default selection will always be instruction item
    default = cell2struct(default,prompt,1);
    prompt = repmat(prompt,1,2);
    options.AlignControls = 'on';
    choice = inputsdlg(prompt,title,format,default,options);

    % Unless the user closes the pop up window, clear all existing light
    if choice.Light ~= 1
        existingLight = findobj(ha,'type','light');
        if ~isempty(existingLight)
            delete(existingLight);
        end
        % Apply selected light to all figure axes
        if choice.Light == 3
            for ia = 1:size(ha,1)
                axes(ha(ia)); % Set the current axis
                h(ia) = camlight('left');
            end
        elseif choice.Light == 4
            for ia = 1:size(ha,1)
                axes(ha(ia));
                h(ia) = camlight('right');
            end
        elseif choice.Light == 5
            for ia = 1:size(ha,1)
                axes(ha(ia));
                h(ia) = camlight('headlight');
            end
        end
    end
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