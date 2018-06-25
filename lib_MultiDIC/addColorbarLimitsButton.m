function []=addColorbarLimitsButton(hf)
%% WHAT THIS FUNCTION..


%%
%get icon
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

%% colorbar function  colorbarFunc

function colorbarFunc(~,~,inputCell)

hf=inputCell{1};

prompt = {'Colorbar minimum value:','Colorbar maximum value:'};
dlg_title = 'Set colorbar limits';

allAxes=findall(hf,'Type','axes');
currentLimits=caxis(allAxes(1));
defaultOptions = {num2str(currentLimits(1)),num2str(currentLimits(2))};
s=25+max([cellfun(@numel,prompt) cellfun(@numel,defaultOptions)]);

Q = inputdlg(prompt,dlg_title,[1 s],defaultOptions);
if ~isempty(Q)
    minC=str2double(Q{1});
    maxC=str2double(Q{2});
    
    
    for ic=1:size(allAxes,1)
        caxis(allAxes(ic),[minC maxC]);
    end
end

end
