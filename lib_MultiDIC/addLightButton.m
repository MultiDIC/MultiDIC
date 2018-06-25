function []=addLightButton(hf)

%get icon
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
uipushtool(hb(1),'TooltipString','light','CData',S,'Tag','colorbar_button','ClickedCallback',{@LightFunc,{hf}});

end

%% colorbar function  colorbarFunc

function LightFunc(~,~,inputCell)
drawnow

hf=inputCell{1};

ha=findobj(hf,'type','axes');

for ia=1:size(ha,1)
    existingLight=findobj(ha(ia),'type','light');
    if ~isempty(existingLight)
        delete(existingLight);
    end
    hl(ia)=camlight(ha(ia)); lighting flat
    camlight(hl(ia),'headlight');
end

drawnow



% prompt = {'Set quiver scale factor'};
% dlg_title = 'Set quiver scale factor';

% allQuiver=findall(hf,'Type','quiver');
% if ~ isempty(allQuiver)
%     for ic=1:size(allQuiver,1)
%         allQuiver(ic).AutoScale='on';
%     end
%     currentFactor=allQuiver(1).AutoScaleFactor;
%     defaultOptions = {num2str(currentFactor)};
%     s=25+max([cellfun(@numel,prompt) cellfun(@numel,defaultOptions)]);
%     
%     Q = inputdlg(prompt,dlg_title,[1 s],defaultOptions);
%     if ~isempty(Q)
%         newFactor=str2double(Q{1});
%         
%         for ic=1:size(allQuiver,1)
%             allQuiver(ic).AutoScaleFactor=newFactor;
%         end
%     end
%     
% 
% else
%     msgbox('There is no quiver in the figure');
% end

end
