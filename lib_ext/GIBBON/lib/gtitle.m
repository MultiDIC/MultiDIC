function [varargout]=gtitle(varargin)

% function hText=gtitle(titleString,fontSize,hf)

%% Parse input

switch nargin
    case 1
        titleString=varargin{1};
        fontSize=[];
        hf=[];        
    case 2
        titleString=varargin{1};
        fontSize=varargin{2};
        hf=[];
    case 3
        titleString=varargin{1};
        fontSize=varargin{2};
        hf=varargin{3};
end

%%

if isempty(fontSize)
    fontSize=11;
end

if isempty(hf)
    hf=gcf;
end

backGroundColor=hf.Color;

if strcmp(backGroundColor,'k')
    textColor='w';
elseif size(backGroundColor,2)==3
    grayLevels=linspace(0,1,100);
    [~,indMax]=max(abs(grayLevels-mean(double(backGroundColor))));    
    textColor=grayLevels(indMax)*ones(1,3);
else
    textColor='k';    
end

hText = uicontrol(hf,'Style','text','String',titleString,'BackgroundColor',backGroundColor,'HorizontalAlignment','Center','FontSize',fontSize,'FontWeight','bold','Units','Points','ForegroundColor',textColor);

% Set the control to be non-opaque and repaint it
% j_hText = findjobj(hText);
% j_hText.setOpaque(false);
% j_hText.repaint();

hText.Units = 'Points';

figResize([],[],{hf,hText});


hFunc=get(hf,'ResizeFcn');

if iscell(hFunc)
    warning('gtitle replaced the ResizeFcn function. Specify your ResizeFcn in the form @(h,e)figResize(h,e,c) to avoid this behavior');    
    set(hf,'ResizeFcn',@(a,b)figResize(a,b,{hf,hText}));
else
    if isempty(hFunc)
        set(hf,'ResizeFcn',@(a,b)figResize(a,b,{hf,hText}));
    else        
        set(hf,'ResizeFcn',@(a,b)(cellfun(@(x)feval(x,a,b),{hFunc,@(a,b)figResize(a,b,{hf,hText})})));
    end
end


%% Collect output

if nargout>0
    varargout{1}=hText;
end

end

function figResize(~,~,inputCell)
hf=inputCell{1};
hTextInfo=inputCell{2};
unitsNow=hf.Units;
hf.Units='Points';

figPosition=hf.Position;
textBoxWidth=figPosition(3);
textBoxHeight=hTextInfo.Extent(4).*ceil(hTextInfo.Extent(3)./figPosition(3));
textPosition=[0 figPosition(4)-textBoxHeight textBoxWidth textBoxHeight];
hTextInfo.Position = textPosition;

hf.Units=unitsNow;
end

