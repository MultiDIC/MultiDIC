function efw(varargin)

% function efw(hf)
% ------------------------------------------------------------------------
% The efw function, the export figure widget, adds a push button to a
% figure toolbar to link with the export_fig function. Press the button
% to start exporting a figure. Users can specify file names, formats,
% resolution and also additional export_fig options. 
%
% Activate the widget (add button) using efw; e.g.:
% figure; surf(peaks(25)); axis equal; axis tight; efw;
%
% Pressing the efw button will open a basic inputdlg allowing users to
% specify all the usual export_fig options. Hints are given in brackets
% behind the input labels. The default entries are altered according to the
% previous usage within the current figure. 
% 
% The Export Figure Widget requires the external code export_fig created by
% Oliver Woodford and Yair Altman. It can be obtained from the Mathworks
% Central File Exchange: 
% <http://www.mathworks.com/matlabcentral/fileexchange/23629-export-fig>
%
% Kevin Mattheus Moerman
% gibbon.toolbox@gmail.com
% 2015/04/29 Created
% 2015/04/29 Updated to support export of multiple image formats e.g. pdf,
% eps
% 2016/12/09 Changed to rely on figure UserData rather than toolbar,
% simplified calllback function use.
%------------------------------------------------------------------------

%% Parse input arguments
switch nargin
    case 0
        hf = gcf;
    case 1
        hf=varargin{1};
    otherwise
        error('Wrong number of input arguments');
end

if ~ishandle(hf)
    hf = gcf;
end

%% Initialise button
hb = findall(hf,'Type','uitoolbar');

%Check for presence of a efw button
hp = findobj(hb,'Tag','efw_button');
if isempty(hp) %If efw button is not present create one 
    
    % Build icon
    s=[NaN,NaN,0.02,0.64,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,0.64,0.01,NaN,NaN;...
       NaN,NaN,0.4,0.7,0.15,0.16,0.15,0.15,0.16,0.16,0.16,0.15,0.71,0.38,NaN,NaN;...
       NaN,NaN,0.45,0.56,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0.58,0.43,NaN,NaN;...
       NaN,NaN,0.45,0.56,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0.58,0.42,NaN,NaN;...
       NaN,NaN,0.45,0.58,0.03,0.05,0.05,0.05,0.05,0.05,0.05,0.03,0.6,0.44,NaN,NaN;...
       0.33,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,0.33;...
       0.78,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,0.48,0.16,1.0,0.78;...
       1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,0.38,0.01,0.79,1.0;...
       1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0;...
       1.0,1.0,0.78,0.62,0.16,0.17,0.17,0.17,0.17,0.17,0.17,0.16,0.63,0.79,1.0,1.0;...
       0.79,1.0,0.78,0.56,NaN,NaN,1,1,1,1,NaN,NaN,0.58,0.78,1.0,0.79;...
       0.51,1.0,0.78,0.56,NaN,NaN,1,NaN,NaN,1,NaN,NaN,0.58,0.79,1.0,0.51;...
       NaN,0.18,0.54,0.56,NaN,NaN,1,1,1,1,NaN,NaN,0.58,0.53,0.18,NaN;...
       NaN,NaN,0.45,0.56,NaN,NaN,1,NaN,NaN,NaN,NaN,NaN,0.58,0.43,NaN,NaN;...
       NaN,NaN,0.4,0.7,0.15,0.15,1,1,1,1,0.16,0.15,0.71,0.38,NaN,NaN;...
       NaN,NaN,0.02,0.64,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,0.63,0.02,NaN,NaN];
    
    S=zeros(16,16,3);
    S(:,:,1)=0.7.*s;
    S(:,:,2)=0.25.*s;
    S(:,:,3)=0.05.*s;
    
    % Create a uipushtool in the toolbar
    hp=uipushtool(hb,'TooltipString','Export Figure Widget','CData',S,'Tag','efw_button');
    
    if ~isfield(hf.UserData,'efw')
        hf.UserData.efw.defaultPath=fullfile(cd,'efw');
        hf.UserData.efw.imName=['figure',num2str(get(hf,'Number'))];
        hf.UserData.efw.imExt='png';
        hf.UserData.efw.imRes='100';
        hf.UserData.efw.exportFigOpt='-transparent';
    else        
        if ~isfield(hf.UserData.efw('defaultPath'))
            hf.UserData.efw.defaultPath=fullfile(cd,'efw');
        end
        
        if ~isfield(hf.UserData.efw('imName'))
            hf.UserData.efw.imName=['figure',num2str(get(hf,'Number'))];
        end
        
        if ~isfield(hf.UserData.efw('imExt'))
            hf.UserData.efw.imExt='png';
        end
        
        if ~isfield(hf.UserData.efw('imRes'))
            hf.UserData.efw.imRes='100';
        end
        
        if ~isfield(hf.UserData.efw('exportFigOpt'))
            hf.UserData.efw.exportFigOpt='-transparent';
        end
    end
    set(hp,'ClickedCallback',{@start_efw,{hf}});      
end

end

function start_efw(~,~,inputCell)
hf=inputCell{1};
figure(hf);
defStruct=hf.UserData.efw; 
prompt = {'Save path (leave empty to browse to desired folder instead):','Image name:','Image extension (i.e. png,jpg, pdf, eps, bmp, tif, fig, all):','Image resolution (e.g. 120. Ignored for fig):','Extra export_fig options (comma seperated, no spaces e.g. -nocrop,-transparent,-painters):'};
dlg_title = 'Export Figure Widget (see: help efw and help export_fig)';
defaultOptions = {defStruct.defaultPath,defStruct.imName,defStruct.imExt,defStruct.imRes,defStruct.exportFigOpt};

s=25+max([cellfun(@numel,prompt) cellfun(@numel,defaultOptions)]);

Q = inputdlg(prompt,dlg_title,[1 s],defaultOptions);

if ~isempty(Q)
    if isempty(Q{1})
        Q{1}=uigetdir(defStruct.defaultPath,'Select save path');
        if Q{1}==0
            return; 
        end
    end    
    
    if ~exist(Q{1},'dir') %create output folder if it does not exist already
        mkdir(Q{1});
    end
    
    if all(~cellfun(@isempty,Q(1:end-1)))        
        
        fileName=fullfile(Q{1},Q{2});        
        inputCell{1,1}=fullfile(Q{1},Q{2});
                
        if strcmp(Q{3},'all')
            formatAll={'-png','-jpg','-tiff','-bmp','-eps','-pdf'};
            for q=1:1:numel(formatAll)
                inputCell{1,end+1}=formatAll{q};
            end
            savefig(hf,fileName); %Save figure in .fig file
        elseif strcmp(Q{3},'fig') %Just figure
            savefig(hf,fileName); %Save figure in .fig file
        else
            stringSet=Q{3}; %The set of potentially multiple image formats
            stringSetSep = strsplit(stringSet,','); %Split into seperate cell components using commas
            for q=1:1:numel(stringSetSep) 
                stringNoSpaces=regexprep(stringSetSep{q},'[^\w'']',''); %Remove potential extra spaces
                if ~strcmp(stringNoSpaces(1),'-') %If first character is not '-'                
                    stringNoSpaces=['-',stringNoSpaces]; %Add '-' to start, e.g. 'jpg' becomes '-jpg'
                end
                inputCell{1,end+1}=stringNoSpaces; %Add to input list
            end                                    
        end
        
        figRes=['-r',Q{4}];
        inputCell{1,end+1}=figRes;
        
        if ~isempty(Q{5})
            stringSet=Q{5}; %The set of potentially multiple options
            stringSetSep = strsplit(stringSet,',');
            for q=1:1:numel(stringSetSep)
                stringNoSpaces=regexprep(stringSetSep{q},'[^\w'']',''); %Remove potential extra spaces
                if ~strcmp(stringNoSpaces(1),'-') %If first character is not '-'
                    stringNoSpaces=['-',stringNoSpaces]; %Add '-' to start, e.g. 'jpg' becomes '-jpg'
                end
                inputCell{1,end+1}=stringNoSpaces; %Add to input list
            end
        end
        
        logicFig=strcmp(inputCell,'-fig');
        if any(logicFig)
            savefig(hf,fileName); %Save figure in .fig file
            inputCell=inputCell(~logicFig);
        end
        
        export_fig(inputCell{:});
        
        %Override defaults
        defStruct.defaultPath=Q{1};
        defStruct.imName=Q{2};
        defStruct.imExt=Q{3};
        defStruct.imRes=Q{4};
        defStruct.exportFigOpt=Q{5};
        
        hf.UserData.efw=defStruct;
        
    else
        return
    end
end

end
 
%% 
% _*GIBBON footer text*_ 
% 
% License: <https://github.com/gibbonCode/GIBBON/blob/master/LICENSE>
% 
% GIBBON: The Geometry and Image-based Bioengineering add-On. A toolbox for
% image segmentation, image-based modeling, meshing, and finite element
% analysis.
% 
% Copyright (C) 2017  Kevin Mattheus Moerman
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
