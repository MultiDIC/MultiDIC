function exportGif(optStruct)

%% Parse input

%Check for files names
if ~isfield(optStruct,'FileNames')    
    error('Please provide FileNames in the input structure');    
end

%Force optStruct.FileNames to be cell
if ~isa(optStruct.FileNames,'cell')
    optStruct.FileNames={optStruct.FileNames};
end

%Check if input is a single folder or list of files
if numel(optStruct.FileNames)==1    
    checkName=optStruct.FileNames{1};    
    if exist(checkName,'dir')==7 
        %It is a folder so attempt to convert all files in folder to a gif.
        %File names are sorted. 
        FileNames=dir(fullfile(checkName));
        FileNames={FileNames(1:end).name};
        FileNames=sort(FileNames(:));
        c=1;
        for q=1:1:numel(FileNames)
            fullfile(checkName,FileNames{q})
            if exist(fullfile(checkName,FileNames{q}),'file')==2
                optStruct.FileNames{c}=fullfile(checkName,FileNames{q});
                c=c+1;
            end
        end
    end
end

%Use defaults if other inputs are missing    
if ~isfield(optStruct,'DelayTime')
    optStruct.DelayTime=0.5; 
end

if ~isfield(optStruct,'LoopCount')
    optStruct.LoopCount=inf; 
end

if ~isfield(optStruct,'FileNameGif')
    optStruct.FileNameGif='exportGif'; 
end

%Check if save location exists
[Savepath,~,fileExt] = fileparts(optStruct.FileNameGif);
if exist(Savepath,'dir')~=7 %create output folder if it does not exist already
    mkdir(Savepath);
end

%Check if extension is given
if isempty(fileExt)    
    optStruct.FileNameGif=[optStruct.FileNameGif,'.gif'];
% elseif ~strcmp(fileExt,'.gif')    
end

%%

numFiles=numel(optStruct.FileNames);
hw = waitbar(0,'Exporting .gif animation');
for q=1:1:numFiles   
    D = imread(optStruct.FileNames{q});    
    [A,map] = rgb2ind(D,256);    
    if q == 1
        imwrite(A,map,optStruct.FileNameGif,'gif','LoopCount',optStruct.LoopCount,'DelayTime',optStruct.DelayTime);
    else
        imwrite(A,map,optStruct.FileNameGif,'gif','WriteMode','append','DelayTime',optStruct.DelayTime);
    end    
    waitbar(q/numFiles,hw,['Exporting .gif animation ',num2str(round(100*q/numFiles)),'%']);
end
close(hw);

 
%% 
% _*GIBBON footer text*_ 
% 
% License: <https://github.com/gibbonCode/GIBBON/blob/master/LICENSE>
% 
% GIBBON: The Geometry and Image-based Bioengineering add-On. A toolbox for
% image segmentation, image-based modeling, meshing, and finite element
% analysis.
% 
% Copyright (C) 2018  Kevin Mattheus Moerman
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
