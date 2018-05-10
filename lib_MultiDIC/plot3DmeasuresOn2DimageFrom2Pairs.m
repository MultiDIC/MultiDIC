function [] = plot3DmeasuresOn2DimageFrom2Pairs(varargin)
%% fucntion for plotting 3D measures on the original 2D images
%plot3DmeasuresOn2Dimages
%plot3DmeasuresOn2Dimages(DIC2DpairResults,DIC3DpairResults)
%%
switch nargin
    case 0
        PathInitial=pwd;
%         select a 2D pair result
        [fileName1, filePath1] = uigetfile(PathInitial,'Select a 2D pair Left result structure');
        [fileName2, filePath2] = uigetfile(filePath1,'Select a 2D pair Right result structure');
%         select the matching 3D pair result
        [fileName3, filePath3] = uigetfile(filePath1,'Select the 3D pair Left result structure');
        [fileName4, filePath4] = uigetfile(filePath1,'Select the 3D pair Right result structure');
                
        DIC2DpairResultsL=load([filePath1 fileName1]);
        DIC2DpairResultsL=DIC2DpairResultsL.DIC2DpairResults;
        DIC2DpairResultsR=load([filePath2 fileName2]);
        DIC2DpairResultsR=DIC2DpairResultsR.DIC2DpairResults;
        
        DIC3DpairResultsL=load([filePath3 fileName3]);
        DIC3DpairResultsL=DIC3DpairResultsL.DIC3DpairResults;
        DIC3DpairResultsR=load([filePath4 fileName4]);
        DIC3DpairResultsR=DIC3DpairResultsR.DIC3DpairResults;
        
        DIC2DpairResults={DIC2DpairResultsL,DIC2DpairResultsR};
        DIC3DpairResults={DIC3DpairResultsL,DIC3DpairResultsR};

    case 2
        DIC2DpairResults=varargin{1};
        DIC3DpairResults=varargin{2};

        
    otherwise
        error('wrong number of input arguments');
end

%% select what to plot
Prompt={'\bf{Select which parameters to plot}';... % 1
    'Points with color as DispMgn (Displacement magnitude)';... % 2
    'Points with color as DispX (X Displacement)';... % 3
    'Points with color as DispY (Y Displacement)';... % 4
    'Points with color as DispZ (Z Displacement)';... % 5
    'Surfaces with color as J (surface area change)';... % 6
    'Surfaces with color as Lambda1 (1st principal stretch)';... % 7
    'Surfaces with color as Lamda2 (2nd principal stretch)';... % 8
    'Surfaces with color as Epc1 (1st principal Lagrangian strain)';... % 9
    'Surfaces with color as Epc2 (2nd principal Lagrangian strain)';... % 10
    'Surfaces with color as epc1 (1st principal Almansi strain)';... % 11
    'Surfaces with color as epc2 (2nd principal Almansi strain)';... % 12
    'Surfaces with color as Emgn (Lagrangian strain tensor magnitude)';... % 13
    'Surfaces with color as emgn (Almansi strain tensor magnitude)';... % 14
    'Select All';}; % 15

Title='Select which parameters to plot';

Formats=struct;
Formats(1,1).type='text';
Formats(2,1).type='check';
Formats(3,1).type='check';
Formats(4,1).type='check';
Formats(5,1).type='check';
Formats(6,1).type='check';
Formats(7,1).type='check';
Formats(8,1).type='check';
Formats(9,1).type='check';
Formats(10,1).type='check';
Formats(11,1).type='check';
Formats(12,1).type='check';
Formats(13,1).type='check';
Formats(14,1).type='check';
Formats(15,1).type='check';

DefAns=cell(numel(Prompt),1);
DefAns{1}=[];
for ii=2:numel(Prompt)
    DefAns{ii}=false;
end

Options.Resize='on';
Options.FontSize=10;

[Answer,Canceled] = inputsdlg(Prompt, Title, Formats, DefAns, Options);
if Answer{15}
    for ii=2:14
        Answer{ii}=true;
    end
end
if Canceled
    return
end

%% Select plotting options
optStruct=struct;
if ~Canceled && sum(cell2mat(Answer(2:end)))>0
    answer = inputdlg({'Enter maximum correlation coefficient to keep points (leave blank for keeping all points)'},'Input',[1,50]);
    optStruct.CorCoeffCutOff=str2double(answer{1}); % maximal correlation coefficient for display (use [] for default which is max)
else
    return
end

%% load images
nImages=DIC2DpairResultsL.nImages;
ImPaths=DIC2DpairResultsL.ImPaths;
ImSet=cell(nImages,1);
for ii=1:nImages
    ImSet{ii}=imread(ImPaths{ii+nImages});
end

%% plot animated images with the 3D measures on them
if Answer{2}
    anim8_DIC_image_3Dmeasure_points_2n(ImSet,DIC2DpairResults,DIC3DpairResults,'DispMgn',optStruct);
end
if Answer{3}
    anim8_DIC_image_3Dmeasure_points_2n(ImSet,DIC2DpairResults,DIC3DpairResults,'DispX',optStruct);
end
if Answer{4}
    anim8_DIC_image_3Dmeasure_points_2n(ImSet,DIC2DpairResults,DIC3DpairResults,'DispY',optStruct);
end
if Answer{5}
    anim8_DIC_image_3Dmeasure_points_2n(ImSet,DIC2DpairResults,DIC3DpairResults,'DispZ',optStruct);
end
if Answer{6}
    anim8_DIC_image_3Dmeasure_faces_2n(ImSet,DIC2DpairResults,DIC3DpairResults,'J',optStruct);
end
if Answer{7}
    anim8_DIC_image_3Dmeasure_faces_2n(ImSet,DIC2DpairResults,DIC3DpairResults,'Lamda1',optStruct);
end
if Answer{8}
    anim8_DIC_image_3Dmeasure_faces_2n(ImSet,DIC2DpairResults,DIC3DpairResults,'Lamda2',optStruct);
end
if Answer{9}
    anim8_DIC_image_3Dmeasure_faces_2n(ImSet,DIC2DpairResults,DIC3DpairResults,'Epc1',optStruct);
end
if Answer{10}
    anim8_DIC_image_3Dmeasure_faces_2n(ImSet,DIC2DpairResults,DIC3DpairResults,'Epc2',optStruct);
end
if Answer{11}
    anim8_DIC_image_3Dmeasure_faces_2n(ImSet,DIC2DpairResults,DIC3DpairResults,'epc1',optStruct);
end
if Answer{12}
    anim8_DIC_image_3Dmeasure_faces_2n(ImSet,DIC2DpairResults,DIC3DpairResults,'epc2',optStruct);
end
if Answer{13}
    anim8_DIC_image_3Dmeasure_faces_2n(ImSet,DIC2DpairResults,DIC3DpairResults,'Emgn',optStruct);
end
if Answer{14}
    anim8_DIC_image_3Dmeasure_faces_2n(ImSet,DIC2DpairResults,DIC3DpairResults,'emgn',optStruct);
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