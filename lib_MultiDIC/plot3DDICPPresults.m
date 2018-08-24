function []=plot3DDICPPresults(varargin)
%% function for plotting 3D-DIC results in STEP3.
% Plotting the 3D reconstruction of points correlated with Ncorr
% The function opens a selection window for all the possible measures to plot
% After selection, the animation figures are plotted
%
% INPUT options:
% DIC3DPPresults


%%
switch nargin
    case 0 % in case no results were entered
        % ask user to load results from 1 or more camera pairs and turn
        % into a cell array
        [file,path] = uigetfile(pwd,'Select a DIC3DPPresults structure');
        result=load([path file]);
        DIC3DPPresults=result.DIC3DPPresults;
        optStruct=struct;
    case 1
        % use given struct
        DIC3DPPresults=varargin{1};
        optStruct=struct;
    case 2
        % use given struct
        DIC3DPPresults=varargin{1};
        optStruct=varargin{2};
    otherwise
        error('wrong number of input arguments');
end

%% select what to plot
Prompt={'\bf{Select which parameters to plot}';... % 1
    'Surfaces with color as pair index';... % 2
    'Points with color as pair index';... % 3
    'Surfaces with color as combined correlation coefficient';... % 4
    'Points with color as correlation coefficient';... % 5
    'Surfaces with color as dispMgn (Displacement magnitude)';... % 6
    'Points with color as dispMgn (Displacement magnitude)';... % 7
    'Surfaces with color as dispX (X Displacement)';... % 8
    'Points with color as dispX (X Displacement)';... % 9
    'Surfaces with color as dispY (Y Displacement)';... % 10
    'Points with color as dispY (Y Displacement)';... % 11
    'Surfaces with color as dispZ (Z Displacement)';... % 12
    'Points with color as dispZ (Z Displacement)';... % 13
    'Surfaces with color as J (surface area change)';... % 14
    'Surfaces with color as lambda1 (1st principal stretch)';... % 15
    'Surfaces with color as Lambda2 (2nd principal stretch)';... % 16
    'Surfaces with color as Epc1 (1st principal Lagrangian strain)';... % 17
    'Surfaces with color as Epc2 (2nd principal Lagrangian strain)';... % 18
    'Surfaces with color as epc1 (1st principal Almansi strain)';... % 19
    'Surfaces with color as epc2 (2nd principal Almansi strain)';... % 20
    'Surfaces with color as Emgn (Lagrangian strain tensor magnitude)';... % 21
    'Surfaces with color as emgn (Almansi strain tensor magnitude)';... % 22
    'Surfaces with color as Lagrangian equivalent strain';... % 23
    'Surfaces with color as Eulerian equivalent strain';... % 24
    'Surfaces with color as Max Lagrangian shear strain';... % 25
    'Surfaces with color as Max Eulerian shear strain';... % 26
    'Surfaces with color as Lamda1+direction (1st principal stretch value and direction)';... % 27
    'Surfaces with color as Lamda2+direction (2nd principal stretch value and direction)';... % 28
    'Surfaces with color as Epc1+direction (1st principal Lagrangian strain value and direction)';... % 29
    'Surfaces with color as Epc2+direction (2nd principal Lagrangian strain value and direction)';... % 30
    'Surfaces with color as epc1+direction (1st principal Almansi strain value and direction)';... % 31
    'Surfaces with color as epc2+direction (2nd principal Almansi strain value and direction)';... % 32
    'Surfaces with color as Epc1+Epc2+direction (1st and 2nd principal Lagrangian strain values and directions)';... % 33
    'Surfaces with color as epc1+epc2+direction (1st and 2nd principal Almansi strain values and directions)'; % 34
    'Surfaces with color as Max Lagrangian shear strain + directions';... % 35
    'Surfaces with color as Max Eulerian shear strain + directions';... % 36
    'Surfaces with color as Lamda1+Lamda2+direction (1st and 2nd principal stretch values and directions)';... % 37
    'Surfaces with color as FaceIsoInd (triangular face isotropy index)';... % 38
    'Select to remove rigid body motion';... % 39
    'Surfaces with color as FaceColors (grayscale from images)';... % 40
    'Select All'; }; % 41

Title='Select which parameters to plot';

Formats=struct;
Formats(1,1).type='text'; %1
Formats(1,2).type='none'; 
Formats(2,1).type='check'; %2
Formats(2,2).type='check'; %3
Formats(3,1).type='check'; %4
Formats(3,2).type='check'; %5
Formats(4,1).type='check'; %6
Formats(4,2).type='check'; %7
Formats(5,1).type='check'; %8
Formats(5,2).type='check'; %9
Formats(6,1).type='check'; %10
Formats(6,2).type='check'; %11
Formats(7,1).type='check'; %12
Formats(7,2).type='check'; %13
Formats(8,1).type='check'; %14
Formats(8,2).type='none'; 
Formats(9,1).type='check'; %15
Formats(9,2).type='check'; %16
Formats(10,1).type='check'; %17
Formats(10,2).type='check'; %18
Formats(11,1).type='check'; %19
Formats(11,2).type='check'; %20
Formats(12,1).type='check'; %21
Formats(12,2).type='check'; %22
Formats(13,1).type='check'; %23
Formats(13,2).type='check'; %24
Formats(14,1).type='check'; %25
Formats(14,2).type='check'; %26
Formats(15,1).type='check'; %27
Formats(15,2).type='check'; %28
Formats(16,1).type='check'; %29
Formats(16,2).type='check'; %30
Formats(17,1).type='check'; %31
Formats(17,2).type='check'; %32
Formats(18,1).type='check'; %33
Formats(18,2).type='check'; %34
Formats(19,1).type='check'; %35
Formats(19,2).type='check'; %36
Formats(20,1).type='check'; %37
Formats(20,2).type='none'; 
Formats(21,1).type='check'; %38
Formats(21,2).type='check'; %39
Formats(22,1).type='check'; %40
Formats(22,2).type='check'; %41
 
DefAns=cell(numel(Prompt),1);
DefAns{1}=[];
for ii=2:numel(Prompt)
    DefAns{ii}=false;
end

Options.Resize='on';
Options.FontSize=10;

[Answer,Canceled] = inputsdlg(Prompt, Title, Formats, DefAns, Options);
if Answer{41}
    for ii=[2:38 40]
        Answer{ii}=true;
    end
end
if Answer{39}
    RBMlogic=true;
else
    RBMlogic=false;
end
if Canceled
    return
end
%% Select plotting options
if ~Canceled && sum(cell2mat(Answer(2:end)))>0
    answer = inputdlg({'Enter maximum correlation coefficient to keep points (leave blank for keeping all points)'},'Input',[1,50]);
    CorCoeffCutOff=str2double(answer{1}); % maximal correlation coefficient for display (use [] for default which is max)
    if isnan(CorCoeffCutOff)
        CorCoeffCutOff=[];
    end
else
    return
end

%% create option struct for plotting
% complete the struct fields
if ~isfield(optStruct,'zDirection')
    optStruct.zDirection=1;
end
optStruct.maxCorrCoeff=CorCoeffCutOff;

% optStruct.lineColor='k';
% optStruct.smoothLogic=0;
% optStruct.colorMap=cMap;
% optStruct.supTitleString=pointMeasureString;

%% plot according to answer

if Canceled
    return
end

if Answer{2}
    anim8_DIC3DPP_faceMeasure(DIC3DPPresults,'FacePairInds',RBMlogic,optStruct);
end
if Answer{3}
    anim8_DIC3DPP_pointMeasure(DIC3DPPresults,'pairInd',RBMlogic,optStruct);
end
if Answer{4}
    anim8_DIC3DPP_faceMeasure(DIC3DPPresults,'FaceCorrComb',RBMlogic,optStruct);
end
if Answer{5}
    anim8_DIC3DPP_pointMeasure(DIC3DPPresults,'corrComb',RBMlogic,optStruct);
end
if Answer{6}
    anim8_DIC3DPP_faceMeasure(DIC3DPPresults,'DispMgn',RBMlogic,optStruct);
end
if Answer{7}
    anim8_DIC3DPP_pointMeasure(DIC3DPPresults,'DispMgn',RBMlogic,optStruct);
end
if Answer{8}
    anim8_DIC3DPP_faceMeasure(DIC3DPPresults,'DispX',RBMlogic,optStruct);
end
if Answer{9}
    anim8_DIC3DPP_pointMeasure(DIC3DPPresults,'DispX',RBMlogic,optStruct);
end
if Answer{10}
    anim8_DIC3DPP_faceMeasure(DIC3DPPresults,'DispY',RBMlogic,optStruct);
end
if Answer{11}
    anim8_DIC3DPP_pointMeasure(DIC3DPPresults,'DispY',RBMlogic,optStruct);
end
if Answer{12}
    anim8_DIC3DPP_faceMeasure(DIC3DPPresults,'DispZ',RBMlogic,optStruct);
end
if Answer{13}
    anim8_DIC3DPP_pointMeasure(DIC3DPPresults,'DispZ',RBMlogic,optStruct);
end
if Answer{14}
    anim8_DIC3DPP_faceMeasure(DIC3DPPresults,'J',RBMlogic,optStruct);
end
if Answer{15}
    anim8_DIC3DPP_faceMeasure(DIC3DPPresults,'Lamda1',RBMlogic,optStruct);
end
if Answer{16}
    anim8_DIC3DPP_faceMeasure(DIC3DPPresults,'Lamda2',RBMlogic,optStruct);
end
if Answer{17}
    anim8_DIC3DPP_faceMeasure(DIC3DPPresults,'Epc1',RBMlogic,optStruct);
end
if Answer{18}
    anim8_DIC3DPP_faceMeasure(DIC3DPPresults,'Epc2',RBMlogic,optStruct);
end
if Answer{19}
    anim8_DIC3DPP_faceMeasure(DIC3DPPresults,'epc1',RBMlogic,optStruct);
end
if Answer{20}
    anim8_DIC3DPP_faceMeasure(DIC3DPPresults,'epc2',RBMlogic,optStruct);
end
if Answer{21}
    anim8_DIC3DPP_faceMeasure(DIC3DPPresults,'Emgn',RBMlogic,optStruct);
end
if Answer{22}
    anim8_DIC3DPP_faceMeasure(DIC3DPPresults,'emgn',RBMlogic,optStruct);
end
if Answer{23}
    anim8_DIC3DPP_faceMeasure(DIC3DPPresults,'Eeq',RBMlogic,optStruct);
end
if Answer{24}
    anim8_DIC3DPP_faceMeasure(DIC3DPPresults,'eeq',RBMlogic,optStruct);
end
if Answer{25}
    anim8_DIC3DPP_faceMeasure(DIC3DPPresults,'EShearMax',RBMlogic,optStruct);
end
if Answer{26}
    anim8_DIC3DPP_faceMeasure(DIC3DPPresults,'eShearMax',RBMlogic,optStruct);
end
if Answer{27}
    anim8_DIC3DPP_faceMeasureDirection(DIC3DPPresults,'Lamda1',RBMlogic,optStruct);
end
if Answer{28}
    anim8_DIC3DPP_faceMeasureDirection(DIC3DPPresults,'Lamda2',RBMlogic,optStruct);
end
if Answer{29}
    anim8_DIC3DPP_faceMeasureDirection(DIC3DPPresults,'Epc1',RBMlogic,optStruct);
end
if Answer{30}
    anim8_DIC3DPP_faceMeasureDirection(DIC3DPPresults,'Epc2',RBMlogic,optStruct);
end
if Answer{31}
    anim8_DIC3DPP_faceMeasureDirection(DIC3DPPresults,'epc1',RBMlogic,optStruct);
end
if Answer{32}
    anim8_DIC3DPP_faceMeasureDirection(DIC3DPPresults,'epc2',RBMlogic,optStruct);
end
if Answer{33}
    anim8_DIC3DPP_faceMeasureDirection(DIC3DPPresults,{'Epc1','Epc2'},RBMlogic,optStruct);
end
if Answer{34}
    anim8_DIC3DPP_faceMeasureDirection(DIC3DPPresults,{'epc1','epc2'},RBMlogic,optStruct);
end
if Answer{35}
    % develop shear direction plot
    msgbox('Shear directions cannot be plotted. It is under construction...');
end
if Answer{36}
    % develop shear direction plot
    msgbox('Shear directions cannot be plotted. It is under construction...');
end
if Answer{37}
    anim8_DIC3DPP_faceMeasureDirection(DIC3DPPresults,{'Lamda1','Lamda2'},RBMlogic,optStruct);
end
if Answer{38}
    anim8_DIC3DPP_faceMeasure(DIC3DPPresults,'FaceIsoInd',RBMlogic,optStruct);
end
if Answer{40}
    anim8_DIC3DPP_faceMeasure(DIC3DPPresults,'FaceColors',RBMlogic,optStruct);
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