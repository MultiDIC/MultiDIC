function [Points,CorCoeffVec,F,CF] = extractNcorrResults(handles_ncorr,IMref)
%% function for Extracting the results from Ncorr and calculate correlated image points, correlation coefficients, faces and face colors, in step 2
%[Points,CorCoeffVec,F,CF] = extractNcorrResults(handles_ncorr,IM);
%
% INPUT:
% * handles_ncorr - the Ncorr handle, after running the analysis (up until "Displacement", no need for "Strains"). 
% * IM - the reference image, for obtaining the grayscale values for the
% faces
%
% OUTPUT:
% * Points - a 2nX1 cell containing 2D image points of the correlated
% points between the 1st image (reference) and all the other images. Points
% which could not be correlated get NaN instead of a value, but the number
% of points is equal for all images, to preserve correspondance.
% * CorCoeffVec - a 2nX1 cell containing vectors the same length as Points,
% with the correlation coefficient for each point (see Ncorr documentation
% for understanding the meaning of the correlation coefficient).
% * F - triangular faces defined on the reference image (regular grid)
% based on Delaunay
% * CF - grayscale values for each face in F, based on the reference image

%%
DICresults = handles_ncorr.data_dic;

% Extract results
Disp=DICresults.displacements;
DispInfo=DICresults.dispinfo;
nCur=size(Disp,2);

Factor=DispInfo.spacing+1;
ROI_DIC=cell(nCur,1);
CorCoeff=cell(nCur,1); CorCoeffVec=cell(nCur,1); Points=cell(nCur,1);
Uvec=cell(nCur,1); Vvec=cell(nCur,1);

for ii=1:nCur
    
    ROI_DIC{ii}=Disp(ii).roi_dic.mask;
    
    CorCoeff{ii}=Disp(ii).plot_corrcoef_dic;
    Uref=Disp(ii).plot_u_ref_formatted;
    Vref=Disp(ii).plot_v_ref_formatted;
    
    [YrefROIVec,XrefROIVec] = find(ROI_DIC{1});
    [YcurROIVec,XcurROIVec] = find(ROI_DIC{ii});
    
    PtempRef=[XrefROIVec,YrefROIVec];
    PtempRef=(PtempRef-1)*Factor+1; % switch from sapcing to pixels
    Pref=PtempRef;
    
    CorCoeffVec{ii}=CorCoeff{ii}(ROI_DIC{1});
    CorCoeffVec{ii}(CorCoeffVec{ii}==0)=NaN;
    
    % displacements from ref to cur
    UrefROIVec=Uref(ROI_DIC{1});
    UrefROIVec(UrefROIVec==0)=NaN;
    VrefROIVec=Vref(ROI_DIC{1});
    VrefROIVec(VrefROIVec==0)=NaN;
    
    Uvec{ii}=UrefROIVec;
    Vvec{ii}=UrefROIVec;
    % current points
    Points{ii}=[Pref(:,1)+UrefROIVec,Pref(:,2)+VrefROIVec];
    
    % save face colors for further 3D analysis
    if ii==1
        % pixel colors
        IMrefSmall=IMref(1:Factor:end,1:Factor:end);
        IMrefSmallMasked=IMrefSmall;
        IMrefSmallMasked(~ROI_DIC{1})=[];
        ColorRef=IMrefSmallMasked(:);
    end
    
end

% Create triangulation
DT = delaunayTriangulation(Pref);
F = DT.ConnectivityList;
V = Pref;

% remove irregular triangles
EdgeLengths = patchEdgeLengths(F,V);
EdgeLengths = [EdgeLengths(1:3:length(EdgeLengths)) EdgeLengths(2:3:length(EdgeLengths)) EdgeLengths(3:3:length(EdgeLengths))];
EdgeLengthsMax=max(EdgeLengths,[],2);

F(EdgeLengthsMax>1.1*sqrt(2)*Factor,:)=[];

% flip direction to have normals pointing out
F=F(:,[1 3 2]);

% face colors (average node colors)
CF=mean(ColorRef(F),2);


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