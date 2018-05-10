function DIC_3Djoined_results = joinPairs(DIC3DAllPairsResults)
%% function for joining multiple 3D-DIC pair results into one structure (plotting in step3)
%
% Inputs:
% * DIC3DAllPairsResults
% 
% Outputs:
% * DIC_3Djoined_results
% 

%%
nPairs=numel(DIC3DAllPairsResults);
nFrames=size(DIC3DAllPairsResults{1}.Points3D,1);

% create cells with all the variables for joining
Fi=cell(nPairs,1);
FCi=cell(nPairs,1);
Fisoi=cell(nPairs,1);
Vi=cell(nPairs,nFrames);
VTi=cell(nPairs,nFrames);
corrCombi=cell(nPairs,nFrames);
FaceCorrCombi=cell(nPairs,nFrames);
DispVi=cell(nPairs,nFrames);
DispVTi=cell(nPairs,nFrames);
DispMi=cell(nPairs,nFrames);
DispMTi=cell(nPairs,nFrames);
VFci=cell(nPairs,nFrames);
VFcTi=cell(nPairs,nFrames);
Ji=cell(nPairs,nFrames);
Lamda1i=cell(nPairs,nFrames);
Lamda2i=cell(nPairs,nFrames);
Emgni=cell(nPairs,nFrames);
emgni=cell(nPairs,nFrames);
Epc1i=cell(nPairs,nFrames);
Epc2i=cell(nPairs,nFrames);
Epc1veci=cell(nPairs,nFrames);
Epc1vecCuri=cell(nPairs,nFrames);
Epc2veci=cell(nPairs,nFrames);
Epc2vecCuri=cell(nPairs,nFrames);
epc1i=cell(nPairs,nFrames);
epc2i=cell(nPairs,nFrames);
epc1veci=cell(nPairs,nFrames);
epc2veci=cell(nPairs,nFrames);


for ip=1:nPairs
    DIC_3Dpair_struct=DIC3DAllPairsResults{ip};
    
    Fi{ip}=DIC_3Dpair_struct.Faces;
    FCi{ip}=DIC_3Dpair_struct.FaceColors;
      
    for it=1:nFrames
        Vi{ip,it}=DIC_3Dpair_struct.Points3D{it};
        VTi{ip,it}=DIC_3Dpair_struct.Points3Dtransformed{it};
        corrCombi{ip,it}=DIC_3Dpair_struct.corrComb{it};
        FaceCorrCombi{ip,it}=DIC_3Dpair_struct.FaceCorrComb{it};
        DispVi{ip,it}=DIC_3Dpair_struct.Disp.DispVec{it};
        DispVTi{ip,it}=DIC_3Dpair_struct.Disp.DispVecTransformed{it};
        DispMi{ip,it}=DIC_3Dpair_struct.Disp.DispMgn{it};
        DispMTi{ip,it}=DIC_3Dpair_struct.Disp.DispMgnTransformed{it};
        VFci{ip,it}=DIC_3Dpair_struct.FaceCentroids{it};
        VFcTi{ip,it}=DIC_3Dpair_struct.FaceCentroidsTransformed{it};
        Ji{ip,it}=DIC_3Dpair_struct.Deform.J{it};
        Lamda1i{ip,it}=DIC_3Dpair_struct.Deform.Lamda1{it};
        Lamda2i{ip,it}=DIC_3Dpair_struct.Deform.Lamda2{it};
        Emgni{ip,it}=DIC_3Dpair_struct.Deform.Emgn{it};
        emgni{ip,it}=DIC_3Dpair_struct.Deform.emgn{it};
        Epc1i{ip,it}=DIC_3Dpair_struct.Deform.Epc1{it};
        Epc2i{ip,it}=DIC_3Dpair_struct.Deform.Epc2{it};
        Epc1veci{ip,it}=DIC_3Dpair_struct.Deform.Epc1vec{it};
        Epc1vecCuri{ip,it}=DIC_3Dpair_struct.Deform.Epc1vecCur{it};
        Epc2veci{ip,it}=DIC_3Dpair_struct.Deform.Epc2vec{it};
        Epc2vecCuri{ip,it}=DIC_3Dpair_struct.Deform.Epc2vecCur{it};
        epc1i{ip,it}=DIC_3Dpair_struct.Deform.epc1{it};
        epc2i{ip,it}=DIC_3Dpair_struct.Deform.epc2{it};
        epc1veci{ip,it}=DIC_3Dpair_struct.Deform.epc1vec{it};
        epc2veci{ip,it}=DIC_3Dpair_struct.Deform.epc2vec{it};
        Fisoi{ip,it}=DIC_3Dpair_struct.FaceIsoInd{it};
    end

end

% join variables that are time invariant
[Fj,~,~]=joinElementSets(Fi,Vi(:,1),[]);
[~,~,FCj]=joinElementSets(Fi,Vi(:,1),FCi);


for it=1:nFrames
    % join variables that are time variant and scalar
    
    [~,Vj{1,it},DispMj{1,it}]=joinElementSets(Fi,Vi(:,it),DispMi(:,it));
    [~,VTj{1,it},DispMTj]=joinElementSets(Fi,VTi(:,it),DispMTi(:,it));
    [~,~,corrCombj{1,it}]=joinElementSets(Fi,Vi(:,it),corrCombi(:,it));
    [~,~,FaceCorrCombj{1,it}]=joinElementSets(Fi,Vi(:,it),FaceCorrCombi(:,it));
    [~,~,Jj{1,it}]=joinElementSets(Fi,Vi(:,it),Ji(:,it));
    [~,~,Lamda1j{1,it}]=joinElementSets(Fi,Vi(:,it),Lamda1i(:,it));
    [~,~,Lamda2j{1,it}]=joinElementSets(Fi,Vi(:,it),Lamda2i(:,it));
    [~,~,Emgnj{1,it}]=joinElementSets(Fi,Vi(:,it),Emgni(:,it));
    [~,~,emgnj{1,it}]=joinElementSets(Fi,Vi(:,it),emgni(:,it));
    [~,~,Epc1j{1,it}]=joinElementSets(Fi,Vi(:,it),Epc1i(:,it));
    [~,~,Epc2j{1,it}]=joinElementSets(Fi,Vi(:,it),Epc2i(:,it));
    [~,~,epc1j{1,it}]=joinElementSets(Fi,Vi(:,it),epc1i(:,it));
    [~,~,epc2j{1,it}]=joinElementSets(Fi,Vi(:,it),epc2i(:,it));
    [~,~,Fisoj{1,it}]=joinElementSets(Fi,Vi(:,it),Fisoi(:,it));

    % join variables that are time variant and vectors
    DispVj{1,it}=DispVi{1,it};
    DispVTj{1,it}=DispVTi{1,it};
    VFcj{1,it}=VFci{1,it};
    VFcTj{1,it}=VFcTi{1,it};
    Epc1vecj{1,it}=Epc1veci{1,it};
    Epc1vecCurj{1,it}=Epc1vecCuri{1,it};
    Epc2vecj{1,it}=Epc2veci{1,it};
    Epc2vecCurj{1,it}=Epc2vecCuri{1,it};
    epc1vecj{1,it}=epc1veci{1,it};
    epc2vecj{1,it}=epc2veci{1,it};
    for ip=2:nPairs
        DispVj{1,it}=[DispVi{1,it}; DispVi{ip,it}];
        DispVTj{1,it}=[DispVTi{1,it}; DispVTi{ip,it}];
        VFcj{1,it}=[VFcj{1,it}; VFci{ip,it}];
        VFcTj{1,it}=[VFcTj{1,it}; VFcTi{ip,it}];
        Epc1vecj{1,it}=[Epc1vecj{1,it}; Epc1veci{ip,it}];
        Epc1vecCurj{1,it}=[Epc1vecCurj{1,it}; Epc1vecCuri{ip,it}];
        Epc2vecj{1,it}=[Epc2vecj{1,it}; Epc2veci{ip,it}];
        Epc2vecCurj{1,it}=[Epc2vecCurj{1,it}; Epc2vecCuri{ip,it}];
        epc1vecj{1,it}=[epc1vecj{1,it}; epc1veci{ip,it}];
        epc2vecj{1,it}=[epc2vecj{1,it}; epc2veci{ip,it}];
    end
end

DIC_3Djoined_results=struct;
DIC_3Djoined_results.Faces=Fj;
DIC_3Djoined_results.FaceColors=FCj;
DIC_3Djoined_results.FaceIsoInd=Fisoj;
DIC_3Djoined_results.FaceCentroids=VFcj;
DIC_3Djoined_results.FaceCentroidsTransformed=VFcTj;
DIC_3Djoined_results.Points3D=Vj;
DIC_3Djoined_results.Points3Dtransformed=VTj;
DIC_3Djoined_results.corrComb=corrCombj;
DIC_3Djoined_results.FaceCorrComb=FaceCorrCombj;

DIC_3Djoined_results.Disp.DispVec=DispVj;
DIC_3Djoined_results.Disp.DispMgn=DispMj;
DIC_3Djoined_results.Disp.DispVecTransformed=DispVTj;
DIC_3Djoined_results.Disp.DispMgnTransformed=DispMTj;

DIC_3Djoined_results.Deform.J=Jj;
DIC_3Djoined_results.Deform.Lamda1=Lamda1j;
DIC_3Djoined_results.Deform.Lamda2=Lamda2j;
DIC_3Djoined_results.Deform.Emgn=Emgnj;
DIC_3Djoined_results.Deform.emgn=emgnj;
DIC_3Djoined_results.Deform.Epc1=Epc1j;
DIC_3Djoined_results.Deform.Epc2=Epc2j;
DIC_3Djoined_results.Deform.Epc1vec=Epc1vecj;
DIC_3Djoined_results.Deform.Epc1vecCur=Epc1vecCurj;
DIC_3Djoined_results.Deform.Epc2vec=Epc2vecj;
DIC_3Djoined_results.Deform.Epc2vecCur=Epc2vecCurj;
DIC_3Djoined_results.Deform.epc1=epc1j;
DIC_3Djoined_results.Deform.epc2=epc2j;
DIC_3Djoined_results.Deform.epc1vec=epc1vecj;
DIC_3Djoined_results.Deform.epc2vec=epc2vecj;

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