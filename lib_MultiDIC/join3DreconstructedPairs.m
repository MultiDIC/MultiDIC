function DIC_3Djoined_results = join3DreconstructedPairs(DIC3DAllPairsResults)
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
Vi=cell(nPairs,nFrames);
corrCombi=cell(nPairs,nFrames);
FaceCorrCombi=cell(nPairs,nFrames);
DispVi=cell(nPairs,nFrames);
DispMi=cell(nPairs,nFrames);
VFci=cell(nPairs,nFrames);



for ip=1:nPairs
    DIC_3Dpair_struct=DIC3DAllPairsResults{ip};
    
    Fi{ip}=DIC_3Dpair_struct.Faces;
    FCi{ip}=DIC_3Dpair_struct.FaceColors;
      
    for it=1:nFrames
        Vi{ip,it}=DIC_3Dpair_struct.Points3D{it};
        corrCombi{ip,it}=DIC_3Dpair_struct.corrComb{it};
        FaceCorrCombi{ip,it}=DIC_3Dpair_struct.FaceCorrComb{it};
        DispVi{ip,it}=DIC_3Dpair_struct.Disp.DispVec{it};
        DispMi{ip,it}=DIC_3Dpair_struct.Disp.DispMgn{it};
        VFci{ip,it}=DIC_3Dpair_struct.FaceCentroids{it};
    end

end

% join variables that are time invariant
[Fj,~,~]=joinElementSets(Fi,Vi(:,1),[]);
[~,~,FCj]=joinElementSets(Fi,Vi(:,1),FCi);


for it=1:nFrames
    % join variables that are time variant and scalar
    
    [~,Vj{1,it},DispMj{1,it}]=joinElementSets(Fi,Vi(:,it),DispMi(:,it));
    [~,~,corrCombj{1,it}]=joinElementSets(Fi,Vi(:,it),corrCombi(:,it));
    [~,~,FaceCorrCombj{1,it}]=joinElementSets(Fi,Vi(:,it),FaceCorrCombi(:,it));

    % join variables that are time variant and vectors
    DispVj{1,it}=DispVi{1,it};
    VFcj{1,it}=VFci{1,it};

    for ip=2:nPairs
        DispVj{1,it}=[DispVi{1,it}; DispVi{ip,it}];
        VFcj{1,it}=[VFcj{1,it}; VFci{ip,it}];
    end
end

DIC_3Djoined_results=struct;
DIC_3Djoined_results.Faces=Fj;
DIC_3Djoined_results.FaceColors=FCj;
DIC_3Djoined_results.FaceCentroids=VFcj;
DIC_3Djoined_results.Points3D=Vj;
DIC_3Djoined_results.corrComb=corrCombj;
DIC_3Djoined_results.FaceCorrComb=FaceCorrCombj;

DIC_3Djoined_results.Disp.DispVec=DispVj;
DIC_3Djoined_results.Disp.DispMgn=DispMj;

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