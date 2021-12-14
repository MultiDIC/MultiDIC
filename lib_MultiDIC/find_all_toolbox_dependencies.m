% find dependencies

% Naviagate to the MultiDIC main folder
clear;

folderPaths={"lib_MultiDIC";...
    "main_scripts";...
    fullfile("lib_ext","arrow3");...
    fullfile("lib_ext","export_fig");...
    fullfile("lib_ext","findjobj");...
    fullfile("lib_ext","GIBBON","lib");...
    fullfile("lib_ext","ncorr_2D_matlab-master");...
    fullfile("lib_ext","numsubplots");...
    fullfile("lib_ext","selectdata");...
    fullfile("lib_ext","subtightplot");...
    fullfile("lib_ext","uipickfiles")};

kk=1;
for ii=1:length(folderPaths)
    listing{ii,1} = dir(fullfile(folderPaths{ii,1},'*.m'));
    for jj=1:length(listing{ii,1})
        filePath=fullfile(listing{ii,1}(jj).folder,listing{ii,1}(jj).name);
        filesPaths{kk,1}=fullfile(listing{ii,1}(jj).folder,listing{ii,1}(jj).name);
        kk=kk+1;
    end
end

names=cell(length(filesPaths),7);
for ii=1:length(filesPaths)
    names_temp=dependencies.toolboxDependencyAnalysis(filesPaths{ii,1});
    names(ii,1)= filesPaths(ii,1);
    names(ii,2:length(names_temp)+1)= names_temp;
end

names_all=dependencies.toolboxDependencyAnalysis(filesPaths);
% [fList,plist] = matlab.codetools.requiredFilesAndProducts(filesPaths);

return
