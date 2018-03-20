function installMultiDIC
%% Add MultiDIC library path so functions are known to use here

mPath=fileparts(mfilename('fullpath')); %Get the current path
addpath(genpath(fullfile(mPath,'lib_MultiDIC'))); % add libraries
addpath(genpath(fullfile(mPath,'lib_ext'))); % add external library

h=msgbox('MultiDIC installed successfully');
hp=get(h, 'position');
set(h, 'position', [hp(1) hp(2) 150 50]); %makes box bigger
ah = get( h, 'CurrentAxes' );
ch = get( ah, 'Children' );
set(ch, 'FontSize',10);

end