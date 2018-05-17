function [numFreeBytes]=freeMemory

try
    if ispc        
        %Can be slow on windows
        [~,mem_stat] = memory;
        numFreeBytes = mem_stat.PhysicalMemory.Available; % numFreeBytes = mem_stat.MaxPossibleArrayBytes; %Alternative
%          [~,S]=system('wmic OS get FreePhysicalMemory');
    elseif isunix && ~ismac
        % Output format of running 'free -b | grep Mem'
        %  total       used       free     shared    buffers     cached
        
        %         [~,S] = unix('free -b | grep Mem'); % Excute free command and collect output in strin S
        %         mem_stat = str2double(regexp(S, '[0-9]*', 'match')); % Get the numbers
        %         numFreeBytes = mem_stat(3) + mem_stat(end) ;
        
        [~,numFreeBytesStr]=unix('free -b | awk ''/Mem/{print $3} /Mem/{print $6}''');
        numFreeBytes=sum(str2double(numFreeBytesStr));
    elseif ismac 
        %UNTESTED!
        [~,S] = unix('vm_stat | grep free'); % Excute vm_stat
        mem_stat = strfind(S,' '); % Detect spaces
        numFreeBytes = str2double(S(mem_stat(end):end))*4096; %Take last value and convert pages to bytes
    end
catch 
    error('Could not determine free memory');
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
