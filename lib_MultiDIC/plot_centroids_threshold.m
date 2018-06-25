function [Cfinal] = plot_centroids_threshold(IM,Np)
%% function to compute the centroids of black dots on white background from an image in step 1
% the user can tune the threshold level and click finish when done
%
% INPUT:
% * IM: grayscale masked image with black representing the dots
% * Np: number of points to be identified on the image
%
% OUTPUT:
% * Cfinal: centroids of dots

%%
Ctemp=[];

% find gray level seperating white from black (gray from black)
level = multithresh(IM,2);
level = double(level(1));

% turn gray image to black&white
IMbw = imbinarize(uint8(IM),level/255);

% find centroids
[Ctemp,J]=findCentroids(IMbw,Np);

% plot BW images over original with centroids on figures
fh=figure; hold all
fh.Units='normalized'; fh.Position=[.02 .1 .9 .8];
axis off
ax1 = axes;
imagesc(IM); hold on
ax2 = axes;
hb=imagesc(J); hold on
linkaxes([ax1,ax2]); % Link them together

% Hide the top and bottom axes
ax1.Visible = 'off';
pbaspect(ax1,[size(IM,2) size(IM,1) 1])
ax2.Visible = 'off';
pbaspect(ax2,[size(IM,2) size(IM,1) 1])
colormap(ax1,'gray');
colormap(ax2,'parula'); % Give each one its own colormap

%set transparacy
alphaValue=0.5;
alpha_data=double(J);
alpha_data(alpha_data~=0)=alphaValue;
set(hb, 'AlphaData', alpha_data);
set(ax2,'Position',ax1.Position);

% plot centroids
hp=plot(Ctemp(:,1),Ctemp(:,2),'b+','LineWidth',1,'MarkerSize',5);
ht=gtitle(['Threshold level = ' num2str(level)],20);

% set axes limits
[r,c]=find(IM~=255);
xlim([min(c) max(c)]);
ylim([min(r) max(r)]);

% getpixelposition(fh)

% Add the slider and slider label text to the figure
uicontrol('Style','slider','Units','Normalized','Position',[.1 .3 .02 .4],'value',level(1),'min',0,'max',255,'Callback',@updateCentroids);
uicontrol('Parent',fh,'Style','text','Units','Normalized','Position',[.1 .25 .04 .04],'String','0','HorizontalAlignment','left','fontsize',12,'BackgroundColor',fh.Color);
uicontrol('Parent',fh,'Style','text','Units','Normalized','Position',[.1 .72 .06 .04],'String','255','fontsize',12,'HorizontalAlignment','left','BackgroundColor',fh.Color);
uicontrol('Parent',fh,'Style','text','Units','Normalized','Position',[.08 .85 .1 .1],'String','gray level threshold','fontsize',12,'HorizontalAlignment','left','BackgroundColor',fh.Color);
uicontrol('Style','pushbutton','Units','Normalized','Position',[.1 .2 .08 .04],'String','Finish','fontsize',12,'Callback','uiresume');

% wait until user clicks on 'Finish' for resuming
uiwait(gcf);

Cfinal=Ctemp;

%%
    function updateCentroids(source,~)

        levelTemp = source.Value;
        % turn gray image to black&white
        IMbw = imbinarize(uint8(IM),levelTemp/255);
        
        % find centroids
        [Ctemp,J]=findCentroids(IMbw,Np);
        
        %update plot
        hp.XData=Ctemp(:,1);
        hp.YData=Ctemp(:,2);
        hb.CData=J;
        alpha_data=double(J);
        alpha_data(alpha_data~=0)=alphaValue;
        set(hb, 'AlphaData', alpha_data);
        set(ht, 'String', ['Threshold level = ' num2str(round(levelTemp))]);


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