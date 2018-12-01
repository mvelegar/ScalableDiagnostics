% Plot results
%% Prelims
clear variables; close all; clc;

% Read the landmap info in
% Landmap information
ncfile_map='landmap.nc';
landmap=ncread(ncfile_map,'LANDMAP');
landmap=permute(landmap,[2 1]);

% Spatial variables
x=ncread(ncfile_map,'lon');% Longitude(-180:5:175)
y=ncread(ncfile_map,'lat');% Latitude(-89,-86:4:86,89)
% No dynamics above elevation 30, so limit analysis till this elevation
nlev=30; z=1:nlev;
nLon=length(x); nLat=length(y);
nCells=nLon*nLat*nlev;

% Limit lat and lev
% Using 32 latitudes, since Ben and I discovered yesterday the weird
% cliff at the south pole in the dominant spatial mode
lat_vec_lim=[-58 62]; lat_vec_ind_lim=NaN(size(lat_vec_lim));
for i=1:length(lat_vec_lim)
    lat_vec_ind_lim(i)=find(y==lat_vec_lim(i));
end
nlat=lat_vec_ind_lim(2)-lat_vec_ind_lim(1)+1;
latLim=y(lat_vec_ind_lim(1):lat_vec_ind_lim(2));
landmapLim=landmap(lat_vec_ind_lim(1):lat_vec_ind_lim(2),:);

% The 6 chemical species of interest
chem_species=cellstr(...
    ['NO  ';
    'NO2 ';
    'O3  ';
    'OH  ';
    'ISOP';
    'CO  ';]);

% For plotting
[X3d,Y3d,Z3d]=meshgrid(x,y,z);
[X2d,Y2d]=meshgrid(x,y);

% Time info
aDay = 72;
% 30 days in July,September,November,April,June
% 31 days in August,October,December,January,March,May
% 28 days in February
nDays = 31*6+30*5+28;
nSnapsMonth=[30 31 30 31 30 31 31 28 31 30 31 30]'*aDay;
nSnaps=load('/Volumes/BRIO/preprocessed/nSnaps.mat');
nSnaps=nSnaps.nSnaps;
nSnapsTotal=sum(nSnaps); months=length(nSnaps);
addpath('cmocean_v1.4/cmocean/');

nFig=1;
%% Load the singular values, temporal modes, spatial modes
inFileNameSVDPre='results/svd_';
inFileNameNMFPre='results/nmf_';
inFileNameSPCAPre='results/spca_';
inFileNameSuf='.h5';

% I have 50 modes for SVD, 8 modes for NMF, 25 modes for SPCA
nModesAll=[50 8 25];
% Take the first 8 modes for now, to plot the U and Vt
nModes=8; iChem = 3;

WfromNMF=NaN(nCells,nModes);
HfromNMF=NaN(nSnapsTotal,nModes);

ZfromSPCA=load('results/Z.mat'); ZfromSPCA=ZfromSPCA.Z;
froNorm=load('results/froNorm.mat');
froNorm=froNorm.froNorm;

inFile=[inFileNameSVDPre,chem_species{iChem},inFileNameSuf];
temp=h5read(inFile,'/U'); temp=temp';
UfromSVD=temp(:,1:nModes);
temp=h5read(inFile,'/Vt');
VfromSVD=temp(:,1:nModes);
SfromSVD=h5read(inFile,'/S');
cumEngfromSVD=h5read(inFile,'/cum');

inFile=[inFileNameNMFPre,chem_species{iChem},inFileNameSuf];
% The ordered W and H
temp=h5read(inFile,'/W'); temp=temp'; W=temp;
temp=h5read(inFile,'/H'); H=temp;
% S
SfromNMF=h5read(inFile,'/S');
% Normalize such that 2-norm is 1
C=sum(W.^2,1); R=sum(H.^2,1);
for iMode=1:nModes
    HfromNMF=H(:,iMode)./sqrt(R(iMode));
    WfromNMF=W(:,iMode)./sqrt(C(iMode));
end
cumEngfromNMF=h5read(inFile,'/cum');


inFile=[inFileNameSPCAPre,chem_species{iChem},inFileNameSuf];
temp=h5read(inFile,'/A'); temp=temp';
AfromSPCA=temp(:,1:nModes);
eValsfromSPCA=h5read(inFile,'/eigvals');
cumEngfromSPCA=h5read(inFile,'/cum');

inFile=[inFileNameSPCAPre,chem_species{iChem},inFileNameMidSparse,...
    inFileNameSuf];
temp=h5read(inFile,'/A'); temp=temp';
AfromSPCASparse=temp(:,1:nModes);
eValsfromSPCASparse=h5read(inFile,'/eigvals');
cumEngfromSPCASparse=h5read(inFile,'/cum');

%% Cumulative energy spectrum and Temporal modes from SVD
% Plot first 4 temporal modes
nDom=4;
colorSetNdom=distinguishable_colors(nDom);
fontSize=28;
nSnapsTotalPlot=nSnapsTotal/aDay; nDaysMonth=[30,31,30,31,30,31,31,28,31,30,31,30];
figure(nFig);
S=cumEngfromSVD*100;
ha = tight_subplot(2,1,[.15 .01],[.2 .1],[.15 .1]); axes(ha(1));
plot(S,'k','Linewidth',4);
grid on; xlim([1 nModesAll(1)-100]);
set(gca,'LineWidth',4,'ColorOrder',colorSetNdom,...
    'XTick',[1,50:50:150],'XTickLabel',{'1','50','100','150'},...
    'Fontsize',fontSize);

ax = gca;
ax.YAxis.TickLabelFormat = '%,.1f';
xlabel('$\mathrm{Modes}$','Interpreter','Latex','Fontsize',fontSize);
ylabel('$\%\;\mathrm{energy}$','Interpreter','Latex',...
    'Fontsize',fontSize);


axes('Position',[.73 .675 .15 .15])
box on
plot(S,'ko','Linewidth',10,'Markersize',16);
hold on;
for iDom=1:nDom
    semilogy(iDom,S(iDom),'o','Linewidth',10,'Markersize',16);
end
hold off; grid on;
set(gca,'LineWidth',4,'ColorOrder',colorSetNdom,...
    'YTick',linspace(S(1),S(10),4),...
    'Fontsize',24);
ax = gca;
ax.YAxis.TickLabelFormat = '%,.1f';
xlim([1 10]); ylim([S(1) S(10)]);
V=VfromSVD(aDay/2:aDay:end,1:nDom);
for iDom=1:nDom
    V(:,iDom)=V(:,iDom)/norm(V(:,iDom),inf);
end

axes(ha(2));
plot(V,'Linewidth',4); xlim([0 192]); ylim([-1.5 1.5]);
set(gca,'LineWidth',4,'ColorOrder',colorSetNdom,...
    'XTick',cumsum(nDaysMonth)-15,'XTickLabel',...
    {'Jul','Aug','Sep','Oct','Nov','Dec','Jan','Feb','Mar','Apr',...
    'May','Jun'},'XTickLabelRotation',45,...
    'YTick',-1:1,'YTickLabel',{'-1','0','1'},'Fontsize',fontSize);
grid on; xlim([1 nSnapsPlot]);
xlabel('$\mathrm{Time}$','Interpreter','Latex','Fontsize',fontSize);
ylabel('$\mathrm{Temporal\;Modes}$','Interpreter','Latex',...
    'Fontsize',fontSize);
drawnow;

%% Cumulative energy spectrum and Temporal modes from NMF
nFig=nFig+1;
S=cumEngfromNMF*100;
figure(nFig);
ha = tight_subplot(2,1,[.15 .01],[.2 .1],[.15 .1]); axes(ha(1));
plot(S,'k','Linewidth',4);
set(gca,'LineWidth',4,'ColorOrder',colorSetNdom,...
    'Fontsize',fontSize);
ax = gca;
ax.YAxis.TickLabelFormat = '%,.0f';
grid on; xlim([1 nModesAll(2)]); ylim([10 35]);
xlabel('$\mathrm{Modes}$','Interpreter','Latex','Fontsize',fontSize);
ylabel('$\%\;\mathrm{energy}$','Interpreter','Latex',...
    'Fontsize',fontSize);


axes('Position',[.73 .66 .15 .15])
box on
plot(S,'ko','Linewidth',10,'Markersize',16);
hold on;
for iDom=1:nDom
    semilogy(iDom,S(iDom),'o','Linewidth',10,'Markersize',16);
end
hold off; grid on;
set(gca,'LineWidth',4,'ColorOrder',colorSetNdom,...
    'YTick',linspace(S(1),S(5),4),...
    'Fontsize',24);
ax = gca;
ax.YAxis.TickLabelFormat = '%,.1f';
xlim([1 5]); ylim([S(1) S(5)]);

axes(ha(2));
H=HfromNMF(aDay/2:aDay:end,1:nDom);
for iDom=1:nDom
    H(:,iDom)=H(:,iDom)/norm(H(:,iDom),inf);
end
plot(H,'Linewidth',4); xlim([0 192]); ylim([-0.2 1.2]);
set(gca,'LineWidth',4,'ColorOrder',colorSetNdom,...
    'XTick',cumsum(nDaysMonth)-15,'XTickLabel',...
    {'Jul','Aug','Sep','Oct','Nov','Dec','Jan','Feb','Mar','Apr',...
    'May','Jun'},'XTickLabelRotation',45,...
    'YTick',0:0.5:1,'YTickLabel',{'0','0.5','1'},'Fontsize',fontSize);
grid on; xlim([1 nSnapsPlot]);
xlabel('$\mathrm{Time}$','Interpreter','Latex','Fontsize',fontSize);
ylabel('$\mathrm{Temporal\;Modes}$','Interpreter','Latex',...
    'Fontsize',fontSize)

%% Cumulative energy spectrum and Temporal modes from SPCA
% Plot first 4 temporal modes
nDom=4;
colorSetNdom=distinguishable_colors(nDom);
close all; nFig=nFig+1; fontSize=28;
nSnapsPlot=nSnapsTotal/aDay;
nDaysMonth=[30,31,30,31,30,31,31,28,31,30,31,30];
figure(nFig);

S=cumEngfromSPCA*100;
ha = tight_subplot(2,1,[.15 .01],[.2 .1],[.15 .1]); axes(ha(1));
plot(S,'k','Linewidth',4);
grid on; xlim([1 nModesAll(3)]); ylim([30 75]);
set(gca,'LineWidth',4,'ColorOrder',colorSetNdom,...
    'XTick',[1,5:5:25],'XTickLabel',{'1','5','10','15','20','25'},...
    'Fontsize',fontSize);
ax = gca;
ax.YAxis.TickLabelFormat = '%,.0f';
xlabel('$\mathrm{Modes}$','Interpreter','Latex','Fontsize',fontSize);
ylabel('$\%\;\mathrm{energy}$','Interpreter','Latex',...
    'Fontsize',fontSize);


axes('Position',[.73 .66 .15 .15])
box on
plot(S,'ko','Linewidth',10,'Markersize',16);
hold on;
for iDom=1:nDom
    semilogy(iDom,S(iDom),'o','Linewidth',10,'Markersize',16);
end
hold off; grid on;
set(gca,'LineWidth',4,'ColorOrder',colorSetNdom,...
    'YTick',linspace(S(1),S(10),4),...
    'Fontsize',24);
ax = gca;
ax.YAxis.TickLabelFormat = '%,.1f';
xlim([1 10]); ylim([S(1) S(10)]);

V=ZfromSPCA(aDay/2:aDay:end,1:nDom);
for iDom=1:nDom
    V(:,iDom)=V(:,iDom)/norm(V(:,iDom),inf);
end

axes(ha(2));
plot(V,'Linewidth',4); xlim([0 192]); ylim([-1.5 1.5]);
set(gca,'LineWidth',4,'ColorOrder',colorSetNdom,...
    'XTick',cumsum(nDaysMonth)-15,'XTickLabel',...
    {'Jul','Aug','Sep','Oct','Nov','Dec','Jan','Feb','Mar','Apr',...
    'May','Jun'},'XTickLabelRotation',45,...
    'YTick',-1:1,'YTickLabel',{'-1','0','1'},'Fontsize',fontSize);
grid on; xlim([1 nSnapsPlot]);
xlabel('$\mathrm{Time}$','Interpreter','Latex','Fontsize',fontSize);
ylabel('$\mathrm{Temporal\;Modes}$','Interpreter','Latex',...
    'Fontsize',fontSize);
drawnow;

%% The 2D Plots
% For Figure 4: Plot the spatial mode 1 for O3 for 4 elevations
nFig=1; close all;
levs=[1,10:10:30]; iMode=1;
% addpath('kakearney-cptcmap-pkg-845bf83/cptcmap/');
% cmap='GMT_jet';
addpath('cmocean_v1.4/cmocean/');
cmap = cmocean('balance');
U=UfromSVD(:,iMode);
U=reshape(U,nLon,nLat,nlev); U=squeeze(U); U=U(:,:,levs);
U=permute(U,[2 1 3]);
figure(nFig);
ha = tight_subplot(2,2,[.01 .01],[.15 0.01],[.15 .15]);
cmin=min(U(:)); cmax=max(U(:));

for iLev=1:length(levs)
    axes(ha(iLev));
    pcolor(X2d,Y2d,U(:,:,iLev));
    if iLev == 1 || iLev == 3
        set(gca,...
            'ytick',[-90,-62,-30,2,30,62,90],...
            'LineWidth',4,...
            'FontSize',fontSize,...
            'YTickLabel',{'-90','-62','-30','2','30','62','90'});
    end
    if iLev == 3 || iLev==4
        set(gca,...
            'xtick',[-90,0,90],'LineWidth',4,...
            'FontSize',fontSize,'XTickLabel',...
            ({'-90','0','90'}));
    end
    if iLev == 2
        set(gca,'LineWidth',4);
    end
    if iLev == 1 || iLev == 3
        ylabel('$\mathrm{Lat}$','Interpreter','Latex','Fontsize',fontSize);
    end
    if iLev == 3 || iLev == 4
        xlabel('$\mathrm{Lon}$','Interpreter','Latex','Fontsize',fontSize);
    end
    % center the colormap
    cLimCenter=max(abs([cmin cmax]));
    caxis([-cLimCenter cLimCenter]); 
    colormap(cmap);
    shading interp;
    hold on;
    [~,h]=contour(X2d,Y2d,landmap,[0.9999 0.9999]);
    h.LineColor='k'; h.LineWidth=2; hold off;
    
end

set(ha(1:2),'XTickLabel',''); set(ha(2:2:4),'YTickLabel','');
hp1 = get(ha(1),'Position');
hp2 = get(ha(2),'Position');
hp3 = get(ha(3),'Position');
hp4 = get(ha(4),'Position');

cb_x = hp4(1) + hp4(3) + 0.005;
cb_y = hp4(2);
cb_w = 0.025;
cb_h = hp2(2) + hp2(4) - hp4(2);
colorbar('Position', [cb_x cb_y cb_w cb_h]);
ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0
    1],'Box','off','Visible','off','Units','normalized',...
    'clipping' , 'off');

drawnow; pause(0.05); clear U;
nFig=nFig+1;

%======================================================================
U=AfromSPCA(:,iMode);
U=reshape(U,nLon,nLat,nlev); U=squeeze(U); U=U(:,:,levs);
U=permute(U,[2 1 3]);
figure(nFig);
ha = tight_subplot(2,2,[.01 .01],[.15 0.01],[.15 .15]);
cmin=min(U(:)); cmax=max(U(:));
for iLev=1:length(levs)
    axes(ha(iLev));
    pcolor(X2d,Y2d,U(:,:,iLev));
    if iLev == 1 || iLev == 3
        set(gca,...
            'ytick',[-90,-62,-30,2,30,62,90],...
            'LineWidth',4,...
            'FontSize',fontSize,...
            'YTickLabel',{'-90','-62','-30','2','30','62','90'});
    end
    if iLev == 3 || iLev==4
        set(gca,...
            'xtick',[-90,0,90],'LineWidth',4,...
            'FontSize',fontSize,'XTickLabel',...
            ({'-90','0','90'}));
    end
    if iLev == 2
        set(gca,'LineWidth',4);
    end
    if iLev == 1 || iLev == 3
        ylabel('$\mathrm{Lat}$','Interpreter','Latex','Fontsize',fontSize);
    end
    if iLev == 3 || iLev == 4
        xlabel('$\mathrm{Lon}$','Interpreter','Latex','Fontsize',fontSize);
    end
    % center the colormap
    cLimCenter=max(abs([cmin cmax]));
    caxis([-cLimCenter cLimCenter]);
    colormap(cmap);
    shading interp;
    hold on;
    [~,h]=contour(X2d,Y2d,landmap,[0.9999 0.9999]);
    h.LineColor='k'; h.LineWidth=2; hold off;
end

set(ha(1:2),'XTickLabel',''); set(ha(2:2:4),'YTickLabel','');
hp1 = get(ha(1),'Position');
hp2 = get(ha(2),'Position');
hp3 = get(ha(3),'Position');
hp4 = get(ha(4),'Position');

cb_x = hp4(1) + hp4(3) + 0.005;
cb_y = hp4(2);
cb_w = 0.025;
cb_h = hp2(2) + hp2(4) - hp4(2);
colorbar('Position', [cb_x cb_y cb_w cb_h]);
ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0
    1],'Box','off','Visible','off','Units','normalized',...
    'clipping' , 'off');
drawnow; pause(0.05); clear U;
nFig=nFig+1;
%======================================================================
U=WfromNMF(:,iMode);
U=reshape(U,nLon,nLat,nlev); U=squeeze(U); U=U(:,:,levs);
U=permute(U,[2 1 3]);
figure(nFig);
ha = tight_subplot(2,2,[.01 .01],[.15 0.01],[.15 .15]);
cmin=min(U(:)); cmax=max(U(:));
iDom=1;
for iLev=1:length(levs)
    axes(ha(iLev));
    pcolor(X2d,Y2d,U(:,:,iLev));
    if iLev == 1 || iLev == 3
        set(gca,...
            'ytick',[-90,-62,-30,2,30,62,90],...
            'LineWidth',4,...
            'FontSize',fontSize,...
            'YTickLabel',{'-90','-62','-30','2','30','62','90'});
    end
    if iLev == 3 || iLev==4
        set(gca,...
            'xtick',[-90,0,90],'LineWidth',4,...
            'FontSize',fontSize,'XTickLabel',...
            ({'-90','0','90'}));
    end
    if iLev == 2
        set(gca,'LineWidth',4);
    end
    if iLev == 1 || iLev == 3
        ylabel('$\mathrm{Lat}$','Interpreter','Latex','Fontsize',fontSize);
    end
    if iLev == 3 || iLev == 4
        xlabel('$\mathrm{Lon}$','Interpreter','Latex','Fontsize',fontSize);
    end
    caxis([cmin cmax]); 
    cmap=cmocean('amp');
    colormap(cmap);
    shading interp;
    hold on;
    [~,h]=contour(X2d,Y2d,landmap,[0.9999 0.9999]);
    h.LineColor='k'; h.LineWidth=2; hold off;
end

set(ha(1:2),'XTickLabel',''); set(ha(2:2:4),'YTickLabel','');
hp1 = get(ha(1),'Position');
hp2 = get(ha(2),'Position');
hp3 = get(ha(3),'Position');
hp4 = get(ha(4),'Position');

cb_x = hp4(1) + hp4(3) + 0.005;
cb_y = hp4(2);
cb_w = 0.025;
cb_h = hp2(2) + hp2(4) - hp4(2);
colorbar('Position', [cb_x cb_y cb_w cb_h]);
ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0
    1],'Box','off','Visible','off','Units','normalized',...
    'clipping' , 'off');
drawnow; pause(0.05); clear U;
nFig=nFig+1;

%% For Figure 5: Plot the spatial modes 1-8 for O3 at surface
addpath('cmocean_v1.4/cmocean/');

nDom=8;

U=UfromSVD(:,1:nDom);
U=reshape(U,nLon,nLat,nlev,nDom); U=U(:,:,1,:); U=squeeze(U);
U=permute(U,[2 1 3]);
figure(nFig);
ha = tight_subplot(4,2,[.01 .01],[.15 0.01],[.15 .15]);
cmin=min(U(:)); cmax=max(U(:));
for iDom=1:nDom
    axes(ha(iDom));
    pcolor(X2d,Y2d,U(:,:,iDom));
    if mod(iDom,2) == 1
        set(gca,...
            'ytick',[-62,2,62],...
            'LineWidth',4,...
            'FontSize',fontSize,...
            'YTickLabel',{'-62','2','62'});
        ylabel('$\mathrm{Lat}$','Interpreter','Latex','Fontsize',fontSize);
    end
    if iDom > 6
        set(gca,...
            'xtick',[-90,0,90],'LineWidth',4,...
            'FontSize',fontSize,'XTickLabel',...
            ({'-90','0','90'}));
        xlabel('$\mathrm{Lon}$','Interpreter','Latex','Fontsize',fontSize);
    end
    if mod(iDom,2) == 0
        set(gca,'LineWidth',4);
    end
    % center the colormap
    cLimCenter=max(abs([cmin cmax]));
    caxis([-cLimCenter cLimCenter]); 
    cmap = cmocean('balance');
    colormap(cmap);
    shading interp;
    hold on;
    [~,h]=contour(X2d,Y2d,landmap,[0.9999 0.9999]);
    h.LineColor='k'; h.LineWidth=2; hold off;
end

set(ha(1:6),'XTickLabel',''); set(ha(2:2:8),'YTickLabel','');
hp2 = get(ha(2),'Position');
hp8 = get(ha(8),'Position');

cb_x = hp8(1) + hp8(3) + 0.005;
cb_y = hp8(2);
cb_w = 0.025;
cb_h = hp2(2) + hp2(4) - hp8(2);
colorbar('Position', [cb_x cb_y cb_w cb_h]);
ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0
    1],'Box','off','Visible','off','Units','normalized',...
    'clipping' , 'off');
drawnow; pause(0.05); clear U;
nFig=nFig+1;

%======================================================================
U=AfromSPCA(:,1:nDom);
U=reshape(U,nLon,nLat,nlev,nDom); U=U(:,:,1,:); U=squeeze(U);
U=permute(U,[2 1 3]);
figure(nFig);
ha = tight_subplot(4,2,[.01 .01],[.15 0.01],[.15 .15]);
cmin=min(U(:)); cmax=max(U(:));
for iDom=1:nDom
    axes(ha(iDom));
    pcolor(X2d,Y2d,U(:,:,iDom));
    if mod(iDom,2) == 1
        set(gca,...
            'ytick',[-62,2,62],...
            'LineWidth',4,...
            'FontSize',fontSize,...
            'YTickLabel',{'-62','2','62'});
        ylabel('$\mathrm{Lat}$','Interpreter','Latex','Fontsize',fontSize);
    end
    if iDom > 6
        set(gca,...
            'xtick',[-90,0,90],'LineWidth',4,...
            'FontSize',fontSize,'XTickLabel',...
            ({'-90','0','90'}));
        xlabel('$\mathrm{Lon}$','Interpreter','Latex','Fontsize',fontSize);
    end
    if mod(iDom,2) == 0
        set(gca,'LineWidth',4);
    end
    % center the colormap
    cLimCenter=max(abs([cmin cmax]));
    caxis([-cLimCenter cLimCenter]);  
    cmap = cmocean('balance');
    colormap(cmap);
    shading interp;
    hold on;
    [~,h]=contour(X2d,Y2d,landmap,[0.9999 0.9999]);
    h.LineColor='k'; h.LineWidth=2; hold off;
end

set(ha(1:6),'XTickLabel',''); set(ha(2:2:8),'YTickLabel','');
hp2 = get(ha(2),'Position');
hp8 = get(ha(8),'Position');

cb_x = hp8(1) + hp8(3) + 0.005;
cb_y = hp8(2);
cb_w = 0.025;
cb_h = hp2(2) + hp2(4) - hp8(2);
colorbar('Position', [cb_x cb_y cb_w cb_h]);
ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0
    1],'Box','off','Visible','off','Units','normalized',...
    'clipping' , 'off');
drawnow; pause(0.05); clear U;
nFig=nFig+1;
%======================================================================
U=WfromNMF(:,1:nDom);
U=reshape(U,nLon,nLat,nlev,nDom); U=U(:,:,1,:); U=squeeze(U);
U=permute(U,[2 1 3]);
figure(nFig);
ha = tight_subplot(4,2,[.01 .01],[.15 0.01],[.15 .15]);
cmin=min(U(:)); cmax=max(U(:));
for iDom=1:nDom
    axes(ha(iDom));
    pcolor(X2d,Y2d,U(:,:,iDom));
    if mod(iDom,2) == 1
        set(gca,...
            'ytick',[-62,2,62],...
            'LineWidth',4,...
            'FontSize',fontSize,...
            'YTickLabel',{'-62','2','62'});
        ylabel('$\mathrm{Lat}$','Interpreter','Latex','Fontsize',fontSize);
    end
    if iDom > 6
        set(gca,...
            'xtick',[-90,0,90],'LineWidth',4,...
            'FontSize',fontSize,'XTickLabel',...
            ({'-90','0','90'}));
        xlabel('$\mathrm{Lon}$','Interpreter','Latex','Fontsize',fontSize);
    end
    if mod(iDom,2) == 0
        set(gca,'LineWidth',4);
    end
    caxis([cmin cmax]); 
    cmap=cmocean('amp');
    colormap(cmap);
    shading interp;
    hold on;
    [~,h]=contour(X2d,Y2d,landmap,[0.9999 0.9999]);
    h.LineColor='k'; h.LineWidth=2; hold off;
end

set(ha(1:6),'XTickLabel',''); set(ha(2:2:8),'YTickLabel','');
hp2 = get(ha(2),'Position');
hp8 = get(ha(8),'Position');

cb_x = hp8(1) + hp8(3) + 0.005;
cb_y = hp8(2);
cb_w = 0.025;
cb_h = hp2(2) + hp2(4) - hp8(2);
colorbar('Position', [cb_x cb_y cb_w cb_h]);
ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0
    1],'Box','off','Visible','off','Units','normalized',...
    'clipping' , 'off');

drawnow; pause(0.05); clear U;
nFig=nFig+1;


