% depthPitch_LRAUV.m
% Last modified Jan 12, 2015
% Ben Raanan

% INSTRUCTIONS - depthPitch_LRAUV
%--------------------------------------------------------------------------
%
% 1) Un-zip folder and add to MATLAB path including sub-folders. If you?re
% asked about file encryption when un-zipping say yes to all (MAC vs PC bug).
%
% 2) Execute depthPitch_LRAUV from MATLAB command line (edit vehicle, year
% and search term in depthPitch_LRAUV.m as needed).
%
% 3) The first time you run the script you will be prompted to select the
% server path folder for: smb://atlas.shore.mbari.org/LRAUV/
% (e.g., \\atlas\LRAUV\).
%
% 4) After pointed to the server, matFilePaths_LRAUV.m crawls the server to
% create a map of the .mat files present on it. This will likely take a few
% (to several) minutes to run? The map (path inventory) is saved locally as
% matFilePaths_LRAUV.mat and is used later as a search database (see help
% findmat_LRAUV for more on how to use keywords).
% You will be be prompted refresh matFilePaths_LRAUV.mat if the file is
% older than 30 days.


% EDIT THESE VARAIBLES TO SEARCH FOR DIFFERENT DATASETS
%--------------------------------------------------------------------------
vehicle = 'Daphne';                 % 'Tethys', 'Daphne', 'Makai'
year=2015;                          % 2010 - 2015
searchterm = '201501302128_201501310254';     % see: help findmat_LRAUV 
                                              % Beach-party  dataset: '201309121813_201309140344' ; '20130911T162528'
                                              % Shark-attack dataset: '201309301141_201310070703'
                                              % Intresting anom (2014): '201410011829_201410020520' 201502032113_201502051041
%--------------------------------------------------------------------------

clearvars -except vehicle year searchterm
clear global
close all

% Get work directory paths
%--------------------------------------------------------------------------
fname=which('depthPitch_LRAUV.m');

workd='~/Documents/MATLAB/MBARI/';          % working directory
figd ='~/Documents/MBARI/';                 % folder for saving figures
outmatfolder=[workd 'mat/workver/'];

% workd=fname(1:end-18);                      % working directory
% figd=workd;                                 % folder for saving figures
% outmatfolder=[workd 'scripts' filesep 'mat' filesep];    % folder for saving mat files

%--------------------------------------------------------------------------

% save plot? (Y=1)
sv = 0;


% locate, fix and load dataset
%--------------------------------------------------------------------------

% syslogs
load syslog_AllVehicles_comp

% get mat file path(s)
[matpath, matname] = findmat_LRAUV( vehicle, year, searchterm, outmatfolder );

% read-in syslog and get CRITICAL error msgs and time-stamps
[ syslog ] = getsyslog( matpath );    % requires connection to LRAUV server
[ ccomp, critical_timestamp] = islogcritical( vehicle, year, matpath );

% get parameters of intrest and fix time-series
[ filename ] = processLargeMAT_LRAUV( vehicle, year, searchterm, outmatfolder );

% load reduced .mat file (workable format)
[ time, time_step, xstruct, names, controls ] = initialize_LRAUV_SIM( filename{:} );


% calc fin angle of attack (critical angle of attack - |alpha|>15)
[ attackAngle, stall, colNames ] = isstall( controls, xstruct );



% unpack and reduce data set
%--------------------------------------------------------------------------
% focus on time window
%{
t1 = time(1);   % datenum(2013,09,12,18,13,37);
t2 = time(end); % datenum(2013,09,12,23,20,00);

tInd = find(time>=t1 & time<=t2); % 1:length(time);
%}

% cut time vector for reshape bin avging (~5 sec bins)
int = round(3/time_step); timemod=time; c=0;
while mod(length(timemod),int)~=0 && c<200
    c=c+1;
    timemod = time(1:end-c);
end; clear timemod

t1 = time(1);
t2 = time(end-c);
tInd = find(time>=t1 & time<=t2); % 1:length(time);



% extract data
%-----------------------------------------------------------------
time   =  time(tInd);
z      =  xstruct.z(tInd);

xstruct.theta = kickout(xstruct.theta);
if isdeg(xstruct.theta)
    pitch  =  xstruct.theta(tInd);
else
    pitch  =  xstruct.theta(tInd)*(180/pi);
end
if isdeg(xstruct.Cmd.pitchCmd)
    pcmd   =  xstruct.Cmd.pitchCmd(tInd);
else
    pcmd   =  xstruct.Cmd.pitchCmd(tInd)*(180/pi);
end

if isdeg(controls(tInd,1))
    ele    =  controls(tInd,1); eleu = 'rad';
else
    ele    =  controls(tInd,1)*(180/pi); eleu = 'deg';
end
%}
dzcmd       =  xstruct.Cmd.depthRateCmd(tInd);
speed       =  xstruct.u(tInd);
speedCmd    =  xstruct.Cmd.speedCmd(tInd);
depthCmd    =  xstruct.Cmd.depthCmd(tInd);
vmode       =  xstruct.Cmd.verticalMode(tInd);
mass_p      =  xstruct.mass_p(tInd);
buoyancy_p  =  xstruct.buoyancy_p(tInd);
stall       =  stall(tInd,:);
attackAngle =  attackAngle(tInd,:)*(180/pi);
dzvh   = xstruct.z_rate_vh(tInd);
% dzdiff = xstruct.z_rate(tInd);



% bin avg. depth to smooth depth rate dz/dt
%-----------------------------------------------------------------
% f=find(mod(length(time),1:100)==0)';
% f = [f , time_step*f]

tbin = median(reshape(time, int, length(time)/int),1);
pbin = mean(reshape(pitch, int, length(z)/int),1);
% zbin = mean(reshape(z, int, length(z)/int),1);
dzvhbin = mean(reshape(dzvh, int, length(z)/int),1);
% dzdiffbin = mean(reshape(dzdiff, int, length(z)/int),1);

% dzbin = zeros(size(zbin));
% dzbin(2:end) = diff(zbin)./(time_step*int);

% interp smoothed signal to match dataset time grid
dzbi     = -interp1(tbin,dzvhbin,time);
pbi      = interp1(tbin,pbin,time);
dz       = -xstruct.z_rate(tInd); % (un-smoothed depth-rate)
% dzvhbi   = -interp1(tbin,dzvhbin,time);
% dzdiffbi = -interp1(tbin,dzdiffbin,time);



% compare
%{
figure;
plot(time, dzbi,'linewidth',2);         hold on;
plot(time, dzdiffbi,'linewidth',2);
plot(time, dzvhbi,'linewidth',2);       hold off;
legend('dzbi','dzdiffbi','dzvhbi')
set(gca,'layer','top','fontWeight','bold','fontsize',14)
axis tight; grid on; box on;
dynamicDateTicks(gca)
%}


% filtering criteria: index 1 m/s profiles
%-----------------------------------------------------------------
ind=(speedCmd>=0.7 & z>2 & time>=t1 & time<=t2);

% in = 1:5:length(ind); % reduce resoution
% ind = ind(in);

% extract 1 m/s data
timeInd  = time;                timeInd(~ind)   = NaN;
zInd     = z;                   zInd(~ind)      = NaN;
speedInd = speed;               speedInd(~ind)  = NaN;
dzInd    = dz;                  dzInd(~ind)     = NaN;
dzbiInd  = dzbi;                dzbiInd(~ind)   = NaN;
pitchInd = pitch;               pitchInd(~ind)  = NaN;
pbiInd   = pbi;                 pbiInd(~ind)    = NaN;
eleInd   = ele;                 eleInd(~ind)    = NaN;
mass_pInd = mass_p;             mass_pInd(~ind) = NaN;
buoyancy_pInd = buoyancy_p;     buoyancy_pInd(~ind) = NaN;
speedCmdInd = speedCmd;         speedCmdInd(~ind)   = NaN;


% Get mission names
%--------------------------------------------------------------------------
[ mission ] = getmission( syslog, time, z );


% Flag CRITICALs (if known) markes +/- 3 sec around CRITICAL error message
%--------------------------------------------------------------------------
%
if ~isempty(critical_timestamp)
    mals=zeros(size(critical_timestamp,1)+1,length(time));
    for k=1:size(critical_timestamp,1)
        mals(k+1,:) = (time>=critical_timestamp(k,2) & time<=critical_timestamp(k,3));
    end
    mal = logical(max(mals));
else
    mal = false(size(time));
end


maltime  = time;     maltime(~mal)   = NaN;
malz     = z;        malz(~mal)      = NaN;
maldzbi  = dzbi;     maldzbi(~mal)   = NaN;
malpitch = pitch;    malpitch(~mal)  = NaN;

% beach party malfunction:
% mal = (time>=datenum(2013,09,12,22,13,00) & time<=datenum(2013,09,12,22,24,00) & speedCmd >= 0.8 |...
%     time>=datenum(2013,09,12,21,08,02) & time<=datenum(2013,09,12,21,08,18));

%}

%% GUI plot
%--------------------------------------------------------------------------
% prep for plot
%-----------------------------------------------------------------

% Set filtering conditions
global eleCangle dzC qual
eleCangle = 10;      % Elevator angle threshold
dzC       = 0.1;    % Depth-rate threshold
qual = (ele<-eleCangle & dzbi'>-dzC) | (ele>eleCangle & dzbi'<dzC) |...
    (ele<0 & dzbi'<-0.4) | (pitch'<-30 & ele<0) | (pitch'<-35);


% plot
%-----------------------------------------------------------------
Xcode=0;
close all;
clear fig s1 s2 s3 s4 p11 p12 p13 p14 sc21 sc22 sc23 p41
global s4 p11 p12 p13 p14 sc21 sc22 sc23 p41
fig = figure(1);
set(gcf,'Units','normalized','Position',[0 0 1 1],...
    'PaperPositionMode','auto') %,'visible','off');

if verLessThan('matlab','8.4')
    set(gcf,'Renderer','zbuffer');
end



% subplot 1
%-----------------------------------------------------------------
s1 = subplot('position',[0.05 0.4 0.55 0.55]); % [left bottom width height]
hold on
p11=plot(time,z,'color', [0.301    0.7450    0.9330],'linewidth',0.8);
p12=plot(time,zInd,'color', [0    0.4470    0.7410],'linewidth',3);
p13=plot(maltime,malz,'.','color', [0.8500    0.3250    0.0980],'markersize',10);
p14=plot(timeInd(qual),zInd(qual),'o','color', [0.4940    0.1840    0.5560],'markersize',6);

if Xcode==1
    % Get Xcode EFC outputs
    %-----------------------------------------------------------------
    EFC = importdata('~/Desktop/seed/EFC_output.txt');
    indEFC = closest(EFC.data(:,1),time);
    plot(timeInd(indEFC),zInd(indEFC),'ro','markersize',6,'markerfacecolor','r');
end

hold off
lg1 = legend([p11,p12,p13,p14],'Platform speed < 0.8','Platform speed > 0.8',...
    'Log CRITICAL (\pm5sec)', 'Classified','location','best');
set([p13, p14],'Visible','off')
ylabel('Depth (m)','fontweight','bold','fontsize',20);
title([vehicle ' log: ' searchterm],...
    'fontweight','bold','fontsize',18, 'Interpreter', 'none');
set(gca,'YDir','reverse','layer','top','fontWeight','bold','fontsize',14)
axis tight; grid on; box on;
dynamicDateTicks(gca)
ps1 = get(s1, 'position'); % [left bottom width height]




% subplot 2
%-----------------------------------------------------------------
%
s2 = subplot('position',[0.65 0.4 0.3 0.55]); % [left bottom width height]
S=ones(size(ind))*9;
range = linspace(-15,15,1028); % max(peakf(:,2))
[cmap,col,colx] = scatcolor(ele',range);

hold on;
sc21 = scatter( dzbiInd, pitchInd,S,colx,'markerfacecolor','flat');
sc22 = scatter( dzbiInd(mal), pitchInd(mal),'MarkerEdgeColor',[0.8500    0.3250    0.0980]);
sc23 = scatter( dzbiInd(qual), pitchInd(qual),'MarkerEdgeColor',[0.4940    0.1840    0.5560]);
lg2 = legend([sc21, sc22, sc23],'Platform speed > 0.8',...
    'Log CRITICAL (\pm5sec)', 'Classified','location','best');
set([sc22, sc23],'Visible','off');
hold off;
if max(abs(dzbiInd))>0.5 || max(abs(pitchInd))>36
    axis tight; set(gca,'xtick',-1:0.1:1);
else
    xlim([-0.5 0.5]); ylim([-35 35]); set(gca,'xtick',-1:0.1:1);
end
title(['Bin avreged data (' num2str(int*time_step,2) ' sec)'],...
    'fontweight','bold','fontsize',18);
ylabel('Platform pitch (deg)','fontweight','bold','fontsize',20);
xlabel('Depth rate (m/s; negative=down)','fontweight','bold' ,...
    'fontsize',20)
set(gca,'layer','top','fontWeight','bold','fontsize',14);
box on; grid on;
ps2= get(s2, 'position'); % [left bottom width height]
hold off;



% subplot 3
%-----------------------------------------------------------------
s3 = subplot('position',[ps2(1) ps2(2)-0.1 ps2(3) 0.025]); % [left bottom width height]
scatcolorbar(range,cmap,'Elevator angle (deg)')
set(gca,'layer','top','fontWeight','bold','fontsize',14,'xgrid','on');
box on;
ps3=get(s3, 'Position');
%}


% subplot 4
%-----------------------------------------------------------------
s4 = subplot('position',[0.05 0.1 0.55 0.25]); % [left bottom width height]
p41=plot(time,dzbi,'color',[    0    0.4470    0.7410],'linewidth',2);
ylabel('Bin avg. depth rate (m/s)','fontweight','bold','fontsize',16);
xlabel('Hour (UTC)','fontweight','bold','fontsize',20);
set(gca,'layer','top','fontWeight','bold','fontsize',14);
box on; grid on;
linkaxes([s1,s4],'x')
dynamicDateTicks([s1,s4],'linked')




% GUI
%--------------------------------------------------------------------------

% prep GUI vars
%-----------------------------------------------------------------
clear guivars br b11 b12 b2 b3
global guivars br b11 b12 b2 b3

guivars.names = {'Bin avg. depth rate (m/s)','Platform depth rate (m/s)',...
    'Commanded depth rate (m/s)',...
    'Platform pitch angle (deg)','Commanded pitch (deg)',...
    'Bin avg. pitch angle (deg)',...
    'Elevator angle (deg)','Ele angle of attack (deg)',...
    'Buoyancy position (m)','Mass psition (m)',...
    'Platform velocity (m/s)','Commanded velocity (m/s)',...
    'Vertical Mode'};  % 'Commanded depth (m)',

guivars.time = time ;
guivars.data(1,:) = dzbi ;       guivars.data(2,:) = dz ;
guivars.data(3,:) = dzcmd ;
guivars.data(4,:) = pitch ;      guivars.data(5,:) = pcmd ;
guivars.data(6,:) = pbi;
guivars.data(7,:) = ele ;        guivars.data(8,:) = attackAngle(:,2) ;
guivars.data(9,:) = buoyancy_p ; guivars.data(10,:) = mass_p ;
guivars.data(11,:) = speed ;     guivars.data(12,:) = speedCmd ;
guivars.data(13,:) = vmode;      % guivars.data(13,:) = depthCmd;

% guivars.mission = mission;

guivars.slider.dzbi = dzbiInd ;  guivars.slider.pitchInd = pitchInd;
guivars.slider.timeInd = timeInd ;     guivars.slider.zInd = zInd;

% Toggle button: select data (brush)
br  = brush;
b11 = false(size(time)); b12 = false(size(time));
b2  = false(size(time)); b3  = false(size(time));


% launch GUI
%-----------------------------------------------------------------
GUI_depthPitch_LRAUV


% Script source water-mark
%-----------------------------------------------------------------
uicontrol('Style', 'text','String', [fname ' - ' datestr(clock)],...
    'Units','normalized',...
    'Position', [0 0 1 0.025]);
set(gcf,'visible','on')

% save
%--------------------------------------------------------------------------
%
if sv==1
    pslg1 = get(lg1,'position');
    set(lg1,'location','none','position',pslg1)
    matnamefig = matname{:};
    print(figure(1),[ figd 'binAvg_' vehicle '_' matnamefig(1:end-4) '' ],'-r300','-dpng');
    close(fig); clear pslg1
end
%}




%% plot missions
%--------------------------------------------------------------------------
%
if ~isempty(mission.name)
    figure(2);
    set(gcf,'Units','normalized','Position',[0 0.3 1 0.5],... % [left bottom width height]
        'PaperPositionMode','auto')
    plot(time , mission.z,'linewidth',2)
    set(gca,'YDir','reverse')
    lg3=legend(mission.namelist,'Interpreter','none','location','best');
    title([vehicle ' missions (start to start) in log: ' searchterm ],...
        'fontweight','bold','fontsize',18, 'Interpreter', 'none');
    ylabel('Depth (m)','fontweight','bold','fontsize',20);
    xlabel('Hour (UTC)','fontweight','bold','fontsize',20);
    set(gca,'layer','top','fontWeight','bold','fontsize',14);
    dynamicDateTicks(gca)
    grid on; box on;
    
    if sv==1
        pslg3 = get(lg3,'position');
        set(lg3,'location','none','position',pslg3)
        matnamefig = matname{:};
        print(figure(2),[figd 'mission_' vehicle '_' matnamefig(1:end-4) ],'-r300','-dpng');
        close; clear pslg3
    end
end
%}
ccomp
%% plot parameter subplots
%{
% Set condition
%--------------------------------------------------------------------------
cond=(stall(:,2)==1); % speed<=0.9 | pitch'<-25 & ele<-2.5 | dzbi<-0.35


% plot
%--------------------------------------------------------------------------
figure;
set(gcf,'Units','normalized','Position',[0 0 1 1],...
    'PaperPositionMode','auto');

ax1=subplot(4,1,1);
plot(time,z); axis tight;
set(gca,'yDir','reverse')
datetick('x','HH:MM')
hold on;
plot(time(cond),z(cond),'.');
ylabel('Depth (m)','fontweight','bold','fontsize',14);
ps1=get(ax1, 'Position');
%{
lg=legend('> 0.9 m/s',' < 0.9 m/s',...
    'location','eastoutside');
%}
lg=legend('~Qualify',' Qualify',...
    'location','eastoutside');
set(lg,'fontsize',14)
set(gca,'layer','top','fontWeight','bold')
axis tight; grid on;
set(ax1,'position',[ps1(1), ps1(2), ps1(3), ps1(4)])

ax2=subplot(4,1,2);
plot(time,pitch); axis tight;
datetick('x','HH:MM')
hold on;
plot(time(cond),pitch(cond),'.');
ylabel(['Platform pitch (deg)'],'fontweight','bold','fontsize',14);
set(gca,'layer','top','fontWeight','bold')
axis tight; grid on;

ax3=subplot(4,1,3);
plot(time,ele); axis tight;
datetick('x','HH:MM')
hold on;
plot(time(cond),ele(cond),'.');
ylabel(['Elevator angle (deg)'],'fontweight','bold','fontsize',14);
set(gca,'layer','top','fontWeight','bold')
axis tight; grid on;

ax4=subplot(4,1,4);
plot(time,dzbi); axis tight;
datetick('x','HH:MM')
hold on;
plot(time(cond),dzbi(cond),'.');
ylabel(['Depth rate' 10 '(m/s; negative=down)'],'fontweight','bold','fontsize',14);
set(gca,'layer','top','fontWeight','bold')
axis tight; grid on;

linkaxes([ax1,ax2,ax3,ax4],'x')
dynamicDateTicks([ax1,ax2,ax3,ax4], 'linked')
%}

%% plot angle of attack/stall
%--------------------------------------------------------------------------
%{
zp = z;
zp(stall(:,2)==0 | isnan(stall(:,1)))=NaN;
attackStall = attackAngle;
attackStall(stall==0)=NaN;

figure;
set(gcf,'Units','normalized','Position',[0 0.1 1 0.8],... % [left bottom width height]
    'PaperPositionMode','auto');
s1=subplot(211);
plot(time, z,'linewidth',2); hold on;
plot(time, zp,'linewidth',2)
set(gca,'YDir','reverse')
grid on; box on;
set(gca,'fontweight','bold', 'fontsize',11)

s2=subplot(212);
plot(time, attackAngle(:,2)); hold on;
plot(time, attackStall(:,2))
ylim([-15 15]);
grid on; box on;
set(gca,'fontweight','bold', 'fontsize',11)
linkaxes([s1,s2],'x')
dynamicDateTicks([s1,s2],'linked')

%
figure;
set(gcf,'Units','normalized','Position',[0 0.3 1 0.6],... % [left bottom width height]
    'PaperPositionMode','auto');
subplot('position',[0.05 0.1 0.9 0.8])
plot(time, attackAngle(:,2)); hold on;
plot(time, attackStall(:,2))
ylim([-15 15]);
grid on; box on;
set(gca,'fontweight','bold', 'fontsize',11)
dynamicDateTicks(gca)

%}

%% depth-rate plots
%--------------------------------------------------------------------------
%{
figure;
set(gcf,'Units','normalized','Position',[0 0.4 1 0.6],...
    'PaperPositionMode','auto');
s1=subplot(211);
plot(tbin,zbin)
set(gca,'YDir','reverse')
ylabel('Depth (m)','fontweight','bold', 'fontsize',14)
set(gca,'fontweight','bold', 'fontsize',11)
grid on

s2=subplot(212);
plot(tbin,-dzbin); hold on
plot(time,-dzcmd+0.05)
ylabel('Depth rate (m/s)','fontweight','bold', 'fontsize',14)
xlabel('12-13 Sep-2013','fontweight','bold', 'fontsize',14)
lg=legend('5sec bin avg.','Commanded');
lg.FontSize=14;
set(gca,'fontweight','bold', 'fontsize',11)
grid on
linkaxes([s1,s2],'x')
dynamicDateTicks([s1,s2],'linked')

figure;
set(gcf,'Units','normalized','Position',[0 0.4 1 0.6],...
    'PaperPositionMode','auto');
plot(time,dz); hold on
plot(time,dzbi,'linewidth',2)
plot(time,-dzcmd+0.05,'linewidth',2)
axis tight; grid on
ylabel('Depth rate (m/s)','fontweight','bold', 'fontsize',18)
xlabel('12-13 Sep-2013','fontweight','bold', 'fontsize',18)
lg=legend('Depth rate','5sec bin avg.','Commanded');
lg.FontSize=14;
set(gca,'fontweight','bold', 'fontsize',12)
dynamicDateTicks(gca)
%}

%% pitch plots
%--------------------------------------------------------------------------
%{
figure;
set(gcf,'Units','normalized','Position',[0 0.4 1 0.6],...
    'PaperPositionMode','auto');
s1=subplot(211);
plot(tbin,zbin)
set(gca,'YDir','reverse')
ylabel('Depth (m)','fontweight','bold', 'fontsize',14)
set(gca,'fontweight','bold', 'fontsize',11)
grid on

s2=subplot(212);
plot(tbin,pbin); hold on
plot(time,pcmd)
ylabel('Pitch angle (deg)','fontweight','bold', 'fontsize',14)
xlabel('12-13 Sep-2013','fontweight','bold', 'fontsize',14)
lg=legend('5sec bin avg.','Commanded');
lg.FontSize=14;
set(gca,'fontweight','bold', 'fontsize',11)
grid on
linkaxes([s1,s2],'x')
dynamicDateTicks([s1,s2],'linked')


figure;
set(gcf,'Units','normalized','Position',[0 0.4 1 0.6],...
    'PaperPositionMode','auto');
plot(time,pitch); hold on
plot(time,pbi,'linewidth',2)
plot(time,pcmd,'linewidth',2)
axis tight; grid on
ylabel('Depth rate (m/s)','fontweight','bold', 'fontsize',18)
xlabel('12-13 Sep-2013','fontweight','bold', 'fontsize',18)
lg=legend('Depth rate','5sec bin avg.','Commanded');
lg.FontSize=14;
set(gca,'fontweight','bold', 'fontsize',12)
dynamicDateTicks(gca)
%}

%% Sample seed for Xcode 

tgrab = ginput
grabi = find( time>=tgrab(1,1) & time<=tgrab(2,1) );

[ ElevOffset ] = ElevRunningAvg_LRAUV( eleInd(grabi), eleu, 1 );
title([vehicle ' log: ' searchterm],...
    'fontweight','bold','fontsize',18, 'Interpreter', 'none');

matnamefig = matname{:};
print(figure(4),[figd 'ElevRunningAvg_' vehicle '_' matnamefig(1:end-4)],'-r300','-dpng');
% close

grabt = time';  %(grabi)';
grabz = z';     %(grabi)';
grabP = pitch'; %(grabi)';
grabE = ele(grabi); 
grabdz = -dz';
grabSC = speedCmd';

tbl = table(grabt, grabz, grabdz, grabSC, deg2rad(grabP), deg2rad(grabE),...
    'VariableNames',{'Time','Depth','depthRate','speedCmd','pitchAngle','elevAngle'});
writetable(tbl, '~/Desktop/seed/grab.txt', 'Delimiter', ' ')
% dlmwrite('~/Desktop/grab.txt' ,tbl)

%% Xcode import and plot

% from Xcode
% EFC = importdata('~/Desktop/seed/EFC_output.txt');
% indEFC = closest(EFC.data(:,1),time);

% from vehicle
% c=load(matpath{:});
EFC = c.CBIT.empericalClassifierFaultDetected;
indEFC = closest(EFC.time(EFC.value==1),time);


figure;
set(gcf,'Units','normalized','Position',[0 0.3 1 0.5],... % [left bottom width height]
        'PaperPositionMode','auto')
set(gca,'YDir','reverse','layer','top','fontWeight','bold','fontsize',14)
hold on;  grid on; box on;
plot(time,z,'color', [0.301    0.7450    0.9330],'linewidth',0.8);
plot(time,zInd,'color', [0    0.4470    0.7410],'linewidth',3);
plot(maltime,malz,'.','color', [0.8500    0.3250    0.0980],'markersize',10);
plot(timeInd(indEFC),zInd(indEFC),'ro','markersize',8,'markerfacecolor','r');
plot(timeInd(qual),zInd(qual),'o','color', [0.4940    0.1840    0.5560],...
    'markerfacecolor',  [0.4940    0.1840    0.5560], 'markersize',6);
% plot(EFC.data(:,1),EFC.data(:,2),'ro','markersize',6);
ylabel('Depth (m)','fontweight','bold','fontsize',20);
title([vehicle ' log: ' searchterm],...
    'fontweight','bold','fontsize',18, 'Interpreter', 'none');
lg4 = legend('Platform speed < 0.8','Platform speed > 0.8',...
    'Log CRITICAL (\pm5sec)', 'Classified on-board', 'Classified post-process',...
    'location','best');
axis tight;
dynamicDateTicks(gca)


if sv==1
        pslg4 = get(lg4,'position');
        set(lg4,'location','none','position',pslg4)
        matnamefig = matname{:};
        print(gcf,[figd 'EFC_' vehicle '_' matnamefig(1:end-4) ],'-r300','-dpng');
        close; clear pslg4
end
    