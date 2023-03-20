function [Glob_Amp,Low_Indx,High_Indx,High_RBC2Mem,Low_RBC2Mem,Alt_Bin,Alt_RBC2Mem,Alt_Min2Min,status] = newbin_wiggles(disfid,gasfid,RBC2Mem,TR,write_path);

%Sanity Check figure subplot dimensions
sc_x = 3;
sc_y = 2;

%%
[Dis_Shift,~] = Tools.SinglePointDixon_FID(disfid,RBC2Mem,gasfid);

[~,rm_dis] = rmoutliers(abs(Dis_Shift(1,:)),'movmedian',8);

NProj = length(disfid);


Dis_Shift(:,rm_dis) = [];
gasfid(:,rm_dis) = [];
disfid(:,rm_dis) = [];

Dis_Shift = Dis_Shift(1,:);

Time_Axis = (0:(NProj-1))*TR/1000;
Time_Axis(rm_dis) = [];
SampF = 1/(TR/1000);
Freq_Axis = SampF*((-(NProj/2)+1):(NProj/2))/NProj;
%If we scale the dissolved k0 by the gas k0, that pretty well detrends data
myfit = fit(Time_Axis',abs(Dis_Shift'),'exp2');
Dis_Shift = Dis_Shift./(myfit(Time_Axis))';

RBC_Fid = real(Dis_Shift);
Bar_Fid = imag(Dis_Shift);

NProj = length(Dis_Shift);

if(mean(RBC_Fid(1,:))<0)
    RBC_Fid = -RBC_Fid;
end
if(mean(Bar_Fid(1,:))<0)
    Bar_Fid = -Bar_Fid;
end

%Create a smoothing window
Sm_window = floor(1/(TR/1000)/5);

%Isolate k0
RBCDet = RBC_Fid(1,:);
BarDet = Bar_Fid(1,:);

%% First sanity check:
WigBinFig = figure('Name','Sanity_Check')
subplot(sc_x,sc_y,1);
plot(Time_Axis,abs(Dis_Shift));
title('Detrended dissolved')
subplot(sc_x,sc_y,2);
plot(Time_Axis,RBCDet,Time_Axis,BarDet);
title('RBC and Membrane Signals')

%% Smooth and fit:
RBCSm = smooth(RBCDet,Sm_window);
%To adjust for any residual wiggle, fit a series of 8 sines to the wiggle -
%Won't be able to get a good global amplitude assessment, but should pick
%out the peaks pretty well.
RBCfit = fit(Time_Axis',RBCSm-mean(RBCSm),'sin8');

% Get global amplitude only from the time between 2 and 6 seconds.
GlobAmpData = RBCSm;
GlobAmpData(Time_Axis>6) = [];
GlobAmpData(Time_Axis<2) = [];
GlobTime = Time_Axis;
GlobTime(Time_Axis>6) = [];
GlobTime(Time_Axis<2) = [];


RBCAmpfit = fit(GlobTime',GlobAmpData-mean(GlobAmpData),'sin1');
%% Sanity Check 2 - Check Fit!
subplot(sc_x,sc_y-1,2)
plot(Time_Axis,RBCDet,'.k');hold on;plot(Time_Axis,RBCfit(Time_Axis)+mean(RBCSm),'r','LineWidth',2);
plot(GlobTime,RBCAmpfit(GlobTime)+mean(GlobAmpData),'b','LineWidth',2);
Glob_Amp = RBCAmpfit.a1*2/mean(GlobAmpData)*100;
title(['RBC Fit: Global Amp = ' num2str(Glob_Amp,4)]);

%% Old School Binning
[~,maxlocsAll] = findpeaks(RBCfit(Time_Axis),Time_Axis,'MinPeakDistance',0.5);
[~,minlocsAll] = findpeaks(-RBCfit(Time_Axis),Time_Axis,'MinPeakDistance',0.5);

NPts = 30;
High_Indx = [];
for i = 1:length(maxlocsAll)
    idx(i) = find(Time_Axis==maxlocsAll(i));
    tmp_indx = (idx(i)-NPts/2):(idx(i)+NPts/2);
    High_Indx = [High_Indx tmp_indx];
end
h_idx = idx;

Low_Indx = [];
for i = 1:length(minlocsAll)
    idx(i) = find(Time_Axis==minlocsAll(i));
    tmp_indx = (idx(i)-NPts/2):(idx(i)+NPts/2);
    Low_Indx = [Low_Indx tmp_indx];
end
l_idx = idx;

High_Indx(High_Indx<1) = [];
Low_Indx(Low_Indx<1) = [];
High_Indx(High_Indx>NProj) = [];
Low_Indx(Low_Indx>NProj) = [];
%Let's force these to have the same number of points - pretty easy to do in
%this case - pull points from the end since those will be low SNR anyway
if length(High_Indx)>length(Low_Indx)
    High_Indx = High_Indx(1:length(Low_Indx));
elseif length(High_Indx)<length(Low_Indx)
    Low_Indx = Low_Indx(1:length(High_Indx));
end
%We also need the "High" and "Low" RBC/Barrier ratio for separation of
%images later
High_RBC2Mem = mean(RBCDet(High_Indx))/mean(BarDet(High_Indx));
Low_RBC2Mem = mean(RBCDet(Low_Indx))/mean(BarDet(Low_Indx));

subplot(sc_x,sc_y-1,3)
plot(Time_Axis,RBCDet,'.k');
hold on
plot(Time_Axis(High_Indx),RBCDet(High_Indx),'.g','MarkerSize',10);
plot(Time_Axis(Low_Indx),RBCDet(Low_Indx),'.m','MarkerSize',10);
mean_RBC_amp = (mean(RBCDet(High_Indx))-mean(RBCDet(Low_Indx)))/mean(RBCDet)*100;
title(['Binning: Mean Osc Amp = ' num2str(mean_RBC_amp,4)]);
%% Alternative Binning (getting a full cycle)

[Alt_Bin,Alt_Min2Min] = Tools.tenbin_wiggles(RBCDet,Time_Axis,h_idx,l_idx);
for i = 1:size(Alt_Bin,1)
    Alt_RBC2Mem(i) = mean(RBCDet(Alt_Bin(i,:)))/mean(BarDet(Alt_Bin(i,:)));
end
close;
%% Save data for revisiting in the future
save(fullfile(write_path,'NewWiggle_Binning.mat'),'Low_Indx','High_Indx','High_RBC2Mem','Low_RBC2Mem','Alt_Bin','Alt_RBC2Mem','Alt_Min2Min')

saveas(WigBinFig,fullfile(write_path,'WiggleBinningDebug.png'))

good = questdlg('Is binning appropriate?');

switch good
    case 'Yes'
        status = 0;
    case 'No'
        status = 1;
    case 'Cancel'
        status = 1;
    case ''
        status = 1;
end
