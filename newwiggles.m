function [status,mean_amps] = newwiggles(Dis_Fid,Gas_Fid,Dis_Traj,Gas_Traj,H1_Image_Dis,LoRes_Gas_Image,Proton_Mask,VentBinMask,RBC_Mask,RBC2Mem,TR,ImSize,scanDateStr,write_path)

sc_x = 2;
sc_y = 2;

%% Prep work
Proton_Mask = logical(Proton_Mask);
VentBinMask = logical(VentBinMask);
RBC_Mask = logical(RBC_Mask);

parent_path = which('newwiggles');
idcs = strfind(parent_path,filesep);%determine location of file separators
parent_path = parent_path(1:idcs(end-1)-1);%remove file

if ~isfolder(fullfile(write_path,'Wiggle_figs'))
    mkdir(fullfile(write_path,'Wiggle_figs'));
end

if(exist(fullfile(parent_path,'AncillaryFiles','HealthyCohort.mat'),'file') == 2) %if Healthy cohort data exists, use
    HealthyFile = dir(fullfile(parent_path,'AncillaryFiles','HealthyCohort.mat'));%get file info
    HealthyData = load(fullfile(HealthyFile.folder,HealthyFile.name), '-regexp','^(?!Comb|Ind)...');%import all variable but figures
    %RBCOscThresh = HealthyData.thresholds.RBCOsc;
    %First couple of times, healthy distribution won't be around
    HealthyDistPresent = 0;
else 
    RBCOscThresh = [8.9596-2*10.5608,8.9596-1*10.5608,8.9596,8.9596+1*10.5608,8.9596+2*10.5608,8.9596+3*10.5608,8.9596+4*10.5608];%8
    HealthyDistPresent = 0;
end
%For now
RBCOscThresh = [8.9596-2*10.5608,8.9596-1*10.5608,8.9596,8.9596+1*10.5608,8.9596+2*10.5608,8.9596+3*10.5608,8.9596+4*10.5608];%8
RBCOscThresh2 = [17.7060   21.6670   25.6280   29.5890   33.5500   37.5109   41.4719];%8

%% Rotate images

Proton_Mask = Tools.canonical2matlab(Proton_Mask);
RBC_Mask = Tools.canonical2matlab(RBC_Mask);

LoRes_Gas_Image = Tools.canonical2matlab(LoRes_Gas_Image);
H1_Image_Dis = Tools.canonical2matlab(H1_Image_Dis);

%% Binning:
[Glob_Amp,Low_Indx,High_Indx,High_RBC2Mem,Low_RBC2Mem,Alt_Bin,Alt_RBC2Mem,Alt_Min2Min,status] = newbin_wiggles(Dis_Fid,Gas_Fid,RBC2Mem,TR,write_path);
if status ~= 0
    return
end

%% Detrend Dissolved Phase data
DisDet = Dis_Fid;

disk0 = abs(Dis_Fid(1,:));
mean_dis = mean(disk0);
t = 0:TR:((length(disk0)-1)*TR);
myfit = fit(t',disk0','exp2');

for i = 1:size(Dis_Fid,2)
    DisDet(:,i) = Dis_Fid(:,i)/myfit(t(i));
end

%DisDet is now centered around 1. I want to preserve the intensity so as to
%be able to report oscillation in units of RBC/Gas (don't know if that will
%be useful, but want to have a look)

DisDet = DisDet*mean_dis;
figure('Name','Sanity Checks')

subplot(sc_x,sc_y,1)
plot(t,disk0,t,myfit(t));
title('Detrending Data')
subplot(sc_x,sc_y,2)
plot(t,abs(DisDet(1,:)))
title('Check that Detrending was successful')

%% Skip the keyhole and just reconstruct
HighOnly = DisDet(:,High_Indx);
HighOnly_Traj = Dis_Traj(:,:,High_Indx);
LowOnly = DisDet(:,Low_Indx);
LowOnly_Traj = Dis_Traj(:,:,Low_Indx);

%% Keyholing
%Hardcode to key radius of 9 per Junlan's paper
Key_Rad = 9;

%Empty Keyhole
Keyhole = DisDet;
Keyhole(1:Key_Rad,:) = NaN;

High_Key = Keyhole;
Low_Key = Keyhole;

KeyholeG = Gas_Fid;
KeyholeG(1:Key_Rad,:) = NaN;

High_KeyG = KeyholeG;
Low_KeyG = KeyholeG;

%Scale keyhole to match intensity of key
mean_high = mean(abs(DisDet(1,High_Indx)));
mean_low = mean(abs(DisDet(1,Low_Indx)));

mean_highG = mean(abs(Gas_Fid(1,High_Indx)));
mean_lowG = mean(abs(Gas_Fid(1,Low_Indx)));

for i = 1:length(Keyhole)
    High_Key(:,i) = High_Key(:,i)*mean_high/abs(DisDet(1,i));
    Low_Key(:,i) = Low_Key(:,i)*mean_low/abs(DisDet(1,i));
    
    High_KeyG(:,i) = High_KeyG(:,i)*mean_highG/abs(Gas_Fid(1,i));
    Low_KeyG(:,i) = Low_KeyG(:,i)*mean_lowG/abs(Gas_Fid(1,i));
end

%Put high and low keys in keyhole
High_Key(1:Key_Rad,High_Indx) = DisDet(1:Key_Rad,High_Indx);
Low_Key(1:Key_Rad,Low_Indx) = DisDet(1:Key_Rad,Low_Indx);

High_KeyG(1:Key_Rad,High_Indx) = Gas_Fid(1:Key_Rad,High_Indx);
Low_KeyG(1:Key_Rad,Low_Indx) = Gas_Fid(1:Key_Rad,Low_Indx);

if ~isfolder(fullfile(write_path,'cs_raw'))
    mkdir(fullfile(write_path,'cs_raw'))
end
save_for_cs(High_Key,Dis_Traj,fullfile(write_path,'cs_raw','High_Key_Raw.mat'));
save_for_cs(Low_Key,Dis_Traj,fullfile(write_path,'cs_raw','Low_Key_Raw.mat'));
save_for_cs(Dis_Fid,Dis_Traj,fullfile(write_path,'cs_raw','Dis_Raw.mat'));
save_for_cs(Gas_Fid,Gas_Traj,fullfile(write_path,'cs_raw','Gas_Raw.mat'));

%Sanity Checks
subplot(sc_x,sc_y,3)
imagesc(abs(High_Key))
title('High Key')
subplot(sc_x,sc_y,4)
imagesc(abs(Low_Key))
title('Low Key')

%% Keyhole Alternative Binning

Key_Rad = 9;%radpts(Key_Rad_Pts)/(size(disfid,1)/2)*0.5;

Keyhole = DisDet;
Keyhole(1:Key_Rad,:) = NaN;

%I should scale the keyhole!
Alt_Keys = zeros([size(DisDet) size(Alt_Bin,1)]);
for i = 1:size(Alt_Bin,1)
    Alt_Keys(:,:,i) = Keyhole;
    mean_key = mean(abs(DisDet(1,Alt_Bin(i,:))));
    for j = 1:length(Keyhole)
        Alt_Keys(:,j,i) = Alt_Keys(:,j,i)*mean_key/abs(DisDet(1,j));
    end
    Alt_Keys(1:Key_Rad,Alt_Bin(i,:),i) = DisDet(1:Key_Rad,Alt_Bin(i,:));
    
    save_for_cs(squeeze(Alt_Keys(:,:,i)),Dis_Traj,fullfile(write_path,'cs_raw',['Alt_bin_key' num2str(i) '.mat']));
end

%% Now, need to reconstruct High, Low, and Unaltered Data
trajx = reshape(Dis_Traj(1,:,:),1,[])';
trajy = reshape(Dis_Traj(2,:,:),1,[])';
trajz = reshape(Dis_Traj(3,:,:),1,[])';

traj_r = [trajx trajy trajz];

gastrajx = reshape(Gas_Traj(1,:,:),1,[])';
gastrajy = reshape(Gas_Traj(2,:,:),1,[])';
gastrajz = reshape(Gas_Traj(3,:,:),1,[])';

gastraj_r = [gastrajx gastrajy gastrajz]*1.5;
%orientation issue!
gastraj_hold = gastraj_r;
gastraj_r(:,1) = gastraj_hold(:,1);
gastraj_r(:,2) = gastraj_hold(:,3);
gastraj_r(:,3) = gastraj_hold(:,2);

%Reshape to column vectors
dis_data_r = reshape(DisDet,1,[])';
high_data_r = reshape(High_Key,1,[])';
low_data_r = reshape(Low_Key,1,[])';
gas_data_r = reshape(Gas_Fid,1,[])';

high_gas_r = reshape(High_KeyG,1,[])';
low_gas_r = reshape(Low_KeyG,1,[])';

%Remove NaNs
all_nan = isnan(dis_data_r);
high_nan = isnan(high_data_r);
low_nan = isnan(low_data_r);
dis_data_r(all_nan) = [];
high_data_r(high_nan) = [];
low_data_r(low_nan) = [];
high_traj_r = traj_r;
high_traj_r(high_nan,:) = [];
low_traj_r = traj_r;
low_traj_r(low_nan,:) = [];
dis_traj_r = traj_r;
dis_traj_r(all_nan,:) = [];

high_nanG = isnan(high_gas_r);
low_nanG = isnan(low_gas_r);

high_gas_r(high_nanG) = [];
low_gas_r(low_nanG) = [];
high_gtraj_r = gastraj_r;
high_gtraj_r(high_nanG,:) = [];
rad = sqrt(high_gtraj_r(:,1).^2+high_gtraj_r(:,2).^2+high_gtraj_r(:,3).^2);
toobig = find(rad>0.5);
high_gas_r(toobig) = [];
high_gtraj_r(toobig,:) = [];
low_gtraj_r = gastraj_r;
low_gtraj_r(low_nanG,:) = [];
rad = sqrt(low_gtraj_r(:,1).^2+low_gtraj_r(:,2).^2+low_gtraj_r(:,3).^2);
toobig = find(rad>0.5);
low_gas_r(toobig) = [];
low_gtraj_r(toobig,:) = [];
%Reconstruct
All_Dis = Reconstruction.Dissolved_Phase_LowResRecon(ImSize,dis_data_r,dis_traj_r);
High_Dis = Reconstruction.Dissolved_Phase_LowResRecon(ImSize,high_data_r,high_traj_r);
Low_Dis = Reconstruction.Dissolved_Phase_LowResRecon(ImSize,low_data_r,low_traj_r);
High_Gas = Reconstruction.Dissolved_Phase_LowResRecon(ImSize,high_gas_r,high_gtraj_r);
Low_Gas = Reconstruction.Dissolved_Phase_LowResRecon(ImSize,low_gas_r,low_gtraj_r);

HighOnly_Dis = Reconstruction.Dissolved_Phase_LowResRecon(ImSize,reshape(HighOnly,1,[])',[reshape(HighOnly_Traj(1,:,:),1,[])',reshape(HighOnly_Traj(2,:,:),1,[])',reshape(HighOnly_Traj(3,:,:),1,[])']);
LowOnly_Dis = Reconstruction.Dissolved_Phase_LowResRecon(ImSize,reshape(LowOnly,1,[])',[reshape(LowOnly_Traj(1,:,:),1,[])',reshape(LowOnly_Traj(2,:,:),1,[])',reshape(LowOnly_Traj(3,:,:),1,[])']);


High_Gas_Scaled = High_Gas*sind(14.6)/sind(1);
Low_Gas_Scaled = Low_Gas*sind(14.6)/sind(1);
%% What about alternative binning:

Alt_Images = zeros(ImSize,ImSize,ImSize,size(Alt_Bin,1));
for i = 1:size(Alt_Bin,1)
    thiskey = squeeze(Alt_Keys(:,:,i));
    thistraj = traj_r;
    thiskey_r = reshape(thiskey,1,[])';
    thisnan = isnan(thiskey_r);
    thiskey_r(thisnan) = [];
    thistraj(thisnan,:) = [];
    Alt_Images(:,:,:,i) = Reconstruction.Dissolved_Phase_LowResRecon(ImSize,thiskey_r,thistraj);
end

%% Separate RBC and Barrier Images
[All_Mem,All_RBC, ~, ~] = Tools.SinglePointDixon_V2(All_Dis,-RBC2Mem,LoRes_Gas_Image,Proton_Mask);
[High_Mem,High_RBC, ~, ~] = Tools.SinglePointDixon_V2(High_Dis,-High_RBC2Mem,LoRes_Gas_Image,Proton_Mask);
[Low_Mem,Low_RBC, ~, ~] = Tools.SinglePointDixon_V2(Low_Dis,-Low_RBC2Mem,LoRes_Gas_Image,Proton_Mask);

[HighOnly_Mem,HighOnly_RBC, ~, ~] = Tools.SinglePointDixon_V2(High_Dis,-High_RBC2Mem,LoRes_Gas_Image,Proton_Mask);
[LowOnly_Mem,LowOnly_RBC, ~, ~] = Tools.SinglePointDixon_V2(Low_Dis,-Low_RBC2Mem,LoRes_Gas_Image,Proton_Mask);

Alt_RBC = zeros(ImSize,ImSize,ImSize,size(Alt_Bin,1));
Alt_Bar = zeros(ImSize,ImSize,ImSize,size(Alt_Bin,1));
for i = 1:size(Alt_Bin,1)
    [Alt_Bar(:,:,:,i),Alt_RBC(:,:,:,i),~,~] = Tools.SinglePointDixon_V2(squeeze(Alt_Images(:,:,:,i)),-Alt_RBC2Mem(i),LoRes_Gas_Image,Proton_Mask);
end
Alt_RBC = abs(Alt_RBC);
All_RBC = abs(All_RBC);
% if ~MaskPres
%     %Just get RBC points with SNR > 2.5
%     [~,RBC_Mask] = Tools.erode_dilate(All_RBC,1,2.5);
%     RBC_Mask = logical(RBC_Mask.*Vent_Mask);
% end
sc_x = 4;sc_y = 3;
figure('Name','Images Sanity Check');
[centerslice,~,~] = AllinOne_Tools.getimcenter(Proton_Mask);
subplot(sc_x,sc_y,1)
imagesc(squeeze(abs(LoRes_Gas_Image(:,:,centerslice))));axis square;axis off;colormap(gray);title('All Gas')
subplot(sc_x,sc_y,2)
imagesc(squeeze(abs(High_Gas(:,:,centerslice))));axis square;axis off;colormap(gray);title('High Gas')
subplot(sc_x,sc_y,3)
imagesc(squeeze(abs(Low_Gas(:,:,centerslice))));axis square;axis off;colormap(gray);title('Low Gas')
subplot(sc_x,sc_y,4)
imagesc(squeeze(abs(All_Dis(:,:,centerslice))));axis square;axis off;colormap(gray);title('All Dissolved')
subplot(sc_x,sc_y,5)
imagesc(squeeze(abs(High_Dis(:,:,centerslice))));axis square;axis off;colormap(gray);title('High Dissolved')
subplot(sc_x,sc_y,6)
imagesc(squeeze(abs(Low_Dis(:,:,centerslice))));axis square;axis off;colormap(gray);title('Low Dissolved')
subplot(sc_x,sc_y,7)
imagesc(squeeze(abs(All_Mem(:,:,centerslice))));axis square;axis off;colormap(gray);title('All Membrane')
subplot(sc_x,sc_y,8)
imagesc(squeeze(abs(High_Mem(:,:,centerslice))));axis square;axis off;colormap(gray);title('High Membrane')
subplot(sc_x,sc_y,9)
imagesc(squeeze(abs(Low_Mem(:,:,centerslice))));axis square;axis off;colormap(gray);title('Low Membrane')
subplot(sc_x,sc_y,10)
imagesc(squeeze(abs(All_RBC(:,:,centerslice))));axis square;axis off;colormap(gray);title('All RBC')
subplot(sc_x,sc_y,11)
imagesc(squeeze(abs(High_RBC(:,:,centerslice))));axis square;axis off;colormap(gray);title('High RBC')
subplot(sc_x,sc_y,12)
imagesc(squeeze(abs(Low_RBC(:,:,centerslice))));axis square;axis off;colormap(gray);title('Low RBC')
%% Test to see if the many key approach even remotely works:
mean_RBC = zeros(1,size(Alt_Bin,1));
mean_Bar = zeros(1,size(Alt_Bin,1));
for i = 1:size(Alt_Bin,1)
    thisRBC = squeeze(Alt_RBC(:,:,:,i));
    thisBar = squeeze(Alt_Bar(:,:,:,i));
    mean_RBC(i) = mean(thisRBC(Proton_Mask(:)));
    mean_Bar(i) = mean(thisBar(Proton_Mask(:)));
end

%% Scale images such that 1 is the overall mean RBC signal
allmean = mean(mean_RBC);
Alt_RBC = Alt_RBC/allmean;
All_RBC2 = All_RBC/allmean;

%% Now, loop through points in mask and get maximum and minimum. Also, get overall phase shift
Phase_Step = 360 / size(Alt_Bin,1);

%Let's find the "phase" of the overall mean RBC:
if Alt_Min2Min
    [~,extr_ind] = max(mean_RBC);
    mean_Phase = extr_ind*Phase_Step;
else
    [~,extr_ind] = min(mean_RBC);
    mean_Phase = extr_ind*Phase_Step;
end

%Amp = zeros(ImSize,ImSize,ImSize);
%imslice(abs(Alt_RBC));
Amp = (max(Alt_RBC,[],4) - min(Alt_RBC,[],4))/mean(All_RBC2(RBC_Mask==1));

if Alt_Min2Min 
    [~,extr_ind] = max(Alt_RBC,[],4);
else
    [~,extr_ind] = min(Alt_RBC,[],4);
end
Phase = mean_Phase* ones(size(RBC_Mask)) - Phase_Step*extr_ind;
Phase(RBC_Mask==0) = -361;
Amp(RBC_Mask==0) = -361;

Amp = Amp*100;

%% Oscillations from high and low keys

% Regional oscillation amplitude scaled by the mean RBC
RBC_Osc_bymean = (High_RBC - Low_RBC)./mean(All_RBC(RBC_Mask(:)))*100.*RBC_Mask;
% Regional oscillation amplitude scaled voxel-by-voxel by RBC
RBC_Osc_byvoxel = (High_RBC - Low_RBC)./All_RBC*100.*RBC_Mask;
% Regional oscillation amplitude using "absolute" RBC/Gas units
RBC_Osc_bygas = (High_RBC./abs(High_Gas_Scaled) - Low_RBC./abs(Low_Gas_Scaled)).*RBC_Mask*100;
% Regional oscillations scaled voxel-by-voxel by a blurred RBC image
RBC_Osc_bysmvoxel = (High_RBC - Low_RBC)./imgaussfilt3(All_RBC, 2)*100.*RBC_Mask;
% Regional oscillation amplitude using "absolute" RBC/Membrane units
RBC_Osc_bymem = (High_RBC./High_Mem - Low_RBC./Low_Mem).*RBC_Mask;

RBC_Osc_bymean(RBC_Mask~=1) = nan;
RBC_Osc_byvoxel(RBC_Mask~=1) = nan;
RBC_Osc_bygas(RBC_Mask~=1) = nan;
RBC_Osc_bysmvoxel(RBC_Mask~=1) = nan;
RBC_Osc_bymem(RBC_Mask~=1) = nan;

mean_amps.global = Glob_Amp;
mean_amps.bymean = mean(RBC_Osc_bymean(RBC_Mask==1));
mean_amps.byvoxel = mean(RBC_Osc_byvoxel(RBC_Mask==1));
mean_amps.bysmvoxel = mean(RBC_Osc_bysmvoxel(RBC_Mask==1));
mean_amps.bygas = mean(RBC_Osc_bygas(RBC_Mask==1));
mean_amps.bymem = mean(RBC_Osc_bymem(RBC_Mask==1));
mean_amps.newway = mean(Amp(RBC_Mask==1));

mean_amps.bymeanstd = std(RBC_Osc_bymean(RBC_Mask==1));
mean_amps.byvoxelstd = std(RBC_Osc_byvoxel(RBC_Mask==1));
mean_amps.bysmvoxelstd = std(RBC_Osc_bysmvoxel(RBC_Mask==1));
mean_amps.bygasstd = std(RBC_Osc_bygas(RBC_Mask==1));
mean_amps.bymemstd = std(RBC_Osc_bymem(RBC_Mask==1));
mean_amps.newwaystd = std(Amp(RBC_Mask==1));

%% Displays
bymean_thresh = [-0.780046424797197,1.968126803008373,5.718859369280803,10.611222901322080,16.780770854075527,24.360064104267494,33.479067048585776];
byvoxel_thresh = [-1.495168592287216,1.917144117304416,6.317796381463472,11.898545644108060,18.872482232086725,27.474859854495627,37.963905525139560];
bysmvoxel_thresh = [-1.815136634737343,2.131440574369190,6.714843996672521,12.012178137002987,18.106981712062048,25.089537935849638,33.057189516365880];
bygas_thresh = [-0.001016994342597,0.006472889493039,0.016205972214379,0.028950956224809,0.045775398806456,0.068176637522764,0.098276994987147];
bymem_thresh = [-7.266246332076620e-04,0.019123803643213,0.041495536308714,0.066579610047285,0.094572584889091,0.125676466379108,0.160098631924971];
newway_thresh = [4.924758605973880,7.157254425571484,10.217395767848377,14.428884694139139,20.248820019481293,28.325552712779235,39.582921954070656];

bymean_edges = [-100 linspace(-10, bymean_thresh(7) + 10,100) 100];
byvoxel_edges = [-100 linspace(-10, byvoxel_thresh(7) + 10,100) 100];
bysmvoxel_edges = [-100 linspace(-10, bysmvoxel_thresh(7) + 10,100) 100];
bygas_edges = [-100 linspace(-0.01, bygas_thresh(7) + 0.02,100) 100];
bymem_edges = [-100 linspace(-0.01, bymem_thresh(7) + 0.04,100) 100];
newway_edges = [-0.5 linspace(0, newway_thresh(7) + 20,100) 100];

EightBinMap = [1 0 0; 1 0.7143 0; 0.4 0.7 0.4; 0 1 0; 184/255 226/255 145/255; 243/255 205/255 213/255; 225/255 129/255 162/255; 197/255 27/255 125/255]; %Used for barrier

osc_hist_fig = figure('Name','Histograms of Different Oscillation Amplitude Methods');
subplot(2,3,1)
AllinOne_Tools.CalculateDissolvedHistogram_subplot(RBC_Osc_bymean(RBC_Mask==1),bymean_edges,bymean_thresh,EightBinMap,[],[]);
%histogram(RBC_Osc_bymean(RBC_Mask==1));
title('High-Low/mean')

subplot(2,3,2)
AllinOne_Tools.CalculateDissolvedHistogram_subplot(RBC_Osc_byvoxel(RBC_Mask==1),byvoxel_edges,byvoxel_thresh,EightBinMap,[],[]);
%histogram(RBC_Osc_byvoxel(RBC_Mask==1));
title('High-Low/voxelwise')

subplot(2,3,3)
AllinOne_Tools.CalculateDissolvedHistogram_subplot(RBC_Osc_bysmvoxel(RBC_Mask==1),bysmvoxel_edges,bysmvoxel_thresh,EightBinMap,[],[]);
%histogram(RBC_Osc_bysmvoxel(RBC_Mask==1));
title('High-Low/smoothed voxelwise')

subplot(2,3,4)
AllinOne_Tools.CalculateDissolvedHistogram_subplot(RBC_Osc_bygas(RBC_Mask==1),bygas_edges,bygas_thresh,EightBinMap,[],[]);
%histogram(RBC_Osc_bygas(RBC_Mask==1));
title('High-Low/gas')

subplot(2,3,5)
AllinOne_Tools.CalculateDissolvedHistogram_subplot(RBC_Osc_bymem(RBC_Mask==1),bymem_edges,bymem_thresh,EightBinMap,[],[]);
%histogram(RBC_Osc_bymem(RBC_Mask==1));
title('High-Low/membrane')

subplot(2,3,6)
AllinOne_Tools.CalculateDissolvedHistogram_subplot(Amp(RBC_Mask==1),newway_edges,newway_thresh,EightBinMap,[],[]);
%histogram(Amp(RBC_Mask==1));
title('Alternative Binning')
set(osc_hist_fig,'Position',[281 100 1558 847])
saveas(osc_hist_fig,fullfile(write_path,'Oscillation_Histograms.png'))

abs_osc_fig = figure('Name','Slices of Different Oscillation Amplitude Methods');
subplot(2,3,1)
imagesc(squeeze(RBC_Osc_bymean(:,:,centerslice)));axis square;axis off;caxis([-5 20]);
title('High-Low/mean')

subplot(2,3,2)
imagesc(squeeze(RBC_Osc_byvoxel(:,:,centerslice)));axis square;axis off;caxis([-5 20]);
title('High-Low/voxelwise')

subplot(2,3,3)
imagesc(squeeze(RBC_Osc_bysmvoxel(:,:,centerslice)));axis square;axis off;caxis([-5 20]);
title('High-Low/smoothed voxelwise')

subplot(2,3,4)
imagesc(squeeze(RBC_Osc_bygas(:,:,centerslice)));axis square;axis off;caxis([0 0.12]);
title('High-Low/gas')

subplot(2,3,5)
imagesc(squeeze(RBC_Osc_bymem(:,:,centerslice)));axis square;axis off;caxis([0 0.2]);
title('High-Low/membrane')

subplot(2,3,6)
imagesc(squeeze(Amp(:,:,centerslice)));axis square;axis off;caxis([0 30]);
title('Alternative Binning')

%% Binning: Hardcode thresholds for now:


bymean_bin = AllinOne_Tools.BinImages(RBC_Osc_bymean,bymean_thresh).*RBC_Mask;
byvoxel_bin = AllinOne_Tools.BinImages(RBC_Osc_byvoxel,byvoxel_thresh).*RBC_Mask;
bysmvoxel_bin = AllinOne_Tools.BinImages(RBC_Osc_bysmvoxel,bysmvoxel_thresh).*RBC_Mask;
bygas_bin = AllinOne_Tools.BinImages(RBC_Osc_bygas,bygas_thresh).*RBC_Mask;
bymem_bin = AllinOne_Tools.BinImages(RBC_Osc_bymem,bymem_thresh).*RBC_Mask;
Amp_bin = AllinOne_Tools.BinImages(Amp,newway_thresh).*RBC_Mask;

EightBinMap = [0 0 0; 1 0 0; 1 0.7143 0; 0.4 0.7 0.4; 0 1 0; 184/255 226/255 145/255; 243/255 205/255 213/255; 225/255 129/255 162/255; 197/255 27/255 125/255]; %Used for barrier

binned_osc_fig = figure('Name','Slices of Different Oscillation Amplitude Methods Binned');
subplot(2,3,1)
imagesc(squeeze(bymean_bin(:,:,centerslice)));axis square;axis off;colormap(EightBinMap);caxis([0 8]);
title('High-Low/mean')

subplot(2,3,2)
imagesc(squeeze(byvoxel_bin(:,:,centerslice)));axis square;axis off;colormap(EightBinMap);caxis([0 8]);
title('High-Low/voxelwise')

subplot(2,3,3)
imagesc(squeeze(bysmvoxel_bin(:,:,centerslice)));axis square;axis off;colormap(EightBinMap);caxis([0 8]);
title('High-Low/smoothed voxelwise')

subplot(2,3,4)
imagesc(squeeze(bygas_bin(:,:,centerslice)));axis square;axis off;colormap(EightBinMap);caxis([0 8]);
title('High-Low/gas')

subplot(2,3,5)
imagesc(squeeze(bymem_bin(:,:,centerslice)));axis square;axis off;colormap(EightBinMap);caxis([0 8]);
title('High-Low/membrane')

subplot(2,3,6)
imagesc(squeeze(Amp_bin(:,:,centerslice)));axis square;axis off;colormap(EightBinMap);caxis([0 8]);
title('Alternative Binning')

%% Saving Results
save(fullfile(write_path,'Oscillation_Amps_Multiple_Ways.mat'),'RBC_Osc_bymean','RBC_Osc_byvoxel','RBC_Osc_bygas','RBC_Osc_bysmvoxel','RBC_Osc_bymem','Amp','Glob_Amp','bymean_bin','byvoxel_bin','bysmvoxel_bin','bygas_bin','bymem_bin','Amp_bin')
