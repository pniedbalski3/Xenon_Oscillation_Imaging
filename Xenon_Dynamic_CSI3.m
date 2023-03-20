
%% load data
%[myfile,mypath] = uigetfile;
% 
%dicom_obj = dicom_open(fullfile(mypath,myfile));

raw = load_spectra_stack([8 8]);

twix_obj = mapVBVD;
%% reshape raw so we have pts x read x phase x dynamics
%raw = double(permute(raw,[1 2 4 3]));

%% Get parameters from twix
dwell_time = twix_obj.hdr.MeasYaps.sRXSPEC.alDwellTime{1,1};  % appears most robust and generic
freq = twix_obj.hdr.Dicom.lFrequency;  % seems to work in all sequences
tr = twix_obj.hdr.Config.TR;

dwell_time=2*dwell_time*1E-9; % convert to seconds for calculations

t = double((0:(size(raw,1)-1))*dwell_time');

%% Quick fft 

raw2 = zeros([size(raw,1)*2 size(raw,2,3,4)]);
raw2(1:size(raw,1),:,:,:) = raw;
Images = zeros(size(raw2));

for i = 1:size(raw2,4)
    for j = 1:size(raw2,5)
        Images(:,:,:,i,j) = fftshift(ifftn(ifftshift(raw2(:,:,:,i,j))));
    end
end

%% Pre-allocate memory
spec_pts = size(raw,1);

Gas = zeros(size(raw,2),size(raw,3),size(raw,4));
RBC = zeros(size(raw,2),size(raw,3),size(raw,4));
Mem = zeros(size(raw,2),size(raw,3),size(raw,4));
Mem_Shift = zeros(size(raw,2),size(raw,3),size(raw,4));
RBC_Shift = zeros(size(raw,2),size(raw,3),size(raw,4));
Mem_FWHM = zeros(size(raw,2),size(raw,3),size(raw,4));
RBC_FWHM = zeros(size(raw,2),size(raw,3),size(raw,4));
Mem_Phase = zeros(size(raw,2),size(raw,3),size(raw,4));
RBC_Phase = zeros(size(raw,2),size(raw,3),size(raw,4));
Spectra = zeros(spec_pts,size(raw,2),size(raw,3),size(raw,4));
RBCfit_only = zeros(spec_pts,size(raw,2),size(raw,3),size(raw,4));
Memfit_only = zeros(spec_pts,size(raw,2),size(raw,3),size(raw,4));
Gasfit_only = zeros(spec_pts,size(raw,2),size(raw,3),size(raw,4));
%% Fit everything!
tic
for i = 1:size(raw,4) %Failed at dynamic 18... try starting at dynamic 19 for now. Make sure to switch back.
    for j = 1:size(raw,2)
        for k = 1:size(raw,3)
           % tic
           try
               % Since receiver is at gas frequency, want gas in position 1
               % rather than 3. Hopefully this will give a better fit?
               disfitObj = NMR_TimeFit_v(raw(1:spec_pts,j,k,i),t(1:spec_pts),[100 1 0.5],[0 -6700  -7430],[30 200 250],[0 200 0],[0 0 0],2,length(t)*2); % first widths lorenzian, 2nd are gauss
               disfitObj = disfitObj.fitTimeDomainSignal();
               %disfitObj.plotTimeAndSpectralFit;
               AppendedDissolvedFit = dwell_time*fftshift(fft(disfitObj.calcComponentTimeDomainSignal(t(1:spec_pts)),[],1),1);
               Spectra(:,j,k,i) = disfitObj.spectralDomainSignal;
               RBCfit_only(:,j,k,i) = AppendedDissolvedFit(:,3);
               Memfit_only(:,j,k,i) = AppendedDissolvedFit(:,2);
               Gasfit_only(:,j,k,i) = AppendedDissolvedFit(:,1);
               Gas(j,k,i) = disfitObj.area(1);
               RBC(j,k,i) = disfitObj.area(3);
               Mem(j,k,i) = disfitObj.area(2);
               Mem_Shift(j,k,i) = disfitObj.freq(2);
               RBC_Shift(j,k,i) = disfitObj.freq(3);
               Mem_FWHM(j,k,i) = disfitObj.fwhm(2);
               RBC_FWHM(j,k,i) = disfitObj.fwhm(3);
               Mem_Phase(j,k,i) = disfitObj.phase(2);
               RBC_Phase(j,k,i) = disfitObj.phase(3);
           catch
               disp(['Row ' num2str(j) ', Column ' num2str(k) ', Dynamic ' num2str(i) ' failed']);
           end
           % toc
           
        end
    end
    disp(['Dynamic ' num2str(i) ' of ' num2str(size(raw,4)) ' complete']); 
end
toc

freq = disfitObj.f;

%% Displaying data
Coil = 1;
Rep = 40;
%% Easy Display
figure;
for i = 1:size(raw,2)
    for j = 1:size(raw,3)
        subplot(size(raw,2),size(raw,3),j+size(raw,3)*(i-1)); 
        plot(freq,squeeze(abs(Spectra(:,i,j,Rep,Coil))));
        xlim([-10000 -3000])
        axis off
        %set(gca,'xdir','reverse');
        
    end
end

%% More involved display

DTR = .11;
figure;
for i = 1:size(raw,2)
    for j = 1:size(raw,3)
        subplot(size(raw,2),size(raw,3),j+size(raw,3)*(i-1)); 
        plot(disfitObj.f,squeeze(abs(Spectra(:,i,j,Rep,Coil))));
        %set(gca,'xdir','reverse');
        axis off
    end
end

figure('Name','Gas');
rept = 0:DTR:((size(raw,4)-1)*DTR); 
for i = 1:size(raw,2)
    for j = 1:size(raw,3)
        subplot(size(raw,2),size(raw,3),j+size(raw,3)*(i-1)); 
        plot(rept,squeeze(Gas(i,j,:)));
        %set(gca,'xdir','reverse');
        axis off
    end
end

figure('Name','Mem');
rept = 0:DTR:((size(raw,4)-1)*DTR); 
for i = 1:size(raw,2)
    for j = 1:size(raw,3)
        subplot(size(raw,2),size(raw,3),j+size(raw,3)*(i-1)); 
        plot(rept,squeeze(Mem(i,j,:)));
        %set(gca,'xdir','reverse');
        axis off
    end
end
%%
figure('Name','RBC');
rept = 0:DTR:((size(raw,4)-1)*DTR); 
my_tile = tiledlayout(size(raw,2),size(raw,3));
my_tile.TileSpacing = 'compact';
my_tile.Padding = 'compact';
for i = 1:size(raw,2)
    for j = 1:size(raw,3)
        %subplot(size(raw,2),size(raw,3),j+size(raw,3)*(i-1)); 
        nexttile
        plot(rept(10:end),squeeze(RBC(i,j,10:end)),'g');
        %set(gca,'xdir','reverse');
        axis off
        axis square
        
    end
end
%%
figure('Name','RBC/Gas');
rept = 0:DTR:((size(raw,4)-1)*DTR); 
for i = 1:size(raw,2)
    for j = 1:size(raw,3)
        subplot(size(raw,2),size(raw,3),j+size(raw,3)*(i-1)); 
        plot(rept,squeeze(RBC(i,j,:)./Gas(i,j,:)));
        %set(gca,'xdir','reverse');
        axis off
    end
end

figure('Name','RBC/Mem');
rept = 0:DTR:((size(raw,4)-1)*DTR); 
for i = 1:size(raw,2)
    for j = 1:size(raw,3)
        subplot(size(raw,2),size(raw,3),j+size(raw,3)*(i-1)); 
        plot(rept,squeeze(RBC(i,j,:)./Mem(i,j,:)));
        %set(gca,'xdir','reverse');
        axis off
    end
end
%%
figure('Name','RBC Shift');
rept = 0:DTR:((size(raw,4)-1)*DTR); 
for i = 1:size(raw,2)
    for j = 1:size(raw,3)
        subplot(size(raw,2),size(raw,3),j+size(raw,3)*(i-1)); 
        plot(rept,squeeze(RBC_Shift(i,j,:)));
        %set(gca,'xdir','reverse');
        %axis off
    end
end

%% Fitting, etc
Amp = zeros(size(raw,2),size(raw,3));
Phase = zeros(size(raw,2),size(raw,3));
Rate = zeros(size(raw,2),size(raw,3));
myfig = figure('Name','Check Detrend and Fit');
for i = 1:size(raw,2)
    for j = 1:size(raw,3)
        [Amp(i,j),Phase(i,j),Rate(i,j)] = detrend_fit_wiggles_2(rept(10:(end-25))',squeeze(RBC(i,j,(10:(end-25)))),squeeze(Gas(i,j,(10:(end-25)))),myfig);
    end
end