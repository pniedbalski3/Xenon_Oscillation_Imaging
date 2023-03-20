clear;clc;close all;

%%
folders = {'P:\IRB_STUDY00146119_Xe_Imaging\20221028_Xe-047',...
'P:\IRB_STUDY00146119_Xe_Imaging\20221115_Xe-049',...
'P:\IRB_STUDY00146119_Xe_Imaging\20221018_Xe-045',...
'P:\IRB_STUDY00146616_COVID_Xe_MRI\20221026_CXe-034_02',...
'P:\IRB_STUDY00148159_MRI-SSc-ILD\QPI01_006',...
'P:\IRB_STUDY00148159_MRI-SSc-ILD\20230112_QPI01-001_V2',...
'P:\IRB_STUDY00146119_Xe_Imaging\20230109_Xe-054',...
'P:\IRB_STUDY00146119_Xe_Imaging\20230111_Xe-052',...
'P:\IRB_STUDY00146119_Xe_Imaging\20230112_Xe-051',...
'P:\IRB_STUDY00146119_Xe_Imaging\20230118_Xe-048',...
'P:\IRB_STUDY00146119_Xe_Imaging\20230206_Xe-055',...
'P:\IRB_STUDY00146119_Xe_Imaging\20230216_Xe-056',...
'P:\IRB_STUDY00146119_Xe_Imaging\20230308_Xe-057',...
};

folders = {'P:\IRB_STUDY00148159_MRI-SSc-ILD\20221118_QPI01-006'};

search_name = 'fid_xe_calibration';

Gas = zeros(140,length(folders));
RBC = zeros(140,length(folders));
Mem = zeros(140,length(folders));
RBC_Shift = zeros(140,length(folders));
Mem_Shift = zeros(140,length(folders));
figure('Name','RBC Wiggles')

Amp = zeros(length(folders),1);
Phase = zeros(length(folders),1);
Rate = zeros(length(folders),1);
Subject{length(folders),1} = ' ';

%%

headers = {'Subject','RBC_Osc','RBCMem_Osc'};
AllWiggles = cell2table(cell(0,size(headers,2)));
AllWiggles.Properties.VariableNames = headers;

for i = 1:length(folders)
    tmp_fold = folders{i};
    if strcmp(tmp_fold((end-1):end),'02')
        Subject{i} = tmp_fold((end-9):end);
    else
        Subject{i} = tmp_fold((end-6):end);
        if strcmp(Subject{i}(1),'_')
            Subject{i} = tmp_fold((end-5):end);
        end
    end
    files = struct2cell(dir(fullfile(tmp_fold,'Raw')));
    names = files(1,:);
    
    try
        myfile = names{find(contains(names,search_name),1,'last')};

        twix = mapVBVD(fullfile(tmp_fold,'Raw',myfile));
        dwell_time = twix.hdr.MeasYaps.sRXSPEC.alDwellTime{1,1};
        dwell_time=dwell_time*1E-9; 
        nFids = twix.hdr.Config.NRepMeas;
        nPts = twix.hdr.Config.VectorSize;
        nDis = twix.hdr.Phoenix.sWipMemBlock.alFree{3};

        theFID = squeeze(double(twix.image()));
        te = twix.hdr.Phoenix.alTE{1}; 

        disData = theFID(:,1:nDis);

        %Kill the first and last 100 points (gets us 4.5s of data and avoids descent to steady state and low snr at end)?)
        disData(:,1:100) = [];
        disData(:,(end-100):end) = [];
        t = double((0:(length(disData)-1))*dwell_time');

        for j = 1:size(disData,2)
            disfitObj = NMR_TimeFit_v(disData(:,j),t,[1 1 1],[0 -700  -7400],[250 200 30],[0 200 0],[0 0 0],2,length(t)*2); % first widths lorenzian, 2nd are gauss
            disfitObj = disfitObj.fitTimeDomainSignal();
            %disfitObj.plotTimeAndSpectralFit;
            Gas(j,i) = disfitObj.area(3);
            RBC(j,i) = disfitObj.area(1);
            Mem(j,i) = disfitObj.area(2);
            Mem_Shift(j,i) = disfitObj.freq(2);
            RBC_Shift(j,i) = disfitObj.freq(1);
        end
        save(fullfile(tmp_fold,'Spectra_fittings.mat'),'Gas','RBC','Mem','Mem_Shift','RBC_Shift');
        
        xdata = 0:0.015:(0.015*(size(RBC,1)-1));
        myfit = fit(xdata',smooth(RBC(:,i)),'exp2');

        det_RBC = RBC(:,i)./myfit(xdata);
        det_RBC = det_RBC-1;

        sinfit = fit(xdata',det_RBC,'sin1');

        Amp(i) = 2*sinfit.a1*100;
        Phase(i) = sinfit.c1;
        Rate = sinfit.b1/2/pi*60;
        
        myfig = figure('Name','Spectra Fitting');
        %figure('Name','How Much Data');
        %plot(abs(disData(1,:)));
        subplot(1,2,1);plot(xdata,det_RBC,xdata,sinfit(xdata));
        title(['Osc Amp = ',num2str(Amp(i),3)]);
        
        R2M = RBC(:,i)./Mem(:,i);
        
        R2M = R2M - mean(R2M);
        sinfit2 = fit(xdata',smooth(R2M),'sin1');
        R2MAmp(i) = 2*sinfit2.a1;
        subplot(1,2,2)
        plot(xdata,R2M,xdata,sinfit2(xdata));
        title(['R2M Osc Amp = ',num2str(R2MAmp(i),3)]);
        
        saveas(myfig,fullfile(tmp_fold,'Wiggle_Spectra_Figs.png'))
        
        newline = {Subject{i},Amp(i),R2MAmp(i)};
        AllWiggles = [AllWiggles;newline];
    catch
        Amp(i) = nan;
        Phase(i) = nan;
        Rate(i) = nan;
    end
end

%%
% figure;
% for i = 1:length(folders)
%     subplot(1,3,i);plot(smooth(RBC(:,i)));
% end
% 


