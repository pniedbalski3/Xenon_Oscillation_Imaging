function [Amp,Phase,Rate] = detrend_fit_wiggles_2(time_axis,RBC,Gas,fighandle)
%% Detrending
%decay_fit = fit(time_axis,Membrane,'exp2');

%decay_fit = fit(time_axis,Membrane,'exp2');
%Det_RBC = RBC./Gas;
if nargin < 4
    myfig = figure('Name','Check Detrend and Fit');
else
    myfig = fighandle;
end


%Getting some outliers... kill those here.
[new_RBC,rmpts] = rmoutliers(RBC);
time_axis(rmpts) = [];

decay_fit = fit(time_axis,new_RBC,'exp2');
Det_RBC = new_RBC./(decay_fit(time_axis));

subplot(2,3,1)
plot(time_axis,new_RBC,'-*r')
title('Raw RBC')
hold on 
plot(time_axis,decay_fit(time_axis),'k');
hold off

subplot(2,3,2)
plot(time_axis,Det_RBC,'-*r')
title('Detrended RBC');

Det_RBC = smooth(Det_RBC/mean(Det_RBC),3);
subplot(2,3,3)
plot(time_axis,Det_RBC,'-*r')
title('Fully Detrended RBC')

%% BandPass Filter
TR = time_axis(2)-time_axis(1);
SampF = 1/(TR);
% with removing outliers, that causes issues... let's see if just a simple
% smoothing is good enough.
% try
%     BP_data = bandpass((Det_RBC-1),[0.4 2.0],SampF)+1;
% catch
%     BP_data = Det_RBC;%bandpass((Det_RBC-1),[0.4 2.5],SampF);
% end
%BP_data = lowpass((Det_RBC-1),2.5,SampF);

BP_data = Det_RBC;
subplot(2,3,4)
plot(time_axis,Det_RBC,'r')
hold on
%plot(time_axis,BP_data+1,'b')
plot(time_axis,BP_data,'b')
title('BandPass RBC')
hold off

%% Sine fit
try
    [sin_fit,gof] = fit(time_axis,BP_data-1,'sin1');
    subplot(2,3,5)
    plot(time_axis,Det_RBC,'r')
    hold on
    plot(time_axis,BP_data,'b')
    plot(time_axis,sin_fit(time_axis)+1,'k')
    title('Sine Fit')
    legend('Detrended','Bandpass','Sine Fit')
    hold off

%% Save my params
    if gof.rsquare > 0.3
        Amp = 2*sin_fit.a1*100;
        Phase = sin_fit.c1;
        Rate = sin_fit.b1/2/pi*60;
    else
        Amp = nan;
        Phase = nan;
        Rate = nan;
    end
catch
    Amp = nan;
    Phase = nan;
    Rate = nan;
end