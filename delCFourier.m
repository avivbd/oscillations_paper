function delCFourier()
%%
clc
clear
close all


%get the data 
[t_full,delC_full,gr,t,delC,t_interp, delC_interp, delC_sinterp] = getData;

%make Lomb-Scargle periodogram
makeLomb(t,delC,t_interp,delC_interp);

% flag = 'cum';%['dev';'cum'] deviation or cumulative decomposition
flag = 'dev';

%make the FT and draw the spectra
[f,Pxx,locs1_1, locs2_1, locs3_1,locs4_1,ifftsig_snip1, ifftsig_snip2, ifftsig_snip3,ifftsig_snip4] = ...
    makeFT(t_interp,delC_interp,t_full,delC_full,flag);

%if you want to plot only one part of the spectrum (the filtered signal)
 [delC_filtered,F_filtered] = plotFiltered(f,locs1_1, locs2_1, locs3_1,locs4_1,...
      ifftsig_snip1, ifftsig_snip2, ifftsig_snip3,ifftsig_snip4,t_interp,t,delC);

%make a spectrogram with the filtered data
makeSpect(t_interp,delC_filtered,F_filtered)

%make a spectrogram with the all the data
% makeSpect(t_interp,delC_interp,f)

%make a windowed PSD with the non-interpolated data
WindowPSD(t,delC,f,Pxx)

%make a windowed PSD with the interpolated data
% WindowPSD(t_interp,delC_interp,f,Pxx)


end



function [f,Pxx,locs1_1, locs2_1, locs3_1,locs4_1,...
    ifftsig_snip1, ifftsig_snip2, ifftsig_snip3,ifftsig_snip4] = ...
    makeFT(t_interp,delC_interp,t_full,delC_full,flag)        



scrsz = get(0,'ScreenSize');
% 
h.fig_spect = figure('Position',[1 0 scrsz(3) scrsz(3)*11/8.5]);


%do the FT
L = length(t_interp);
dt = t_interp(2) - t_interp(1);
Fs = 1/dt;

NFFT = L;% or set to 2^nextpow2(L); if you want zero padding
Y = fft(delC_interp,NFFT);
Pxx2 = abs(Y).^2/Fs/L;
Pxx = [Pxx2(1); 2*Pxx2(2:NFFT/2)];
f = 0 : Fs/(NFFT-1) : Fs/2;

tsfit0 = zeros(NFFT,1); tsfit0(1) = Y(1); 
ifftsig_snip0 = ifft(tsfit0,'Symmetric'); 

%get the first few points (long term trend)
locs1_1 = 2:2; %cannot include the first point, added below
locs1_2 = [1, locs1_1, NFFT-locs1_1+2]; %get symmetrically from both sides

tsfit1 = zeros(NFFT,1);
tsfit1(locs1_2) = Y(locs1_2);

% tsfit = Y; %for testing - do you get what you put in?

ifftsig_snip1 = ifft(tsfit1,'Symmetric');
ifftsig_snip1 = ifftsig_snip1(1:length(t_interp)); %cut off the padding


switch flag
    case 'cum'
        subplot(2,4,5); 
        hl = plot(t_full,delC_full,'.','Color',[0.729 0.832 0.957]);
%         alpha(hl,0.5)
        %add reverse transformed f(t) segments to plot of data
        line(t_interp,ifftsig_snip1,'Color','r')
        
        %uncomment to add average to plot        
%         line(t_interp,ifftsig_snip0,'Color','b','LineStyle','-.')
        
        
        axis square
        grid on
        box on
        xlabel('Time M.y.')
        ylabel('\delta^{13}C')
        set(gca,'XDir','reverse')
        xlim([0 541])
        

        subplot(2,4,1);
        loglog(f,Pxx,'o-','Color',[0.729 0.832 0.957])
        %plot transformed segments on periodogram
        %add first datapoint 
        
        
        line(f([1 locs1_1]),[(abs(Y(1))).^2/Fs/L ; 2*(abs(Y(locs1_1))).^2/Fs/L],...
            'Color','r','Marker','o')
        axis square
        xlabel('Cycles/M.y.')
        ylabel('Power')
        set(gca,'XTick',[10^-2 10^-1 10^0])
        
        line([10^-2 10^-2],[10^-6 10^4],'Color','k','LineStyle',':')
        line([10^-1 10^-1],[10^-6 10^4],'Color','k','LineStyle',':')
        line([10^0 10^0],[10^-6 10^4],'Color','k','LineStyle',':')
        line([10^-4 10^2],[10^-4 10^-4],'Color','k','LineStyle',':')
        line([10^-4 10^2],[10^-2 10^-2],'Color','k','LineStyle',':')
        line([10^-4 10^2],[10^0 10^0],'Color','k','LineStyle',':')
        line([10^-4 10^2],[10^2 10^2],'Color','k','LineStyle',':')
        
        box on
        xlim([0.001 10])
        ylim([0 1000])
        
        
    case 'dev'
        subplot(2,4,5); 
%         plot(t_full,delC_full,'.')
        hl = plot(t_full,delC_full,'.','Color',[0.729 0.832 0.957]);
        %add reverse transformed f(t) segments to plot of data
        line(t_interp,ifftsig_snip1,'Color','r')
        %uncomment to add average to plot        
%         line(t_interp,ifftsig_snip0,'Color','b','LineStyle','-.')        
        ylim([-6 10])
        axis square
        grid on
        box on
        xlabel('Time M.y.')
        ylabel('\delta^{13}C')
        set(gca,'XDir','reverse')
        xlim([0 541])
        

        subplot(2,4,1);
        loglog(f,Pxx,'o-','Color',[0.729 0.832 0.957])
        %plot transformed segments on periodogram
        line(f([1 locs1_1]),[(abs(Y(1))).^2/Fs/L ; 2*(abs(Y(locs1_1))).^2/Fs/L],...
            'Color','r','Marker','o')
        axis square
        box on
        
        line([10^-2 10^-2],[10^-6 10^4],'Color','k','LineStyle',':')
        line([10^-1 10^-1],[10^-6 10^4],'Color','k','LineStyle',':')
        line([10^0 10^0],[10^-6 10^4],'Color','k','LineStyle',':')
        line([10^-4 10^2],[10^-4 10^-4],'Color','k','LineStyle',':')
        line([10^-4 10^2],[10^-2 10^-2],'Color','k','LineStyle',':')
        line([10^-4 10^2],[10^0 10^0],'Color','k','LineStyle',':')
        line([10^-4 10^2],[10^2 10^2],'Color','k','LineStyle',':')
        
        xlabel('Cycles/M.y.')
        ylabel('Power')
        xlim([0.001 10])
        ylim([0 1000])
        set(gca,'XTick',[10^-2 10^-1 10^0])
end






% the next few
locs2_1 = locs1_1(end)+1:67;%50
locs2_2 = [locs2_1, NFFT-locs2_1+2]; %the second block of  coeffs

tsfit2 = zeros(NFFT,1);
tsfit2(locs2_2) = Y(locs2_2);

% tsfit = Y; %for testing - do you get what you put in?

%add back those frequencies from round 1
%inverse transform the signal including the lower frequencies
ifftsig_snip2 = ifft(tsfit2,'Symmetric');
ifftsig_snip2 = ifftsig_snip2(1:length(t_interp)); %cut off the padding


switch flag
    case 'cum'
        subplot(2,4,2);
        loglog(f,Pxx,'o-','Color',[0.729 0.832 0.957])
        %plot transformed segments on periodogram
        line(f([1 locs1_1 locs2_1]),[(abs(Y(1))).^2/Fs/L ; 2*(abs(Y(locs1_1))).^2/Fs/L ; 2*(abs(Y(locs2_1))).^2/Fs/L],'Color','k','Marker','o')
        axis square
        box on
        line([10^-2 10^-2],[10^-6 10^4],'Color','k','LineStyle',':')
        line([10^-1 10^-1],[10^-6 10^4],'Color','k','LineStyle',':')
        line([10^0 10^0],[10^-6 10^4],'Color','k','LineStyle',':')
        line([10^-4 10^2],[10^-4 10^-4],'Color','k','LineStyle',':')
        line([10^-4 10^2],[10^-2 10^-2],'Color','k','LineStyle',':')
        line([10^-4 10^2],[10^0 10^0],'Color','k','LineStyle',':')
        line([10^-4 10^2],[10^2 10^2],'Color','k','LineStyle',':')
        xlabel('Cycles/M.y.')
        ylabel('Power')
        xlim([0.001 10])
        ylim([0 1000])
        set(gca,'XTick',[10^-2 10^-1 10^0])
        
        subplot(2,4,6); 
        plot(t_full,delC_full,'.','Color',[0.729 0.832 0.957])
        %add reverse transformed f(t) segments to plot of data
        line(t_interp,ifftsig_snip1 + ifftsig_snip2,'Color','k')
        axis square
        grid on
        box on
        xlabel('Time [M.y.]')
        ylabel('\delta^{13}C')
        set(gca,'XDir','reverse')
        xlim([0 541])
    
    case 'dev'
        subplot(2,4,2);
        loglog(f,Pxx,'o-','Color',[0.729 0.832 0.957])
        %plot transformed segments on periodogram
        line(f(locs2_1),...
            2*(abs(Y(locs2_1))).^2/Fs/L,...
            'Color','k','Marker','o')
        axis square
        box on
        line([10^-2 10^-2],[10^-6 10^4],'Color','k','LineStyle',':')
        line([10^-1 10^-1],[10^-6 10^4],'Color','k','LineStyle',':')
        line([10^0 10^0],[10^-6 10^4],'Color','k','LineStyle',':')
        line([10^-4 10^2],[10^-4 10^-4],'Color','k','LineStyle',':')
        line([10^-4 10^2],[10^-2 10^-2],'Color','k','LineStyle',':')
        line([10^-4 10^2],[10^0 10^0],'Color','k','LineStyle',':')
        line([10^-4 10^2],[10^2 10^2],'Color','k','LineStyle',':')
        xlabel('Cycles/M.y.')
        ylabel('Power')
        xlim([0.001 10])
        ylim([0 1000])
        set(gca,'XTick',[10^-2 10^-1 10^0])

        subplot(2,4,6); 
%         plot(t_full,delC_full,'.')
        %add reverse transformed f(t) segments to plot of data
        line(t_interp,ifftsig_snip2,'Color','k')
        axis square
        ylim([-6 10])
        grid on
        box on
        xlabel('Time [M.y.]')
        ylabel('\delta^{13}C')
        set(gca,'XDir','reverse')
        xlim([0 541])
end




% the next few
locs3_1 = locs2_1(end)+1:953;%the third block of coeffs
locs3_2 = [locs3_1, NFFT-locs3_1+2]; %get the conjugates

tsfit3 = zeros(NFFT,1);
tsfit3(locs3_2) = Y(locs3_2);

% tsfit = Y; %for testing - do you get what you put in?

%add back those frequencies from round 1
%inverse transform the signal including the lower frequencies
ifftsig_snip3 = ifft(tsfit3,'Symmetric');
ifftsig_snip3 = ifftsig_snip3(1:length(t_interp)); %cut off the padding

switch flag
    case 'cum'
        subplot(2,4,3);
        loglog(f,Pxx,'o-','Color',[0.729 0.832 0.957])
        %plot transformed segments on periodogram
        line(f([1 locs1_1 locs2_1 locs3_1]),...
            [(abs(Y(1))).^2/Fs/L ; 2*(abs(Y(locs1_1))).^2/Fs/L ; ...
            2*(abs(Y(locs2_1))).^2/Fs/L ; 2*(abs(Y(locs3_1))).^2/Fs/L],...
            'Color','g','Marker','o')
        axis square
        box on
        line([10^-2 10^-2],[10^-6 10^4],'Color','k','LineStyle',':')
        line([10^-1 10^-1],[10^-6 10^4],'Color','k','LineStyle',':')
        line([10^0 10^0],[10^-6 10^4],'Color','k','LineStyle',':')
        line([10^-4 10^2],[10^-4 10^-4],'Color','k','LineStyle',':')
        line([10^-4 10^2],[10^-2 10^-2],'Color','k','LineStyle',':')
        line([10^-4 10^2],[10^0 10^0],'Color','k','LineStyle',':')
        line([10^-4 10^2],[10^2 10^2],'Color','k','LineStyle',':')
        xlabel('Cycles/M.y.')
        ylabel('Power')
        xlim([0.001 10])
        ylim([0 1000])
        set(gca,'XTick',[10^-2 10^-1 10^0])
        
        subplot(2,4,7); 
        plot(t_full,delC_full,'.','Color',[0.729 0.832 0.957])
        %add reverse transformed f(t) segments to plot of data
        line(t_interp,ifftsig_snip1 + ifftsig_snip2 + ifftsig_snip3,'Color','g')
        axis square
        grid on
        box on
        xlabel('Time [M.y.]')
        ylabel('\delta^{13}C')
        set(gca,'XDir','reverse')
        xlim([0 541])
        
    case 'dev'
        subplot(2,4,3);
        loglog(f,Pxx,'o-','Color',[0.729 0.832 0.957])
        %plot transformed segments on periodogram
        line(f(locs3_1),...
             2*(abs(Y(locs3_1))).^2/Fs/L,...
            'Color','g','Marker','o')
        axis square
        box on
        line([10^-2 10^-2],[10^-6 10^4],'Color','k','LineStyle',':')
        line([10^-1 10^-1],[10^-6 10^4],'Color','k','LineStyle',':')
        line([10^0 10^0],[10^-6 10^4],'Color','k','LineStyle',':')
        line([10^-4 10^2],[10^-4 10^-4],'Color','k','LineStyle',':')
        line([10^-4 10^2],[10^-2 10^-2],'Color','k','LineStyle',':')
        line([10^-4 10^2],[10^0 10^0],'Color','k','LineStyle',':')
        line([10^-4 10^2],[10^2 10^2],'Color','k','LineStyle',':')
        xlabel('Cycles/M.y.')
        ylabel('Power')
        xlim([0.001 10])
        ylim([0 1000])
        set(gca,'XTick',[10^-2 10^-1 10^0])
        
        
        subplot(2,4,7); 
%         plot(t_full,delC_full,'.')
        %add reverse transformed f(t) segments to plot of data
        line(t_interp,ifftsig_snip3,'Color','g')
        axis square
        ylim([-6 10])
        grid on
        box on
        xlabel('Time [M.y.]')
        ylabel('\delta^{13}C')
        set(gca,'XDir','reverse')
        xlim([0 541])
        
end



% the next few
locs4_1 = (locs3_1(end)+1):NFFT/2;%all the other  coeffs
locs4_2 = [locs4_1, NFFT-locs4_1+2]; %get the conjugates

tsfit4 = zeros(NFFT,1);
tsfit4(locs4_2) = Y(locs4_2);

%add back those frequencies from round 1
%inverse transform the signal including the lower frequencies
ifftsig_snip4 = ifft(tsfit4,'Symmetric');
ifftsig_snip4 = ifftsig_snip4(1:length(t_interp)); %cut off the padding

switch flag
    case 'cum'

        subplot(2,4,4);
        loglog(f,Pxx,'o-','Color',[0.729 0.832 0.957])
        %plot transformed segments on periodogram
        line(f([1 locs1_1 locs2_1 locs3_1 locs4_1]),...
            [(abs(Y(1))).^2/Fs/L ; 2*(abs(Y(locs1_1))).^2/Fs/L ; ...
            2*(abs(Y(locs2_1))).^2/Fs/L ; 2*(abs(Y(locs3_1))).^2/Fs/L; ...
            2*(abs(Y(locs4_1))).^2/Fs/L],...
            'Color','c','Marker','o')
        axis square
        box on
        line([10^-2 10^-2],[10^-6 10^4],'Color','k','LineStyle',':')
        line([10^-1 10^-1],[10^-6 10^4],'Color','k','LineStyle',':')
        line([10^0 10^0],[10^-6 10^4],'Color','k','LineStyle',':')
        line([10^-4 10^2],[10^-4 10^-4],'Color','k','LineStyle',':')
        line([10^-4 10^2],[10^-2 10^-2],'Color','k','LineStyle',':')
        line([10^-4 10^2],[10^0 10^0],'Color','k','LineStyle',':')
        line([10^-4 10^2],[10^2 10^2],'Color','k','LineStyle',':')
        % tsfit = Y; %for testing - do you get what you put in?
        xlabel('Cycles/M.y.')
        ylabel('Power')
        xlim([0.001 10])
        ylim([0 1000])
        set(gca,'XTick',[10^-2 10^-1 10^0])

        subplot(2,4,8); 
        plot(t_full,delC_full,'.','Color',[0.729 0.832 0.957])
        %add reverse transformed f(t) segments to plot of data
        line(t_interp,ifftsig_snip1 + ifftsig_snip2 + ifftsig_snip3 + ifftsig_snip4,...
            'Color','c')

        axis square
        grid on
        box on
        xlabel('Time [M.y.]')
        ylabel('\delta^{13}C')
        set(gca,'XDir','reverse')
        xlim([0 541])
        
    case 'dev'
        subplot(2,4,4);
        loglog(f,Pxx,'o-','Color',[0.729 0.832 0.957])
        %plot transformed segments on periodogram
        line(f(locs4_1),...
            2*(abs(Y(locs4_1))).^2/Fs/L,...
            'Color','c','Marker','o')
        axis square
        box on
        % tsfit = Y; %for testing - do you get what you put in?
        line([10^-2 10^-2],[10^-6 10^4],'Color','k','LineStyle',':')
        line([10^-1 10^-1],[10^-6 10^4],'Color','k','LineStyle',':')
        line([10^0 10^0],[10^-6 10^4],'Color','k','LineStyle',':')
        line([10^-4 10^2],[10^-4 10^-4],'Color','k','LineStyle',':')
        line([10^-4 10^2],[10^-2 10^-2],'Color','k','LineStyle',':')
        line([10^-4 10^2],[10^0 10^0],'Color','k','LineStyle',':')
        line([10^-4 10^2],[10^2 10^2],'Color','k','LineStyle',':')
        xlabel('Cycles/M.y.')
        ylabel('Power')
        xlim([0.001 10])
        ylim([0 1000])
        set(gca,'XTick',[10^-2 10^-1 10^0])

        subplot(2,4,8); 
%         plot(t_full,delC_full,'.')
        %add reverse transformed f(t) segments to plot of data
        line(t_interp,ifftsig_snip4,...
            'Color','c')
        ylim([-6 10])
        axis square
        grid on
        box on
        xlabel('Time [M.y.]')
        ylabel('\delta^{13}C')
        set(gca,'XDir','reverse')
        xlim([0 541])
        
end

end



function [delC_filtered,F_filtered] = plotFiltered(f,locs1_1, locs2_1, ...
    locs3_1,locs4_1, ifftsig_snip1, ifftsig_snip2, ifftsig_snip3,ifftsig_snip4 , t_interp,t,delC)


% F_filtered = f([ locs2_1 locs3_1  ]);
% delC_filtered = ifftsig_snip2 + ifftsig_snip3 ;

F_filtered = f([  locs3_1  ]);
delC_filtered =  ifftsig_snip3 ;


figure
plot(t_interp, delC_filtered)

 axis tight 
set(gca,'XDir','reverse')
xlabel('Time M.y.')
ylabel('\delta^{13}C')

end



function makeSpect(t_interp,delC_filtered,F_filtered)
%%
% spectrogram(delC_filtered,[],[],[],[]);
hf = figure('Position',[437    15   923   790]);
subplot 121
plot(t_interp, delC_filtered)
grid on
 axis tight square
set(gca,'XDir','reverse')
xlabel('Time M.y.')
ylabel('\delta^{13}C')

spectax = subplot(122);
n_win = 50;
window = round(length(delC_filtered)/n_win);%
noverlap = round(0.9*window);
F = F_filtered';
Fs = length(t_interp)/max(t_interp);%3.099250/541*1e4;
%window duration is given by
window_duration = 1/Fs*window;

[S,F,T,P] = spectrogram(delC_filtered,window,noverlap,F,Fs);
him = image(T,F(2:end),log10(P(2:end,:)),'CDataMapping','scaled');
axis tight;
ylabel('Frequency')
axis square
set(gca,'XDir','reverse','YScale','log','YDir','normal')
xlabel('Time M.y.')


 mycmap = MYCMAP();
 colormap(gca,mycmap)

 colorbar('location','NorthOutside')

%%


end


function [pxx,f] = makeLomb(t_unique,delC_unique,t_interp,delC_interp)
%% Lomb-Scargle


%get a lomb-scargle psd
[pxx,f] = plomb(delC_unique,t_unique,'psd');

%plot the lomb-scargle psd
figure
subplot 121
semilogx(f,pxx)
hold on
hax = gca;
axis square
set(hax,'YGrid','on')
xlabel('Frequency [cycles/Ma]')
ylabel('Power')
axis tight
set(hax,'XTick',10.^(-4:1))

%convert some part of the signal back into the time domain

loc = 1:6; %which indexes

delC_regression(:,1) = delC_unique;

%preallocate

for i = 1:length(loc)

    f0 = f(loc(i));

    %Now know the frequency. Find the amplitude and phase via linear regression
    ft = 2*pi*f0*t_unique;

    ABC = [ones(size(ft)) cos(ft) sin(ft)] \ delC_regression(:,i);

    fx = 2*pi*f0*t_unique;

    y(:,i) = [ones(size(fx)) cos(fx) sin(fx)] * ABC;
    
    delC_regression(:,i+1) = delC_regression(:,i) - y(:,i);

end

ycumsum = cumsum(y,2);

%highlight the locations that were converted back to the time domain
semilogx(f(1:6),pxx(1:6),'o-r')

%plot the components that you want back in the time domain (1,2,...)
subplot 122
hax2 = gca;
plot(t_interp,delC_interp,'Color','k')
set(gca,'XDir','reverse')
hold on
plot(t_unique,ycumsum(:,end),'r','Parent',hax2)

xlabel('Time (My)')
ylabel('\delta^{13}C')

legend(hax2,'data', 'L-S')

axis square

end


function WindowPSD(t,x,f,Pxx)
% n = linspace(000,400,3)
% n = 200;
n=0;
% NW = 4;
NW = [0.75,12];
% figure('Position', [ 0     0   1400   800]);
figure1 = figure();
axes1 = axes('Parent',figure1,...
    'Position',[0.13 0.11 0.775 0.773333333333333]);

% loglog(f,Pxx,'o')
% hold on



for i=1:length(NW)
    myFun(t,x,NW(i))
end

nn=180;
loglog(f(nn:end), .1*(f(nn:end).^-2),'g')
% loglog(f(nn:end), .1*(f(nn:end)).^-1.8,'g')
% loglog(f(1:nn), pi/10 ./(f(1:nn)).^1,'g')




%% 
    function myFun(t,x,NW)
        Fs = 1/mean(diff(t));
        [psdmtm,Fxx] = pmtm(x,NW,length(x),Fs,...
            'Droplasttaper',false);        
        h = loglog(Fxx(round(NW)+2:end),psdmtm(round(NW)+2:end),'o');
        hold on
        xlim([1e-2,10])
        ylim([10^-5, 10^2])
    end


legend(['NW = ' num2str(NW(1))],['NW = ' num2str(NW(2))],'y = 0.1*f^{-2}')


xlabel('Freq [M.yr.^{-1}]')
pos = get(gca,'Position');
tcks = sort((1./(get(gca,'XTick'))));
xLim = sort(10.^(-log10(get(gca,'XLim'))));
yLim = sort(((get(gca,'YLim'))));
box off
set(gca,'YGrid','on')

%vertical grid lines
vert = -2:1:0;
for i = 1:length(vert)
    line([10^vert(i) 10^vert(i)],yLim,'Color','k','LineStyle',':')
end

 
ha2 = axes('Position',pos,'XAxisLocation','top','YAxisLocation','right',...
     'Color','none','YTick',[],'XDir','reverse','xscale','log',...
     'XLim',xLim);

xlabel('Period [M.yr.]')
%  uistack(ha2,'bottom')
 
title('Thompson Multitaper PSD Estimate');


end


function mycmap = MYCMAP()

mycmap = [
         0         0    0.5625
         0         0    0.5747
         0         0    0.5868
         0         0    0.5990
         0         0    0.6111
         0         0    0.6233
         0         0    0.6354
         0         0    0.6476
         0         0    0.6597
         0         0    0.6719
         0         0    0.6840
         0         0    0.6962
         0         0    0.7083
         0         0    0.7205
         0         0    0.7326
         0         0    0.7448
         0         0    0.7569
         0         0    0.7691
         0         0    0.7812
         0         0    0.7934
         0         0    0.8056
         0         0    0.8177
         0         0    0.8299
         0         0    0.8420
         0         0    0.8542
         0         0    0.8663
         0         0    0.8785
         0         0    0.8906
         0         0    0.9028
         0         0    0.9149
         0         0    0.9271
         0         0    0.9392
         0         0    0.9514
         0         0    0.9635
         0         0    0.9757
         0         0    0.9878
         0         0    1.0000
         0    0.2000    1.0000
         0    0.4000    1.0000
         0    0.6000    1.0000
         0    0.8000    1.0000
         0    1.0000    1.0000
    0.2500    1.0000    0.7500
    0.5000    1.0000    0.5000
    0.7500    1.0000    0.2500
    1.0000    1.0000         0
    1.0000    0.8000         0
    1.0000    0.6000         0
    1.0000    0.4000         0
    1.0000    0.2000         0
    1.0000         0         0
    0.9848         0         0
    0.9697         0         0
    0.9545         0         0
    0.9091         0         0
    0.8636         0         0
    0.8182         0         0
    0.7727         0         0
    0.7273         0         0
    0.6818         0         0
    0.6364         0         0
    0.5909         0         0
    0.5455         0         0
    0.5000         0         0
    ];

end