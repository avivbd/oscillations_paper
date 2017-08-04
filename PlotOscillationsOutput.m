function PlotOscillationsOutput
 close all
 clear
 clc

load('Output_Final') %put file in same dir or change path

% which figures to make?
% fig 1: amp as a function of f
% fig 2: delta as a func of t 
% fig 3: like fig 1 but only Fv and kbp
% fig 4: like fig 3 but only Fv, kbp, and kborg
% fig 5: like fig 2 but Fv and kbp for 1, 10, 100 m.y.

flag.whicfig = 3;


%which fields to plot? {'Fv','kws','Fws','kwp','kbp','kborg','Fwp','Fborg','Fbp'}
F = {'Fv','kws','Fws','kwp','kbp','kborg','Fwp','Fborg','Fbp'};


%save the resulting plots ['yes','no']
flag.save_fig1 = 'no'; %fig 1: amp as a function of f

flag.save_fig2 = 'no'; %fig 2: delta as a func of t for f = 1 m.y




%%%%%%%%%%%%%%Fig 1%%%%%%%%%%%%%%%%%%%%%%

switch flag.whicfig
    case 1

        %create the figure for delC amp as a function of f
        hfig1 = figure('PaperType','usletter','PaperPositionMode','auto',...
         'Position',[50 50 1200 770]); %letter paper portrait

        %loop through the variables and make subplots
        for i = 1:length(F)
            plotfun1(i,(Output.(F{i}).amp)',Output.(F{i}).f,(Output.(F{i}).a),...
                Output.(F{i}).delta,Output.(F{i}).t,F{i},flag)
        end

    case 2
    
%%%%%%%%%%%%%%Fig 2%%%%%%%%%%%%%%%%%%%%%%    
% Period of 1 m.y.   

         hfig2 = figure('PaperType','usletter','PaperPositionMode','auto',...
                'Position',[50 50 1200 770]); %letter paper portrait

         for i = 1:length(F)
             plotfun2(i,Output.(F{i}).amp,Output.(F{i}).f,(Output.(F{i}).a)',...
                 Output.(F{i}).delta,Output.(F{i}).t,F{i},flag)
         end
         
     case 3
    
%%%%%%%%%%%%%%Fig 3%%%%%%%%%%%%%%%%%%%%%%    
  

         hfig3 = figure;
         j = [1,5];
         for i = 1:length(j)
             plotfun3(j,i,(Output.(F{j(i)}).amp)',Output.(F{j(i)}).f,(Output.(F{j(i)}).a),...
                 Output.(F{j(i)}).delta,Output.(F{j(i)}).t,F{j(i)},flag)
         end   
         
     case 4
          hfig4 = figure('Position',[1 5 1055 693]);
            j = [1,5,6];
         for i = 1:length(j)
             plotfun4(j,i,(Output.(F{j(i)}).amp)',Output.(F{j(i)}).f,(Output.(F{j(i)}).a),...
                 Output.(F{j(i)}).delta,Output.(F{j(i)}).t,F{j(i)},flag)
         end
    
    case 5
          hfig5 = figure;
          F = {'Fv','kbp'};
          m = [16, 12]; %size of amp matrix 
          f_ind = [19, 16, 12];
          count = 1;
          for j = 1:length(f_ind)
             for i = 1:length(F)
                 plotfun5(count,m(i),i,j,Output.(F{i}).amp,Output.(F{i}).f,(Output.(F{i}).a)',...
                     Output.(F{i}).delta,Output.(F{i}).t,F{i},flag,f_ind(j))
                 count = count + 1;
               
                 
             end
          end
           
          annotation(hfig5,'line',[0.5 0.5],...
    [0.974083056478405 0.0564784053156147]);

end
    
end


function plotfun1(i,amp,f,a,delta,t,name,flag)

subplot(3,3,i)

[m,n] = size(amp);%size of amp matrix 
%which is length(a) rows and length(f) columns
ColorSet = flipud(jet(m));%for line colors

%plot amp as a function of frequency and period
[Ax,hLine1,~] = ...
plotyy(f,amp,1./f,amp,'semilogx','semilogx');

% set(Ax,'Parent', hfig)
set(hLine1,'LineStyle','-');%,'Marker','o'
box off
set(Ax(2),'XAxisLocation','top','XDir','Reverse')
delete(get(Ax(2),'Children'))
set(Ax(2),'YTick',[])
% ylim(Ax(1), [0 18])
  xlim(Ax(1), [10^(-3),10^2])
  xlim(Ax(2), [10^(-2),10^3])
hlines = (get(Ax(1),'Children'));

for ii = 1:m
set(hlines(ii),'Color', ColorSet(ii,:))
end

set(Ax(1), 'YTickMode','auto')
set(Ax(1), 'XTick',[0.001 0.01 0.1 1 10 100 1000])
set(Ax(2), 'XTick',[0.001 0.01 0.1 1 10 100 1000])
 
%organize the axes labels so that they are only on the margins
 switch i
     case {7 }
         xlabel(Ax(1),'Frequency [m.y.^{-1}]')
         ylabel(Ax(1),'Amplitude of \delta^{13}C Oscillations [O]')
     case{8 9}
        xlabel(Ax(1),'Frequency [m.y.^{-1}]')
      case {1}
          ylabel(Ax(1),'Amplitude of \delta^{13}C Oscillations [O]')
          xlabel(Ax(2),'Period [m.y.]')
     case {2 3}
        xlabel(Ax(2),'Period [m.y.]')
     case {4}
        ylabel(Ax(1),'Amplitude of \delta^{13}C Oscillations [O]')
 end

%put title inside plots
ht = text('Units','normalized','Position',[0.2 0.92],'String',...
     ['Sinusoidal Modulation of ' name],'FontSize',13);

%put in a legend for the a's
M = flipud([repmat('a = ',length(a),1) num2str(round(a'*100)/100)]);
hleg = legend(hlines, M);
set(hleg,'Position',[0.927726296880477 0.32125 0.0583793375394322 0.40375],...
    'units','normalized')

switch flag.save_fig1
    case 'yes'
        set(gcf,'PaperPositionMode','auto')
        print(hfig,'OscillationsOutputAll', '-depsc2','-painters','-r300')
end






end

function plotfun2(i,amp,f,a,delta,t,name,flag)

subplot(3,3,i)
hold on

m = length(a); %size of amp matrix 
ColorSet = (jet(m));%for line colors

T = 1./f';
%choose the period
% f_ind = find(f == 0.01);
f_ind = 19;
T(f_ind)



for ii = 1:m
plot(t{ii,f_ind}(1:end)/1e6,delta{ii,f_ind},'Color',ColorSet(ii,:))
end

xlabel('Time [m.y.]')
ylabel('\delta^{13}C ')
box on
% axis tight

ht = text('Units','normalized','Position',[0.2 0.92],'String',...
     ['Sinusoidal Modulation of ' name],'FontSize',13);








end

function plotfun3(j,i,amp,f,a,delta,t,name,flag)

subplot(1,2,i)
switch i
    case 1
        amp = amp(:,1:16);
        
    case 2
        amp = amp(:,1:12);
end

[m,n] = size(amp);%size of amp matrix 
%which is length(a) rows and length(f) columns
ColorSet = flipud(jet(n));%for line colors

%plot amp as a function of frequency and period
[Ax,hLine1,~] = ...
plotyy(f,amp,1./f,amp,'semilogx','semilogx');

% set(Ax,'Parent', hfig)
set(hLine1,'LineStyle','-');%,'Marker','o'
box off
set(Ax(2),'XAxisLocation','top','XDir','Reverse')
delete(get(Ax(2),'Children'))
set(Ax(2),'YTick',[])
% ylim(Ax(1), [0 18])
  xlim(Ax(1), [10^(-3),10^2])
  xlim(Ax(2), [10^(-2),10^3])
hlines = (get(Ax(1),'Children'));

for ii = 1:n
set(hlines(ii),'Color', ColorSet(ii,:))
end

set(Ax(1), 'YTickMode','auto')
set(Ax(1), 'XTick',[0.001 0.01 0.1 1 10 100 1000])
set(Ax(2), 'XTick',[0.001 0.01 0.1 1 10 100 1000])
 
xlabel(Ax(1),'Frequency [m.y.^{-1}]')
xlabel(Ax(2),'Period [m.y.]')
ylabel(Ax(1),'Amplitude of \delta^{13}C Oscillations')
axis(Ax,'square')
% switch i
%     case 1
%         caxis([0 0.6])
%         colorbar('location','WestOutside')
%     case 2
%         caxis([0 0.42])
%         colorbar('location','eastOutside')
% end


%put title inside plots
ht = text('Units','normalized','Position',[0.1 0.92],'String',...
     ['Sinusoidal Modulation of ' name],'FontSize',13);

% %put in a legend for the a's
% M = flipud([repmat('a = ',length(a),1) num2str(round(a'*100)/100)]);
% hleg = legend(hlines, M);
% set(hleg,'Position',[0.927726296880477 0.32125 0.0583793375394322 0.40375],...
%     'units','normalized')
% switch i
%     case 2
%         colorbar
% end

end

function plotfun4(j,i,amp,f,a,delta,t,name,flag)

subplot(1,3,i)
switch i
    case 1
        amp = amp(:,1:16);
        
    case 2
        amp = amp(:,1:12);
        
    case 3
        amp = amp(:,1:20);    
end

[m,n] = size(amp);%size of amp matrix 
%which is length(a) rows and length(f) columns
ColorSet = flipud(jet(n));%for line colors

%plot amp as a function of frequency and period
[Ax,hLine1,~] = ...
plotyy(f,amp,1./f,amp,'semilogx','semilogx');

% set(Ax,'Parent', hfig)
set(hLine1,'LineStyle','-');%,'Marker','o'
box off
set(Ax(2),'XAxisLocation','top','XDir','Reverse')
delete(get(Ax(2),'Children'))
set(Ax(2),'YTick',[])
% ylim(Ax(1), [0 18])
  xlim(Ax(1), [10^(-3),10^2])
  xlim(Ax(2), [10^(-2),10^3])
hlines = (get(Ax(1),'Children'));

for ii = 1:n
set(hlines(ii),'Color', ColorSet(ii,:))
end

set(Ax(1), 'YTickMode','auto')
set(Ax(1), 'XTick',[0.001 0.01 0.1 1 10 100 1000])
set(Ax(2), 'XTick',[0.001 0.01 0.1 1 10 100 1000])
 
xlabel(Ax(1),'Frequency [m.y.^{-1}]')
xlabel(Ax(2),'Period [m.y.]')
ylabel(Ax(1),'Amplitude of \delta^{13}C Oscillations')
axis(Ax,'square')
% switch i
%     case 1
%         caxis([0 0.6])
%         colorbar('location','WestOutside')
%     case 2
%         caxis([0 0.42])
%         colorbar('location','eastOutside')
% end


%put title inside plots
ht = text('Units','normalized','Position',[0.1 0.92],'String',...
     ['Sinusoidal Modulation of ' name],'FontSize',13);

% %put in a legend for the a's
% M = flipud([repmat('a = ',length(a),1) num2str(round(a'*100)/100)]);
% hleg = legend(hlines, M);
% set(hleg,'Position',[0.927726296880477 0.32125 0.0583793375394322 0.40375],...
%     'units','normalized')
% switch i
%     case 2
%         colorbar
% end



end

function plotfun5(count,m,i,j,amp,f,a,delta,t,name,flag,f_ind)

subplot(3,2,count)
hold on


ColorSet = (jet(m));%for line colors

for ii = 1:m
    switch i
        case 1
            plot(t{ii,f_ind}/1e6,delta{ii,f_ind},'Color',ColorSet(ii,:))%(350:end)
        case 2
            plot(t{ii,f_ind}/1e6,delta{ii,f_ind},'Color',ColorSet(ii,:))
    end

end

axis tight
box on

ylabel('\delta^{13}C')
xlabel('Time [m.y.]')

switch i
    case 1
        ylim([-1.5 3])
    case 2
        ylim([-3.5 8])
end
         

switch count
       case 1
           title('Sinusoidal Modulation of Fv');  
       case 2
           title('Sinusoidal Modulation of kwp');                 
end



end
