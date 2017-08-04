function [  ] = OscillationsSense(  )

tic
clc
close all


makeplots = 1; %make plots? Can also just save output. 
save_file = 'no'; %save yes or no. 
flag_testing = 'testing'; %testing or full_run. Full run will take several hours. 

%also note that code is currently parallelized (using parfor). If you want
%to run it serially then change the parfor in line 47 to a regular for. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





switch flag_testing
    case 'full_run' %warning this will take a long time (a few hours)
        flag = struct('which_os', {'Fv', 'kws', 'Fws', 'kwp', 'Fwp', 'kborg', 'Fborg', 'kbp', 'Fbp'});
        f = logspace(-3,2,26);
        a  = linspace(0,1,26);
    case 'testing' %for testing
        flag = struct('which_os', {'Fv'});
        f = logspace(-3,2,10);
        a = linspace(0,1,5);
end

z = (1:length(flag))';
Period = 1./f;

[Tg,ag,zg] =  meshgrid(Period,a,z);
% TT = Tg(:); aa = ag(:); zz = zg(:);
[m,n,p] = size(zg); %a x f x Flag

%preallocate mem
 delta = cell(m,n,p);
 t = cell(m,n,p);
 MP = cell(m,n,p);
 MC = cell(m,n,p);
 amp = NaN(m,n,p);
 
%call the model for different combinations of amplitude a and period T 
% (in a*sin(2*pi*t/T)) 
parfor i = 1:m*n*p
    [t{i}, MP{i}, MC{i}, delta{i}] = OscillationsModel(Tg(i)*1e6,ag(i),flag(zg(i)));         
end

%find and remove cells of runs that did not integrate
[I,J] = find(cellfun(@length,delta)<1000);

for i = 1:length(I)
    MP{I(i),J(i)} = []; 
    MC{I(i),J(i)} = []; 
    delta{I(i),J(i)} = []; 
    t{I(i),J(i)} = []; 
end

%find the cells that did integrate and calculate the amplitude of
%oscillations in the final two cycles
Ind = find(cellfun(@length,delta)==1000);
for i = 1:length(Ind)
    amp(Ind(i)) = (( max(delta{Ind(i)}(800:1000)) - min(delta{Ind(i)}(800:1000)) )/2);
end


switch makeplots
     case 1
         for i = 1:length(flag)
         plotfunc(f,(amp(:,:,i)),a,flag(i).which_os)
         end
end

switch save_file
    case 'yes'
%         filname =  ['Os_output_' date];
%         for i = 1:length(flag)
       for i = 1:length(flag)     
        Output.(flag(i).which_os).amp = amp(:,:,i);
        Output.(flag(i).which_os).f = f;
        Output.(flag(i).which_os).a = a;
        Output.(flag(i).which_os).delta = delta(:,:,i);
        Output.(flag(i).which_os).t = t(:,:,i);
        Output.(flag(i).which_os).MP = MP(:,:,i);
        Output.(flag(i).which_os).MC = MC(:,:,i);
        
       end
       save('Output','Output','-append')
        
end



function plotfunc(f,amp,a,F)

colorset = flipud(jet(length(a)));
hfig2 = figure('Position',[664   359   560   420]);

[Ax,hLine1,~] = plotyy(f,amp,1./f,amp,'semilogx','semilogx');

set(hLine1,'LineStyle','-','Marker','o')
hlines = get(Ax(1),'Children');

 for i = 1:length(a)%:-1:1
    set(hlines(i),'Color', colorset(i,:))
 end
 
box off
set(Ax(2),'XAxisLocation','top','XDir','Reverse')
delete(get(Ax(2),'Children'))
set(Ax(2),'YTick',[])

% ylim(Ax(1), [0 8])

set(Ax(1), 'YTickMode','auto')
xlabel(Ax(1),'Frequency [m.y.^{-1}]')
xlabel(Ax(2),'Period [m.y.]')
ylabel(Ax(1),'Amplitude of \delta^{13}C Oscillations [O]')

hleg = legend([repmat('a = ',length(a),1) num2str(a(:))]);
M = [repmat('a = ',length(a),1) num2str(round(a'*100)/100) ];
legend(hlines,flipud(M))

 ht = text('Units','normalized','Position',[0.329 0.92],'String',...
     ['Sinusoidal Modulation of ' F],'FontSize',13);


 
 