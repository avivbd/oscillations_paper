function [T,MP,MC, delta] = OscillationsModel(Period,a,flag)
%Function containing carbon cycle model. Carries out time integration. 
% Called by OscillatingFborgFreq


switch nargin
    case 0
        close all
        clc
        clear
        tic
        Period = 1e6;  
        a = 0.5;
        flag.Fworg_os = 0;
        flag.which_os = 'kwp';

end


tmax = 10*Period;
t_interp = linspace(0,tmax,1000);

%Fluxes, masses, and coefficents.
F.v = 4e12;
MCss = 3.8e18;
MPss = 2e15;
Fbo = 12e12;
F.wo = 12e12;
Fws = 4e12;
Fwp = 3.6e10;
Fbp = Fwp;
k.ws = Fws/MCss;
k.bo = Fbo/MPss;
k.bp = Fbp/MPss;
k.wp = Fwp/MCss;
eps = 28;


del.v = -5;
del.wo = -25;
pCO20 = 500;

MC0 = MCss*1;
delss = 1;
del0 = (delss*MCss + (MC0 - MCss)*del.v)/(MC0);



flag.Fworg_os = 0; %oscillations in organic carbon weathering 1= yes

y0 = [MC0,MPss];

switch flag.which_os
    case {'Fv','Fws','Fborg','Fwp','Fbp','kbp','kwp','kws','kborg','kbp_kborg'}

        sol_opts = odeset('NonNegative',[1 2],'MaxStep',0.01*Period);  %'RelTol',1e-14
        [T, M] = ...
                ode23t(@(t,y) myodeM(t,y,k,F,a,Period,flag,MPss,MCss),...
                           t_interp ,y0,sol_opts); % Solve ODE for M

            sol2_opts = odeset('MaxStep',0.01*Period);  %'RelTol',1e-14
            [~, delta] = ...
                ode23t(@(t,y) myode_delta(t,y,k,F,a,Period,M(:,1),M(:,2),eps,del,flag,MPss,T),...
                           T, del0,sol2_opts); % Solve ODE for delta
    
   

        
end


%reconstruct the fluxes
MP = M(:,2);
MC = M(:,1);


switch nargin
    case 0
        plotfun(T,M,delta, a,Period, MCss,MPss,eps,pCO20)
end


function dy = myodeM(t,y,k,F,a,Period,flag,MPss,MCss)

sint =  a*sin(2*pi*t/Period); 

switch flag.which_os
    
    case 'kborg'
        switch flag.Fworg_os
            case 1
                dy(1) = F.v + F.wo*(1 - sint ) - k.ws*y(1) - k.bo*(1 + sint )*y(2);
            otherwise
                dy(1) = F.v + F.wo - k.ws*y(1) - k.bo*(1 + sint )*y(2);
        end
        
    case 'Fv'
        dy(1) = F.v*(1 + sint) + F.wo - k.ws*y(1) - k.bo*y(2);
    
    case 'Fborg'
        dy(1) = F.v + F.wo - k.ws*y(1) - k.bo*(1 + sint )*MPss;
        
    case 'kws'
       dy(1) = F.v + F.wo - k.ws*(1 + sint )*y(1) - k.bo*y(2); 
   
    case 'Fws'
       dy(1) = F.v + F.wo - k.ws*(1 + sint )*MCss - k.bo*y(2); 
       
    case {'kwp','kbp','Fwp','Fbp'}
       dy(1) = F.v + F.wo - k.ws*y(1) - k.bo*y(2); 
       
    case 'kbp_kborg'

       dy(1) = F.v + F.wo - k.ws*y(1) - k.bo*(1 + sint )*y(2);
    
       
end


switch flag.which_os
    case 'kwp'
        dy(2) = k.wp*(1 + sint )*y(1) - k.bp*y(2);   
    case 'kbp'
        dy(2) = k.wp*y(1) - k.bp*(1 + sint )*y(2); 
    case 'Fwp'
        dy(2) = k.wp*(1 + sint )*MCss - k.bp*y(2);   
    case 'Fbp'
        dy(2) = k.wp*y(1) - k.bp*(1 + sint )*MPss; 
    case 'kbp_kborg'
        dy(2) = k.wp*y(1) - k.bp*(1 + sint )*y(2);   
    otherwise
        dy(2) = k.wp*y(1) - k.bp*y(2);   
end


dy = dy(:);


function dy = myode_delta(t,y,k,F,a,Period,MCinterp,MPinterp,eps,del,flag,MPss,T)
sint = a*sin(2*pi*t/Period); 
MC = interp1(T,MCinterp,t); 
MP = interp1(T,MPinterp,t); 

switch flag.which_os
    case 'kborg'
        switch flag.Fworg_os
            case 1
                dy = (F.v*(del.v - y) + F.wo*(1 - sint )*(del.wo - y) + k.bo*(1 + sint )*MP*eps)/MC;
            otherwise
                dy = (F.v*(del.v - y) + F.wo*(del.wo - y) + k.bo*(1 + sint )*MP*eps)/MC;
        end
    case 'Fv'
        dy = (F.v*(1 + sint )*(del.v - y) + F.wo*(del.wo - y) + k.bo*MP*eps)/MC;
    
    case 'Fborg'
        dy = (F.v*(del.v - y) + F.wo*(del.wo - y) + k.bo*(1 + sint )*MPss*eps)/MC;  
        
    
    case {'kws','kwp','kbp','Fbp','Fwp','Fws'}
        dy = (F.v*(del.v - y) + F.wo*(del.wo - y) + k.bo*MP*eps)/MC;
        
    case 'kbp_kborg'
        dy = (F.v*(del.v - y) + F.wo*(del.wo - y) + k.bo*(1 + sint )*MP*eps)/MC;
        
       
        
end


function plotfun(T,Y,delta, a,Period,MCss,MPss,eps,pCO20)

figure('Position',[418     5   660   800])

        subplot 221
        plot(T,a*sin(2*pi*T/Period))
        xlabel('Time (y)')
        ylabel('Sinusoidal forcing')
        grid on
        axis square

        subplot 222
        [Ax,~,~] = plotyy(T,Y(:,1),T,Y(:,2));
        xlabel('Time (y)')
        ylabel(Ax(1),'M_C') % left y-axis
        ylabel(Ax(2),'M_P') % right y-axis
        grid on
        axis(Ax,'square')

        subplot 223
        plot(T,pCO20*(Y(:,1)/MCss).^2)
        xlabel('Time (y)')
        ylabel('pCO_2')
        grid on
        axis square


        subplot 224
        plot(T,delta)
%         hold on 
%         plot(T,-5 + 0.21425*(1+a*sin(2*pi*T/Period))*eps,'g')
        xlabel('Time (y)')
        ylabel('\delta^{13}C')
        grid on
        axis square
%         legend('Model results','steady-state result')
        
        toc





