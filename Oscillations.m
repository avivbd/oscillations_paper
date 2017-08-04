function []=Oscillations()
% Produces Figures XXX in Bachan 201X. 
% Model sketch:
%                                   
% 
%                       _______
%                      |       |
%      kwp * MC  -->   |   M1  |   -->  kbp * MP
%                      |       | 
%                      |_______|
% 
%                       
%                       
%                    
%                                          
%                       _______
%                      |       |
%          Fin   -->   |   M2  |   --> kbc * MC
%                      |       |   --> kbo * MP
%                      |_______|
%                      
% 
% 
% dM1/dt = F1in - aM1 + bM2
% dM2/dt = F2in - dM2 + cM1
%%   
   clc
   clear
   close all
   rng(1)
  

   


   
        %have five perturbing factors: kbo kbc kwp kbp kwo
        %choose the amounts by which you'd like to change them
        %make another row for a different run
        %           kbo  kbc kwp kbp kwo  
        factor = [  
                  1       1   1    1  0
                  1      .1   1   .1  0
                  3      .1   1   .1  0
                 ];   
             
          forcing = [4];  
                % 0 none
                 % 2 gaussian shaped
                 % 3 pulse
                 % 4 instantenous carbon mass addition
                 % 5 intsantenous P addition
                 % 6 oscillating forcing in C
                 % 7 stochastic forcing (not implemented yet)
                 % 8 step 
                 % 9 exponential decay C
                 % 11 sigmoidal decay of C addition
        
        [m,n] = size(factor);%get dims

        output = cell(m,1);%preallocate mem
    for j = 1:length(forcing)
           for i = 1:m %number of model iterations with different factors
                %call model function
                [output{i}] = modelfun(factor(i,:),forcing(j));

                %uncomment to plot all results
                output{i}.h = plotfun1(output{i},i);
                
                
           end
    end

    
plotfun2(output)
     

assignin('base', 'Output',output)
   
end


function h = plotfun1(Output,i)

ColOrd = get(gca,'ColorOrder');
[Cm,Cn] = size(ColOrd);
ColRow = rem(i,Cm);

    if ColRow == 0

      ColRow = Cm;

    end

    % Get the color
    col = ColOrd(ColRow,:);



    
h.fig = figure(1);
set(gcf,'Position',[440   100   560   698])


%m rows n columns of subplots
m = 3; n = 2;
    
    h.ax1 = subplot(m,n,1);
    hold on
    plot(Output.t, Output.delC,'Color',col );
%     hold on
%     plot(Output.t, Output.delC_linz+ IC.delC,'Color','b');
    ylabel('\delta^{13}C')
    grid on
    box on


    h.ax2 = subplot(m,n,2); 
    hold on
%     [AX2, h2(1), h2(2)] = plotyy(t, delC,t,delC_alt);
%     set(AX2(1),'YLim',get(AX2(2),'YLim'))
    plot(Output.t, Output.pCO2,'Color',col );
    ylabel('pCO_2')
    grid on
    box on
    
     
    h.ax3 = subplot(m,n,3);
    hold on
    plot(Output.t, Output.MP ,'Color',col);
    ylabel('MP')
    grid on
    box on
     
    h.ax4 = subplot(m,n,4);
    hold on
    plot(Output.t, Output.MC ,'Color',col);
    ylabel('MC')
    grid on
    box on
    
    h.ax5 = subplot(m,n,5);
    h.l5(1) = plot(Output.t, Output.ALKP_mean,'Color',col);
    hold(h.ax5(1),'on')
    ylabel(h.ax5(1),'ALK:P')
    grid on
    box on
    
    h.ax6 = subplot(m,n,6);
     h.l6 = plot(Output.t,Output.CP_mean,'Color',col);
%     h.l6 = plot(Output.t,Output.F_extra_out,'Color',col);
    hold(h.ax6,'on')
     ylabel(h.ax6,'C:P')
    grid(h.ax6,'on')
    box on


end

function plotfun2(output)


figure('Position',[360 5 560 693]);
set(gcf,'PaperPositionMode','auto')
ColOrd = get(gca,'ColorOrder');


Ax1 = subplot(3,1,1);
plot(output{1}.t,output{1}.delC,'Color',ColOrd(1,:));
grid on

daspect([2e6 10 1])


Ax2 = subplot(3,1,2);
plot(output{2}.t,output{2}.delC,'Color',ColOrd(1,:));
grid on
ylim([-4 6])
daspect([2e6 10 1])

Ax3 = subplot(3,1,3);
plot(output{3}.t,output{3}.delC,'Color',ColOrd(1,:));
grid on
ylim([-8 8])
daspect([2e6 10 1])



end

function [Output] = modelfun(factor,forcing) 
% This is where the model is run
   
   %choose type of perturbation
   switches.per_type = forcing;
    
    switches.stoch = 'n'; %stochastic component? ['y','n']
                 
    switches.kbp_switch = 'P'; %choose P or O dependency on P burial ['O','P']
    
    switches.Fbcarb_switch = 'sil'; %choose carbonate burial parameterization 
%   or silicate burial in which Fwcarb + Fwsil = Fbcarb. Main difference is that 
%   Fworg appears in isotope balance if using 'sil' ['carb','sil']
    
    switches.vars_switch = 'dev'; %Use deviation vars ['dev'] or regular ones ['reg'] (use 'dev' for oscillations)
    
    %run length [yrs] 
    t_run_start = -1e6;
    t_run_end = 6e6;
    
    t_perturbation_start = 0e6;
    t_perturbation_end = 5e6;
    
    
    %initial mass
    M0.P = 2e15;
    M0.C = 3.8e18;
    M0.O = 3.78e19;

    %for export
    Output.M0.P = M0.P;
    Output.M0.C = M0.C;
    Output.M0.O = M0.O;
    
    
    %isotpic parameters
    p.delC_in = -5;
    p.epsilon = 25;
    
    %weathering fluxes
    F0.in = 50e12;
    p.forg = 0.2;
    F0.volc = 5e12;
    F0.volc_red = p.forg*F0.volc;
    F0.volc_ox = (1 - p.forg)*F0.volc;
    F0.worg = p.forg*(F0.in - F0.volc);
    F0.wcarb = (1- p.forg)*(F0.in - F0.volc);
    F0.wsil  = F0.volc_ox;
    F0.wp = 36e9;
    F0.bp = 36e9;
    
    %isotopic compositions of weathering fluxes
    p.delwcarb = p.delC_in + (F0.worg + F0.volc_red)/F0.in*p.epsilon;
    p.delworg = p.delwcarb - p.epsilon;
    
    %burial fluxes
    F0.bcarb = (1- p.forg)*F0.in;
    F0.borg = p.forg*F0.in;
    
    %isotopic composition of carbonate burial flux
    p.delC = p.delC_in + p.forg*p.epsilon;
%     Fworg_0 = 14e12;
    
    %flux coefficients
    k0.wo = F0.worg/M0.O;
    k0.bc = F0.bcarb/M0.C;
    k0.ws = F0.wsil/M0.C;
    k0.bo = F0.borg/M0.P;
    k0.wp = F0.wp/M0.C;  
    
    switch switches.kbp_switch
        case 'P'
            k0.bp = F0.bp/M0.P;
        case 'O'
            k0.bp = F0.bp/M0.O;
    end
    
    
    %save for export to base workspace
    Output.p.delC = p.delC;
    Output.k0.bo = k0.bo;
    Output.k0.bc = k0.bc;
    Output.k0.wp = k0.wp;
    Output.k0.bp = k0.bp;
    Output.k0.wo = k0.wo;
    Output.k0.ws = k0.ws;



    %modify initial values by certain amounts set in factor above
    k.bo = factor(1)*k0.bo;
    k.bc = factor(2)*k0.bc;
    k.ws = factor(2)*k0.ws;
    k.wp = factor(3)*k0.wp;
    k.bp = factor(4)*k0.bp;
    k.wo = factor(5)*k0.wo;

    %save for export to base workspace
    Output.k.bo = k.bo;
    Output.k.wo = k.wo;
    Output.k.bc = k.bc;
    Output.k.ws = k.ws;
    Output.k.wp = k.wp;
    Output.k.bp = k.bp;
    
    Output.kappa = (4*k.wp*k.bo)./(k.bp - k.bc).^2;
    disp(Output.kappa)
     
%Uncomment to disply condition for oscillations
%     disp({'For oscillations',(kbp - kbc)^2 4*kwp*kbo})
    
     
    Output.CP_0 = k0.bo/k0.bp;
    Output.ALKP_0 = k0.bc/k0.wp;
%Uncomment to display initial C:P and ALK:P     
%     disp({'C:P_0', Output.CP_0, 'ALK:P_0', Output.ALKP_0})

     
    Output.CP = k.bo/k.bp;
    Output.ALKP = k.bc/k.wp;
%Uncomment to display perturbation C:P and ALK:P 
%     disp({'C:P', Output.CP, 'ALK:P', Output.ALKP})

    A = [ -k.bp  k.wp; -k.bo -k.bc];
    eigA = eig(A);
    Output.Eig = eigA;
%Uncomment to display calculated eigenvalues of coeff matrix
%     disp({'eig A' num2str(eigA(1)) num2str(eigA(2))})
    
%to compare analytical and computed eigenvalues    
%     lambda1 =( - (kbp + kbc) + sqrt( (kbp + kbc)^2 - 4*(kbp*kbc + kwp*kbo) ) )/2;
%     lambda2 =( - (kbp + kbc) - sqrt( (kbp + kbc)^2 - 4*(kbp*kbc + kwp*kbo) ) )/2;
%     disp(['lambda_1:  ' num2str(lambda1); 'lambda_2:  ' num2str(lambda2)])

    

    %% 
    %perturbation
    n=1000;     %number of points in the interpolation    
    
    switch switches.per_type 
        
        case 1
            M0.extra = 0; %Total mass of stuff to be released over t_perturbation
    
            p.tinterp=[t_run_start  ...
            linspace(t_perturbation_start,t_perturbation_end,n+2) ...
            t_run_end];
    
            F0.extra_interp=[0 0 M0.extra*normpdf(linspace(-4,4,n),0,1)/...
            max(normpdf(linspace(-4,4,n),0,1)) 0 0];    % gaussian shaped excursion with an area equal to Mextra
        
        case 2
            M0.extra = 5e12; %Total mass of stuff to be released over t_perturbation
    
            p.tinterp=[t_run_start  ...
            linspace(t_perturbation_start,t_perturbation_end,n+2) ...
            t_run_end];
    
            F0.extra_interp=[0 0 M0.extra*normpdf(linspace(-4,4,n),0,1)/...
            max(normpdf(linspace(-4,4,n),0,1)) 0 0];    % gaussian shaped excursion with an area equal to Mextra
        
        case 3 
            
            M0.extra = 1e13;
            
            p.tinterp=[t_run_start  ...
            linspace(t_perturbation_start,t_perturbation_end,n+2) ...
            t_run_end];
    
            F0.extra_interp=[0 0 M0.extra*ones(1,n) 0 0];    
         
        case 4 
            M0.extra = 1.9e18;           
            
        case 5 
            M0.extra = 1e15;   
            
        case 6
%             M0.extra = 1e13;
            M0.extra = 5e12;
            
            p.tinterp=[t_run_start  ...
            linspace(t_perturbation_start,t_perturbation_end,n+2) ...
            t_run_end];
    
            F0.extra_interp=[0 0 M0.extra*sin(2*pi*p.tinterp(1:end-4)/1e6) 0 0 ];    
         
        case 7
            M0.extra = 1e13;
%             Mextra = 5e12;
            
            p.tinterp=[t_run_start  ...
            linspace(t_perturbation_start,t_perturbation_end,n+2) ...
            t_run_end];
    
            F0.extra_interp=[0 0 
                M0.extra*sin(2*pi*p.tinterp(1:end-4)/1e6) 0 0 ];    
        case 8
            
            Fextra = 5e13;
            M0.extra = 0;
            
            p.tinterp=[t_run_start  ...
            linspace(t_perturbation_start,t_perturbation_end,n+2) ...
            t_run_end];
    
            F0.extra_interp=[0 0 Fextra*ones(1,n+2) ]; 
         
        case 9
            
            Fextra = 5e13;
            M0.extra = 0;
            
            tlinspace = linspace(t_perturbation_start,t_perturbation_end,n+2);
            p.tinterp=[t_run_start  tlinspace   t_run_end];
    
            F0.extra_interp=[0 0 Fextra*exp(-1/1e6*tlinspace) ]; 
            
        case 10
            M0.extra = 1.9e18;  
            M0.extra_P = 1e15;    
        
        case 11
            
            Fextra = 5e13;
            M0.extra = 0;
            
            tlinspace = linspace(t_perturbation_start,t_perturbation_end,n+2);
            p.tinterp=[t_run_start  tlinspace   t_run_end];
    
            F0.extra_interp = [0 0 1*Fextra - Fextra./(1 + exp(-1/1e6*tlinspace+5)) ];     
            
    end

    
    

    
%%    
    
%initial conditions 
switch switches.vars_switch
    case 'dev'
        ic = [0 0 p.delC 0 ]; %initial conditions
    case 'reg'
        ic = [M0.P M0.C p.delC M0.O];
end

% a = [];

%call the solver
switches.solvernum = 1;
%the solver up to the pertubation.
sol1 = ode15s(@(t,y) odefun(t,y,k0,switches,F0,M0,p),[t_run_start, t_perturbation_start],ic); 
%the start of pertubation to end.



switch switches.vars_switch
    case 'dev'
        switch switches.per_type 
            case 4
                MC_pre = sol1.y(2,end) + M0.C;
                delC_pre = sol1.y(3,end) ;
                delC_post = (MC_pre*delC_pre + M0.extra*p.delC_in)/(MC_pre + M0.extra);            
                ic = [0; M0.extra; delC_post; 0]; 
        
        end
    case 'reg'
        switch switches.per_type 
        case 4
            delC_pre = sol1.y(3,end) ;
            MC_pre = sol1.y(2,end);
            delC_post = (MC_pre*delC_pre + M0.extra*p.delC_in)/(MC_pre + M0.extra);            
            ic = [0; MC_pre + M0.extra; delC_post; 0]; 
        end
        
        
end


switch switches.per_type 
%         case 4
%             delC_pre = sol1.y(3,end) ;
%             delC_post = (MC_pre*delC_pre + M0.extra*p.delC_in)/(MC_pre + M0.extra);   
%             switch switches.vars_switch
%                 case 'dev'
% %             MC_post = MC_pre + M0.extra;
%                 case 'reg'
%             end
%             ic = [0; M0.extra; delC_post; 0]; 
            
        case 5
            ic = sol1.y(:,end) + [M0.extra; 0; 0; 0]; 
            
        case 10
            delC_pre = sol1.y(3,end) ;
            delC_post = (MC_pre*delC_pre + M0.extra*p.delC_in)/(MC_pre + M0.extra);            
            ic = sol1.y(:,end) + [M0.extra_P; M0.extra; delC_post; delC_post]; 
            
%         otherwise
%             ic = sol1.y(:,end);
end

switches.solvernum = 2;
sol2 = ode15s(@(t,y) odefun(t,y,k,switches,F0,M0,p) ,[t_perturbation_start, t_run_end], ic); 


%%
%Output

    Output.t1    = sol1.x;
    Output.MP1   = sol1.y(1,:);
    Output.MC1   = sol1.y(2,:);
    Output.delC1 = sol1.y(3,:);
    Output.MO1   = sol1.y(4,:);
    
    Output.t2    = sol2.x;
    Output.MP2   = sol2.y(1,:);
    Output.MC2   = sol2.y(2,:);
    Output.delC2 = sol2.y(3,:);
    Output.MO2   = sol2.y(4,:);

    Output.t  = [sol1.x sol2.x];
    Output.MP  = [sol1.y(1,:) sol2.y(1,:)] ;
    Output.MC  = [sol1.y(2,:) sol2.y(2,:)] ;
    Output.delC = [sol1.y(3,:) sol2.y(3,:)] ;
    Output.MO = [sol1.y(4,:) sol2.y(4,:)];
    
    Output.Fborg_pert1 = k0.bo.*Output.MP1;
    Output.Fborg_pert2 = k.bo.*Output.MP2;
    Output.Fborg_pert  = [Output.Fborg_pert1 Output.Fborg_pert2];
    
    switch switches.Fbcarb_switch
        case 'carb'
            Output.Fbcarb_pert1 = k0.bc.*Output.MC1;
            Output.Fbcarb_pert2 = k.bc.*Output.MC2;
            Output.Fbcarb_pert  = [Output.Fbcarb_pert1 Output.Fbcarb_pert2];
        case 'sil'
            Output.Fwsil_pert1 = k0.ws.*Output.MC1;
            Output.Fwsil_pert2 = k.ws.*Output.MC2;
            Output.Fwsil_pert  = [Output.Fwsil_pert1 Output.Fwsil_pert2];
            
            Output.Fbcarb_pert1 = Output.Fwsil_pert1 + F0.bcarb;
            Output.Fbcarb_pert2 = Output.Fwsil_pert2 + F0.bcarb;
            Output.Fbcarb_pert  = [Output.Fbcarb_pert1 Output.Fbcarb_pert2];

    end
    
    switch switches.kbp_switch
        case 'O'
            Output.Fbp_pert1 = k0.bp.*Output.MO1;
            Output.Fbp_pert2 = k.bp.*Output.MO2;
            Output.Fbp_pert  = [Output.Fbp_pert1 Output.Fbp_pert2];
        case 'P'
            Output.Fbp_pert1 = k0.bp.*Output.MP1;
            Output.Fbp_pert2 = k.bp.*Output.MP2;
            Output.Fbp_pert  = [Output.Fbp_pert1 Output.Fbp_pert2];
    end
         
    
    Output.Fwp_pert1  = k0.wp.*Output.MP1;
    Output.Fwp_pert2  = k.wp.*Output.MP2;
    Output.Fwp_pert   = [Output.Fwp_pert1 Output.Fwp_pert2];
    
    
    Output.Fworg_pert1  = k0.wo.*Output.MO1;
    Output.Fworg_pert2  = k.wo.*Output.MO2;
    Output.Fworg_pert   = [Output.Fwp_pert1 Output.Fwp_pert2];
    
    switch switches.vars_switch
        case 'dev'
            Output.Fborg  = Output.Fborg_pert + F0.borg;
            Output.Fworg  = Output.Fworg_pert + F0.worg;
            Output.Fbcarb = Output.Fbcarb_pert + F0.bcarb;
            Output.Fbp    = Output.Fbp_pert + F0.bp;
            Output.Fwp    = Output.Fwp_pert + F0.wp;
        case 'reg'
            Output.Fborg  = Output.Fborg_pert;
            Output.Fworg  = Output.Fworg_pert;
            Output.Fbcarb = Output.Fbcarb_pert;
            Output.Fbp    = Output.Fbp_pert;
            Output.Fwp    = Output.Fwp_pert;
    end

    if min(Output.Fborg)<0;
        warning('Negative Fborg flux!')
    elseif min(Output.Fbcarb)<0;
        warning('Negative Fbcarb flux!')
    elseif min(Output.Fbp)<0
        warning('Negative Fbp flux!')
    elseif min(Output.Fwp)<0
        warning('Negative Fwp flux!')
    end
    
    switch switches.per_type 
        case {0,4,5,10}
            Output.F_extra_out = zeros(length(Output.t),1);
        case {1,2,3,6,8,9,11}
            %recalculate perturbation
            Output.F_extra_out = interp1(p.tinterp, F0.extra_interp,Output.t,'pchip');
            
    end
    
    Output.CP_mean = Output.Fborg./Output.Fbp;
    Output.ALKP_mean = Output.Fbcarb./Output.Fwp;
    Output.pCO2 = 560*((Output.MC + M0.C)/M0.C).^2;

    switch switches.vars_switch
        case 'dev'
            Output.MP = Output.MP + M0.P;
            Output.MC = Output.MC + M0.C;
            Output.MO = Output.MO + M0.O;
    end
    

    end

   
%This is where you put your diff eq
function dy=odefun(t,y,k,switches,F0,M0,p)

    MP = y(1);
    MC = y(2);
    delC = y(3);
    MO = y(4);
    
    Fborg = k.bo*MP;
    Fworg = k.wo*MO;
    Fwp = k.wp*MC;
    
    switch switches.stoch
        case 'y'
            w = 1e12*randn(1);
    end
    
    switch switches.Fbcarb_switch
        case 'carb'
            Fbcarb = k.bc*MC;
        case 'sil'
            Fwsil = k.ws*MC;

    end
  
    switch switches.kbp_switch
        case 'O'
            Fbp = k.bp*MO;
        case 'P'
            Fbp = k.bp*MP;
    end
        
    
    switch switches.per_type 
        case {0,4,5,10}
            F_extra = 0;
        case {1,2,3,6,8,9,11}
            F_extra = interp1(p.tinterp, F0.extra_interp,t,'pchip');
    end

    %diff equations
    
    dy(1) = Fwp  - Fbp; 
    
    switch switches.vars_switch
        case 'dev'
            switch switches.stoch
                case 'y'
                    switch switches.Fbcarb_switch
                    case 'carb'
                            dy(2) =  F_extra - Fbcarb - Fborg + w;  
%                         case 'sil'
%                             dy(2) = F_extra +  Fworg - Fwsil - Fborg ; 
                    end
                case 'n'
                    switch switches.Fbcarb_switch
                        case 'carb'
                            dy(2) =  F_extra - Fbcarb - Fborg ;  
                        case 'sil'
                            dy(2) = F_extra +  Fworg - Fwsil - Fborg ;  
                    end
            end
        case 'reg'
            switch switches.stoch
                case 'y'
                    switch switches.Fbcarb_switch
                        case 'carb'
                                dy(2) =  F0.in + F_extra - Fbcarb - Fborg + w;  
    %                         case 'sil'
    %                             dy(2) =  F0.volc + F_extra +  Fworg - Fwsil - Fborg + w ;  
                    end
                case 'n'
                    switch switches.Fbcarb_switch
                        case 'carb'
                            dy(2) =  F0.in + F_extra - Fbcarb - Fborg ;  
                        case 'sil'
                            dy(2) =  F0.volc + F_extra +  Fworg - Fwsil - Fborg ;  
                    end
            end
    end
    
    switch switches.vars_switch
        case 'dev'
            switch switches.Fbcarb_switch
                case 'carb'
                    dy(3) = ( ...
                     (F_extra + F0.in )*(p.delC_in - delC)...
                     + (F0.borg + Fborg)*p.epsilon ...
                     ) /(MC + M0.C);    
                case 'sil'
                    dy(3) = ( ...
                     + (F_extra )*(p.delC_in - delC)...
                     + (F_extra + F0.volc_ox)*(p.delwcarb - delC)...
                     + F0.volc_red*(p.delworg - delC)...
                     + F0.wcarb*(p.delwcarb - delC)...
                     + (F0.worg + Fworg)*(p.delworg - delC) ...
                     + (F0.borg + Fborg)*p.epsilon ...
                     ) /(MC + M0.C); 
            end
        case 'reg'
            switch switches.Fbcarb_switch
                case 'carb'
                    dy(3) = ( ...
                     (F_extra + F0.in )*(p.delC_in - delC)...
                     + (Fborg)*p.epsilon ...
                     ) /(MC);    
                case 'sil'
                    dy(3) = ( ...
                     + (F_extra )*(p.delC_in - delC)...
                     + (F_extra + F0.volc_ox)*(p.delwcarb - delC)...
                     + F0.volc_red*(p.delworg - delC)...
                     + F0.wcarb*(p.delwcarb - delC)...
                     + (Fworg)*(p.delworg - delC) ...
                     + (Fborg)*p.epsilon ...
                     ) /(MC); 
            end
    end
%     

    switch switches.vars_switch
        case 'dev'
            dy(4) =  Fborg - Fworg ;
        case 'reg'
            dy(4) =  Fborg - Fworg - F0.volc_red;
    end
    
    
    
    dy = dy(:);
    
    
end    
    
    
    

 
 

