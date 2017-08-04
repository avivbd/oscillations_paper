function Figure1Plot
%% 

clear
close all
clc
%%
%Updated Jun 15th 2016

%option 1 four panels two above two below
%option 2 three panels below
flag_panels = 1;


%Save and overwrite the current copy of the figure? ['yes',~]
print_delC_fig = '~'; 

%get the data
[t, delC, gr] = getData();



%% Plot d13C data


switch flag_panels
    
    case 1
         h.fig = figure('Position', [440   190   560   500]);
        set(h.fig,'PaperPositionMode','auto');

%         h.lines = gscatter(t(1:8144),delC(1:8144),gr(1:8144),[],'o');
        h.lines = gscatter(t,delC,gr,[],'o');

        set(h.lines, 'MarkerSize',2)

        hold on

        h.ax = gca;
        set(h.ax,'XDir','reverse')
        
         xlabel('')

        ylabel('\delta^{13}C')

        legend off

        box on

        grid on

        ylim([-7 +10])
        xlim([0 541])
        
        
        
        fontsize = 13; %font size for axes notation
        linstyle = '-'; %linestyle for data boxes on main plot
        lincolor = 'b'; %Line color for data boxes on main plot
        
        %set position of main panel with all data
        MainPos = [0.1 0.375 0.85 0.3];
        set(h.ax,'Position',MainPos,'FontSize',fontsize);
        h.axp = axes('Position',get(h.ax,'Position'),...
            'Xtick',[],'YAxisLocation','right',...
            'YLim',get(h.ax,'YLim'),...
            'YTick',get(h.ax,'YTick'),...
            'Color','none','FontSize',fontsize);
        
        %Location of the four panels
        UpperLeft  = [0.1 0.75 0.4 0.2];
        UpperRight = [0.55 0.75 0.4 0.2];
        LowerRight  = [0.55 0.1 0.4 0.2];
        LowerLeft = [0.1 0.1 0.4 0.2];
        
        %Assoicated annotation A B C D
        annotation(h.fig,'textbox',...
        [UpperLeft(1)-0.01...
        UpperLeft(2)+UpperLeft(4)-0.01... 
        0.063 0.048],'String',{'A.'},...
        'LineStyle','none');
        
        annotation(h.fig,'textbox',...
        [UpperRight(1)-0.01...
        UpperRight(2)+UpperRight(4)-0.01... 
        0.063 0.048],'String',{'B.'},...
        'LineStyle','none');
    
        annotation(h.fig,'textbox',...
        [LowerLeft(1)-0.01...
        LowerLeft(2)+LowerLeft(4)-0.01... 
        0.063 0.048],'String',{'C.'},...
        'LineStyle','none');
    
        annotation(h.fig,'textbox',...
        [LowerRight(1)-0.01...
        LowerRight(2)+LowerRight(4)-0.01... 
        0.063 0.048],'String',{'D.'},...
        'LineStyle','none');
    
        
        
        
        %%%%%%%%
        
        %Panel with cambrian data
        h.ax1 = axes('Parent',h.fig,...
        'Position', UpperLeft);
    
        upln = 7.75;
        lowln = -6.5;
        leftln = 491;
        rightln = 535;
    
        lb = min(find(leftln<t & t< rightln));
        ub = max(find(leftln<t & t< rightln));

        h.l1 = plot(t(lb:ub),delC(lb:ub),'o');

        set(h.l1,'MarkerSize',2)
        set(h.ax1,'YTick',[-5 0 5 10],...
            'FontSize',fontsize,'XDir','reverse')
        
        xlim([leftln rightln])
        ylim([lowln upln])
        grid on
        ylabel('\delta^{13}C')
        
        %box
        line([leftln rightln ],[upln upln],'LineStyle',linstyle,'Color',...
            lincolor,'Parent',h.ax);
        line([leftln rightln ],[lowln lowln],'LineStyle',linstyle,'Color',...
            lincolor,'Parent',h.ax)
        line([leftln leftln ],[lowln upln],'LineStyle',linstyle,'Color',...
            lincolor,'Parent',h.ax)
        line([rightln rightln ],[lowln upln],'LineStyle',linstyle,'Color',...
            lincolor,'Parent',h.ax)
        
        %Anotation on main plot 
        text(leftln + -0.7,upln - 0.7,'A.','Parent',h.ax);
        
        %%%%%
        
        %panel with Early Paleozoic data
        h.ax2 = axes('Parent',h.fig,...
        'Position', UpperRight);
 
        upln = 9.5;
        lowln = -4.5;
        leftln = 400;
        rightln = 470;
    
        lb = min(find(leftln<t & t< rightln));
        ub = max(find(leftln<t & t< rightln));

        h.l2 = plot(t(lb:ub),delC(lb:ub),'o');
        set(h.l2,'MarkerSize',2)
        set(h.ax2,'YTick',[-5 0 5 10],'FontSize',fontsize,...
            'XDir','reverse','YAxisLocation','right')
        xlim([leftln rightln])
        ylim([lowln upln])
        grid on
        
%         xlabel('Time [Ma]')
        
        %box
        line([leftln rightln ],[upln upln],'LineStyle',linstyle,'Color',...
            lincolor,'Parent',h.ax);
        line([leftln rightln ],[lowln lowln],'LineStyle',linstyle,'Color',...
            lincolor,'Parent',h.ax)
        line([leftln leftln ],[lowln upln],'LineStyle',linstyle,'Color',...
            lincolor,'Parent',h.ax)
        line([rightln rightln ],[lowln upln],'LineStyle',linstyle,'Color',...
            lincolor,'Parent',h.ax)
        
        text(leftln -2 ,upln - 0.8,'B.','Parent',h.ax);
        
        
        
        %panel Early Triassic data
        h.ax3 = axes('Parent',h.fig,...
        'Position', LowerLeft);
 
        upln = 8.5;
        lowln = -3;
        leftln = 244;
        rightln = 253;
    
        lb = min(find(leftln<t & t< rightln));
        ub = max(find(leftln<t & t< rightln));

        h.l3 = plot(t(lb:ub),delC(lb:ub),'o');
        set(h.l3,'MarkerSize',2)
        set(h.ax3,'YTick',[-5 0 5 10],...
            'FontSize',fontsize,'XDir','reverse')
        xlim([leftln rightln])
        ylim([lowln upln])
        grid on
        xlabel('Time [Ma]')
        
        %box
         line([leftln rightln ],[upln upln],'LineStyle',linstyle,'Color',...
             lincolor,'Parent',h.ax);
         line([leftln rightln ],[lowln lowln],'LineStyle',linstyle,'Color',...
             lincolor,'Parent',h.ax)
        line([leftln leftln ],[lowln upln],'LineStyle',linstyle,'Color',...
            lincolor,'Parent',h.ax)
         line([rightln rightln ],[lowln upln],'LineStyle',linstyle,'Color',...
             lincolor,'Parent',h.ax)
        
        text(leftln + -0.7 ,upln - 0.35,'C.','Parent',h.ax);
        ylabel('\delta^{13}C')
        
        %Add arrow for end-Permian extinction
        annotation(h.fig,'arrow',[0.135 0.135],...
                                 [0.251 0.196],'Color',[1 0 0]);


         %%%%%%%%
         
        %make panel with Early Jurassic data
        h.ax4 = axes('Parent',h.fig,...
        'Position',LowerRight);        

        upln = 8;
        lowln = -6;
        leftln = 199;
        rightln = 202;
    
        lb = min(find(leftln<t & t< rightln));
        ub = max(find(leftln<t & t< rightln));
        
        h.l4 = plot(t(lb:ub),delC(lb:ub),'o');
        set(h.l4,'MarkerSize',2)
%         set(h.l4,'Color',get(h.lines(6),'Color'))
        set(h.ax4,'YTick',[-5 0 5 10],'FontSize',fontsize,...
            'XDir','reverse','YAxisLocation','right')
        xlim([leftln rightln])
        ylim([lowln upln])
        grid on
        xlabel('Time [Ma]')
        
        
        %Add arrow for end-Triassic extinction
        annotation(h.fig,'arrow',[0.603 0.603],...
                                 [0.299 0.254]...
                                    ,'Color',[1 0 0]);

        
        %box
        line([leftln rightln ],[upln upln],'LineStyle',linstyle,'Color',...
            lincolor,'Parent',h.ax);
        line([leftln rightln ],[lowln lowln],'LineStyle',linstyle,'Color',...
            lincolor,'Parent',h.ax)
        line([leftln leftln ],[lowln upln],'LineStyle',linstyle,'Color',...
            lincolor,'Parent',h.ax)
        line([rightln rightln ],[lowln upln],'LineStyle',linstyle,'Color',...
            lincolor,'Parent',h.ax)
        
        
        text(leftln - 1 ,upln - 0.6,'D.','Parent',h.ax);
        
        set(findall(h.fig,'Type','text'),'FontSize',13)

        
        
        
        
set(h.ax,'XTick',0:100:500)
    
    
    case 2
        % [left, bottom, width, height]
        % h.fig = figure('Position', [440   190   560   500]);
        %560 by 500 results in a figure which is 197.526 mm by 176.183 mm
        %want a figure that is 180 mm by 120mm so 510.3 by 340.5
        h.fig = figure('Position', [440   190   510.3   340.5]);

        set(h.fig,'PaperPositionMode','auto');

        h.lines = gscatter(t,delC,gr,[],'o');
        hold on 
        % plot(t_cramer,delC_cramer,'*k')

        set(h.lines, 'MarkerSize',2)

        hold on

        h.ax = gca;
        set(h.ax,'XDir','reverse')

        xlabel('Time [Ma]')

        ylabel('\delta^{13}C  [permil]')

        legend off

        box on

        grid on

        ylim([-7 +10])
        xlim([0 541])

        %make an inset figure for the Late Permian Early Triassic (approx 245--260)

        fontsize = 9; %font size for axes notation
        linstyle = '-'; %linestyle for data boxes on main plot
        lincolor = 'b'; %Line color for data boxes on main plot

        % [left, bottom, width, height]
        %set position of main panel with all data
        MainPos = [0.1 0.375 0.85 0.6];

        set(h.ax,'Position',MainPos,'FontSize',fontsize);

        % h.axp = axes('Position',get(h.ax,'Position'),...
        %     'Xtick',[],'YAxisLocation','right',...
        %     'YLim',get(h.ax,'YLim'),...
        %     'YTick',get(h.ax,'YTick'),...
        %     'Color','none','FontSize',fontsize);

        %Location of the three panels
        % [left, bottom, width, height]
        Left = [0.1 0.1 0.25 0.17];
        Center = [0.4 0.1 0.25 0.17];
        Right  = [0.7 0.1 0.25 0.17];








        %Assoicated annotation A B C D
        annotation(h.fig,'textbox',...
        [Left(1)-0.01...
        Left(2)+Left(4)... 
        0.063 0.048],'String',{'A.'},...
        'LineStyle','none');

        annotation(h.fig,'textbox',...
        [Center(1)-0.01...
        Center(2)+Center(4)... 
        0.063 0.048],'String',{'B.'},...
        'LineStyle','none');

        annotation(h.fig,'textbox',...
        [Right(1)-0.01...
        Right(2)+Right(4)... 
        0.063 0.048],'String',{'C.'},...
        'LineStyle','none');






        %Panel with cambrian data
        h.ax1 = axes('Parent',h.fig,...
        'Position', Left);

        upln = 8;
        lowln = -6.5;
        leftln = 520;
        rightln = 541;

        lb = min(find(leftln<t & t< rightln));
        ub = max(find(leftln<t & t< rightln));

        h.l1 = plot(t(lb:ub),delC(lb:ub),'o');
        %         h.l1 = scatter(t(lb:ub),delC(lb:ub),5,gr(lb:ub));
        set(h.l1,'MarkerSize',2)
        %         set(h.l1,'Color',get(h.lines(12),'Color'))
        %         set(h.ax1,'YAxisLocation','right','YTick',[-5 0 5 10],...
        %             'FontSize',fontsize,'XDir','reverse')
        set(h.ax1,'YTick',[-4 0 4 8],...
            'FontSize',fontsize,'XDir','reverse')

        xlim([leftln rightln])
        ylim([lowln upln])
        grid on
        % ylabel('\delta^{13}C  [permil]')

        %make a box on the Cambrian Data        
        line([leftln rightln ],[upln upln],'LineStyle',linstyle,'Color',...
            lincolor,'Parent',h.ax);
        line([leftln rightln ],[lowln lowln],'LineStyle',linstyle,'Color',...
            lincolor,'Parent',h.ax)
        line([leftln leftln ],[lowln upln],'LineStyle',linstyle,'Color',...
            lincolor,'Parent',h.ax)
        line([rightln rightln ],[lowln upln],'LineStyle',linstyle,'Color',...
            lincolor,'Parent',h.ax)

        %         xlc =  MainPos(1) + leftln/(541 - 0) * (MainPos(3));
        %         ylc =  MainPos(2) + upln/(10 + 7) * (MainPos(4));
        %         ylc = (7.75- (-7))/(10-(-7))*0.3 + MainPos(2);

        %          [0.847 0.62 0.06 0.048]        
        %         %Anotation on main plot 
        text(leftln + -0.7,upln - 0.7,'A.','Parent',h.ax);
        %         


        %         xlabel('Time [Ma]')

        %%%%%

        %panel with Early Paleozoic data
        h.ax2 = axes('Parent',h.fig,...
        'Position', Center);

        upln = 8.0;
        lowln = -4.5;
        leftln = 440;
        rightln = 465;

        lb = min(find(leftln<t & t< rightln));
        ub = max(find(leftln<t & t< rightln));

        h.l2 = plot(t(lb:ub),delC(lb:ub),'o');
        set(h.l2,'MarkerSize',2)
        set(h.ax2,'YTick',[-4 0 4 8],'FontSize',fontsize,...
            'XDir','reverse');%,'YAxisLocation','right')
        xlim([leftln rightln])
        ylim([lowln upln])
        grid on

        %         xlabel('Time [Ma]')

        %make a box on the Early Paleozoic Data
        line([leftln rightln ],[upln upln],'LineStyle',linstyle,'Color',...
            lincolor,'Parent',h.ax);
        line([leftln rightln ],[lowln lowln],'LineStyle',linstyle,'Color',...
            lincolor,'Parent',h.ax)
        line([leftln leftln ],[lowln upln],'LineStyle',linstyle,'Color',...
            lincolor,'Parent',h.ax)
        line([rightln rightln ],[lowln upln],'LineStyle',linstyle,'Color',...
            lincolor,'Parent',h.ax)

        text(leftln -2 ,upln - 0.8,'B.','Parent',h.ax);

        %         ylabel('\delta^{13}C  [permil]')


        % %panel with Mid Paleozoic data
        % h.ax3 = axes('Parent',h.fig,...
        % 'Position', LowerLeft);
        % 
        % upln = 5.2;
        % lowln = -4.0;
        % leftln = 370;
        % rightln = 385;
        % 
        % lb = min(find(leftln<t & t< rightln));
        % ub = max(find(leftln<t & t< rightln));
        % 
        % h.l3 = plot(t(lb:ub),delC(lb:ub),'o');
        % set(h.l3,'MarkerSize',2)
        % set(h.ax3,'YTick',[-5 0 5 10],...
        %     'FontSize',fontsize,'XDir','reverse')
        % xlim([leftln rightln])
        % ylim([lowln upln])
        % grid on
        % xlabel('Time [Ma]')
        % 
        % %make a box on the Mid Paleozoic Data
        % line([leftln rightln ],[upln upln],'LineStyle',linstyle,'Color',...
        %     lincolor,'Parent',h.ax);
        % line([leftln rightln ],[lowln lowln],'LineStyle',linstyle,'Color',...
        %     lincolor,'Parent',h.ax)
        % line([leftln leftln ],[lowln upln],'LineStyle',linstyle,'Color',...
        %     lincolor,'Parent',h.ax)
        % line([rightln rightln ],[lowln upln],'LineStyle',linstyle,'Color',...
        %     lincolor,'Parent',h.ax)
        % 
        % text(leftln + -0.7 ,upln - 0.35,'C.','Parent',h.ax);
        % ylabel('\delta^{13}C  [permil]')

         %make panel with Early Triassic data
        h.ax4 = axes('Parent',h.fig,...
        'Position',Right);        

        upln = 8.0;
        lowln = -1.5;
        leftln = 242;
        rightln = 254;

        lb = min(find(leftln<t & t< rightln));
        ub = max(find(leftln<t & t< rightln));

        h.l4 = plot(t(lb:ub),delC(lb:ub),'o');
        set(h.l4,'MarkerSize',2)
        %         set(h.l4,'Color',get(h.lines(6),'Color'))
        set(h.ax4,'YTick',[-4 0 4 8],'XTick',[242 246 250 254],'FontSize',fontsize,...
            'XDir','reverse');%,'YAxisLocation','right')
        xlim([leftln rightln])
        ylim([lowln upln])
        grid on
        % xlabel('Time [Ma]')


        %make a box on the Early Triassic Data
        line([leftln rightln ],[upln upln],'LineStyle',linstyle,'Color',...
            lincolor,'Parent',h.ax);
        line([leftln rightln ],[lowln lowln],'LineStyle',linstyle,'Color',...
            lincolor,'Parent',h.ax)
        line([leftln leftln ],[lowln upln],'LineStyle',linstyle,'Color',...
            lincolor,'Parent',h.ax)
        line([rightln rightln ],[lowln upln],'LineStyle',linstyle,'Color',...
            lincolor,'Parent',h.ax)

        text(leftln - 1 ,upln - 0.6,'C.','Parent',h.ax);

        %make the fonts bigger
        set(findall(h.fig,'Type','text'),'FontSize',fontsize)

        %get the ticks right
        set(h.ax,'XTick',0:100:500)

        %Add arrow for end-Permian extinction
        annotation(h.fig,'arrow',[0.737254901960783 0.736029411764704],...
            [0.249266862170088 0.194134897360704],'Color',[1 0 0]);
end
  


switch print_delC_fig
    case 'yes'
        print -dpdf -noui -loose delCPhanPanels2

end

 
            
        
end
    
    





