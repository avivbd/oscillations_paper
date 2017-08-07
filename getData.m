
function [t_full,delC_full,gr,t,delC,t_interp, delC_interp, delC_sinterp] = getData() 

%do you want some summary statistics and plots?
flag_dataPlots = 'no';

%% Load Data
%load previously saved data - put in same directory as this file 
load('PhanData.mat')

load('Hettangian_LBC_d13C.mat')


%% Concat data
delC_ST = delC; %rename the Saltzman and Thomas data
t_ST = t;

%remove some overlap between the the Cramer data and the Saltzman and Thomas datasets
t_cramer(1:680) = [];
delC_cramer(1:680) = [];

%the Hettangian is represented by only 15 datapoints in the S+T
%dataset. Remove them
ind_Het = find(and( t_ST<201.8,t_ST>199.03));
t_ST(ind_Het) = [];
delC_ST(ind_Het) = [];


%append the Cramer Cenozoic data and the LBC data to the Saltzman and Thomas data
delC_full = [delC_ST; delC_cramer; LBC_d13C];
t_full = [t_ST ; t_cramer; LBC_t];

[t_full,Ind] = sort(t_full);
delC_full = delC_full(Ind);

%make grouping variables for plotting by age/period
gr = MakeGroups(t_full);


%% Data Cleaning


%remove some egregious outliers from the S+T data
IndOutL = [1789 2020 2501 2502 4593 5267 5396 6389 6506 7301 7318 7319 7395 ...
          7396  7526 7281 7298 7299 5654 7375 7376 7506 7539];
      
t_ST(IndOutL) = [];
delC_ST(IndOutL) = [];


%clean up the Cramer data a bit by removing a few egregious outliers
outL = [90 451 469 717 900 914 3157 3676 4369 5125 11032 12102 12667 ...
    13332 13500 14906 16662 16771 219 220 262 379 401 791 920 1062 ...
    1075 1720 1741 2561 2810 2843];

t_cramer(outL) = [];
delC_cramer(outL) = [];

%remove some outliers from the LBC data
IndOutL = [ 33    62   293   315   419   433   440   460   504   569 415 ...
                417    468   572   642   753   786   788   791   802  832];
LBC_t(IndOutL) = [];
LBC_d13C(IndOutL) = [];


%append the cleaned data 
delC = [delC_ST; delC_cramer; LBC_d13C];
t = [t_ST ; t_cramer; LBC_t];


%remove NaN datapoints
nanind = isnan(delC);
delC(nanind) = [];
t(nanind) = [];


%replace duplicates with their average
[t_unique,~,idx] = unique(t);
delC_unique = accumarray(idx,delC,[],@mean);

%sort
[B, I] = sort(t_unique);
t = t_unique(I);
delC = delC_unique(I);


%Rates of change. No more than 1 permil per 10,000 years. Anything more
%than that is considered an artifact
n_std = 1e4;% 
  
diff_delC = diff(delC);
diff_t = diff(t);
dydx = diff_delC./diff_t; %change in permil per m.y.


% get locations of outliers
ind_diff_delC_outliers = find(or(dydx>n_std,dydx<-n_std));

% make a dataset that does not contain them
delC_no_outliers = delC;
t_no_outliers = t;
delC_no_outliers(ind_diff_delC_outliers) = [];
t_no_outliers(ind_diff_delC_outliers) = [];

delC = delC_no_outliers;
t = t_no_outliers;



%Preallocate variables for an interpolated dataset
t_pre_interp = t;
delC_pre_interp = delC;

%prior to interpolation must downsample the cenozoic which has about 
% 20x the number of points of the other time bins
IndCeno = 1:max(find(t<66.14));
t_ceno_ds = downsample(t(IndCeno),20);
delC_ceno_ds = downsample(delC(IndCeno),20);

%remove the cenozoic and add back the downsampled version
t_pre_interp(IndCeno) = [];
delC_pre_interp(IndCeno) = [];

t_pre_interp = [t_ceno_ds ; t_pre_interp];
delC_pre_interp = [delC_ceno_ds ; delC_pre_interp];

t_interp = (linspace(min(t_pre_interp),max(t_pre_interp),length(t_pre_interp) ))';
delC_interp = interp1(t,delC,t_interp);

% smoothing
[fitresult, gof] = createFit(t_pre_interp, delC_pre_interp);
delC_sinterp = fitresult(t_interp);


%%

switch flag_dataPlots
    case 'yes'
        % Plot fit with data.
        figure( 'Name', 'Smoothing Spline' );
        p = fitoptions(fitresult);
        excludedPoints = logical(p.Exclude)';
        h = plot(fitresult);
        hold on
        plot(t_full, delC_full,'.')
%         plot(t_pre_interp(excludedPoints), delC_pre_interp(excludedPoints),'x')
        legend('Smoothing spline', 'delC full', 'Location', 'NorthEast' );
        % Label axes
        xlabel t
        ylabel delC
        grid on

        
end


end


function gr = MakeGroups(t)

%assign groups for plot colors. Based on GTS 2012.
Cenozoic        = [66.14,   0];
Cretaceous      = [146.39,  66.14];
Jurassic		= [201.30,  146.39];
Triassic        = [252.16,  201.30]; 
Permian         = [298.88,  252.16];
Carboniferous   = [358.94,  298.88];
Devonian        = [419.20,  358.94];
Silurian        = [443.83,  419.20];
Ordovician      = [485.37,  443.83];
Cambrian        = [541,     485.37];


gr = zeros(size(t));


for i = 1:length(t)
    
    
    if and(Cenozoic(1)>=t(i),t(i)>Cenozoic(2))
        gr(i) = 1;
    end 
    
    if and(Cretaceous(1)>=t(i),t(i)>Cretaceous(2))
        gr(i) = 2;
    end 
    
    if and(Jurassic(1)>=t(i),t(i)>Jurassic(2))
        gr(i) = 3;
    end 
    
    if and(Triassic(1)>=t(i),t(i)>Triassic(2))
        gr(i) = 4;
    end
    
    if and(Permian(1)>=t(i),t(i)>Permian(2))
        gr(i) = 5;
    end
    
    if and(Carboniferous(1)>=t(i),t(i)>Carboniferous(2))
        gr(i) = 6;
    end
    
    if and(Devonian(1)>=t(i),t(i)>Devonian(2))
        gr(i) = 7;
    end
    
    if and(Silurian(1)>=t(i),t(i)>Silurian(2))
        gr(i) = 8;
    end
    
    if and(Ordovician(1)>=t(i),t(i)>Ordovician(2))
        gr(i) = 9;
    end
    
    if and(Cambrian(1)>=t(i),t(i)>Cambrian(2))
        gr(i) = 10;
    end
    
    
end


end


function dataPlots(t_full,delC_full,gr,t,delC,t_interp,delC_interp)

gscatter(t_full,delC_full,gr)
hold on
% line(t,delC)
line(t_interp,delC_interp)


end


function [fitresult, gof] = createFit(t_pre_interp, delC_pre_interp)
%% Fit: 'untitled fit 1'.
[xData, yData] = prepareCurveData( t_pre_interp, delC_pre_interp );

ft = fittype( 'smoothingspline' );
opts = fitoptions( 'Method', 'SmoothingSpline' );
opts.SmoothingParam = 0.999999994163153;
opts.Normalize = 'on';

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );


end

