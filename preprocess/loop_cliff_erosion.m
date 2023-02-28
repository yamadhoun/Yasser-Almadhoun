%
global Kc
% global location
%
tic
%
%
% LOOP : CLIFF_EROSION
%
% LOOP : Cliff erosion rate according to various formulations
% close all;
% clear all;
% clc;

%% Cliff erosion model option:
% TWL Kc=0.510; Hackney Kc=0.012; Trenhaile Kc=0.02/2.50; Energy_Flux Kc=0.20/0.10; Zhang Kc=0.035
pars.roh        = 1025   ;          
pars.g          = 9.81   ;      
pars.d          = 2.0    ;                                                 % see pars.d = water.d; below line 70
pars.gamma      = 0.78   ;        

pars.alpha      = 0.100  ;         
pars.Beta       = 0.025  ;         
pars.theta      = 135    ;
pars.Beach      = 50     ;

pars.Ej         = 3.5    ;        
pars.Hc         = 25     ;       
pars.sigma      = 1.0    ;                 
pars.Bc         = 0.1    ;         
pars.fc         = 0.65   ;  

pars.form       = 'TWL'  ;
% pars.Kc         = Kc ;       
% pars.Kz         = 0.035  ;      
pars.Xdecay     = 0.07   ;                                                 % exponential decay coeffiecient 0.03-0.07 equally likely

%% Open and read historical nc file
addpath(genpath('..\Limber_2018_model_input\hindcast_time_series\'))
ncfile = 'socalnearshore_1201to1800.nc';                                   % nc file name  
ncinfo(ncfile)                                                             % To get information about the nc file
ncdisp(ncfile)                                                             % To display nc file

%% To prepar wave and water properties for erosion models

% To read wave and water properties: 'var' existing in nc file
wave.Hs = ncread(ncfile,'Hs') ;                                            % wave height of transects SWH
wave.Tp = ncread(ncfile,'Tp') ;                                            % wave period of transects TP
wave.Dm = ncread(ncfile,'Dm') ;                                            % wave direction of transects MWD
wave.t = ncread(ncfile,'time') ;                                           % recording time of waves conditions

water.stn = ncread(ncfile,'stn') ;                                         % number of transects
water.lat = ncread(ncfile,'latitude') ;                                    % location of transects
water.long = ncread(ncfile,'longitude') ;                                  % location of transects
water.d = ncread(ncfile,'depth') ;                                         % water depths at transects

location.stn = ncread(ncfile,'stn') ;                                      % number of transects
location.lat = ncread(ncfile,'latitude') ;                                 % location of transects
location.long = ncread(ncfile,'longitude') ;                               % location of transects

%% special for hindcast data from nc file socalnearshore_1201to1800.nc
water.long = nonzeros(water.long);
water.lat = nonzeros(water.lat);
water.d = nonzeros(water.d);
water.stn = nonzeros(water.stn);

location.lat = nonzeros(location.lat);
location.long = nonzeros(location.long);
% location.stn = nonzeros(location.stn);
for i=1:length(water.d)
    if water.d(i)~=0 && isnan(water.d(i))~=1
        location.stn_new(i)=location.stn(i);
    end
end

wave.Hs = wave.Hs(:, 1:length(water.d));
wave.Tp = wave.Tp(:, 1:length(water.d));
wave.Dm = wave.Dm(:, 1:length(water.d)); 
% end of above section

[n, m] = size(wave.Hs);
ni = n - rem(n,8);
trh = wave.t;
if ni > 0
    trh(ni+1:n) = [];
end
trh = reshape(trh,8,[]);

timestart = wave.t(1,1);                                                   % the start time of the SLR curvce in [eg year]
timeend = wave.t(end,1);                                                   % the end time of the SLR curvce in [eg year]
SLRstart = 0.0;                                                            % the SLR at the start time of the SLR curvce in [m]
SLRend = 0.05;                                                             % the SLR at the end time of the SLR curvce in [m]
[SLRquad] = SLR_quad(timestart,timeend,SLRstart,SLRend);

SLNAVD88 = readtable('MSL_TIDE_NAVD88_dailyMax.txt');                      % MSL+Tide (wrt NAVD88)
SWL = SLNAVD88.Tide_NAVD88;

%% get Kc calibration values for the cliff erosion model 
% [Kc] = Kc_calibration(stn);
% read Kc_calibration.txt text file
Kc_calibrate = readtable("Kc_calibration_table.txt");
% assign Kc calibration values to location.stn
for i = 1 : length(location.stn)
    for j = 1 : length(Kc_calibrate.Hindcast_stn)
        if location.stn(i) == Kc_calibrate.Hindcast_stn(j)
            switch pars.form
                case 'TWL'
                    Kc(i,1) = Kc_calibrate.Kc_TWL(j);
                case 'Hackney'
                    Kc(i,1) = Kc_calibrate.Kc_Hack(j);
                case 'Trenhaile'
                    Kc(i,1) = Kc_calibrate.Kc_Tren(j);
                case 'Energy_Flux'
                    Kc(i,1) = Kc_calibrate.Kc_EF(j);
                case 'Zhang'
                    Kz(i,1) = Kc_calibrate.Kz_Zhang(j);
            end
        end
    end
end

% % % pars.Kc=0.510;
% % % % pars.Kc=Kc;
%% loop through calculations 
for k = 1 : m
    switch pars.form
        case 'Zhang'
            pars.Kz = Kz(k);
        otherwise
            pars.Kc = Kc(k);
    end
    
    pars.d = water.d(k);
    Hsh = wave.Hs(:,k);
    Tph = wave.Tp(:,k);
    Dmh = wave.Dm(:,k);
    if ni > 0
        Hsh(ni+1:n) = [];
        Tph(ni+1:n) = [];
        Dmh(ni+1:n) = [];
    end
    Hsh = reshape(Hsh,8,[]);
    Tph = reshape(Tph,8,[]);
    Dmh = reshape(Dmh,8,[]);

    for i = 1:size(Hsh,2)
        Hsd(i) = max(Hsh(:,i));
        Tpd(i) = mean(Tph(:,i));
        Dmd(i) = mean(Dmh(:,i));
    end
    Hs = Hsd';
    Tp = Tpd';
    Dm = Dmd';
    SL = SWL(size(Hs,1),1); %+ SLRquad(1:size(Hs,1),2);

    for j=1:size(Hs,1)
        if isnan(Hs(j))==1;
            Hs(j)=min(Hs);
            Tp(j)=min(Tp);
            Dm(j)=min(Dm);
        end
    end
    % Hs(end,1)
    % Tp(end,1)
    [dndt_cliff, Om, TWL] = cliff_erosion(Hs, Tp, SL, pars);

    format shortg
    outcliffday      = [ dndt_cliff, Om, TWL ];
    outcliffdaycum   = [ cumsum(dndt_cliff), Om, TWL ];

    sieze = size(outcliffday,1) - rem(size(outcliffday,1),365) + 1 ;
    siezefull = size(outcliffday,1);
    siezefix = fix(sieze) - rem(fix(sieze), 365);
    siezerem = siezefull - siezefix;
    %% start from here.....................................................
    dndt_cliff(end+1:siezefull+(365-siezerem),:) = zeros;
    dndt_cliff_m = reshape(dndt_cliff,365,[]);

    for i = 1:size(dndt_cliff_m,2)
        dndt_cliff_yr(i) = max(dndt_cliff_m(:,i));
    end
    
    cliff.dndt_yearly_R(:,k) = dndt_cliff_yr';
    cliff.dndt_yearly_VR(:,k) = dndt_cliff_yr' *pars.Hc *1;                         % per 1-m run alongshore
    cliff.dndt_yearly_VS(:,k) = dndt_cliff_yr' *pars.Hc *1 * (1-pars.fc);           % per 1-m run alongshore to ShorelineS
    cliff.dndt_yearly_cumR(:,k) = cumsum(dndt_cliff_yr');
    cliff.dndt_yearly_cumVR(:,k) = cumsum(dndt_cliff_yr') *pars.Hc *1;              % per 1-m run alongshore
    cliff.dndt_yearly_cumVS(:,k) = cumsum(dndt_cliff_yr') *pars.Hc *1 *(1-pars.fc); % per 1-m run alongshore to ShorelineS
end

% GRS 1980 to Cartesian
Raduis = 6378137.0                                              ;          % raduis of earth in [m]
location.x = Raduis * cosd(location.lat) .* cosd(location.long) ;          % x coordinate in [m]
location.y = Raduis * cosd(location.lat) .* sind(location.long) ;          % y coordinate in [m]
location.z = Raduis * sind(location.lat)                        ;          % z coordinate in [m]


%%%%%%%%%%%%%%%%%%%%%%%%% draw coastal stretch segment %%%%%%%%%%%%%%%%%%%%
figure
hold on
box on
plot(location.long, location.lat, "LineWidth",2);
title('Caostline Stretch for A Coastal Segment from Stations 1201-1677')
set(gca,'fontsize',18)
xlabel('Longitude [deg]')
ylabel('Latitude [deg]')
legend('Coastline Stretch')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% sea level rise from start to end with default values
figure
plot(SLRquad(:,1), SLRquad(:,2), 'b', 'linewidth', 2);
legend('SLR','Location','best')
title(['Sea Level Rise Over Time from ' num2str(year(wave.t(1))) ' to ', num2str(year(wave.t(end)))])
xlabel('Time [Years]')
ylabel('SLR [m]')
set(gca,'FontSize',18)
box on
xlim([min(SLRquad(:,1)) max(SLRquad(:,1))])
datetick

%% sea water level over time
figure
plot(SLRquad(1:length(SWL),1), SWL, 'b', 'linewidth', 1.25);
legend('SWL','Location','best')
title(['Sea Water Level Over Time from ' num2str(year(wave.t(1))) ' - ', num2str(year(wave.t(end)))])
xlabel('Time (Hours)')
ylabel('SWL (m)')
set(gca,'FontSize',18)
box on
xlim([min(SLRquad(1:length(SWL),1)) max(SLRquad(1:length(SWL),1))])
datetick

%% nearshore water depth at stations
figure
plot(water.stn(1:length(water.d)), water.d, 'b', 'linewidth', 1.25);
legend('Water Depth Nearshore','Location','best')
title(['Water Depth Nearshore at Stations from ' water.stn(1) ' to ', water.stn(end)])
xlabel('Stations [-]')
ylabel('Water Depth [m]')
set(gca,'FontSize',18)
box on
xlim([min(water.stn(1:length(water.d))) max(water.stn(1:length(water.d)))])


toc
