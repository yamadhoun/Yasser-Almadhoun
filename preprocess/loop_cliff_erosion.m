
global Kc
% global location

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
pars.d          = 2.0    ;       % see pars.d = water.d; below line 70
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
pars.Xdecay     = 0.07   ; % exponential decay coeffiecient 0.03-0.07 equally likely

%% Open and read historical nc file
ncfile = 'D:\Courses 2021-23\Module 15 - MSc research, thesis and defence\thesis_report\presentations\socalnearshore_1201to1800_edited.nc'; % nc file name
ncinfo(ncfile)                                                  % To get information about the nc file
ncdisp(ncfile)                                                  % To display nc file

%% To prepar wave and water properties for erosion models

% To read wave and water properties: 'var' existing in nc file
wave.Hs = ncread(ncfile,'Hs') ;              % wave height of transects SWH
wave.Tp = ncread(ncfile,'Tp') ;              % wave period of transects TP
wave.Dm = ncread(ncfile,'Dm') ;              % wave direction of transects MWD
wave.t = ncread(ncfile,'time') ;             % recording time of waves conditions

water.stn = ncread(ncfile,'stn') ;           % number of transects
water.lat = ncread(ncfile,'latitude') ;      % location of transects
water.long = ncread(ncfile,'longitude') ;    % location of transects
water.d = ncread(ncfile,'depth') ;           % water depths at transects

location.stn = ncread(ncfile,'stn') ;        % number of transects
location.lat = ncread(ncfile,'latitude') ;   % location of transects
location.long = ncread(ncfile,'longitude') ; % location of transects

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

timestart = wave.t(1,1);                      % the start time of the SLR curvce in [eg year]
timeend = wave.t(end,1);                      % the end time of the SLR curvce in [eg year]
SLRstart = 0.0;                               % the SLR at the start time of the SLR curvce in [m]
SLRend = 0.05;                                % the SLR at the end time of the SLR curvce in [m]
[SLRquad] = SLR_quad(timestart,timeend,SLRstart,SLRend);

SLNAVD88 = readtable('MSL_TIDE_NAVD88_dailyMax.txt');  % MSL+Tide (wrt NAVD88)
SWL = SLNAVD88.Tide_NAVD88;
% [ns, ms] = size(SLNAVD88);
% nsi = ns - rem(ns,365);
% if nsi > 0
%     SL(nsi+1:ns) = [];
% end
% SL = reshape(SL,365,[]);
% figure; plot(x.Tide_NAVD88,'color',[0.3010 0.7450 0.9330])
% xlabel('Time (days)'); ylabel('SL (m +NAVD88)');

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
Raduis = 6378137.0                                              ;  % raduis of earth in [m]
location.x = Raduis * cosd(location.lat) .* cosd(location.long) ;  % x coordinate in [m]
location.y = Raduis * cosd(location.lat) .* sind(location.long) ;  % y coordinate in [m]
location.z = Raduis * sind(location.lat)                        ;  % z coordinate in [m]

figure
hold on
box on
plot(location.lat, location.long);
for i = 1 : size(dndt_cliff_m,2) %(sieze-1)/365
    shore = location.long - cliff.dndt_yearly_R(i,:)';
    plot(location.lat,shore);
end

for k = 0 : fix(siezefull/365) %(sieze-1)/365
    lgendnames{k+1} =  sprintf('Year %d\n',k);
end
lgendnames{end+1} = sprintf('Year %6.3f\n',siezefull/365);

legend(lgendnames,'Location','eastoutside')

title(['Cliff retreat over years from ' num2str(year(wave.t(1))) ' - ', num2str(year(wave.t(end)))])
xlabel('Latitude (deg N)')
ylabel('Longitude (deg N)')
set(gca,'FontName','Calibri','FontSize',14)

%%%%%%%%%%%% draw coastal stretch segment %%%%%%%%%%%%
figure
hold on
box on
plot(location.long, location.lat, "LineWidth",2);
title('Caostline Stretch for A Coastal Segment from Stations 1201-1677')
set(gca,'fontsize',18)
xlabel('Longitude [deg]')
ylabel('Latitude [deg]')
legend('Coastline Stretch')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
hold on
box on
plot(location.lat, location.long);
for i = 1 : size(dndt_cliff_m,2) %(sieze-1)/365
    shore = location.long - cliff.dndt_yearly_cumR(i,:)';
    plot(location.lat,shore);
end

legend(lgendnames,'Location','eastoutside')

title(['Cummulaticve cliff retreat over years from ' num2str(year(wave.t(1))) ' - ', num2str(year(wave.t(end)))])
xlabel('Latitude (deg N)')
ylabel('Longitude (deg N)')
set(gca,'FontName','Calibri','FontSize',14)

toc

% %% Shoreline and beach width:----
% % method 1
% dn = 0.0;                                              % offsets in meters
% de = Beach;                                            % offsets in meters
% dLat = dn/Raduis;                                      % coordinate offsets in radians 
% dLong = de/[Raduis *cos(pi*dLat/180)];                 % coordinate offsets in radians 
% shoreloc.lat = location.lat + dLat * 180/pi;           % coordinate offsets in radians 
% shoreloc.long = location.long + dLong * 180/pi;        % coordinate offsets in radians 
% % method 2
% dx=50; dy=0;
% shoreloc.latnew = location.lat + dy/Raduis * 180/pi;           % coordinate offsets in radians 
% shoreloc.longnew = location.long + dx/Raduis * 180/(pi*cosd(location.lat*pi/180));        % coordinate offsets in radians 
% % 
% figure; hold on;
% plot(shoreloc.lat, shoreloc.long,'b');
% plot((shoreloc.lat+shoreloc.latnew)/2, shoreloc.longnew,'m');
% 
% %% latitude/longitude to cartesian
% % ref = [27, 77];         % ref point [lat, lon]
% % p = [27.2046, 77.4977];
% ref = [location.lat(1), location.long(1)];
% p = [location.lat, location.long];
% for i = 1 : size(location.lat, 1)
%     dist = distance(ref, p(i,:));
%     arclen(i) = dist(1);
%     az(i) = dist(2);
% end
% % arclen = 0.4880
% % az = 65.0993
% % x = deg2km(arclen.*cosd(az));
% % y = deg2km(arclen.*sind(az));
% % [x, y]
% % ans = 1×2
% %    22.8478   49.2198
% Raduis = 6378137.0                      ;  % raduis of earth in [m]
% x = Raduis * cosd(arclen') .* cosd(az') ;  % x coordinate in [m]
% y = Raduis * cosd(arclen') .* sind(az') ;  % y coordinate in [m]
% z = Raduis * sind(arclen')              ;  % z coordinate in [m]
% [x y z]
% figure; hold on;
% plot(location.lat, location.long);
% plot(location.lat, location.long);
% plot(shoreloc.lat, shoreloc.long,'b');
% plot((shoreloc.lat+shoreloc.latnew)/2, shoreloc.longnew,'m');
% 
% % toc
% 
% %% Calculate the cliff base line slope
% for i = 1 : length(location.lat) - 1
%     cliff.slopem(i) = ( location.y(i+1) - location.y(i) ) / ( location.x(i+1) - location.x(i) );
% end
% cliff.slopem(end+1) = cliff.slopem(end)';
% cliff.slopetan = atand(slopem)';
% cliff.Rx = cliff.dndt_yearly_R .* sin(slopetan);
% cliff.Ry = cliff.dndt_yearly_R .* cos(slopetan);
% % 
% figure
% hold on
% plot(location.x, location.y)
% plot(location.x-cliff.Rx(30,:), location.y-cliff.Ry(30,:))
% 
% toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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


