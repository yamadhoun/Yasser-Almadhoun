close all; clear all; clc;

%% Open and read historical nc file
addpath(genpath('D:\Courses 2021-23\Module 15 - MSc research, thesis and defence\chapter_3_methodology\test_15dec2022\Shorelines_20221117_test_model\script\Limber_2018_model_input/\'))
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

wave.Hs = wave.Hs(:, 1:length(water.d));
wave.Tp = wave.Tp(:, 1:length(water.d));
wave.Dm = wave.Dm(:, 1:length(water.d)); 

location.stn = nonzeros(location.stn);
location.lat = nonzeros(location.lat);
location.long = nonzeros(location.long);

%% plot on world map
latitude = location.lat;
longitude = location.long;
transect = location.stn;

%% historical retreat
tab = readtable('historical_retreat.txt');
histo.stn = tab.Var1 ;
histo.ret = tab.Var2 ;
histo.long = tab.Var4 ;
histo.lat = tab.Var5 ;
for i = 1:length(histo.stn)
    j = find(histo.stn(i)==transect);
    if ~isempty(j)
        short.Hs(:,i) = wave.Hs(:,j);
        short.Tp(:,i) = wave.Tp(:,j);
        short.Dm(:,i) = wave.Dm(:,j);
        short.d(i) = water.d(j);
        short.stn(:,i) = histo.stn(i);
        short.ret(:,i) = histo.ret(i);
        short.lat(:,i) = histo.lat(i);
        short.long(:,i) = histo.long(i);
    end
end

%% plot transects on world map
% wm = webmap('World Imagery');
% s = geoshape(short.lat,short.long);
% wmline(s,'Color', 'blue', 'Width', 3);
% wmmarker(wm,short.lat,short.long);

%% max Hs and mean Tp daily
[n, m] = size(short.Hs);
ni = n - rem(n,8);
trh = wave.t;
if ni > 0
    trh(ni+1:n) = [];
end
trh = reshape(trh,8,[]);

%% Read paramsx txt file: cliff and coast proporties
fileID = fopen('paramsx.txt');
C = textscan(fileID,'%s %s %f32');
fclose(fileID);
for i=1:size(C{1},1)
    pars.(string(C{1,1}(i))) = C{1,3}(i);
end

%% read SWL
SLNAVD88 = readtable('MSL_TIDE_NAVD88_dailyMax.txt');                      % MSL+Tide (wrt NAVD88)
SWL = SLNAVD88.Tide_NAVD88;

%% per each transect
% k = 24  
count=0;
figure(1);set(gcf, 'Position', get(0, 'Screensize')); hold on; box on;
for k = 1 : m
    pars.d = water.d(k);
    Hsh = short.Hs(:,k);
    Tph = short.Tp(:,k);
    Dmh = short.Dm(:,k);
    
    n = length(Hsh);
    nh = n - rem(n,8);
    if nh > 0
        Hsh(nh+1:n) = [];
        Tph(nh+1:n) = [];
        Dmh(nh+1:n) = [];
    end
    Hsh = reshape(Hsh,8,[]);
    Tph = reshape(Tph,8,[]);
    Dmh = reshape(Dmh,8,[]);

    % max Hs in the day
    for i = 1:size(Hsh,2)
        Hsd(i) = max(Hsh(:,i));
        Tpd(i) = mean(Tph(:,i));
        Dmd(i) = mean(Dmh(:,i));
    end
    Hs = Hsd';
    Tp = Tpd';
    Dm = Dmd';
    td = linspace(wave.t(1),wave.t(1)+length(Hs),length(Hs));
    
    SL = SWL(1:length(Hs));

    for j=1:size(Hs,1)
        if isnan(Hs(j))==1;
            Hs(j)=min(Hs);
            Tp(j)=min(Tp);
            Dm(j)=min(Dm);
        end
    end

    [dndt_cliff,Om,TWL] = erosion_model(Hs,Tp,SL,pars);
    
    %% daily results
    format shortg
    out_cliff_retreat = [ dndt_cliff, Om, TWL, Hs, Tp];
    
    %% from daily to yearly
    ni = length(TWL) - rem(length(TWL),365);
    if ni > 0
        Hs(ni+1:end) = [];
        Tp(ni+1:end) = [];
        TWL(ni+1:end) = [];
        Om(ni+1:end) = [];
        dndt_cliff(ni+1:end) = [];
    end
    Hsd = reshape(Hs,365,[]);
    Tpd = reshape(Tp,365,[]);
    TWLd = reshape(TWL,365,[]);
    Omd = reshape(Om,365,[]);
    dndt_cliffd = reshape(dndt_cliff,365,[]);

    % averaged over the year
    for i = 1:size(TWLd,2)
        Hsy(i) = mean(Hsd(:,i));
        Tpy(i) = mean(Tpd(:,i));
        TWLy(i) = mean(TWLd(:,i));
        Omy(i) = mean(Omd(:,i));
        dndt_cliffy(i) = mean(dndt_cliffd(:,i));
    end
    Hs=[]; Tp=[]; TWL=[]; Om=[]; dndt=[];
    Hs = Hsy';
    Tp = Tpy';
    TWL = TWLy';
    Om = Omy';
    dndt = dndt_cliffy';

    %% plot all transects over 30 years
    ty=1:30;
    figure(1); hold on;
    
    subplot(5,1,1); hold on; box on; grid on;
    plot(ty,Hs,'-s');
    % legend('Hs')
    % xlabel('Time [year]')
    ylabel('Hs [m]')
    % title('Hs over Time Per Transect')
    set(gca,'fontsize',10)

    subplot(5,1,2); hold on; box on; grid on;
    plot(ty,Tp,'-*');
    % legend('Tp')
    % xlabel('Time [year]')
    ylabel('Tp [s]')
    % title('Tp over Time Per Transect')
    set(gca,'fontsize',10)

    subplot(5,1,3); hold on; box on; grid on;
    plot(ty,TWL,'-s');
    % legend('TWL')
    % xlabel('Time [year]')
    ylabel('TWL [m]')
    % title('TWL over Time Per Transect')
    set(gca,'fontsize',10)

    subplot(5,1,4); hold on; box on; grid on;
    plot(ty,Om,'-*');
    % legend('Om')
    % xlabel('Time [year]')
    ylabel('Om [Forcing]')
    % title('Om over Time Per Transect')
    set(gca,'fontsize',10)

    subplot(5,1,5); hold on; box on; grid on;
    plot(ty,dndt,'-o'); 
    % legend('dndt')
    xlabel('Time [year]')
    ylabel('dndt [m/yr]')
    % title('dndt over Time Per Transect')
    set(gca,'fontsize',10)
    Hsh=[]; Hsd=[];
    Tph=[]; Tpd=[];
    
    %% results for mean retreat per year per transect
    result(k,:) = [short.stn(k) mean(TWL) mean(Om) mean(dndt) mean(Hs) mean(Tp)];

    %% count
    count=count+1;

end

%% figure for mean retreat per year per transect: calculated
% averaged over the treansects
figure(2); set(gcf, 'Position', get(0, 'Screensize')); subplot(2,1,1);
hold on; box on; %title('Calc dndt vs Calc Om')
yyaxis right % Om
plot(result(:,1),result(:,3),'o-')
ylabel('Om [Forcing]')
yyaxis left % dndt
plot(result(:,1),result(:,4),'*--')
ylabel('Calc dndt [m/yr]')
xlabel('Transect ID [-]')
set(gca,'FontSize',18)

%% figure for mean retreat per year per transect: observed
subplot(2,1,2);
hold on; box on; %title('Obs dndt vs Calc Om')
yyaxis right % Om
plot(result(:,1),result(:,3),'o-')
ylabel('Om [Forcing]')
yyaxis left % dndt
plot(short.stn,short.ret,'--*k')
ylabel('Obs dndt [m/yr]')
xlabel('Transect ID [-]')
set(gca,'FontSize',18)

%% figure for mean retreat per year per transect: observed
figure(3);hold on; box on; title('Obs dndt vs Calc dndt')
plot(result(:,1),result(:,4),'-*')
plot(short.stn,short.ret,'--*k')
ylabel('dndt [m/yr]')
xlabel('Transect ID [-]')
legend('Calc dndt','Obs dndt')
set(gca,'FontSize',18)

%% calculate Kc
% Kc for TWL/Tren/EF/Hack models
short.Kc = short.ret ./ result(:,3)';
% K for Zhang mdoel
if pars.form == 1
    Fr = pars.sigma .* pars.Bc;
    Fw = Om;
    for i = 1 : length(Fw)
        if Fw(i) >= Fr
            short.K = short.ret / log(Fw(i)/Fr);
        else
            short.K = 0;
        end
    end
end

%% figure Om vs dndt: Kc is the slope
figure(4);hold on; box on; title('Calc Om vs Obs dndt')
plot(result(:,3),short.ret,'o','MarkerFaceColor',[0 0.447 0.741], 'MarkerSize', 5);
ylabel('Obs dndt [m/yr]')
xlabel('Calc Om [Forcing]')
set(gca,'FontSize',18)
% Linear regression through the origin (0,0)
Kc = result(:,3)\short.ret(:);                  
yfit = result(:,3)*Kc;
SStot = sum((short.ret-mean(short.ret)).^2);            % Total Sum-Of-Squares
SSres = sum((short.ret(:)-yfit(:)).^2);                 % Residual Sum-Of-Squares
% Regression coefficient R2
Rsq = 1-SSres/SStot                                     % R^2 or Rsq= 1-var(short.ret-yfit')/var(short.ret);
plot(result(:,3), yfit, '-r')   
legend({'Obs dndt','Fitted dndt through (0,0)'}, 'Location', 'Northwest')
hold off; grid on
axis([0  max(result(:,3))    0  max(short.ret)])

%% stats for Kc
stats_Kc = {'Mean','Min','Max','Range'}
stats_Kc = [mean(short.Kc) min(short.Kc) max(short.Kc) max(short.Kc)-min(short.Kc)]
% stats_Kc = [mean(obs.Kc(11:17)) min(obs.Kc(11:17)) max(obs.Kc(11:17)) max(obs.Kc(11:17))-min(obs.Kc(11:17))]

res.Bias = result(:,4) - short.ret;
res.RMS = rms(result(:,4));
res.RMSE = rmse(result(:,4),short.ret); 
res.error = abs(res.Bias ./ short.ret) * 100;

%% basic fitting: linear regression
coefficients = polyfit(short.ret, result(:,4), 1);
% res.fitted_ret = polyval(coefficients, short.ret);
figure(5); hold on; box on; grid on;
plot(short.ret, result(:,4),'o','MarkerFaceColor',[0 0.447 0.741], 'MarkerSize', 5);
% plot(short.ret, res.fitted_ret, 'rs-', 'LineWidth', 1.25, 'MarkerSize', 3);
% e = std(res.fitted_ret) * ones(size(short.ret));
% errorbar(short.ret,res.fitted_ret,e, 'color', 'k')
% legend('Obs dndt', 'Fitted dndt','location','best');
title('Linear Regression: Obs vs Calc dndt')
xlabel('Obs dndt 1930-2010 [m/yr]')
ylabel(['Calc dndt ' num2str(year(wave.t(1))) '-' num2str(year(wave.t(end))) ' [m/yr]'])
set(gca,'FontSize',18); axis tight; 

