
tic
clear all
close all
clc
% global Kc
% run 'matlab\oetsettings.m'
% 
% Kzz = [0.005 : 0.002 : 0.05];
% kcc = [0.005: 0.005 : 1.000];
% Kcc = [0.05: 0.05 : 1.00];
% for ik = 1 : length(kcc)
%     Kc = kcc(ik);
% calculate cliff erosion from model option
run loop_cliff_erosion.m

% historical_retreat_lookup.m

% Open and read location lat/log from nc file
% ncfile = 'Limber_2018_model_input\hindcast_time_series\socalnearshore_1201to1800.nc'; % nc file name
% ncinfo(ncfile)                                                                         % To get information about the nc file
% ncdisp(ncfile)                                                                         % To display nc file
% location.stn = ncread(ncfile,'stn') ;                                                  % number of transects
% location.lat = ncread(ncfile,'latitude') ;                                             % location of transects
% location.long = ncread(ncfile,'longitude') ;                                           % location of transects

% Open and read historical retreat fro all sections from txt file
tab = readtable('historical_retreat.txt');
histo.stn = tab.Var1 ;
histo.ret = tab.Var2 ;
histo.tag = tab.Var3 ;
histo.lat = tab.Var4 ;
histo.long = tab.Var5;

% lookup and sort historical retreat for certain section from the above
for i = 1: length(location.stn_new)
%     lookstation = location.stn(i);
    for j = 1: length(histo.stn)
        if histo.stn(j) == location.stn_new(i)
            loc.stn(i) = histo.stn(j);
            loc.ret(i) = histo.ret(j);
            loc.lat(i) = histo.lat(j);
            loc.long(i) = histo.long(j);
%             loc.tag(i) = histo.tag(j);
        end
    end
end

% delete zeros data (NA elements)
loc.stn = nonzeros(loc.stn);
loc.ret = nonzeros(loc.ret);
loc.lat = nonzeros(loc.lat);
loc.long = nonzeros(loc.long);
% loc.tag = nonzeros(loc.tag);
%%
% loc.stn = loc.stn(1:length(loc.ret)); % re-do... delete only whats missing vs stn
% loc.lat = loc.lat(1:length(loc.ret));
% loc.long = loc.long(1:length(loc.ret));

% print sorted data to text file
fileID = fopen('historical_retreat_sorted.txt','w');

fprintf(fileID, 'Date and Time (Local): %12s\n\n', datestr(now));
fprintf(fileID,'Historical retreat data for a given section/stations at a given year:\n\n');
fprintf(fileID,'%12s %12s %12s %12s\n','Stations [#]','Retreat [m]','Latitude [deg]','Longitude [deg]');

% loc.ret(end+1) = loc.ret(end);
txtmat = [loc.stn, loc.ret, loc.lat, loc.long];
for i = 1 : size(txtmat,1)
    fprintf(fileID,'%6g\t',txtmat(i,:));
    fprintf(fileID,'\n');
end

% another option: writematrix(txtmat, "historical_retreat_sorted.txt");
fclose(fileID);

% plots retreats
% Historical sea cliff retreat rates derived from 1930s and 2010 cliff edge positions were used to calibrate
% the cliff retreat models (Figure 5) in Limber paper.
figure; semilogy(histo.stn, histo.ret,'o','MarkerFaceColor', ...
                 [0 0.447 0.741], 'MarkerSize', 5);
ylim([10^-4 10^1]);
legend('Historical retreat rates 1930-2010','location','best');
title(' Historical cliff retreat rates within the study area between 1930s and 2010')
xlabel('Transect ID')
ylabel('Historical cliff retreat rate [m/yr]')
set(gca,'FontSize',18);
box on;

% get statistics bet predicted and observed
for i = 1: length(loc.stn)
%     lookstation = loc.stn(i);
    for j = 1: length(location.stn_new)
        if loc.stn(i) == location.stn_new(j)
            loc.cal(i,1) = cliff.dndt_yearly_R(end,j);
            loc.Om(i,1) = outcliffdaycum(j,2);
        end
    end
end
% loc.cal = loc.cal(1:length(loc.ret)); % re-do... delete only whats missing vs stn
% loc.Om = loc.Om(1:length(loc.ret));


loc.Bias = loc.cal - loc.ret;
loc.RMS_ret = rms(loc.ret);
loc.RMS_cal = rms(loc.cal);
loc.RMSE = rmse(loc.ret , loc.cal); % run oetsettings first
loc.error = abs(loc.Bias ./ loc.ret) * 100;
% Willmott 1981 ( 0 <= d <= 1  --> 0 for poor and 1 for excellent)
% loc.d = 1 - sum(loc.ret - loc.cal)^2 / sum((loc.cal - mean(loc.ret)) + (loc.ret - mean(loc.ret)).^2);

% linear regression
coefficients = polyfit(loc.ret, loc.cal, 1);
loc.fitted_ret = polyval(coefficients, loc.ret);
% R2_1 = rsquare(loc.ret,loc.fitted_ret)
% Rsq1 = 1 - sum((y - yCalc1).^2)/sum((y - mean(y)).^2)

loc.Kc_calibration = loc.ret(:)\loc.Om(:);
loc.fitted_slope = loc.ret(:)\loc.cal(:);
loc.fitted_ret_interc = loc.fitted_slope * loc.ret;
% R2_0 = rsquare(loc.ret,loc.fitted_ret_interc)

figure; hold on;
plot(loc.ret, loc.fitted_ret, 'rs-', 'LineWidth', 3, 'MarkerSize', 8);
plot(loc.ret, loc.cal,'o','MarkerFaceColor',[0 0.447 0.741], 'MarkerSize', 8);
e = std(loc.fitted_ret) * ones(size(loc.ret));
errorbar(loc.ret,loc.fitted_ret,e, 'color', 'k')
legend('Observed Retreat', 'Fitted Retreat','location','best');
title('Observed Verses Fitted Retreate Rates: Linear Regression')
xlabel('Observed Retreate Rate 1930-2010 [m/yr]')
ylabel(['Fitted Retreate Rate  ' num2str(year(wave.t(1))) '-' num2str(year(wave.t(end))) ' [m/yr]'])
set(gca,'FontSize',18);
axis tight;
box on
% xlim([0 roundup(max(loc.ret),0.1)]);
% ylim([0 roundup(max(loc.cal),0.1)]);
grid on;
box on;
hold off;

% plot retreat predicted and historical
figure; hold on; 
plot(loc.stn, loc.cal,'or'); 
plot(loc.stn, loc.cal,'-b'); 
plot(loc.stn, loc.ret,'*k'); 
legend(['Predicted ' num2str(year(wave.t(1))) '-' num2str(year(wave.t(end)))],['Predicted ' num2str(year(wave.t(1))) '-' num2str(year(wave.t(end)))],'Historical 1930-2010','location','best');
title('Observed Versus Predicted Retreate Rates for All Transects')
xlabel('Transect ID [-]')
ylabel('Retreat rate [m/yr]')
set(gca,'FontSize',18);
axis tight;
grid on;
box on;
hold off;

figure; hold on; 
plot(loc.stn, loc.cal,'or');  % loc.ret ./ loc.Om .*loc.Om
plot(loc.stn, loc.ret,'*k'); 
legend(['Predicted ' num2str(year(wave.t(1))) '-' num2str(year(wave.t(end)))],'Historical 1930-2010','location','best');
title('Observed Versus Predicted Retreate Rates for All Transects')
xlabel('Transect ID')
ylabel('Retreat rate [m/yr]')
set(gca,'FontSize',18);
% set(gca,'YTickLabel',[0:0.2:1.6]);
axis tight;
grid on;
box on;
hold off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculating the regression coefficient
% Compute the mean of your data
yhat = mean(loc.cal);
xhat = mean(loc.ret);
% Compute regression coefficients in the least-square sense
b = (loc.ret - xhat).*(loc.cal' - yhat)/sum((loc.ret - xhat).^2);   % Regression coefficient
a = yhat - b*xhat;                                                  % Y-intercept
% Plot data
% figure()
% scatter(loc.ret, loc.cal);
% hold on
% plot(loc.ret, a + b.*loc.ret, 'r');
% legend(['b = ' num2str(b)]); % Your regression coefficient is b

% [ik Kc]
% loc.StatsBox(ik,1:length(abs(loc.Bias))) = abs(loc.Bias);

% end


%% plot histograms and pdfs
figure; hist(wave.Hs); col=colorbar; ylabel(col,'Number of Transects'); title('SWH Histogram at All Transects'); xlabel('SWH [m]'); ylabel('Frequency [-]'); set(gca,'FontSize',18);
figure; hist(wave.Hs(:,1)); title('SWH Histogram at Transect 1201'); xlabel('SWH [m]'); ylabel('Frequency [-]'); set(gca,'FontSize',18);
figure; hist(wave.Tp(:,1)); title('WPP Histogram at Transect 1201'); xlabel('WPP [m]'); ylabel('Frequency [-]'); set(gca,'FontSize',18);
figure; hist(wave.Dm(:,1)); title('MWD Histogram at Transect 1201'); xlabel('MWD [m]'); ylabel('Frequency [-]'); set(gca,'FontSize',18);
figure; hist(SWL); title('SWL Histogram for the Coastal Segment 1201-1677'); xlabel('SWL [m]'); ylabel('Frequency [-]'); set(gca,'FontSize',18);
% MSLTIDE=readtable('MSL_TIDE_NAVD88_dailyMax.txt');
% figure; plot(MSLTIDE.Tide_NAVD88,'color',[0.3010 0.7450 0.9330]);

toc

% close all

%% SWH at transect 1201
figure; 
box on
plot(wave.t,wave.Hs(:,1),'color','b','linewidth',1.25)
title('Significant Wave Heights Over time from 1980-2010 at Station 1201')
ylabel('SWH [m]')
xlabel('Time [Years]')
xlim([min(wave.t) max(wave.t)])
set(gca,'fontsize',18)
datetick

%% peak period at transect 1201
figure; 
box on
plot(wave.t,wave.Tp(:,1),'color','b','linewidth',1.25)
title('Wave Peak Period Over time from 1980-2010 at Station 1201')
ylabel('WPP [s]')
xlabel('Time [Years]')
xlim([min(wave.t) max(wave.t)])
axis tight 
ylim([0 28])
set(gca,'fontsize',18)
datetick

%% MEAN DIRECTION at transect 1201
figure; 
box on
plot(wave.t,wave.Dm(:,1),'color','b','linewidth',1.25)
title('Mean Wave Direction Over time from 1980-2010 at Station 1201')
ylabel('MWD [s]')
xlabel('Time [Years]')
xlim([min(wave.t) max(wave.t)])
xticks([0:30:400])
set(gca,'fontsize',18)
axis tight
datetick
