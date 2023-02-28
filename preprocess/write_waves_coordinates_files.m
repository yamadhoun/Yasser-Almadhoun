
function write_waves_coordinates_files()

%% Open and read historical nc file

addpath(genpath('..\Limber_2018_model_input\historical_time_series\'))

ncfile = 'historicalgfdl_1201to1800.nc';     % nc file name  
ncinfo(ncfile)                               % To get information about the nc file
ncdisp(ncfile)                               % To display nc file

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

location.stn = nonzeros(location.stn);
location.lat = nonzeros(location.lat);
location.long = nonzeros(location.long);

wave.Hs = wave.Hs(:, 1:length(water.d));
wave.Tp = wave.Tp(:, 1:length(water.d));
wave.Dm = wave.Dm(:, 1:length(water.d));

for j = 1 :size(wave.Hs,2)
    Hs = wave.Hs(:,j);
    Hs = reshape(Hs,8,[]);
    Tp = wave.Tp(:,j);
    Tp = reshape(Tp,8,[]);
    Dm = wave.Dm(:,j);
    Dm = reshape(Dm,8,[]);
for i = 1 :length(Hs)
    wave.Hsd(i,j) = max(Hs(:,i));
    wave.Tpd(i,j) = max(Tp(:,i));
    wave.Dmd(i,j) = max(Dm(:,i));
end
end

% end of above section

% GRS 1980 to Cartesian
location.Raduis = 6378137.0                                              ;  % raduis of earth in [m]
location.x = location.Raduis * cosd(location.lat) .* cosd(location.long) ;  % x coordinate in [m]
location.y = location.Raduis * cosd(location.lat) .* sind(location.long) ;  % y coordinate in [m]
location.z = location.Raduis * sind(location.lat)                        ;  % z coordinate in [m]

coast.x=location.x;
coast.y=location.y;
xymatcoast = [coast.x, coast.y];

%% write XY coast-line file
fileID = sprintf(['XY_Coast_Line_%d_%d_601_1200.txt'],year(wave.t(1)),year(wave.t(end)));
fileID = fopen(fileID,'w');
for j = 1 : length(xymatcoast)
    fprintf(fileID,' %0.4f \t',xymatcoast(j,:));
    fprintf(fileID,'\n');
end
fclose(fileID);
disp('          ');
disp('          msg:...');
disp('          The coastline coordinates text file has been written...');
disp('          ');


%% waves / water
t = wave.t;
% for daily wave time series
tday(1)=t(1);
for i= 1 :length(t)
    if rem(i,8)==1
        tday(i)=t(i);
    end
end
tday(end)=t(end);
tday=nonzeros(tday);
%
tstart = year(t(1));
tend = year(t(end));

Hs = fillmissing(wave.Hsd,'nearest');
Tp = fillmissing(wave.Tpd,'nearest');
Dm = fillmissing(wave.Dmd,'nearest');

watermat = [water.stn water.d water.lat water.long];

%% write water depth file
waterdepth = [water.stn water.d];
fileID = sprintf(['California%sxyCoordinates' num2str(19762005) '_StN_STATIONS_depth.txt'],'\');
fileID = sprintf(['1_600_Depth' num2str(19762005) '_StN_STATIONS_depth.txt'],'\');                 
fileID = fopen(fileID,'w');
for i = 1 : size(waterdepth,1)
    fprintf(fileID,'%d      %0.4f',waterdepth(i,:));
    fprintf(fileID,'\n');
end
fclose(fileID);
disp('          ');
disp('          msg:...');
fprintf('          The water depth and stations .txt file has been written...\n');
disp('          ');

%% write waves file
t = tday;
for i = 1 : size(Hs,2)
    wavemat = [t Hs(:,i) Tp(:,i) Dm(:,i)];
    fileID = ['wavets_' sprintf('%d', i) 'TS.wvt'];                 
    fileID = fopen(fileID,'w');
    for j = 1 : size(wavemat,1)
        fprintf(fileID,'%0.4f,%0.4f,%0.4f,%0.4f',wavemat(j,:));
        fprintf(fileID,'\n');
    end
    fclose(fileID);
    disp('          ');
    disp('          msg:...');
    fprintf('          The wave climate TS.wvt file number # %d has been written...\n',i);
    disp('          ');
    
% ['wvt\wavets_' sprintf('%d', i) 'TS.wvt     ' sprintf('%d     ', location.x(i)) sprintf('%d\n', location.y(i))]

end

x=location.x;
y=location.y;
fileID = 'wavets.wvt';                                           
fileID = fopen(fileID,'w');
for i = 1 : length(x)
    fprintf(fileID,['wavets_' sprintf('%d', i) 'TS.wvt     ' sprintf('%d     ', location.x(i)) sprintf('%d\n', location.y(i))]);
end
fclose(fileID);
disp('          ');
disp('          msg:...');
fprintf('          The wave climate TS.wvt file number # %d has been written...\n',i);
disp('          ');
% 
% 
disp('          ');
disp('          msg:...');
disp('          The wave climate TS.wvt files has been written...');
disp('          ');

end


