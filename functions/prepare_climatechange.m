function [CC]=prepare_climatechange(S, TIME)
% function [CC]=prepare_climatchange(S, TIME)
%
% Read climate change related corrections and initialize structure 
%
% INPUT:
%   S
%      .ccSLR    change in sea level, constant rate [m/yr] or time series of change since beginning simulation [m]
%      .ccHS     change in wave height constant rate [m/yr] or time series of change since beginning simulation [m]
%      .ccDIR    change in wave direction constant rate [deg/yr] or time series of change since beginning simulation [deg]
%
%
% OUTPUT:
%   CC
%       .timenum     time axis, datenum format
%       .SLR         instantaneous correction rate SLR (m)        
%       .HS          instantaneous correction rate Hs (m)
%       .DIR         instantaneous correction rate PHIw  (radians cartesian)
%
%% Copyright notice
%   --------------------------------------------------------------------
%   Copyright (C) 2022 IHE Delft & Deltares
%
%       Dano Roelvink
%       d.roelvink@un-ihe.org
%       Westvest 7
%       2611AX Delft
%
%       Bas Huisman
%       bas.huisman@deltares.nl
%       Boussinesqweg 1
%       2629HV Delft
%
%   This library is free software: you can redistribute it and/or
%   modify it under the terms of the GNU Lesser General Public
%   License as published by the Free Software Foundation, either
%   version 2.1 of the License, or (at your option) any later version.
%
%   This library is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
%   Lesser General Public License for more details.
%
%   You should have received a copy of the GNU Lesser General Public
%   License along with this library. If not, see <http://www.gnu.org/licenses
%   --------------------------------------------------------------------

    fprintf('  Prepare climate change parameters \n');
    
    %% Initialize
    CC=struct;
    CC.timenum=NaN;    % serves as check whether we need introduce_climatechange or not
    CC.SLR=0.0;
    CC.HS=1.0;
    CC.DIR=0.0;
    
    CC.HScor=1.0;
    CC.PHIcor=0.0;
    CC.SLRo=0.0;
    
    %% 1. Sea level rise
    if isscalar(S.ccSLR)             % constant rate ccSLR=0.002 [m/yr]
       CC.timenum=[];
       CC.SLR=S.ccSLR;
        
    elseif ischar(S.ccSLR) && ~isempty(S.ccSLR)          % we have a timeseries with dMSL wrt initial MSL in time from file: yyyymmdd dMSL [m]
        CCraw = load(S.ccSLR);
        try
            CC.timenum=datenum(num2str(CCraw(:,1)),'yyyymmdd');
        catch
            try
                CC.timenum=datenum(num2str(CCraw(:,1)),'yyyymmddHHMM');
            catch
                CC.timenum=datenum(num2str(CCraw(:,1)),'yyyymmddHHMMSS');
            end
        end
        
        % Check length of input and abort if necessary
        if CC.timenum(1)>TIME.timenum0 || CC.timenum(end)<TIME.tend
           error('Time series for SLR not long enough'); 
        end
        
        CC.SLR=interpNANs(CCraw(:,2));      % correct derivative determined in introduce_climatechange 
        
    end
    
    
    %% 2. Wave characteristics correction
    if isscalar(S.ccHS)
        CC.timenum=[];
        CC.HS=S.ccHS;
        
    elseif ischar(S.ccHS) && ~isempty(S.ccHS)
        CCraw = load(S.ccHS);
        try
            CC.timenum=datenum(num2str(CCraw(:,1)),'yyyymmdd');
        catch
            try
                CC.timenum=datenum(num2str(CCraw(:,1)),'yyyymmddHHMM');
            catch
                CC.timenum=datenum(num2str(CCraw(:,1)),'yyyymmddHHMMSS');
            end
        end
        
        % Check length of input and abort if necessary
        if CC.timenum(1)>TIME.timenum0 || CC.timenum(end)<TIME.tend
           error('Time series for CC wave height change not long enough'); 
        end
        
        CCraw(:,2)=interpNANs(CCraw(:,2));
        CC.HS=[0 diff(CCraw(:,2))./diff(CC.timenum)*365.];    % rate/yr
        
    end
    
    if isscalar(S.ccDIR)
        CC.timenum=[];
        CC.DIR=S.ccDIR;
        
    elseif ischar(S.ccDIR) && ~isempty(S.ccDIR)
        CCraw = load(S.ccDIR);
        try
            CC.timenum=datenum(num2str(CCraw(:,1)),'yyyymmdd');
        catch
            try
                CC.timenum=datenum(num2str(CCraw(:,1)),'yyyymmddHHMM');
            catch
                CC.timenum=datenum(num2str(CCraw(:,1)),'yyyymmddHHMMSS');
            end
        end
        
        % Check length of input and abort if necessary
        if CC.timenum(1)>TIME.timenum0 || CC.timenum(end)<TIME.tend
           error('Time series for CC wave direction change not long enough'); 
        end
        
        CCraw(:,2)=interpNANsDIR(CCraw(:,2));
        CC.DIR=[0 diff(CCraw(:,2))./diff(CC.timenum)*365.];    % rate/yr
        
    end    
end