function [CC]=introduce_climatechange(CC, TIME)
% function [CC]=introduce_climatechange(CC,TIME)
%
% Derive climate change related corrections at considered time instance.
%
% INPUT:
%   CC
%      .timenum time in datenum
%      .SLR     rate of sea level rise [m/yr]
%      .HS      change in wave height since start simulation (m)
%      .DIR     change in wave direction since start of simulation (deg)
%
%   TIME       current moment in time 'tnow' is used (in datnum format [days since 0000-00-01 00:00])
%
% OUTPUT:
%   CC
%              .SLRo            sea level rise at cosnidered time instance (m)
%              .HScor           Significant wave height correction at considered time instance (m)
%              .PHIcor          Wave direction correction at considered time instance (deg)
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

    if isnan(CC.timenum); return; end
    
    %% Constant rate per year
    if isempty(CC.timenum)
        % elapsed time
        dt=(TIME.tnow-TIME.timenum0)/365.;
        CC.SLRo   = CC.SLR;   % goes *dt in calculation dn
        CC.HScor  = CC.HS*dt;
        CC.DIRcor = CC.DIR*dt/180.*pi;
        
    %% Time series    
    else
        dt      = (TIME.tnow-TIME.tprev)/365.;
        
        if length(CC.timenum(:))==length(CC.SLR(:)) && dt>0.0
            slr     = interp1(CC.timenum(:),CC.SLR(:),TIME.tnow);   % SLR amount since simulation start in [m]
    		slrprev = interp1(CC.timenum(:),CC.SLR(:),TIME.tprev);
            if isnan(slr) || isnan(slrprev)
                error('Could not determine SLR rate. Check your input time series');
            end
            CC.SLRo = (slr-slrprev)/dt;   % instantanious rate for coastline_change
        else    
            CC.SLRo=0.0;
        end
        
        
        if length(CC.timenum(:))==length(CC.HS(:))
            hsfactor     = interp1(CC.timenum(:),CC.HS(:),TIME.tnow) / CC.HS(1);
            if isnan(rate)
                error('Could not determine wave height change factor. Check your input time series');
            end
            CC.HScor = hsfactor;
        end
        
        if length(CC.timenum(:))==length(CC.DIR(:))
            sinrate      = interp1(CC.timenum(:),sind(CC.DIR(:)),TIME.tnow);
            cosrate      = interp1(CC.timenum(:),cosd(CC.DIR(:)),TIME.tnow);
            rate         = atan2d(sinrate,cosrate);
            if isnan(rate)
                error('Could not determine wave direction change rate. Check your input time series');
            end
            CC.PHIcor = CC.PHIcor + rate*dt;
        end
    end
    
end  
