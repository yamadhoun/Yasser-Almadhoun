function [TIME, WAVE] = get_coastlineupdate_timestep(TIME, BATHY, WAVE, FORMAT)
% function [TIME, WAVE] = get_coastlineupdate_timestep(TIME, BATHY, WAVE, FORMAT)
% 
% INPUT:
%   TIME, BATHY, WAVE, FORMAT
%
% OUTPUT:
%   TIME, WAVE
%
%% Copyright notice
%   --------------------------------------------------------------------
%   Copyright (C) 2021 IHE Delft & Deltares
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

    try
        indxw=WAVE.WVC(1).indxw;
    catch
        indxw=1;
    end
    
    if (~isempty(WAVE.WVCfile) || ~isempty(WAVE.Wavematfile)) && ~isempty(WAVE.WVC(1).timenum)
       WVCtimenum=WAVE.WVC(1).timenum;
    else 
       WVCtimenum=[];   
    end
    
    %% time step determination for coastline update
    if TIME.tc==0
       % do nothing; timestep fixed at TIME.dt
    else
        %% Cut timestep to output time point
        if ((~isempty(FORMAT.SLplot) && FORMAT.tplot~=0) || (~isempty(BATHY.bathy_update) && BATHY.tupdate~=0))
            if (~isempty(WAVE.WVCfile) || ~isempty(WAVE.Wavematfile)) && ~isempty(WVCtimenum)
                TIME.dt=min(TIME.tc*TIME.adt,(min((WVCtimenum(indxw+1)-TIME.tnow)/365,min([(TIME.tend-TIME.tnow)/365,(FORMAT.tplot-TIME.tnow)/365,(BATHY.tupdate-TIME.tnow)/365]))));
                if TIME.dt==(WVCtimenum(indxw+1)-TIME.tnow)/365
                    WAVE.WVC(1).indxw=indxw+1;
                end
            else
                TIME.dt=min(TIME.tc*TIME.adt,min([(TIME.tend-TIME.tnow)/365,(FORMAT.tplot-TIME.tnow)/365,(BATHY.tupdate-TIME.tnow)/365]));
            end
        else
            %% Cut timestep to wave condition update time point
            if (~isempty(WAVE.WVCfile) || ~isempty(WAVE.Wavematfile)) && ~isempty(WVCtimenum)
                TIME.dt=min(TIME.tc*TIME.adt,(min((WVCtimenum(indxw+1)-TIME.tnow)/365,(TIME.tend-TIME.tnow)/365)));
                if TIME.dt==(WVCtimenum(indxw+1)-TIME.tnow)/365
                    WAVE.WVC(1).indxw=indxw+1;
                end
            else
                TIME.dt=min(TIME.tc*TIME.adt,(TIME.tend-TIME.tnow)/365);
            end
        end
    end
end
