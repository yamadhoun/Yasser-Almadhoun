function [GROYNE] = transport_bypass(TRANSP,WAVE,COAST,STRUC,GROYNE)
%% compute bypass rate at groynes if transport is towards groyne
%
% INPUT:
%    TRANSP
%       .QS          : longshore transport rate
%    WAVE
%       .HStdp       : Hs at toe of dynamic profile
%    COAST
%       .x           : x coordinates of coast section
%       .y           : y coordinates of coast section
%    STRUC
%       .x_hard      : x coordinates of all hard structures
%       .y_hard      : y coordinates of all hard structures
%    GROYNE
%       .x           : x coordinates of groynes
%       .y           : y coordinates of groynes
%       .QS          : bypass transport rates for each groyne, found earlier
%       .phigroyne   : orientation of groyne in °N
% 
% OUTPUT:
%    GROYNE
%       .QS          : updated array of bypass transport at groynes
%
%% Copyright notice
%   --------------------------------------------------------------------
%   Copyright (C) 2020 IHE Delft & Deltares
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

    if ~isempty(GROYNE.QS)

        %% determine cross-shore extent of groynes w.r.t. coastline
        for side=1:2
            idgroyne=find(GROYNE.idcoast(:,side)==COAST.i_mc); 
            if side==1
                % end point (left side of groyne & right side of coastal segment)       
                xbnd=COAST.x(end); 
                ybnd=COAST.y(end); 
                HSbnd=WAVE.HStdp(end); 
                QSbnd=TRANSP.QS(end); 
            else
                % start point (right side of groyne & left side of coastal segment)
                xbnd=COAST.x(1); 
                ybnd=COAST.y(1); 
                HSbnd=WAVE.HStdp(1); 
                QSbnd=TRANSP.QS(1); 
            end
            if (~isempty(idgroyne))
                ist=GROYNE.strucnum(idgroyne); 
                [x_str,y_str,n_str]=get_one_polygon(STRUC.x_hard,STRUC.y_hard,ist); 
                
                % furthest extent of idgroyne relative to shore normal 
                th=GROYNE.phigroyne(idgroyne);                               % orientation of groyne in °N
                proj=(x_str(2:end-1)-xbnd)*sind(th)+(y_str(2:end-1)-ybnd)*cosd(th); 
                [ystr,tipindx0]=max(proj);                              % ystr is the most seaward extent of the groyne on this side
                GROYNE.tipindx(idgroyne,side)=tipindx0+1;               % index of most seaward point on this side
                eps=1e-1; 
                ystr=max(ystr,eps); 
                
                % compute bypass fraction based on proxy for cross-shore extent of longshore transport and dean profile
                if (QSbnd>0 && side==1) || (QSbnd<0 && side==2)  
                    Ap = (1.04+0.086*log(TRANSP.d50))^2;                % sediment scale parameter
                    Ds=Ap*ystr^(2/3);                                   % depth at tip of structure according to dean profile
                    Dlt=TRANSP.Aw./WAVE.gamma*WAVE.HStdp(end);          % most seaward extent of the longshore current
                    BPF=1-Ds/Dlt;                                       % fraction of bypassing sediment
                    if isnan(BPF)                                       % if Aw is zero, switch off bypassing
                       BPF=0; 
                    end 
                    BPF=min(max(BPF,0),1); 
                    GROYNE.QS(idgroyne,side)=BPF*QSbnd; 
                end
            end
        end
    end
end

