function [NOUR]=get_nourishments(TIME,COAST,NOUR)
% function [nour]=get_nourishments(n,NOUR.nourish,t_nour,x_nour,y_nour,n_nour)
%
% INPUT:
%    TIME
%       .tnow
%    COAST
%       .x
%       .y
%       .n           Number of coastline points (length of x and y)
%    NOUR
%       .nourish     Switch for using nourishments (0, 1 or 2) -> 2 is common method used
%       .tstart   Start of nourishments
%       .tend     End of nourishment 
%       .x_nour      x-coordinates of nourishments
%       .y_nour      y-coordinates of nourishments
%       .n_nour      Number of nourishments
%       .rate_m3_per_yr   nourishment rate for each individual measure in [m3/day]
%
% OUTPUT:
%     NOUR
%       .rate_m3_per_m_yr     : actual nourishment rate in [m3/m/day]
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

    NOUR.rate_m3_per_m_yr=zeros(1,COAST.n);
    if NOUR.nourish==1
        for i_n=1:NOUR.n_nour
            if TIME.tnow>=NOUR.tstart(i_n) && TIME.tnow<=NOUR.tend(i_n)
                [ x_n,y_n,~,~,~ ] = get_one_polygon( NOUR.x_nour,NOUR.y_nour,i_n );
                in=inpolygon(COAST.x,COAST.y,x_n,y_n);
                if sum(in)>0
                   cl=(hypot(diff(COAST.x(in)),diff(COAST.y(in))));
                   idnotnan=find(~isnan(cl));
                   clsum=sum(cl(idnotnan));
                   NOUR.rate_m3_per_m_yr=NOUR.rate_m3_per_m_yr+in*NOUR.rate_m3_per_yr(i_n)/clsum;
                end
            end
        end
    elseif NOUR.nourish==2
        for i_n=1:NOUR.n_nour
            if TIME.tnow>=NOUR.tstart(i_n) && TIME.tnow<=NOUR.tend(i_n)
                % find indices of begin and end point of each nourishment
                % compute number of grid cells in-between begin and end point
                dist1=((COAST.x-NOUR.x_nour(i_n,1)).^2+(COAST.y-NOUR.y_nour(i_n,1)).^2).^0.5;
                dist2=((COAST.x-NOUR.x_nour(i_n,2)).^2+(COAST.y-NOUR.y_nour(i_n,2)).^2).^0.5;
                idx1=find(dist1==min(dist1),1);
                idx2=find(dist2==min(dist2),1);
                idx=[min(idx1,idx2):max(idx1,idx2)];
                dxmin=50;
                dxnour=max((diff(NOUR.x_nour(i_n,1:2)).^2 + diff(NOUR.y_nour(i_n,1:2)).^2 ).^0.5,dxmin); 

                % add nourishment 
                % nour [m3/m/year]
                if ~isempty(idx)
                    NOUR.rate_m3_per_m_yr(idx)=NOUR.rate_m3_per_m_yr(idx)+NOUR.rate_m3_per_yr(i_n)/dxnour;
                end
            end
        end
    end
end  
