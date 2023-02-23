function [COAST]=get_coastline_orientation(COAST)
% function [PHIc]=get_coastline_orientation(x,y,n)
% 
% INPUT:
%   x          x-coordinate of coastal segment [m]
%   y          y-coordinate of coastal segment [m]
%   n          number of grid cells of coastal segment
%
% OUTPUT:
%   PHIc       Coastline orientation at each grid cell [°]
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

    % angle in-between grid cell points
    PHIc0=mod(360.0-atan2d(diff(COAST.y),diff(COAST.x)),360);
    
    % coast angles (PHIc) at transport points (only in-between coastline points, not at 'boundary transport points')
    COAST.PHIc=PHIc0;

    % coast angles (PHIc2) at coastline points
    cPHIc=cosd(PHIc0);
    sPHIc=sind(PHIc0);
    cPHIc=[cPHIc(1),(cPHIc(1:end-1)+cPHIc(2:end))/2,cPHIc(end)];
    sPHIc=[sPHIc(1),(sPHIc(1:end-1)+sPHIc(2:end))/2,sPHIc(end)];
    if COAST.cyclic==1
        cPHIc(1)=(cPHIc(1)+cPHIc(end))/2;
        sPHIc(1)=(sPHIc(1)+sPHIc(end))/2;
        cPHIc=cPHIc(1:end-1);
        sPHIc=sPHIc(1:end-1);
    end
    COAST.PHIc2=atan2d(sPHIc,cPHIc);

    % make sure derived coastline angle is smooth
    % COAST.smoothsteps=1;
    % COAST.PHIc=get_smoothdata(PHIc,'angle',COAST.smoothsteps);
end


