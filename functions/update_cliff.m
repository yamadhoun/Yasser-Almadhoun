function [CLIFF]=update_cliff(CLIFF,TIME)
% function [CLIFF]=update_cliff(CLIFF,GROYNE,TIME
% 
% updates the coastline at each time step based on transport gradients,
% nourishment volumes and sea level rise
% modifications to describe coastline behaviour next to groynes
% 
%% Copyright notice
%   --------------------------------------------------------------------
%   Copyright (C) 2023 IHE Delft & Deltares
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

    if isempty(CLIFF.x)
        return; 
    end

    %% Spatial steps
    n=CLIFF.n;
    nq=CLIFF.nq;
    nend=n;

    %% Spatial discretisation
    x0=CLIFF.x_mc;y0=CLIFF.y_mc;
    CLIFF.dx=zeros(1,n);
    CLIFF.dy=zeros(1,n);

    %% COMPUTE GRID SIZE
    for i=1:n
        if CLIFF.cyclic
            im1=mod2(i-1,n);
            ip1=mod2(i+1,n);
            im1q=mod2(i-1,nq);
            ip1q=mod2(i,nq);
        else
            im1=max(i-1,1);
            ip1=min(i+1,n);
            im1q=max(i-1,1);
            ip1q=min(i,nq);
        end

        CLIFF.ds(i)=hypot(CLIFF.x_mc(ip1)-CLIFF.x_mc(im1),CLIFF.y_mc(ip1)-CLIFF.y_mc(im1));

    end

    %% COMPUTE RATE OF COASTLINE CHANGE
    for i=1:n
        if CLIFF.cyclic
            im1=mod2(i-1,n);
            ip1=mod2(i+1,n);
        else
            im1=max(i-1,1);
            ip1=min(i+1,n);
        end
        dn=CLIFF.dndt(i)*TIME.dt;
        CLIFF.dx(i)= dn*(CLIFF.y_mc(ip1)-CLIFF.y_mc(im1))/CLIFF.ds(i);
        CLIFF.dy(i)= dn*(CLIFF.x_mc(ip1)-CLIFF.x_mc(im1))/CLIFF.ds(i);  
        CLIFF.retreat=CLIFF.retreat+dn;                                    % to use as dx in Trenhaile and Energy_Flux models

    end

    %% Update cliffline positions, all points
    for i=1:n
        CLIFF.x_mc(i)=x0(i)+CLIFF.dx(i);
        CLIFF.y_mc(i)=y0(i)+CLIFF.dy(i);
    end
    if CLIFF.cyclic
        CLIFF.x_mc(n)=CLIFF.x_mc(1);
        CLIFF.y_mc(n)=CLIFF.y_mc(1);
    end

end