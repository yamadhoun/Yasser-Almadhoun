function plot_debug(use_debug,COAST,WAVE,TRANSP,im3,ip3)
% function plot_debug(use_debug,COAST.x,COAST.y,WAVE.dPHItdp,TRANSP.debug.QS2,TRANSP.debug.QS0,TRANSP.debug.QS1,im3,ip3)
%
% INPUT :
%    use_debug  : switch for making the debug plot
%    COAST.x            : x-coordinates of coastline section [Nx1]
%    COAST.y            : y-coordinates of coastline section [Nx1]
%    WAVE.dPHItdp       : relative angle of incidence of the waves [Nx1]
%    TRANSP.debug.QS2   : transport rates : transport computation, upwind correction and shadows [Nx1]
%    TRANSP.debug.QS0   : transport rates : transport computation [Nx1]
%    TRANSP.debug.QS1   : transport rates : transport computation and upwind correction [Nx1]
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

    if use_debug==1
        %if length(x)>100
        fc9=gcf;
        figure(112);
        clf;plot_figureproperties(gcf,1220,920,30);
        
        x2=(COAST.x(1:end-1)+COAST.x(2:end))/2;
        y2=(COAST.y(1:end-1)+COAST.y(2:end))/2;
        
        %% dPHI
        subplot(2,2,1);
        plot(COAST.x,COAST.y);hold on;
        title('WAVE.dPHItdp [?]');
        scatter(x2,y2,10,WAVE.dPHItdp);
        set(gca,'clim',[0,360]);
        colorbar;
        try;plot(COAST.x(im3(2:end,1)),COAST.y(im3(2:end,1)),'ko');end
        %xlim(1.0e+04 * [ 7.1911    7.3864 ]);
        %ylim(1.0e+05 * [ 4.5184    4.5359 ]);
        
        %% TRANSP.debug.QS2
        subplot(2,2,2);
        plot(COAST.x,COAST.y);hold on;
        title('TRANSP.debug.QS2');
        scatter(x2,y2,10,TRANSP.debug.QS2);
        colorbar;
        try;plot(COAST.x(im3(2:end,1)),COAST.y(im3(2:end,1)),'ks');end;
        %xlim(1.0e+04 * [ 7.1911    7.3864 ]);
        %ylim(1.0e+05 * [ 4.5184    4.5359 ]);
        %TRANSP.debug.QS2(IDshadow1)=0.;
        %TRANSP.debug.QS2(IDshadow2)=0.;
        
        %% TRANSP.debug.QS0
        subplot(2,2,3);
        plot(COAST.x,COAST.y);hold on;
        title('TRANSP.debug.QS0');
        scatter(x2,y2,10,TRANSP.debug.QS0);
        colorbar;
        try;plot(x2(im3),y2(im3),'ks');end
        %xlim(1.0e+04 * [ 7.1911    7.3864 ]);
        %ylim(1.0e+05 * [ 4.5184    4.5359 ]);
        
        %% TRANSP.debug.QS1
        subplot(2,2,4);
        plot(COAST.x,COAST.y);hold on;
        title('TRANSP.debug.QS1');
        scatter(x2,y2,10,TRANSP.debug.QS1);
        colorbar;
        try;plot(x2(ip3),y2(ip3),'k+');end
        %set(gca,'clim',[0,1]);colorbar;
        %xlim(1.0e+04 * [ 7.1911    7.3864 ]);
        %ylim(1.0e+05 * [ 4.5184    4.5359 ]);
        
        figure(fc9);
        %end
    end
end