function [COAST]=get_foreshore_orientation(COAST)
% function [COAST]=get_foreshore_orientation(COAST)
% 
% INPUT:
%   COAST          Structure with data of the ShorelineS model, of which is used
%    .phif0     user input for the shoreline orientation [°]
%    .x          x-coordinate of coastal segment [m]
%    .y          y-coordinate of coastal segment [m]
%    .n          number of grid cells of coastal segment
%    .PHIc       Coastline orientation at each grid cell [°]
%
% OUTPUT:
%   COAST
%      .PHIf       Lower shoreface orientation at each grid cell for the current coastal segment [°]
%      .PHIf_x     x-coordinate at each grid cell for all of the coastal segments [m]
%      .PHIf_y     y-coordinate at each grid cell for all of the coastal segments [m]
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
            
    % shoreface orientation (PHIf)
    PHIf=[];
    
    if ischar(COAST.PHIf0)   % input from file
        % read table with predefined shoreface orientation (with 3 columns, PHIx, PHIy and PHIf orientation)
        % INPUT EXAMPLE :   S.phi='foreshoreorientation.txt', 3 column text
        % file
        COAST.PHIf0=load(COAST.PHIf0);
    end
    
    if isempty(COAST.PHIf0) && ~isfield(COAST,'PHIf') % input at t=0
        % lower shoreface at original orientation of coastline at t=0
        % NO INPUT      : use coastline orientation (e.g. in case of CERC)
        % this is only done at t=0. Afterwards the assessed PHIf at t=0 is reinterpolated on the grid every time step (i.e. when field 'PHIf' is available). 
        PHIf_x=COAST.xq;
        PHIf_y=COAST.yq;
        PHIf=COAST.PHIc;
        COAST.PHIf0=PHIf;
       
    elseif isscalar(COAST.PHIf0)
        % lower shoreface at fixed orientation for the whole grid
        % INPUT EXAMPLE :   COAST.PHIf0=312;
        PHIf_x=COAST.xq;
        PHIf_y=COAST.yq;
        PHIf=repmat(COAST.PHIf0,[1,COAST.nq]);
    
    elseif iscell(COAST.PHIf0) && length(COAST.PHIf0)==2 % input at t=0
        % lower shoreface angles based on smoothed coastline (only present section in COAST.PHIc)
        % INPUT EXAMPLE :   COAST.PHIf0={'gaussian',7};
        smoothmethod=COAST.PHIf0{1};
        smoothrange=COAST.PHIf0{2};
        nl=length(COAST.PHIc);
        cPHIf=smoothdata(cosd(repmat(COAST.PHIc,[1 3])),smoothmethod,smoothrange);
        sPHIf=smoothdata(sind(repmat(COAST.PHIc,[1 3])),smoothmethod,smoothrange);
        PHIf=atan2d(sPHIf(nl+1:2*nl),cPHIf(nl+1:2*nl));
        PHIf=mod(PHIf,360);
        PHIf_x=COAST.xq;
        PHIf_y=COAST.yq;
        %COAST.PHIf0=PHIf;

    elseif isvector(COAST.PHIf0) && ~iscell(COAST.PHIf0)  % vectors for interpolation at t>=1 
        % lower shoreface interpolation from timesteps after initial t0 onwards
        PHIf_x=COAST.PHIf_x;
        PHIf_y=COAST.PHIf_y;
        PHIf_mc=COAST.PHIf0;
        
        % find the right alongshore location for each of the wave climates
        % re-interpolate PHIf          
        var1=struct;
        var2=struct;
        var2.PHIf_mc=PHIf_mc;
        [~,var2i,~]=get_interpolation_on_grid('weighted_distance',COAST.xq,COAST.yq,PHIf_x,PHIf_y,var1,var2);        
        PHIf=var2i.PHIf_mc;  
        
    elseif isnumeric(COAST.PHIf0)  % array from input structure
        % lower shoreface at predefined orientation in table
        % INPUT EXAMPLE :   PHIf0=[x1,y1,phif1; ... ; xn,yn,phifn]
        PHIf0=COAST.PHIf0;       
        PHIf_x=PHIf0(:,1)';
        PHIf_y=PHIf0(:,2)';
        PHIf_mc=PHIf0(:,3)';
        
        % find the right alongshore location for each of the wave climates
        % re-interpolate PHIf          
        var1=struct;
        var2=struct;
        var2.PHIf_mc=PHIf_mc;
        [~,var2i,~]=get_interpolation_on_grid('weighted_distance',COAST.xq,COAST.yq,PHIf_x,PHIf_y,var1,var2);        
        PHIf=var2i.PHIf_mc;  
        PHIf_x=COAST.xq;
        PHIf_y=COAST.yq;
    else     
        error('get_foreshore_orientation::Invalid foreshore orientation method');
    end
    
    % Output to structure
    COAST.PHIf=PHIf;
    COAST.PHIf_x=PHIf_x;
    COAST.PHIf_y=PHIf_y;
end
