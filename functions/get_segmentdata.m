function [COAST,WAVE,TRANSP,i1,i2]=get_segmentdata(COAST,WAVE,TRANSP,i_mc)
% [COAST,WAVE,TRANSP]=get_segmentdata(COAST,WAVE,TRANSP,i_mc);
%
% INPUT:
%    COAST
%      .x_mc               x-coordinates for all coastal segments
%      .y_mc               y-coordinates for all coastal segments
%    WAVE
%      .PHItdp_mc          incoming nearhsore wave angle for all coastal segments
%      .PHIbr_mc           incoming wave angle at point of breaking for all coastal segments
%    TRANSP
%      .QS_mc              transports for all coastal segments
%      .trform             transport formulation
%    i_mc                  index of segment to be retrieved
%         
% OUTPUT:
%    COAST
%      .x                  x-coordinates of cells for considered coastal segment
%      .y                  y-coordinates of cells for considered coastal segment
%      .n                  number of cells for considered coastal segment
%      .s                  length along the coast for considered coastal segment
%      .cyclic             index whether considered coastal segment is cyclical
%      .n_mc               number of coastal segments
%    WAVE
%      .PHI                incoming wave angle (either nearshore or at point of breaking) for considered coastal segment
%    TRANSP
%      .QS                 transports for considered coastal segment
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

    COAST.i_mc=i_mc; 
    [ COAST.x,COAST.y,COAST.n_mc,i1,i2 ] = get_one_polygon( COAST.x_mc,COAST.y_mc,i_mc );
    [ COAST.s ] = get_one_polygon( COAST.s_mc,COAST.s_mc,i_mc );
    if (strcmpi(TRANSP.trform,'CERC') || strcmpi(TRANSP.trform,'CERC2'))
        [ TRANSP.QS,WAVE.PHI ] = get_one_polygon( TRANSP.QS_mc,WAVE.PHItdp_mc,i_mc );
    else
        [ TRANSP.QS,WAVE.PHI ] = get_one_polygon( TRANSP.QS_mc,WAVE.PHIbr_mc,i_mc );
    end
    if length(COAST.x)<2
        disp('This is too short a coastline segment, good sir.')
        COAST.cyclic=true; % or crash in neumann bnd
    else
        COAST.cyclic = hypot(COAST.x(end)-COAST.x(1),COAST.y(end)-COAST.y(1))<COAST.ds0;
    end
    COAST.n=length(COAST.x);
    COAST.nq=length(TRANSP.QS);
    
end
