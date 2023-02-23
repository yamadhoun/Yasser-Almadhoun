function [TRANSP]=transport_shadow_dune(DUNE,STRUC,WAVE,TRANSP)
% function [TRANSP]=transport_shadow_dune(DUNE,STRUC,WAVE,TRANSP)
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

    if DUNE.used
        x=DUNE.x_dune;
        y=DUNE.y_dune;
        [ ~,~,TRANSP.shadowD ]        = find_shadows_mc( x,y,[DUNE.x_dune],[DUNE.y_dune],WAVE.PHItdp,0 );
        if ~isempty(STRUC.x_hard)&&~isempty(DUNE.x_dune)
            [ ~,~,TRANSP.shadowD_h ]      = find_shadows_mc( x,y,[STRUC.x_hard],[STRUC.y_hard],WAVE.PHItdp,1);
        else
            TRANSP.shadowD_h=[];
        end
    else
        TRANSP.shadowD_h=[];
        TRANSP.shadowD=[];
    end
end
