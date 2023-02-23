function [TRANSP]=prepare_transport(S)
% function [TRANSP]=prepare_transport(S)
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

    fprintf('  Prepare transport \n');
    TRANSP=struct;
    TRANSP.trform=S.trform;                                                            % switch for transport formulation (e.g. S.trform='CERC', 'KAMP', 'MILH' or 'VR14')
    TRANSP.b=S.b;                                                                      % CERC : coeff in simple cerc formula
    TRANSP.qscal=S.qscal;                                                              % Calibration factor of the transport (works for all transport formulas)
    TRANSP.d50=S.d50;                                                                  % KAMP & MILH & VR14 : median grain diameter [m]
    TRANSP.porosity=S.porosity;                                                        % KAMP & MILH & VR14 : TRANSP.porosity (typically 0.4) [-]
    TRANSP.tanbeta=S.tanbeta;                                                          % KAMP & MILH & VR14 : mean bed slope [ratio 1/slope]
    TRANSP.rhos=S.rhos;                                                                % KAMP & MILH & VR14 : density of sand [kg/m3]
    TRANSP.rhow=S.rhow;                                                                % KAMP & MILH & VR14 : density of water [kg/m3]
    TRANSP.g=S.g;                                                                      % KAMP & MILH & VR14 : gravitational acceleration [m2/s]
    TRANSP.alpha=S.alpha;                                                              % KAMP & MILH & VR14 : calibration factor for point of breaking (TRANSP.alpha = 1.8 for Egmond data)
    TRANSP.gamma=S.gamma;                                                              % KAMP & MILH & VR14 : breaking coefficient (Hs/h) with 5% breaking waves
    TRANSP.Pswell=S.Pswell;                                                            % VR14 : Percentage swell (between 0 - 100) [-]
    TRANSP.crit=S.crit;                                                                % stability criterion (not active)
    TRANSP.Aw=S.Aw;                                                                    % factor rep. Hs/actual Hs for determining depth of closure at bypassing groyne (1.27 if time series is used, higher for a representative Hs)
    TRANSP.twopoints=S.twopoints;
    
    % boundary conditions
    TRANSP.boundary_condition_start=S.boundary_condition_start;                        % boundary condition 'CTAN','FIXD','FIXD2','FUNC'
    TRANSP.boundary_condition_end=S.boundary_condition_end;                            % 
    TRANSP.BCfile=S.BCfile;                                                            % Boundary condition function in time (file) , columns (time , QS_start, Qs_end)
    TRANSP.QS_start=S.QS_start;                                                        % [m3/year]
    TRANSP.QS_end=S.QS_end;
    
end
