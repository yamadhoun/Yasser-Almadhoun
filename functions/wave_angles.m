function [WAVE,TRANSP]=wave_angles(COAST,WAVE,TRANSP,STRUC)
% function [WAVE,TRANSP]=wave_angles(COAST,WAVE,TRANSP,STRUC)
%
%
% Computation of critical wave angles and TRANSP.QSmax at dynamic boundary and point of breaking.
% S-Phi curve at point of breaking to obtain 'critical coastline orientation' with respect to offshore conditions ('dPHI')
%
% INPUT:
%    COAST
%       .n           : number of coastline points
%    WAVE
%       .PHItdp      : Incoming wave direction (in degrees North)
%       .HStdp       : Incoming wave height
%       .TP          : Wave period [s]
%       .dnearshore  : depth of dynamic boundary
%       .c1          : coefficient describing the curvature of the QS-Phi curve (only used in case trform='RAY')
%       .c2          : coefficient describing the curvature of the QS-Phi curve (only used in case trform='RAY')
%       .QSoffset    : coefficient describing the curvature of the QS-Phi curve (only used in case trform='RAY')
%    TRANSP 
%       .trform            transport formulation (either 'CERC', 'KAMP', 'MILH', 'CERC3', 'VR14')
%       .gamma             Breaking coefficient [-]
%       .alpha             Calibration factor for point of breaking (e.g. 1.8 for Egmond)
%       .<various transport formulation variables>
% 
% OUTPUT:
%    WAVE
%       .dPHIcrit    : Critical orientation of the coastline with respect to breaking waves ([2] or [Nx2] Radians)
%    TRANSP
%       .QSmax       : Maximum transport foir critical wave angle
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

    if ~strcmpi(TRANSP.trform,'RAY')
        PHIlocal = WAVE.PHItdp;
        HSlocal  = WAVE.HStdp;
        DEPlocal = repmat(WAVE.dnearshore,[1,COAST.nq]);
        if length(WAVE.TP)==1
            WAVE.TP=repmat(WAVE.TP,[1,COAST.nq]);
        end
        if length(DEPlocal)==1
            [kh1,WAVE.ctdp]=wave_GUO2002(WAVE.TP,DEPlocal);
            WAVE.ntdp = 0.5.*(1.+(2.*kh1)./sinh(2.*kh1));
            kh1=repmat(kh1,[1,COAST.nq]);
        else
            WAVE.ctdp=[];
            WAVE.ntdp=[];
            for i=1:COAST.nq
                [kh1(i),WAVE.ctdp(i)]=wave_GUO2002(WAVE.TP(i),DEPlocal(i));
                WAVE.ntdp(i) = 0.5.*(1.+(2.*kh1(i))./sinh(2.*kh1(i)));
            end
        end
    else
        c1 = WAVE.c1;
        c2 = WAVE.c2;
        QSoffset = WAVE.QSoffset;
    end
   
    %% Define critical angles for CERC and CERC2
    % critical angle for CERC1 = 45°    (Fixed, no influence of dynamic boundary)
    % critical angle for CERC2 = 42.39° (Fixed, no influence of dynamic boundary)
    % critical angle for other transport formulations depends on refraction and shoaling on the shoreface
    % these are therefore determined iteratively

    WAVE.dPHIcrit=repmat(45,[1,COAST.nq]);
    WAVE.dPHIcritbr=repmat(45,[1,COAST.nq]);
    WAVE0=WAVE;
    TRANSP.QSmax=zeros(1,COAST.nq);
    
    if strcmpi(TRANSP.trform,'CERC')
        WAVE.dPHIcrit=repmat(45,[1,COAST.nq]);
        WAVE0.dPHItdp=WAVE.dPHIcrit;
        TRANSP=transport(TRANSP,WAVE0,STRUC,'QSmax');
        
    elseif strcmpi(TRANSP.trform,'CERC2')
        WAVE.dPHIcrit=repmat(42.39,[1,COAST.nq]);
        WAVE0.dPHItdp=WAVE.dPHIcrit;
        TRANSP=transport(TRANSP,WAVE0,STRUC,'QSmax');
        
    elseif strcmpi(TRANSP.trform,'RAY')
        WAVE.dPHIcrit=abs(1./(sqrt(2).*c2));
        WAVE0.dPHItdp=WAVE.dPHIcrit;
        TRANSP=transport(TRANSP,WAVE0,STRUC,'QSmax');
        
    else         
        % compute dPHIcrit, dPHIcritbr and QSmax
        [WAVE,TRANSP]=get_Sphimax(WAVE,TRANSP,STRUC);
        
    end
    
    % limit refraction
    refraclimit=90.0;
    WAVE.dPHItdp((WAVE.dPHItdp-WAVE.dPHIo)<-refraclimit)=-refraclimit;
    WAVE.dPHItdp((WAVE.dPHItdp-WAVE.dPHIo)>refraclimit)=refraclimit;
    refraclimit=90.0;
    WAVE.dPHIbr((WAVE.dPHIbr-WAVE.dPHIo)<-refraclimit)=-refraclimit;
    WAVE.dPHIbr((WAVE.dPHIbr-WAVE.dPHIo)>refraclimit)=refraclimit;
    
end
