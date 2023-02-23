function [phiwnd,S]=introduce_wind(S,WndCF,WndCL,tnow)
% function [phiwnd,S]=introduce_wind(S,WndCF,WndCL,tnow)
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

    if S.dune
        if ~isempty(S.WndCfile)
            WndCF.uzi=interp1(WndCF.timenum,WndCF.uz,tnow); 
            WndCF.Diri=mod(atan2d(interp1(WndCF.timenum,sind(WndCF.Dir),tnow),interp1(WndCF.timenum,cosd(WndCF.Dir),tnow)),360);
            S.uz=interpNANs(WndCF.uzi);
            phiwnd=interpNANsDIR(WndCF.Diri);
            
        elseif ~isempty(S.Windclimfile)
            iwc=round((rand*length(WndCL.uz)+0.5));
            phiwnd=WndCL.dir(iwc);
            S.uz=WndCL.uz(iwc);      
        else
            phiwnd=S.phiwnd0;         
        end
    else
        phiwnd=[];
    end
end

%% do we need to include S.z ( elevation of measured data) here? or S.z should be manually inserted and used in dune evolution function?
