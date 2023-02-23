function [Wat_cf,Wat_cl]=prepare_tideconditions(S)
% function [Wat_cf,Wat_cl]=prepare_tideconditions(S)
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
        Wat_cf=struct;
        Wat_cl=struct;
        if ~isempty(S.Watfile)
            WatCraw=load(S.Watfile);warning off;
            Wat_cf.timenum=datenum([num2str(WatCraw(:,1))],'yyyymmddHHMM'); 
            Wat_cf.SWL=interpNANs(WatCraw(:,2));
        elseif ~isempty(S.Watclimfile)
            WCraw=load(S.Watclimfile);
            Wat_cl.SWL=WCraw(:,1);
        end
    else
        Wat_cf=[];
        Wat_cl=[];
    end
    
end
