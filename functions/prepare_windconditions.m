function [WndCF,WndCL]=prepare_windconditions(S)
% function [WndCF,WndCL]=prepare_windconditions(S)
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
        WndCF=struct;
        WndCL=struct;
        if ~isempty(S.WndCfile)
            if S.RWSfiletype
                WndCraw=load(S.WndCfile);warning off;
                WndCF.timenum=datenum([num2str(WndCraw(:,1))],'yyyymmddHHMM'); 
                WndCF.uz=interpNANs(WndCraw(:,2));
                WndCF.Dir=interpNANsDIR(WndCraw(:,3));
            else
                WndCraw=load(S.WndCfile);warning off;
                WndCF.timenum=datenum([num2str(WndCraw(:,1),'%08.0f'),num2str(WndCraw(:,2),'%06.0f')],'yyyymmddHHMMSS');
                WndCF.uz=interpNANs(WndCraw(:,3));
                WndCF.Dir=interpNANsDIR(WndCraw(:,4));
            end
        elseif ~isempty(S.Windclimfile)
            WCraw=load(S.Windclimfile);
            WndCL=struct;
            WndCL.uz=WCraw(:,1);
            WndCL.Dir=WCraw(:,2);
        end
    else
        WndCF=[];
        WndCL=[];
    end
    
end    
    