function [COAST,WAVE,TRANSP]=collect_variables(COAST,WAVE,TRANSP,i_mc)
% function [COAST,WAVE,TRANSP]=collect_variables(COAST,WAVE,TRANSP,i_mc)
%
% INPUT
%     WAVE
%     COAST
%     TRANSP
%     i_mc
%
% OUTPUT
%     WAVE      with updated 'mc'
%     COAST     with updated 'mc'
%     TRANSP    with updated 'mc'
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

    dSds=0;
    fieldnm{1}={'s','PHIc','PHIc2'};
    fieldnm{2}={'PHIo','PHItdp','PHIbr','HSo','HStdp','HSbr','dPHIo','dPHItdp','dPHIbr','hbr'};
    fieldnm{3}={'QS'};

    if i_mc==1
        for kk=1:length(fieldnm{1})
        COAST.([(fieldnm{1}{kk}),'_mc'])=COAST.(fieldnm{1}{kk}); 
        end
        for kk=1:length(fieldnm{2})
        WAVE.([(fieldnm{2}{kk}),'_mc'])=WAVE.(fieldnm{2}{kk}); 
        end
        for kk=1:length(fieldnm{3})
        TRANSP.([(fieldnm{3}{kk}),'_mc'])=TRANSP.(fieldnm{3}{kk}); 
        end
        
    else
        for kk=1:length(fieldnm{1})
        COAST.([(fieldnm{1}{kk}),'_mc'])=[COAST.([(fieldnm{1}{kk}),'_mc']),nan,COAST.(fieldnm{1}{kk})]; 
        end
        for kk=1:length(fieldnm{2})
        WAVE.([(fieldnm{2}{kk}),'_mc'])=[WAVE.([(fieldnm{2}{kk}),'_mc']),nan,WAVE.(fieldnm{2}{kk})]; 
        end
        for kk=1:length(fieldnm{3})
        TRANSP.([(fieldnm{3}{kk}),'_mc'])=[TRANSP.([(fieldnm{3}{kk}),'_mc']),nan,TRANSP.(fieldnm{3}{kk})]; 
        end
        
    end
end
