function [VAR]=get_smoothdata(VAR,type,smoothsteps,movwindow)
% function [VAR]=get_smoothdata(VAR,type,smoothsteps,movwindow)
% 
% INPUT:
%     VAR           column with the variable [1xN] 
%     type          e.g. 'angle'/'vector' or 'scalar' 
%     smoothsteps   number of times the smaoothing is applied 
%     movwindow     (optional) using a moving average instead of 'central averaging' 
% 
% OUTPUT:
%     VAR           smoothed variable [1xN] 
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

    VAR2=VAR;
    for nn=1:smoothsteps 
        if ~isempty(findstr(type,'angle')) ||  ~isempty(findstr(type,'vector'))
            sVAR=sind(VAR);
            cVAR=cosd(VAR);
            sVARavg=[sVAR(1),(sVAR(1:end-1)+sVAR(2:end))/2,sVAR(end)];
            cVARavg=[cVAR(1),(cVAR(1:end-1)+cVAR(2:end))/2,cVAR(end)];
            sVAR2=(sVARavg(1:end-1)+sVARavg(2:end))/2;
            cVAR2=(cVARavg(1:end-1)+cVARavg(2:end))/2;
            VAR=atan2d(sVAR2,cVAR2);
            VAR=mod(VAR,360);   
        elseif ~isempty(findstr(type,'angle')) || ~isempty(findstr(type,'vector')) && nargin==4
            cVAR=smoothdata(cosd(VAR),'movmean',movwindow);
            sVAR=smoothdata(sind(VAR),'movmean',movwindow);
            VAR2=atan2d(sVAR,cVAR);
            VAR2=mod(VAR2,360);
        elseif nargin==4
            VAR2=smoothdata(VAR,'movmean',movwindow);
        else
            VARavg=[VAR(1),(VAR(1:end-1)+VAR(2:end))/2,VAR(end)];
            VAR2=(VARavg(1:end-1)+VARavg(2:end))/2;  
        end
    end  

    factor=0.2;
    VAR=(1-factor)*VAR+factor*VAR2;
    
end
