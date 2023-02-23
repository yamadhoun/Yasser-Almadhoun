function [kd]=wave_diffraction_coeff(omega,kdform)
% function [kd]=wave_diffraction_coeff(omega,kd)
%
% INPUT:
%    omega : angle difference in degrees
%    S.kd  : 
%
% OUTPUT:
%    kd    : 
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

    if kdform=='Roelvink'
        om_x=(omega+90)/180;
        kd=1-exp(-(.5./om_x).^4);
        kd(omega<-90)=1;
        
    elseif kdform=='Kamphuis'
    
        omega=-omega;   % Note: we define omega positive towards shadow zone
        kd=ones(size(omega));
        kd(omega <= 0 & omega >=-90) = abs(0.69+0.008*omega(omega <= 0 & omega >=-90)); % 
        kd(omega >  0 & omega <= 40) = 0.71+0.37*sind(omega(omega >  0 & omega <= 40));
        kd(omega > 40 & omega <= 90) = 0.83+0.17*sind(omega(omega > 40 & omega <= 90));
        kd(omega < -90              ) = 0;
    end
end   
