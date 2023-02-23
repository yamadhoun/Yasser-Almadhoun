function [WAVE,TRANSP]=get_Sphimax(WAVE,TRANSP,STRUC)
% function [dPHImax,QSmax,dPHIbrmax]=get_Sphimax(TRANSP,HStdp,htdp,TP,gamma)
%
% INPUT:
%    WAVE
%        .HStdp  : wave height at nearshore location (depth-of-closure) [1xN]
%        .htdp   : water depth at nearshore location (depth-of-closure)
%        .TP     : wave period at nearshore location (depth-of-closure) 
%    TRANSP      : input structure with TRANSPORT data  
%        .gamma  : wave breaking coefficient
% 
% OUTPUT:
%    dPHImax  : critical angle at nearshore location (depth-of-closure) [1xN]      
%    QSmax    : maximum transport at nearshore location (depth-of-closure) [1xN]        
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


    %% find maximum by fitting a parabola
    eps=1d-4;
    WAVE.dPHIcrit=[];     % at depth-of-closure
    WAVE.dPHIcritbr=[];   % at point of breaking
    TRANSP.QSmax=[];
    nq=length(WAVE.HStdp); % nr of grid cells of the coast
    QSmax=zeros(1,nq);
    for i=1:nq
        % perform iteration for each alongshore grid cell -> map WAVE to WAVE1
        fieldnm=get_fields(WAVE);
        
        for kk=1:length(fieldnm)
            if length(WAVE.(fieldnm{kk}))==nq
                WAVE1.(fieldnm{kk})=WAVE.(fieldnm{kk})(i);
            end
        end
        WAVE1.dnearshore=WAVE.dnearshore;
        WAVE1.gamma=WAVE.gamma;
        WAVE1.diffraction=WAVE.diffraction;
        WAVE2=WAVE1;
        WAVEm=WAVE1;
        
        % computation 1
        WAVE1.dPHItdp=35;
        [WAVE1]=wave_breakingheight(WAVE1,TRANSP);
        [TRANSP]=transport(TRANSP,WAVE1,STRUC,'QSmax');
        dPHI1=WAVE1.dPHItdp;
        dPHI1br=WAVE1.dPHIbr;
        QS1=TRANSP.QSmax; 
        
        % computation 2
        WAVE2.dPHItdp=45;
        [WAVE2]=wave_breakingheight(WAVE2,TRANSP);
        [TRANSP]=transport(TRANSP,WAVE2,STRUC,'QSmax');
        dPHI2=WAVE2.dPHItdp;
        dPHI2br=WAVE2.dPHIbr;
        QS2=TRANSP.QSmax;
        
        % compute avarage dPHI
        dPHIm=.5*(dPHI1+dPHI2);       % at depth-of-closure
        dPHIbrm=.5*(dPHI1br+dPHI2br); % at point of breaking
        err=1e10;iter=0;err0=inf;err=1; 
        while ( err>eps || iter==0 ) && err0>err
            err0=err;
            iter=iter+1;
            
            % computation 3, 4 etc
            WAVEm.dPHItdp=dPHIm;
            [WAVEm]=wave_breakingheight(WAVEm,TRANSP);
            [TRANSP]=transport(TRANSP,WAVEm,STRUC,'QSmax');
            QSm=TRANSP.QSmax;
            
            % compute better estimate for dPHIm using computed dPHI's and QS's
            A=[dPHI1^2,dPHI1,1;dPHI2^2,dPHI2,1;dPHIm^2,dPHIm,1];              % at depth-of-closure
            A2=[dPHI1br^2,dPHI1br,1;dPHI2br^2,dPHI2br,1;dPHIbrm^2,dPHIbrm,1]; % at point of breaking
            B=[QS1;QS2;QSm];
            warning off
            a=A\B;
            a2=A2\B;
            warning on
            dPHImold=dPHIm;
            if QS1>QS2
                QS2=QSm;
                dPHI2=dPHIm;
            else
                QS1=QSm;
                dPHI1=dPHIm;
            end
            dPHIm=-a(2)/(2*a(1));     % at depth-of-closure
            dPHIbrm=-a2(2)/(2*a2(1)); % at point of breaking
            %disp([num2str(iter),' ',num2str(dPHIm,'%8.4f')])
            if dPHIm>50
                dPHIm=999;
            end
            err=abs(dPHIm-dPHImold);
        end
        
        dPHIbrm(dPHIm==999)=nan;
        dPHIm(dPHIm==999)=nan;
        
        WAVE.dPHIcrit(i)=dPHIm;     % at depth-of-closure
        WAVE.dPHIbrcrit(i)=dPHIbrm; % at point of breaking
        QSmax(i)=QSm;
    end
    WAVE.dPHIcrit=interpNANs(WAVE.dPHIcrit);
    WAVE.dPHIbrcrit=interpNANs(WAVE.dPHIbrcrit);
    TRANSP.QSmax=QSmax;
end
