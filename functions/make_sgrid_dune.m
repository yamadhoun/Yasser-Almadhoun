function [x_dune,y_dune,x_dune0,y_dune0,Dfelev0]=make_sgrid_dune(S,x_dune,y_dune,x_dune0,y_dune0,x_mc0,y_mc0,i_mc,x_hard,y_hard,ds0,smoothfac,Dfelev0)
% function [x_dune,y_dune,x_dune0,y_dune0,Dfelev0]=make_sgrid_dune(S,x_dune,y_dune,x_dune0,y_dune0,x_mc0,y_mc0,i_mc,x_hard,y_hard,ds0,smoothfac,Dfelev0)
% 
% smoothing dune foots as same as coasline
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
        %% find number of sections
        nansd=find(isnan(x_dune));
        n_dune=length(nansd)+1;
        for ij=1:n_dune
            [ xc0,yc0,~,~,~ ] = get_one_polygon( x_mc0,y_mc0,i_mc );
            [ xd,yd,~,i1,i2 ] = get_one_polygon( x_dune,y_dune,ij );
            [ xd0,yd0,~,i01,i02 ] = get_one_polygon( x_dune0,y_dune0,ij );
            s=zeros(size(xd));
            s0=zeros(size(xd));
            s(1)=0;
            s0(1)=0;
            for i=2:length(xd)
                s(i)=s(i-1)+hypot(xd(i)-xd(i-1),yd(i)-yd(i-1));
            end
            for ii=2:length(xd0)
                s0(ii)=s0(ii-1)+hypot(xd0(ii)-xd0(ii-1),yd0(ii)-yd0(ii-1));
            end
            snew=s;
            snew0=s0;
            i=2;
            if length(x_hard)>1
                n_st=length(find(isnan(x_hard)))+1;      % number of structures
                 x_hnew=zeros(1,6*n_st);
                 y_hnew=zeros(1,6*n_st);
                 for i_st=1:n_st
                     [ xs,ys,~,~,~ ] = get_one_polygon( x_hard,y_hard,i_st );
                     for kk=2:length(xs)
                         phist=90-atan2d(ys(kk)-ys(kk-1),xs(kk)-xs(kk-1));
                         if phist<0
                             phist=phist+180;
                         elseif phist>180
                             phist=phist-180;
                         end
                         theta_st(kk-1)=phist;
                     end
                     len=5*hypot(max(xs)-min(xs),max(ys)-min(ys));                % Length of arbitrary line
                     o=atan2d((yc0(end)-yc0(1)),(xc0(end)-xc0(1)));
                     if (o >= -45 &&  o<=45) ||  (o < -45 &&  o>-135)                        %% Shoreline horizontal or vertical & land on the down or left side
                         xS1=[xs(1),xs(1)-len*sind(theta_st(1))]; 
                         yS1=[ys(1),ys(1)-len*cosd(theta_st(1))]; 
                         xS2=[xs(end),xs(end)-len*sind(theta_st(end))]; 
                         yS2=[ys(end),ys(end)-len*cosd(theta_st(end))]; 
                     elseif  (abs (o) >= 135) ||  (o > 45 &&  o<135)                  %% Shoreline horizontal or vertical & land on the up or right side 
                         xS1=[xs(1),xs(1)+len*sind(theta_st(1))]; 
                         yS1=[ys(1),ys(1)+len*cosd(theta_st(1))]; 
                         xS2=[xs(end),xs(end)+len*sind(theta_st(end))]; 
                         yS2=[ys(end),ys(end)+len*cosd(theta_st(end))];         
                     end
                     xsnew=[xS1(end),xs,xS2(end)];
                     ysnew=[yS1(end),ys,yS2(end)];
                     [structS]=find_struct_2(xd,yd,xsnew,ysnew);
                     FS= find(structS);
                     FS1=FS(2:2:end);
                     nostruc=length(FS)/2; 
                 end
            else
                FS=[];
                FS1=[];
                nostruc=[];
            end
            while i<=length(snew)
                ds=snew(i)-snew(i-1);
                if ismember(i,FS1)
                    i=i+1;
                    continue
                end
                if ds<ds0/2
                    for is=1:nostruc
                        if i<=FS(2*(is-1)+1)
                            FS(2*(is-1)+1:end) = FS(2*(is-1)+1:end)-1;
                            FS1=FS(2:2:end);
                            if 1                           % shift structures after the grid point only one step back for one iteration
                                break
                            end
                        end
                    end
                    %throw out point i
                    if ismember(i,FS)
                       snew=[snew(1:i-2),snew(i:end)];
                       snew0=[snew0(1:i-2),snew0(i:end)];
                      Dfelev0=[Dfelev0(1:i-2),Dfelev0(i:end)];
                       i=i+1;   
                    else
                       snew=[snew(1:i-1),snew(i+1:end)];
                       snew0=[snew0(1:i-1),snew0(i+1:end)];
                       Dfelev0=[Dfelev0(1:i-1),Dfelev0(i+1:end)];
                       i=i+1;
                    end
                elseif ds>ds0*2
                    for is=1:nostruc
                        if i<=FS(2*(is-1)+1)
                            FS(2*(is-1)+1:end) = FS(2*(is-1)+1:end)+1;
                            FS1=FS(2:2:end);
                            if 1                        % shift structures after the point only one step forward for one iteration
                                break
                            end
                        end
                    end
                    %insert point i
                    snew=[snew(1:i-1),.5*(snew(i-1)+snew(i)),snew(i:end)];
                    snew0=[snew0(1:i-1),.5*(snew0(i-1)+snew0(i)),snew0(i:end)];
                    Dfelev0=[Dfelev0(1:i-1),Dfelev0(i-1),Dfelev0(i:end)];
                    i=i+1;
                else
                    i=i+1;
                end
            end
            snew(2:end-1)=smoothfac*snew(1:end-2)+(1.-2*smoothfac)*snew(2:end-1)+smoothfac*snew(3:end);
            snew0(2:end-1)=smoothfac*snew0(1:end-2)+(1.-2*smoothfac)*snew0(2:end-1)+smoothfac*snew0(3:end);
            xd=interp1(s,xd,snew);
            yd=interp1(s,yd,snew);
            xd0=interp1(s0,xd0,snew0);
            yd0=interp1(s0,yd0,snew0);
              %% insert xd,yd back into x_dune,y_dune
            x_dune=[x_dune(1:i1-1),xd,x_dune(i2+1:end)];
            y_dune=[y_dune(1:i1-1),yd,y_dune(i2+1:end)];
            x_dune0=[x_dune0(1:i01-1),xd0,x_dune0(i02+1:end)];
            y_dune0=[y_dune0(1:i01-1),yd0,y_dune0(i02+1:end)];
        end
    else
        x_dune=[];
        y_dune=[];
        x_dune0=[];
        y_dune0=[];
    end
end
