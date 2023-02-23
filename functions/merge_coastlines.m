function [ COAST,merged ] = merge_coastlines( COAST, i_mc)
% function [ xnew,ynew,merged ] = merge_coastlines( COAST.x,COAST.y )
% 
% UNTITLED Summary of this function goes here
% Detailed explanation goes here
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

global intersectfn

    eps=.1;
    s(1)=0;
    yesplot=0;
    [ COAST.x,COAST.y,COAST.n_mc ] = get_one_polygon( COAST.x_mc,COAST.y_mc,i_mc );
    for i=2:length(COAST.x)
        if isnan(COAST.x(i))
            s(i)=nan;
        elseif isnan(COAST.x(i-1))
            s(i)=0;
        else
            s(i)=s(i-1)+hypot(COAST.x(i)-COAST.x(i-1),COAST.y(i)-COAST.y(i-1));
        end
    end

    [xx,yy]=intersectfn(COAST.x,COAST.y);
    if yesplot
        figure(3);
        plot(COAST.x,COAST.y,'.-b',xx,yy,'ok');
        hold on
        num=[1:length(COAST.x)];
        for i=1:length(COAST.x)
            text(COAST.x(i),COAST.y(i),num2str(num(i)));
        end
    end
    iX=0;
    ind=[];
    for i=1:length(s)-1
        if ~isnan(s(i))&&~isnan(s(i+1))
            for ip=1:length(xx)
                err=abs(s(i+1)-s(i)-hypot(xx(ip)-COAST.x(i)  ,yy(ip)-COAST.y(i)) ...
                    -hypot(xx(ip)-COAST.x(i+1),yy(ip)-COAST.y(i+1)));
                if err<eps
                    iX=iX+1;
                    ind(iX)=i;
                end
            end
        end
    end
    ind;
    if isempty(ind) || length(ind)<2
        xnew=COAST.x;
        ynew=COAST.y;
        merged=0;
    elseif length(ind)<4
        xnew=[COAST.x(1:ind(1)),COAST.x(ind(2)+1:end),nan,COAST.x(ind(1)+1:ind(2)),COAST.x(ind(1)+1)];
        ynew=[COAST.y(1:ind(1)),COAST.y(ind(2)+1:end),nan,COAST.y(ind(1)+1:ind(2)),COAST.y(ind(1)+1)];
        merged=0;
    else
        xnew=[COAST.x(1:ind(1)),COAST.x(ind(4)+1:end),nan,COAST.x(ind(2)+1:ind(3)),COAST.x(ind(2)+1)];
        ynew=[COAST.y(1:ind(1)),COAST.y(ind(4)+1:end),nan,COAST.y(ind(2)+1:ind(3)),COAST.y(ind(2)+1)];
        merged=1;
    end
    if yesplot
        plot(xnew,ynew,'k','linewidth',2)
        hold off
        if length(ind)>1
            ind
            %%pause
        end
    end

    % insert new section in x_mc and y_mc
    [COAST.x_mc,COAST.y_mc]=insert_section(xnew,ynew,COAST.x_mc,COAST.y_mc,i_mc);
    
end
