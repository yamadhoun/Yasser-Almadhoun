function [var1i,var2i,idGRID]=get_interpolation_on_grid(method,xq,yq,xw,yw,var1,var2)
% function [var1i,var2i,idGRID]=get_interpolation_on_grid(method,xq,yq,xw,yw,var1,var2)
%
% find the right alongshore location for each of the wave climates (xq and yq = transport points)
% sort locations alongshore
% choose only the closest wave climate at each grid cell, to make sure that no multiple climates are enforced on a single grid cell.
% Make sure that wave climates are not too far from the grid cells. 
% Otherwise throw out the ones that are at great distance.
%  - remove wave climate points that are further away (in distance) 
%  - from the grid point than the distance to another nearby climate point.
%  - and which are more than 2x the cross-shore distance from the coast than the adjacent climate point
% interpolate the waves at the right alongshore location 
%
%% Copyright notice
%   --------------------------------------------------------------------
%   Copyright (C) 2020 IHE Delft & Deltares
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

    fldnms1=get_fields(var1);
    fldnms2=get_fields(var2);
    var1i=struct;
    var2i=struct;
    idGRID=[];      % these are either the grid locations for each of the applied wave climates 'WVC(kk)' in case of 'alongshore_mapping' or the applied wave cliamte points in sorted-order for each grid cell 'weighted_distance' (with a NaN meaning the wave climate file could not be used, e.g. too far from coastline)
        
    %------------------------------------------------------------------------%
    %% use quadratic distance weighting for each of the wave locations
    %------------------------------------------------------------------------%       
    if strcmpi(method,'weighted_distance') && length(xw)>=1
        
        pwr=2;        % use quadratic distance weighting
        nrpoints=2;   % number of nearest points used (rest at larger distance is omitted)
        for ii=1:length(xq)
            dist=((xq(ii)-xw).^2+(yq(ii)-yw).^2).^0.5;
            [dists,ids]=sort(dist);
            if length(ids)>nrpoints
                dists=dists(1:nrpoints);
                ids=ids(1:nrpoints);
            end

            if length(ids)>=2
                wght0=1-(dists/sum(dists));
                wght=wght0.^pwr/sum(wght0.^pwr);
                %fprintf('ii=%3.0f : ',ii);
                %fprintf('dists=%4.3f, ',dists);
                %fprintf('ids=%4.3f, ',ids);
                %fprintf('wghts=%4.3f, ',wght);
                %fprintf('\n');
                for kk=1:length(fldnms1)
                    var1i.(fldnms1{kk})(ii)=sum(var1.(fldnms1{kk})(ids).*wght);
                end
                for kk=1:length(fldnms2)
                    cphi=sum(cosd(var2.(fldnms2{kk})(ids)).*wght);
                    sphi=sum(sind(var2.(fldnms2{kk})(ids)).*wght);
                    var2i.(fldnms2{kk})(ii)=mod(atan2d(sphi,cphi),360);
                end
                idGRID(ii,:)=ids(:)'; 
            elseif length(ids)==1
                for kk=1:length(fldnms1)
                    var1i.(fldnms1{kk})(ii)=var1.(fldnms1{kk})(ids);
                end
                for kk=1:length(fldnms2)
                    var2i.(fldnms2{kk})(ii)=mod(var2.(fldnms2{kk})(ids),360);
                end
                idGRID(ii,1)=ids; 
            else
                for kk=1:length(fldnms1)
                    var1i.(fldnms1{kk})=0;
                end
                for kk=1:length(fldnms2)
                    var2i.(fldnms2{kk})=0;
                end
            end
        end
    
    %------------------------------------------------------------------------%
    %% Map the wave climate locations on the coast and interpolate alongshore 
    %------------------------------------------------------------------------%
    elseif strcmpi(method,'alongshore_mapping') && length(xw)>=1
        
        % find the right alongshore location for each of the wave climates
        distq = [0,cumsum((diff(xq).^2+diff(yq).^2).^0.5)];
        dist0=[];
        for kk=1:length(xw)
            dist=((xq-xw(kk)).^2+(yq-yw(kk)).^2).^0.5;
            id=find(dist==min(dist),1);
            distxy=distq(id);
            idGRID(kk,1)=id;            
            dist0(kk,1)=distxy;
        end
        idGRID0=idGRID;
        
        % sort locations alongshore
        [WVCsort,IDS]=sort(idGRID);
        [b,IDS2]=sort(IDS);
        dist0=dist0(IDS);
        idGRID=idGRID(IDS);      
        for kk=1:length(fldnms1)
            var1.(fldnms1{kk})=var1.(fldnms1{kk})(IDS);
        end
        for kk=1:length(fldnms2)
            var2.(fldnms2{kk})=var2.(fldnms2{kk})(IDS);
        end
        xw=xw(IDS);
        yw=yw(IDS);
              
        % choose only the closest wave climate at each grid cell, 
        % to make sure that no multiple climates are enforced on a single grid cell.
        dist3=[];
        idDATA=[1:length(idGRID)];
        ids2=[];
        for xx=1:length(idGRID)
            ID=find(idGRID==idGRID(xx));
            dist2=[];
            for dd=1:length(ID)
                dist1=((xq-xw(ID(dd))).^2+(yq-yw(ID(dd))).^2).^0.5;
                dist2(dd)=min(dist1);
            end
            idDATA=setdiff(idDATA,ID(find(dist2~=min(dist2))));
            ids2=[ids2;ID(find(dist2~=min(dist2)))];
            dist3(xx)=min(dist2);
        end
        dist0=dist0(idDATA);
        for kk=1:length(fldnms1)
            var1.(fldnms1{kk})=var1.(fldnms1{kk})(idDATA);
        end
        for kk=1:length(fldnms2)
            var2.(fldnms2{kk})=var2.(fldnms2{kk})(idDATA);
        end
        idGRID=idGRID(idDATA);
        idGRID0(IDS(ids2))=nan;
        
        % Make sure that wave climates are not too far from the grid cells. 
        % Otherwise throw out the ones that are at great distance.
        iter=10;
        IDS0=IDS;
        ids3=[];
        for kk=1:iter 
            dx=[0,cumsum((diff(xq).^2+diff(yq).^2).^0.5)];
            dx0=((xq(end)-xq(1)).^2+(yq(end)-yq(1)).^2).^0.5;
            dxend=dx(end)-dist0(end);
            dx1=dist0(1)+dxend;
            dx=[dx0+dx1,diff(dist0(:))',dx0+dx1];                          % alongshore distance between grid points of 2 consecutive climates
            dcross=[dist3(1)-dist3(end),diff(dist3),dist3(1)-dist3(end)];  % difference in nearest cross-shore distance to wave climates of 2 consecutive climates
            dist3=[dist3(end),dist3,dist3(1)];
            idxx=[1:length(idGRID)];
            for xx=1:length(idGRID)
                % remove wave climate points that are further away (in distance) 
                % from the grid point than the distance to another nearby climate point.
                % and which are more than 2x the cross-shore distance from the coast than the adjacent climate point
                distfactor=2;
                if (dcross(xx)>dx(xx) && dist3(xx+1)>distfactor*dist3(xx)) || ...     % forward check if cross-shore distance to wave output location is not increasing excessively
                   (-dcross(xx+1)>dx(xx+1) && dist3(xx+1)>distfactor*dist3(xx+2))     % backward check if cross-shore distance to wave output location is not increasing excessively
                    idxx=setdiff(idxx,xx);
                    ids3=[ids3;IDS0(xx)];
                end
            end
            dist0=dist0(idxx);
            dist3=dist3(idxx+1);
            idGRID=idGRID(idxx);
            for kk=1:length(fldnms1)
                var1.(fldnms1{kk})=var1.(fldnms1{kk})(idxx);
            end
            for kk=1:length(fldnms2)
                var2.(fldnms2{kk})=var2.(fldnms2{kk})(idxx);
            end
            IDS0=IDS0(idxx);
        end
        
        % interpolate the waves at the right alongshore location 
        if length(dist0)>=2
            for kk=1:length(fldnms1)
                var1i.(fldnms1{kk})=interp1(dist0,var1.(fldnms1{kk}),distq,'linear');
                var1i.(fldnms1{kk})=interpNANs(var1i.(fldnms1{kk}));
            end
            for kk=1:length(fldnms2)
                cphi=interp1(dist0,cosd(var2.(fldnms2{kk})),distq,'linear');
                sphi=interp1(dist0,sind(var2.(fldnms2{kk})),distq,'linear');
                cphi=interpNANs(cphi);
                sphi=interpNANs(sphi);
                var2i.(fldnms2{kk})=mod(atan2d(sphi,cphi),360);
            end
            idGRID(ids3)=nan; 
        elseif length(dist0)==1
            for kk=1:length(fldnms1)
                var1i.(fldnms1{kk})=var1.(fldnms1{kk});
            end
            for kk=1:length(fldnms2)
                var2i.(fldnms2{kk})=mod(var2.(fldnms2{kk}),360);
            end
            idGRID(ids3)=nan; 
        else
            for kk=1:length(fldnms1)
                var1i.(fldnms1{kk})=0;
            end
            for kk=1:length(fldnms2)
                var2i.(fldnms2{kk})=0;
            end
        end
    end
end