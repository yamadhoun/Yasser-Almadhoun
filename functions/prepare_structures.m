function [STRUC]=prepare_structures(S,COAST)
% function [x_hard,y_hard]=prepare_structures(S)
%
% INPUT
%     S
%         .struc
%         .x_hard
%         .y_hard
%         .XYoffset
%         .LDBstructures
%         .diffrac
%         .perm
%         .x_perm
%         .y_perm
%
% OUTPUT
%     STRUC
%         .x_hard
%         .y_hard
%         .xtip
%         .ytip.
%         .hstip
%         .delta0
%         .tip_wet
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

    fprintf('  Prepare structures \n');
    STRUC=struct;
    
    % Wave blocking structure
    STRUC.x_hard=[];
    STRUC.y_hard=[];
    STRUC.s_hard=[];
    
    % Wave transmitting structures
    STRUC.perm=S.perm;                                                                  % switch for using hard structures
    STRUC.LDBperm=S.LDBperm;                                                              % LDB with perm structures ([Nx2] ASCII FILE WITHOUT HEADER) <- leave empty to use interactive mode!
    STRUC.wavetransm=S.wavetransm;
    STRUC.x_perm=[];     % for permeable structures
    STRUC.y_perm=[];     % for permeable structures
    
% Diffraction at a structure
    STRUC.diffraction=S.diffraction;
    STRUC.rotfac=S.rotfac;
    STRUC.kd=S.kd;
    STRUC.xtip=[];       % for diffraction
    STRUC.ytip=[];       % for diffraction
    STRUC.hstip=[];      % for diffraction
    STRUC.delta0=[];     % for diffraction
    STRUC.tip_wet=[];    % for diffraction
    STRUC.xp=[];         % for diffraction
    STRUC.yp=[];         % for diffraction
    STRUC.wetstr_mc=[];  % for diffraction
    
    %% Coastal structures (blocking waves)
    if S.struct
        xy_hard=[];
        if ~isempty(S.LDBstructures)
            xy_hard=load(S.LDBstructures);
        end
        if isfield(S,'x_hard')
            if ~isempty(S.x_hard)
                STRUC.x_hard=S.x_hard;
                STRUC.y_hard=S.y_hard;
            end
        elseif ~isempty(xy_hard)
            xy_hard=load(S.LDBstructures);
            x_hard=xy_hard(:,1)'-S.XYoffset(1);
            y_hard=xy_hard(:,2)'-S.XYoffset(2);
            STRUC.x_hard=x_hard;
            STRUC.y_hard=y_hard;
        elseif isempty(xy_hard) && ~isempty(S.LDBstructures)
            STRUC.x_hard=[];
            STRUC.y_hard=[];
        else
            figure(11);
            axis equal;
            xl=xlim;yl=ylim;
            htxt2=text(xl(1)+0.02*diff(xl),yl(2)-0.01*diff(yl),'Add structure (LMB); Next structure (RMB); Exit (q)');set(htxt2,'HorizontalAlignment','Left','VerticalAlignment','Top','FontWeight','Bold','FontAngle','Italic','Color',[0.1 0.6 0.1]);
            [x_hard,y_hard]=select_multi_polygon('k');
            set(htxt2,'Visible','off');
            STRUC.x_hard=x_hard;
            STRUC.y_hard=y_hard;
        end
    end
    
    % Check whether to include diffraction
    if isempty(STRUC.x_hard)
       STRUC.diffraction=0; 
    end
    
    % tidy up the structures without nans at the start and end
    idnotnan=find(~isnan(STRUC.x_hard));
    if ~isempty(STRUC.x_hard)
        STRUC.x_hard=STRUC.x_hard(idnotnan(1):idnotnan(end));
        STRUC.y_hard=STRUC.y_hard(idnotnan(1):idnotnan(end));
        idnan=find(isnan(STRUC.x_hard));
        iduse=setdiff([1:length(STRUC.x_hard)],idnan(diff(idnan)==1));
        STRUC.x_hard=STRUC.x_hard(iduse);
        STRUC.y_hard=STRUC.y_hard(iduse);
    end
        
    %% Diffraction : initializes variables for diffraction at structures
    %if S.diffrac
    % to be filled in
    %end
    
    %% IDENTIFY THE STRUCTURES WHICH ARE LOCATED IN THE SEA (OR AT LAND)
    nmc=length(COAST.x_mc)-1;
    for ist=1:length(STRUC.x_hard)
        [~,icl]=min(hypot(COAST.x_mc-STRUC.x_hard(ist),COAST.y_mc-STRUC.y_hard(ist)));
        if ~isnan(STRUC.x_hard(ist))
            im1=max(icl-1,1);
            ip1=min(icl+1,nmc+1);
            if isnan(COAST.x_mc(im1)); im1=im1-1;if im1==0;im1=2;end; end
            if isnan(COAST.x_mc(ip1)); ip1=ip1+1;if ip1>nmc+1;ip1=nmc;end; end
            dirm=360-atan2d(COAST.y_mc(ip1)-COAST.y_mc(im1),COAST.x_mc(ip1)-COAST.x_mc(im1)); 
            dirstr=atan2d(STRUC.x_hard(ist)-COAST.x_mc(icl),STRUC.y_hard(ist)-COAST.y_mc(icl));
            STRUC.wetstr_mc(ist)=int8(cosd(dirstr-dirm)>0);  % for octave, logicals cannot be assigned NaNs later on
        end
    end
    STRUC.wetstr_mc(isnan(STRUC.x_hard))=nan;
    
    
    %% Permeable structures
    if S.perm
        if isfield(S,'x_perm')
            STRUC.x_perm=S.x_perm;
            STRUC.y_perm=S.y_perm;
        elseif ~isempty(S.LDBstructures)
            xy_perm=load(S.LDBstructures);
            x_perm=xy_perm(:,1)'-S.XYoffset(1);
            y_perm=xy_perm(:,2)'-S.XYoffset(2);
            STRUC.x_perm=x_perm;
            STRUC.y_perm=y_perm;
        else
        figure(11);
            axis equal;
            xl=xlim;yl=ylim;
            htxt2=text(xl(1)+0.02*diff(xl),yl(2)-0.01*diff(yl),'Add permeable structure (LMB); Next structure (RMB); Exit (q)');set(htxt2,'HorizontalAlignment','Left','VerticalAlignment','Top','FontWeight','Bold','FontAngle','Italic','Color',[0.1 0.6 0.1]);
            [x_perm,y_perm]=select_multi_polygon('k');
            set(htxt2,'Visible','off');
            STRUC.x_perm=x_perm;
            STRUC.y_perm=y_perm;
        end
    end
end
