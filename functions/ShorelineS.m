function [S,O]=ShorelineS(S)
% MODEL : ShorelineS
%
% This model computes shoreline changes as a result of gradients in alongshore
% sediment transport for arbitrary shaped coastlines.
%
% INPUT:
%     S       data structure wiht fields:
%              .
%
% made by:      J.A. Roelvink (2016-present) - IHE Delft
% modified by:  B.J.A. Huisman (2017-present) - Deltares
% modified by:  A.M. Elghandour (2018-present) - IHE Delft
% modified by:  M.E. Ghonim (2019) - IHE Delft
% modified by:  C.M. Mudde (2019) - TU-Delft
% modified by:  J. Reyns (2018-present) - IHE Delft, Deltares
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
%profile on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DEFAULT INPUT PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[S]       = initialize_defaultvalues(S);
[S]       = initialize_randomgenerator(S);
[TIME]    = initialize_time(S);
[COAST]   = prepare_coastline(S);
[DUNE]    = prepare_dunes(S,COAST);
[CLIFF]   = prepare_cliff(S,COAST);                                        % <- new function by YA2022
[CC]      = prepare_climatechange(S,TIME);
[WAVE]    = prepare_waveconditions(S,TIME);
% global WAVE
[STRUC]   = prepare_structures(S,COAST);
[NOUR]    = prepare_nourishment(S,COAST,STRUC);
[TRANSP]  = prepare_transport(S);
[SPIT]    = prepare_spit(S);
[CHANNEL] = prepare_channel(S);
[DELTA]   = prepare_delta(S);
[BATHY]   = initialize_bathyupdate(S);
[FORMAT]  = initialize_plot(S,COAST);    
[O,V]     = initialize_output();
%[xl,yl,phirf,V,iwtw,CLplot,CLplot2,iint,BWplot,BWplot2,innt,ds_cl,qwave,qwind,time,step,int,bermW,x_trans,y_trans,n_trans] = initialize_plot_variables(S,COAST.x_mc,COAST.y_mc,STRUC.x_hard,STRUC.y_hard,x_dune,y_dune,TIME.timenum0); % <- new function by AE&MG
%[WndCF,WndCL]=prepare_windconditions(S);   % <- new function by AE&MG
TIME.adt=1/365/24;
eps=0.1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Loop over time
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('  Loop over time \n'); 
while  TIME.tnow<TIME.tend 
    TIME.it=TIME.it+1; 
    TIME.nt=TIME.it; 
    
    %% Determine climate change effects
    [CC]=introduce_climatechange(CC, TIME);
    
    if S.debug==2
        if ~isoctave, warning off, end 
        save('debug.mat'); 
        if ~isoctave, warning on, end
    end 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% PHASE 0 : GRID                                             %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [COAST,STRUC,GROYNE,TRANSP]=prepare_grid_groyne(COAST,STRUC,TRANSP,S.yesplot);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% PHASE 1 : TRANSPORT                                        %%
    %% loop over coastline sections                               %%
    %% compute sediment transport                                 %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for i_mc=1:COAST.n_mc 
                       
        %% Alongshore coordinate s; regridding where necessary
        [COAST]=make_sgrid_mc(COAST,i_mc);
        
        %% Alongshore coordinate s; regridding where necessary: for cliffs
        [CLIFF]=make_sgrid_mc_cliff(CLIFF,i_mc);                       % <- new function by YA2022

        %% Interpolate wave conditions along the coast
        [WAVE]=introduce_wave(S,WAVE,TIME,COAST,CC);
        
        %% Get coastline orientation PHIc
        [COAST]=get_coastline_orientation(COAST);
        
        %% Get foreshore orientation PHIf
        [COAST]=get_foreshore_orientation(COAST);

        %% Get refracted waves
        % From here on, WAVE.PHItdp and WAVE.HStdp are always given at the toe of the
        % dynamic profile (TDP) and are always with a size of 1 by n
        [WAVE]=wave_refraction(WAVE,COAST,TRANSP);
        
        %% Introduce permeable structures
        [WAVE.HStdp]=introduce_perm_structures(STRUC,WAVE,COAST);
        
        %% Wave diffraction
        if STRUC.diffraction==1
           [STRUC,WAVE]=wave_diffraction(STRUC,COAST,WAVE);           
        end
        
        %% Introduce the wind conditions to the domain
        %[PHIwnd,S]=introduce_wind(S,WndCF,WndCL,TIME.tnow); 
        
        %% Introduce water levels conditions to the domain
        %[SWL]=introduce_tide(S,Wat_cf,Wat_cl,TIME.tnow);    
        
        %% Relative angle of the offshore and nearshore waves 
        WAVE.dPHIo=atan2d(sind(COAST.PHIc-WAVE.PHIo),cosd(COAST.PHIc-WAVE.PHIo));
        WAVE.dPHItdp=atan2d(sind(COAST.PHIc-WAVE.PHItdp),cosd(COAST.PHIc-WAVE.PHItdp));        % limit refraction
        
        %% Wave height in the nearshore (due to refraction and shoaling)
        [WAVE]=wave_breakingheight(WAVE,TRANSP);
        
        %% Nearshore wave direction at point of breaking
        WAVE.PHIbr=mod(COAST.PHIc-WAVE.dPHIbr,360); 
        
        %% Critical angles for transport -> WAVE.dPHIcrit & TRANSP.QSmax
        [WAVE,TRANSP]=wave_angles(COAST,WAVE,TRANSP,STRUC);
        
        %% Longshore Transport
        [TRANSP]=transport(TRANSP,WAVE,STRUC);
        
        %% Shadowing effect on Transport
        [TRANSP,WAVE]=transport_shadow_treat(COAST,STRUC,WAVE,TRANSP);
        %[TRANSP]=transport_shadow_dune(DUNE,STRUC,WAVE,TRANSP);       
        
        %% Upwind correction for high-angle -> using wave angle at 'tdp'
        % TRANSP.QS        : Transport rates in grid cells with upwind correction [1xN] (in [m3/yr] including pores); -> TRANSP.QS & TRANSP.im3 & TRANSP.ip3
        [TRANSP]=get_upwindcorrection(COAST,WAVE,TRANSP,SPIT); % using dPHItdp (based on wave conditions at the nearshore location)        
        
        %% Coastline cells intersected by hard structures
        % STRUC.structS : Indices of structures on the coastline (for xy-points)
        [STRUC,TRANSP]=find_struct(COAST,STRUC,TRANSP);
        
        %% Sand Bypassing and Transmission (GROYNE.QS)
        [GROYNE]=transport_bypass(TRANSP,WAVE,COAST,STRUC,GROYNE);
        
        %% Dune evolution
        %[qs,qw,qss,qww]=dune_evolution(S,x_dune,y_dune,x_dune0,y_dune0,COAST.x_mc,COAST.y_mc,COAST.x_mc0,COAST.y_mc0,i_mc,Dfelev,PHIwnd,SWL,COAST.cyclic,WAVE.HStdp,str_indx,shadowD,shadowD_h,WAVE.dPHIcor); % <- new function by AE&MG
        
        %% Adaptive time step based on transport of separate coastline sections
        [TIME]=get_timestep(TIME,COAST,TRANSP);
        
        %% Boundary condition
        [TRANSP]=transport_boundary_condition(TRANSP,COAST,TIME,GROYNE);
        
        %% Collect QS, s and WAVE.PHItdp in QS_mc, WAVE.PHIo_mc and s_mc -> stores the data of this coastline section in 'mc'
        [COAST,WAVE,TRANSP]=collect_variables(COAST,WAVE,TRANSP,i_mc);

        %% debugging plot for QS and SPHI
        plot_debug(S.debug,COAST,WAVE,TRANSP); % only used when S.debug=1;

        %% cliff retreat rate
        [CLIFF]=cliff_retreat_rate(S,COAST,CLIFF,WAVE,TIME);               % <- new function by YA2022 
        
    end
    TIME.adt_record(TIME.it+1)=TIME.adt;
    
    %% Final timestep for coastline update 
    [TIME,WAVE] = get_coastlineupdate_timestep(TIME, BATHY, WAVE, FORMAT);    
    
    if S.debug==2
        if ~isoctave, warning off, end
        save('debug2.mat');
        if ~isoctave, warning on, end
    end    
    
    %% data assimilation
    % [DA]=get_dataassimilation(DA,TIME,WAVE,'update');

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% PHASE 2 : COASTLINE CHANGE                                 %%
    %% Compute coastline and dune change                          %%
    %% Add nourishments to the coast                              %%
    %% Re-connect coastal sections at groynes                     %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    COAST.x_mc1=COAST.x_mc;
    COAST.y_mc1=COAST.y_mc;
    COAST.n_mc1=COAST.n_mc;
    
    for i_mc=1:COAST.n_mc1

        %% Retrieve data of individual segment (i_mc)
        [COAST,WAVE,TRANSP]=get_segmentdata(COAST,WAVE,TRANSP,i_mc);
        
        %% Nourishment
        % NOUR.nour   : Indices of grid cells with a nourishment & number of nourishments (same length as coastline COAST.x and COAST.y)
        [NOUR]=get_nourishments(TIME,COAST,NOUR);
        
        %% Dune foot change
        %[x_dune,y_dune,S]=dune_foot_change(S,DUNE); % <- new function by AE&MG
        
        %% cliff coordinate update
        [CLIFF]=update_cliff(CLIFF,GROYNE,TIME);

        %% Coastline change  ++cliff retrat volume drift contribution
        [COAST,GROYNE]=coastline_change(S,COAST,TRANSP,STRUC,GROYNE,TIME,NOUR,CC,CLIFF); % <- add cliff Vol same as NOUR
        
        %% insert new section into x_mc and y_mc
        [COAST.x_mc,COAST.y_mc]=insert_section(COAST.x,COAST.y,COAST.x_mc,COAST.y_mc,COAST.i_mc);
        
    end % loop over sections
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% PHASE 3 : OTHER PROCESSES                                  %%
    %% Overwash process                                           %%
    %% Move channels                                              %%
    %% Splitting and Merging coastlines                           %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %% Overwash process
    COAST.n_mc=length(find(isnan(COAST.x_mc)))+1;  
    for i_mc=1:COAST.n_mc
       if COAST.n_mc>=i_mc && isempty(find(isnan(COAST.x),1)) 
            %[OP]=overwash_potential(WAVE.TP,HSo,S.Bheight,S.tanbeta) % added by Ahmed
            [COAST,SPIT] = find_spit_width_mc(i_mc,COAST,WAVE,SPIT,STRUC,TRANSP);
        end
    end
    
    %% reconnect coastlines through groynes
    [COAST]=get_reconnectedgroynes(COAST,GROYNE); 
    
    %% remove spikes from small 'bubble islands' %  & clean up redundant NaNs
    [COAST]=get_spikesremoved(COAST);
    
    %% move channel
    [CHANNEL,COAST]=move_channel2(CHANNEL,COAST,WAVE,TIME);
    
    %% merge multiple coastline sections -> only x_mc and y_mc are currently merged.
    COAST.n_mc=length(find(isnan(COAST.x_mc)))+1;  
    for i_mc=1:COAST.n_mc
        [COAST]=merge_coastlines(COAST,i_mc);
    end
    try
        [COAST]=merge_coastlines_mc(COAST);
    catch
        error('problem in merge_coastlines_mc')
    end
    
    %% clean up redundant NaNs
    [COAST]=cleanup_nans(COAST);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% PHASE 4 : PLOTTING AND STORING COASTLINES                  %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    %% Plotting                                     
    [V,FORMAT,TIME]=plot_coast(CHANNEL,STRUC,COAST,CLIFF,WAVE,TIME,TRANSP,FORMAT,V);
    
    %% Apply shoreline change due to tide           
    [FORMAT,COAST]=update_shoreline(S,COAST,FORMAT);
    
    %% Bathymetry update & Plot                              
    [S,BATHY]=update_bathy(S,BATHY,TIME,COAST);
    
    %% Store data in OUTPUT.MAT (information before coastline update)
    [O,TIME.itout]=save_shorelinesSS(O,S,TIME,COAST,CLIFF,WAVE,TRANSP,STRUC,NOUR,GROYNE);
    
    %% Store all shoreline data
    if ~isoctave, warning off, end
        save(fullfile(pwd,FORMAT.outputdir,'output.mat'),'O','S');
    if ~isoctave, warning on, end
    
    %% Next time step                                    
    TIME.tprev=TIME.tnow;
    TIME.tnow=TIME.tnow + TIME.dt*365;     %calculate the current time after each time step
    
    S.times(TIME.it+1)=TIME.tnow;
    fprintf('   %s \n',datestr(TIME.tnow,'yyyy-mm-dd HH:MM'));
end
    
fprintf('  Post-process results \n');
    
%% Extract shorelines cooridnates /figures at specific dates (x_mc0 and y_mc0)
extract_shoreline(S,STRUC,COAST,FORMAT);

%% for Plot beach berm witdh variation against time
% extract_berm(S,step,qwind,qwave,time,bermW,WBplot2,ds_cl,CSplot2,n_trans);

%% videos
if S.video==1
    make_video(S,V);
end
% Area_BSS_check(S)

%% data assimilation
% [DA]=get_dataassimilation(DA,TIME,WAVE,'save');
%profsave(profile('info'),'profile_results')
%profexplore
end