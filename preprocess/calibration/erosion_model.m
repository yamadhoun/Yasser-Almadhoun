function [dndt_cliff,Om,TWL] = erosion_model(Hs,Tp,SL,pars)

%% OUTPUT:
%      dndt_cliff   : cliff retreat rate in [m/yr]
%      Om           : cliff retreat forcing or accumulation of all wave energy over a given time in [kj/m2]
%      TWL          : total water level in [m]


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%                          Calculations                   %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Wave length and total water level
Lw = 9.81/2/pi *Tp.^2;
R2_percent_Stockdon = 1.1 *0.35 .*pars.Beta .*sqrt(Hs .*Lw) + ...
                      0.5 .*sqrt(Hs .*Lw .*(0.563 *pars.Beta.^2 + 0.004));
eta_setup = 0.18 *pars.Beta .*Tp .*sqrt(Hs *pars.g);
TWL = eta_setup + R2_percent_Stockdon + SL;                                % SL = MSL (above NAVD88) + ...
                                                                           %      tide (above MSL) + setup + R + SLR
%% Cliff erosion models
switch pars.form
    case 1 % 'Zhang'
        C = Lw ./ Tp;
        k = (2*pi)./Lw;
        Cg = C./2 .*(1+(2 .*k .*pars.d) ./sinh(2.*k .*pars.d));
        Fw = sind(pars.theta)*(pars.rho .*pars.g .*Hs.^2 ./8) .*Cg/100000;
        Fr = pars.sigma .* pars.Bc;
        for i = 1 : length(Fw)
            if Fw(i) >= Fr
                dndt_cliff(i) = pars.Kc *log(Fw(i)/Fr);
            else
                dndt_cliff(i) = 0;
            end
        end
        dndt_cliff = dndt_cliff';
        Om = Fw;

    case 2 % 'TWL'
        Om = TWL / pars.Ej;
        dndt_cliff = Om .* pars.Kc;

    case 3 % 'Trenhaile'        
        db = pars.gamma .* Hs;
        Hb = Lw .*(0.142 *tanh((2*pi *db) ./(Lw)));                        % eqn name
        dx = 0; %0.1 *ones(size(Hs));                                      % improve the first approximate
        error = 100;                                                       % per cent or fraction
        nloops = 1;
        while error > 1
        w = 1/pars.alpha .*(db + pars.Ej + (dx .*pars.alpha) - (TWL + SL));
        Om = pars.rho *Hb / (2*pars.gamma) .*exp(-pars.Xdecay .*w) / 1000;
        dndt_cliff = Om .* pars.Kc;
        dx = dx + dndt_cliff;
        if nloops ==1
            error = 100;
            dndt_prev = dndt_cliff;
        else
            error = abs(( dndt_cliff - dndt_prev ) / dndt_prev) * 100;
            dndt_prev = dndt_cliff;
        end
            nloops = 1 + nloops;
        end

    case 4 %'Energy_Flux'
        db = pars.gamma .* Hs;
        Hb = Lw .*(0.142 *tanh((2*pi *db) ./(Lw)));                        % eqn name
        Cg = sqrt(pars.g*db);
        dx = 0; %0.1 *ones(size(Hs));                                      % improve the first approximate
        error = 100;                                                       % per cent or fraction
        nloops = 1;
        while error > 1
        w = 1/pars.alpha .*(db + pars.Ej + (dx .*pars.alpha) - (TWL + SL));
        Om = (1/8 *pars.g *pars.rho *Hb.^2 .*Cg) .*exp(-pars.Xdecay .*w)/10000;
        dndt_cliff = Om .* pars.Kc;
        dx = dx + dndt_cliff;
        if nloops ==1
            error = 100;
            dndt_prev = dndt_cliff;
        else
            error = abs(( dndt_cliff - dndt_prev ) / dndt_prev) * 100;
            dndt_prev = dndt_cliff;
        end
            nloops = 1 + nloops;
        end

    case 5 % 'Hackney'
        Om = 1/8 *pars.rho *pars.g *(Hs + SL).^2 /10000;
        dndt_cliff = Om .* pars.Kc;
end

end
