function dvdt = CalcAllDeriv(t, vcon, p2)
% 
% CalcCoagDeriv calculates the derivates in the population balance
% equations of particles and thorium
%
% Note, use a new, single structure to transfer all the parameters
%
% USAGE:
%
%
% HISTORY:
%  04-05-09: First cut - based heavily on GAJ's code
%
% Adrian Burd, University of Georgia, 2009

% calculate the number of sections. The input vector is made up of
% particles (n_sections) followed by thorium (n_sections + 1);

n_sections = 0.5*(length(vcon) - 1);

% Separate out out mass and thorium from the input vector

vtmp = vcon;
vcon = vtmp(1:n_sections);
vthm = vtmp(n_sections+1:end);

% Now calculate the mass derivatives

vcon_r = vcon';      % vcon passed as column vector - make a row

vcon_shift = [0 vcon_r(1:n_sections-1)];

term1 = vcon_r * p2.b25;
term1 = vcon_r .* term1;

term2 = vcon_r * p2.b1;
term2 = term2 .* vcon_shift;

if t > 300
    term3 = (0.5*p2.growth - p2.sink_loss)* vcon;
elseif t < 300 && t > 290
    growth = (0.5 - 0.15)/(300 - 290)*t + (0.15*300 - 0.5*290)/(300-290); 
    term3 = (growth*p2.growth - p2.sink_loss)* vcon;    
else
    term3 = (0.15*p2.growth - p2.sink_loss)* vcon;
end

dvdt_mass = (term1 + term2)' + term3;

% Now calculate the thorium derivatives. This is just a straightforward
% translation of the FORTRAN code with no attempt at making things pretty
% or vectorized (beyond the obvious) at all

th_gain = zeros(1, n_sections);
th_loss = zeros(1, n_sections);

th_diss = vthm(end);
th_part = vthm(1:end-1);

sum_th_part = sum(th_part);

% Rate of change in the dissolved pool

sum_abs = sum(p2.kabs.*vcon);

d_th_diss = p2.prd/p2.uran/p2.th_bar - p2.lambda_th*th_diss - sum_abs*th_diss + p2.kdabs*sum_th_part;

% Loss from section k

for ksec = 1 : n_sections
    
    sum1 = 0.0;
    sum2 = 0.0;
    
    for isec = ksec+1 : n_sections
        sum1 = sum1 + p2.b5(isec, ksec)*vcon(isec); 
    end
    
    for isec = 1 : ksec-1
        sum2 = sum2 + p2.b3(isec, ksec)*vcon(isec);
    end
    
    th_loss(ksec) = (sum1 + sum2 + p2.b4(ksec,ksec)*vcon(ksec))*th_part(ksec);
    
end

% Gain to section k

for ksec = 1 : n_sections
    
    sum1 = 0.0;
    sum2 = 0.0;
    sum3 = 0.0;
    
    if (ksec > 2)
        for isec = 1 : ksec-2
            sum1 = sum1 + (p2.b1(isec,ksec) - p2.b3(isec,ksec-1))*th_part(isec);
        end
    end
    
    if (ksec > 1)
        sum1 = sum1*vcon(ksec-1);
    end
    
    if (ksec > 2)
        for isec = 1 : ksec-2
            sum2 = sum2 + p2.b3(isec,ksec-1)*vcon(isec);
        end
    end
    
    if (ksec > 1)
        sum2 = sum2*th_part(ksec-1);
    end
    
    if ksec > 1
        for isec = 1 : ksec-1
            sum3 = sum3 + p2.b2(isec,ksec)*vcon(ksec)*th_part(isec);
        end
    end
    
    term1a = 0;
    if ksec > 1
        term1a = p2.b4(ksec-1,ksec-1)*vcon(ksec-1)*th_part(ksec-1);
    end
    
    th_gain(ksec) = sum1 + sum2 + sum3 + term1a;
end

% Now do the settling terms

sett_loss = diag(p2.sink_loss).*th_part;

% Now do the decay, adsorption, and desorption terms

th_part_decay = p2.lambda_th*th_part;
th_part_adsor = th_diss*p2.kabs.*vcon;
th_part_deads = p2.kdabs*th_part;

% Now assemble the thorium derivatives

d_th_part = th_gain' - th_loss' - th_part_decay + th_part_adsor - th_part_deads - sett_loss;

% Now assemble the whoe derivative vector, remembering that the dissolved
% th goes at the end

dvdt = [dvdt_mass; d_th_part; d_th_diss];





