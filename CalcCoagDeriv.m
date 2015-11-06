function dvdt = CalcCoagDeriv(t, vcon, p2)
% 
% CalcCoagDeriv calculates the derivates in the population balance
% equations 
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

n_sections = length(vcon);

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

dvdt = (term1 + term2)' + term3;

