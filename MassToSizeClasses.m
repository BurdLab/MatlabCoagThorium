function x = MassToSizeClasses(p, p2, t, y)
%
% CAlculates the mass transfer by coagulation between the various size
% classes.
%
%

% Find hte number of sections and the number of times

[n_times, n_sections] = size(y);

% Create a matrix of column vectors of mass - Each element of a column has 
% the same value. This will store the mass transfer between size classes. 

mass_within_small = zeros(n_times, 1);
mass_within_inter = zeros(n_times, 1);
mass_within_large = zeros(n_times, 1);

mass_small_2_inter = zeros(n_times, 1);
mass_inter_2_large = zeros(n_times, 1);

mass_small_2_large = zeros(n_times, 1);

% Find out how many size classes there are and in which sections their
% boundaries are

n_size_classes = 2;

% Loop through the times and figure out where the mass goes to by
% aggregation only

f1 = p.sec1_prop;
f2 = p.sec2_prop;

for i_time = 1 : n_times
 
    vol = y(i_time,:);
    vol = vol(ones(n_sections,1),:);
    
% First, do what goes to section k+1. This is made up of two terms: loss
% from k due collisions within k (beta4); loss from k due to collisions
% with all sections smaller than k (sum over beta3).

    term1  = p2.b4 .* vol .* vol;
    term1a = diag(term1);
    term1b = diag(term1a(1:n_sections-1), -1);    

% The matrix term1b is structured such that mass from a column goes to a 
% particular row. Mass that enters a section with size-class boundary needs
% to be divided into 
    
    mass_within_small(i_time) = sum(sum(term1b(1:p.sec1-1, 1:p.sec1-1)));
    mass_within_small(i_time) = mass_within_small(i_time) + term1b(p.sec1,p.sec1-1)*f1;
    
    mass_within_inter(i_time) = sum(sum(term1b(p.sec1+1:p.sec2-1,p.sec1+1:p.sec2-1)));
    mass_within_inter(i_time) = mass_within_inter(i_time) + term1b(p.sec1+1,p.sec1)*(1-f1);
    mass_within_inter(i_time) = mass_within_inter(i_time) + term1b(p.sec2, p.sec2-1)*f2;
    
    mass_small_2_inter(i_time) = term1b(p.sec1,p.sec1-1)*(1-f1) + term1b(p.sec1+1,p.sec1)f1;
    
    
    
    
    
    
    
    
    
    
    
    
    term2 = p2.b3 .* vol' .* vol;
    term2a = sum(term2);
    term2b = diag(term2a(1:n_sections-1), -1);

% Now what goes into section j > k. This is the gain to j from collisions
% with k < j. 

    term3  = p2.b2 .* vol' .* vol;
    term3a = term3';    
    
% Now what goes into section (j+1) > k: This is the loss from k due to
% collisions j>k minus the gain to j.

    term4a = p2.b5 .* vol' .* vol;
    term4b = p2.b2 .* vol' .* vol;
    
    term4c = term4a - term4b';
   
    term4d = [zeros(1,n_sections); term4c(1:n_sections-1,:)];

% Sum these terms to get total "mass to" matrix.    
    
    term5 = term1b + term2b + term3a + term4d;
    
    mass_to(:,:,i_time) = term5;
    
end