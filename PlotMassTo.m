function [F, mass_to, norm_data] = PlotMassTo(p, p2, t, y)
%
% Plots where the mass (volume) goes
%

[n_times, n_sections] = size(y);

% Create a matrix of column vectors of mass - Each element of a column has 
% the same value

mass_to = zeros(n_sections, n_sections, n_times);


% Loop through the times and figure out where the mass goes to by
% aggregation only

for i_time = 1 : n_times
 
    vol = y(i_time,:);
    vol = vol(ones(n_sections,1),:);

    
% First, do what goes to section k+1. This is made up of two terms: loss
% from k due collisions within k (beta4); loss from k due to collisions
% with all sections smaller than k (sum over beta3).

    term1 = p2.b4 .* vol .* vol;
    
    term1a = diag(term1);
    term1b = diag(term1a(1:n_sections-1), -1);    
    
    term2  = p2.b3 .* vol' .* vol;
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

% Plot one frame of the mass to array

tmp1 = reshape(mass_to(:,:,1), n_sections, n_sections);%

norm_data_start = MyPlotMassTo(p, tmp1, n_sections);

tmp2 = reshape(mass_to(:,:,end), n_sections, n_sections);

norm_data = MyPlotMassTo(p, tmp2, n_sections);


% Make a movie of how things change

%F = MakeMovieMassTo(p, mass_to, n_sections, n_times);
F = [];

% How much goes from a section into the different size classes. First, how
% much material goes from the lower sections into the intermediate size
% class.

% mass_to_intermediate1 = sum(tmp1(p.sec1:p.sec2-1, 1:p.sec1))./sum(tmp1(:, 1:p.sec1));
% mass_to_intermediate2 = sum(tmp2(p.sec1:p.sec2-1, 1:p.sec1))./sum(tmp2(:, 1:p.sec1));
% 
% mass_to_large1 = sum(tmp1(p.sec2:end, 1:p.sec1))./sum(tmp1(:, 1:p.sec1));
% mass_to_large2 = sum(tmp2(p.sec2:end, 1:p.sec1))./sum(tmp2(:, 1:p.sec1));


mass_stays1 = sum(tmp1(1:p.section(1)-1, 1:p.section(1)-1))./(sum(tmp1(:, 1:p.section(1)-1)));
mass_stays2 = sum(tmp2(1:p.section(1)-1, 1:p.section(1)-1))./(sum(tmp2(:, 1:p.section(1)-1)));

mass_to_intermediate1 = sum(tmp1(p.section(1):p.section(2)-1, 1:p.section(1)-1))./sum(tmp1(:, 1:p.section(1)-1));
mass_to_intermediate2 = sum(tmp2(p.section(1):p.section(2)-1, 1:p.section(1)-1))./sum(tmp2(:, 1:p.section(1)-1));

mass_to_large1 = sum(tmp1(p.section(2):end, 1:p.section(1)-1))./sum(tmp1(:, 1:p.section(1)-1));
mass_to_large2 = sum(tmp2(p.section(2):end, 1:p.section(1)-1))./sum(tmp2(:, 1:p.section(1)-1));

% figure
% subplot(2,1,1)
% bar(y(1,:)')
% subplot(2,1,2)
% bar(y(end,:)')

% figure
% subplot(2,1,1)
% bar([mass_stays1; mass_to_intermediate1; mass_to_large1]')
% legend('Mass Stays', 'To Intermediate', 'To Large', 'Location', 'EastOutside')
% set(gca, 'YLim', [0 1])
% xlabel('Originating Section')
% ylabel('Fraction of Mass Flow')
% 
% subplot(2,1,2)
% bar([mass_stays2; mass_to_intermediate2; mass_to_large2]')
% legend('Mass Stays', 'To Intermediate', 'To Large', 'Location', 'EastOutside')
% set(gca, 'YLim', [0 1])
% xlabel('Originating Section')
% ylabel('Fraction of Mass Flow')

% figure
% subplot(2,1,1)
% bar([sum(mass_stays1) sum(mass_to_intermediate1) sum(mass_to_large1)])
% subplot(2,1,2)
% bar([sum(mass_stays2) sum(mass_to_intermediate2) sum(mass_to_large2)])



% Plot out where most of the mass is within these size classes

% figure
% subplot(2,1,1)
% pvol1 = [y(1,1:p.section(1)-1)/sum(y(1,1:p.section(1)-1)) y(1,p.section(1):p.section(2)-1)/sum(y(1,p.section(1):p.section(2)-1)) y(1,p.section(2):end)/sum(y(1,p.section(2):end))];
% bar(pvol1')
% hold on
% plot([p.section(1) p.section(1)], [0 1.0], 'r--', [p.section(2) p.section(2)], [0 1.0], 'r--')
% xlabel('Section')
% ylabel('Proportion of particle solid volume')
% %set(gca, 'XTick', [0.5 : p.section(2)+0.5])
% %set(gca, 'XTickLabel', [1 : p.section(2)-1])
% 
% subplot(2,1,2)
% pvol2 = [y(end,1:p.section(1)-1)/sum(y(end,1:p.section(1)-1)) y(end,p.section(1):p.section(2)-1)/sum(y(end,p.section(1):p.section(2)-1)) y(end,p.section(2):end)/sum(y(end,p.section(2):end))];
% bar(pvol2')
% hold on
% plot([p.section(1) p.section(1)], [0 1.0], 'r--', [p.section(2) p.section(2)], [0 1.0], 'r--')
% xlabel('Section')
% ylabel('Proportion of particle solid volume')
% %set(gca, 'XTick', [0.5 : p.section(2)+0.5])
% %set(gca, 'XTickLabel', [1 : p.section(2)-1])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function norm_data = MyPlotMassTo(p, data, n_sections)

% first create the coordinates of the patches. This calculates the
% corrdinates of the vertices of each cell. Or, to think of it another way,
% the coordinates of each intersection of lines. 

xdata = zeros(4, n_sections*n_sections);
ydata = zeros(4, n_sections*n_sections);

for ix = 0 : n_sections-1
    
    xtemp = [ix; ix+1; ix+1; ix];
    
    for jy = 0 : n_sections-1
        
        ytemp = [jy; jy; jy+1; jy+1];
        
        xdata(:,n_sections*ix + (jy+1)) = xtemp;
        ydata(:,n_sections*ix + (jy+1)) = ytemp;
        
    end
    
end

% Now make a color map based on the mass-to data. This is normalized so
% that each column sums to unity. The end result will give the proportion
% of material from a given section that ends up in a different one. Note
% that we have to do a little trick because the last column in data is all
% zeros, so we have to omit it when we do the normalization (or else we
% will get NaN) and then add it back in later before we make the color map.

sum_data = sum(data);
sum_data = sum_data(ones(n_sections, 1), :);

norm_data = data(:, 1:n_sections-1)./sum_data(:, 1:n_sections-1);
norm_data = [norm_data zeros(n_sections, 1)];

cdata = reshape(norm_data, n_sections*n_sections, 1);

% Now plot the figure. First we plot the figure with by coloring each patch
% white. Then we come in afterwards and use the colormap to overwrite the
% color in each patch. 

hf1 = figure;

my_col_map = makeColorMap([1 1 1],[0 0 1],40);

%p = patch(xdata, ydata, norm_data)
pa = patch(xdata, ydata, 'w');
%set(gca, 'CLim', [0.01 1])
colormap('pink')
set(gca, 'FontName', 'Helvetica', 'FontSize', 18)
set(pa, 'FaceColor', 'flat', 'FaceVertexCData', cdata, 'CDataMapping', 'scaled')
colorbar
set(gca, 'FontName', 'Helvetica', 'FontSize', 18)
ylabel('Receiving Section', 'FontName', 'Helvetica', 'FontSize', 18)
xlabel('Originating Section', 'FontName', 'Helvetica', 'FontSize', 18)
%title('Relative transfer from originating to receiving section')
set(gca, 'XTick', [0.5 : 2 : n_sections-0.5])
set(gca, 'YTick', [0.5 : 2 : n_sections-0.5])
set(gca, 'XTickLabel', [1 : 2 : n_sections-1])
set(gca, 'YTickLabel', [1 : 2 : n_sections-1])
hold on
plot([0 p.section(1)-1], [p.section(1)-1, p.section(1)-1], 'w', 'LineWidth', 1)
plot([p.section(1)-1 p.section(1)-1], [0 p.section(1)-1], 'w', 'LineWidth', 1)
plot([0 p.section(2)-1], [p.section(2)-1, p.section(2)-1], 'w', 'LineWidth', 1)
plot([p.section(2)-1 p.section(2)-1], [0 p.section(2)-1], 'w', 'LineWidth', 1)
orient landscape

figure
colormap('gray')
waterfall(norm_data)
set(gca, 'FontName', 'Helvetica', 'FontSize', 18)
view([-66 14])
xlabel('Originating', 'FontName', 'Helvetica', 'FontSize', 18)
ylabel('Receiving', 'FontName', 'Helvetica', 'FontSize', 18)
orient landscape

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function F = MakeMovieMassTo(p, data, n_sections, n_times)

% first create the coordinates of the patches

xdata = zeros(4, n_sections*n_sections);
ydata = zeros(4, n_sections*n_sections);

for ix = 0 : n_sections-1
    
    xtemp = [ix; ix+1; ix+1; ix];
    
    for jy = 0 : n_sections-1
        
        ytemp = [jy; jy; jy+1; jy+1];
        
        xdata(:,n_sections*ix + (jy+1)) = xtemp;
        ydata(:,n_sections*ix + (jy+1)) = ytemp;
        
    end
    
end

move_to_array = zeros(size(data));

for i_time = 1 : n_times
    
    tmp = reshape(data(:,:,i_time), n_sections, n_sections);
    
    sum_data = sum(tmp);
    sum_data = sum_data(ones(n_sections, 1), :);

    norm_data = tmp(:, 1:n_sections-1)./sum_data(:, 1:n_sections-1);
    norm_data = [norm_data zeros(n_sections, 1)];

    move_to_array(:,:,i_time) = norm_data;
end

% Make movie

figure
pa = patch(xdata, ydata, 'w');
axis tight
set(gca, 'XTick', [0.5 : 2 : n_sections-0.5])
set(gca, 'YTick', [0.5 : 2 : n_sections-0.5])
set(gca, 'XTickLabel', [1 : 2 : n_sections-1])
set(gca, 'YTickLabel', [1 : 2 : n_sections-1])
set(gca, 'nextplot', 'replacechildren')

for j = 1 : n_times

    d = reshape(move_to_array(:,:,j), n_sections, n_sections);
    
    cdata = reshape(d, n_sections*n_sections, 1);

    pa = patch(xdata, ydata, 'w');


    colormap('pink')
    set(pa, 'FaceColor', 'flat', 'FaceVertexCData', cdata, 'CDataMapping', 'scaled')
    colorbar
    hold on

    plot([0 p.section(1)], [p.section(1), p.section(1)], 'w')
    plot([p.section(1) p.section(1)], [0 p.section(1)], 'w')
    plot([0 p.section(2)], [p.section(2), p.section(2)], 'w')
    plot([p.section(2) p.section(2)], [0 p.section(2)], 'w')

    F(j) = getframe(gcf);
    pause(0.1)

end  

% h2 = figure;
% movie(h2, F, 4)
