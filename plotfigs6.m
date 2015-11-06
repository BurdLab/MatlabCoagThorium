function x = plotfigs6(norm_data_rect, norm_data_curv)

x = 1;

hf = figure;
subplot(2,1,1)
colormap('gray')
waterfall(norm_data_rect)
set(gca, 'FontName', 'Helvetica', 'FontSize', 18)
view([-66 14])
xlabel('Originating', 'FontName', 'Helvetica', 'FontSize', 18)
ylabel('Receiving', 'FontName', 'Helvetica', 'FontSize', 18)
text(30,5,0.8, 'a.', 'FontName', 'Helvetica', 'FontSize', 18)

subplot(2,1,2)
colormap('gray')
waterfall(norm_data_curv)
set(gca, 'FontName', 'Helvetica', 'FontSize', 18)
view([-66 14])
xlabel('Originating', 'FontName', 'Helvetica', 'FontSize', 18)
ylabel('Receiving', 'FontName', 'Helvetica', 'FontSize', 18)
text(30,5,0.8, 'b.', 'FontName', 'Helvetica', 'FontSize', 18)

orient tall

print -dpdf figure7.pdf