close all

hwhole = histogram(qualitySJ, 'FaceColor', 'red', 'EdgeColor', 'black', 'FaceAlpha', 0.9);
hwhole.BinWidth=2e-2;
tit = title('     SJ values for mesh #02');
set(tit,'FontSize',44);
set(gca,'FontSize',66);
h=gcf;
set(h,'PaperOrientation','landscape');
ylim([0 250000]);