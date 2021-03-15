close all
h=figure;

set(h,'PaperOrientation','landscape');
set(h,'PaperPosition', [1 1 28 19]);

bp=boxplot(normalisedresiduals); 

ylabel('um');
xlabel('Tags');

title('Average varifold error per tag','FontSize',30)
set(gca,'FontSize',30);
set(bp,'LineWidth', 2);
% print(gcf, '-dpdf', 'boxplots_mass_lv.pdf');

