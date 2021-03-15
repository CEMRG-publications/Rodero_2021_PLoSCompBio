close all 
clear 
Readalldata

phenotypes = ["EDV","Myo\_vol","ESV","SV"...
    ,"EF","V1","EF1","ESP","dPdtmax"...
    ,"dPdtmin","PeakP","tpeak","ET","ICT"...
    ,"IRT","tsys","QRS","AT1090","AT"];

fenotypesLVCT = abs(fenotypesLVCT(:,2:end));
fenotypesRVCT = abs(fenotypesRVCT(:,2:end));

var_LV = nanvar(fenotypesLVCT);
var_RV = nanvar(fenotypesRVCT);

mean_LV = nanmean(fenotypesLVCT);
mean_RV = nanmean(fenotypesRVCT);

max_LV = nanmax(fenotypesLVCT);
max_RV = nanmax(fenotypesRVCT);

min_LV = nanmin(fenotypesLVCT);
min_RV = nanmin(fenotypesRVCT);

norm_range_LV = (max_LV - min_LV)./mean_LV;
norm_range_RV = (max_RV - min_RV)./mean_RV;
norm_range_RV = [norm_range_RV NaN]

CV_LV = var_LV./mean_LV;
CV_RV = var_RV./mean_RV;

[norm_range_LV_sorted, norm_range_LV_order] = sort(norm_range_LV);
new_phenotypes = phenotypes(norm_range_LV_order);

norm_range_RV_sorted = norm_range_RV(norm_range_LV_order);


final_indices = [1:15,17,18];


figure('DefaultAxesFontSize',30)
% 
% % Convert color code to 1-by-3 RGB array (0~1 each)
% str = '#b33939';
% color = sscanf(str(2:end),'%2x%2x%2x',[1 3])/255;
str = '#545454';
color = sscanf(str(2:end),'%2x%2x%2x',[1 3])/255;

plot(norm_range_RV_sorted(final_indices),'s','LineWidth',3,'MarkerSize',12,'Color',color,'MarkerFaceColor',color)
hold

plot(norm_range_LV_sorted(final_indices),'.','LineWidth',3,'MarkerSize',40,'Color','k')

% 
% % Convert color code to 1-by-3 RGB array (0~1 each)
% str = '#227093';
% color = sscanf(str(2:end),'%2x%2x%2x',[1 3])/255;
% 
% plot(norm_range_RV,'+','LineWidth',3,'MarkerSize',25,'Color',color)
hold off
xticklabels(new_phenotypes(final_indices))
set(gca,'xtick',1:length(new_phenotypes(final_indices)));
xtickangle(45);
% 
x_line = [0 length(final_indices)+1];
y_line = [0.2 0.2];
% 
% Convert color code to 1-by-3 RGB array (0~1 each)
str = '#ffb142';
color = sscanf(str(2:end),'%2x%2x%2x',[1 3])/255;
 
line(x_line,y_line,'Color','red','LineWidth',3,'LineStyle','--');

legend("RV","LV",'Cutoff');

title({"Range of the phenotypes of the CT","cohort for both ventricles"},...
    'FontSize',50);
ylabel('Normalised range')


% [norm_range_LV_sorted, norm_range_LV_order] = sort(norm_range_LV);
% new_phenotypes_LV = phenotypes(norm_range_LV_order);
% 
% [norm_range_RV_sorted, norm_range_RV_order] = sort(norm_range_RV);
% new_phenotypes_RV = phenotypes(norm_range_RV_order);

% figure('DefaultAxesFontSize',30)

% sp1=subplot(1,2,1);
% plot(norm_range_LV_sorted,'.','LineWidth',3,'MarkerSize',25,'Color','k')
% ylim([0 1.4])
% hold
% xticklabels(new_phenotypes_LV)
% set(gca,'xtick',1:length(new_phenotypes_LV));
% xtickangle(45);
% 
% x_line = [0 20];
% y_line = [0.2 0.2];
% 
% % Convert color code to 1-by-3 RGB array (0~1 each)
% str = '#ffb142';
% color = sscanf(str(2:end),'%2x%2x%2x',[1 3])/255;
%  
% line(x_line,y_line,'Color',color,'LineWidth',3,'LineStyle','--');
% title('LV phenotypes') 
% hold off

% sp2=subplot(1,2,2); 
% plot(norm_range_RV_sorted,'.','LineWidth',3,'MarkerSize',25,'Color','k')
% ylim(sp2,[0 1.4])
% hold
% title('RV phenotypes') 
% xticklabels(new_phenotypes_RV)
% set(gca,'xtick',1:length(new_phenotypes_RV));
% xtickangle(45);
% 
% x_line = [0 20];
% y_line = [0.2 0.2];
% 
% % Convert color code to 1-by-3 RGB array (0~1 each)
% str = '#ffb142';
% color = sscanf(str(2:end),'%2x%2x%2x',[1 3])/255;
%  
% line(x_line,y_line,'Color',color,'LineWidth',3,'LineStyle','--');
% 
% hold off



