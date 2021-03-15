close all 
clear 
Readalldata

phenotypes = ["EDV","Myo\_vol","ESV","SV"...
    ,"EF","V1","EF1","ESP","dPdtmax"...
    ,"dPdtmin","PeakP","tpeak","ET","ICT"...
    ,"IRT","tsys","QRS","AT1090","AT"];

fenotypesLVCT = abs(fenotypesLVCT(:,2:end));
fenotypesRVCT = abs(fenotypesRVCT(:,2:end));

mean_LV = nanmean(fenotypesLVCT);
mean_RV = nanmean(fenotypesRVCT);

max_LV = nanmax(fenotypesLVCT);
max_RV = nanmax(fenotypesRVCT);

min_LV = nanmin(fenotypesLVCT);
min_RV = nanmin(fenotypesRVCT);

norm_range_LV = (max_LV - min_LV)./mean_LV;
norm_range_RV = (max_RV - min_RV)./mean_RV;
norm_range_RV = [norm_range_RV NaN];


[norm_range_LV_sorted, norm_range_LV_order] = sort(norm_range_LV);
new_phenotypes = phenotypes(norm_range_LV_order);

norm_range_RV_sorted = norm_range_RV(norm_range_LV_order);


cutoffs = 0:0.05:1;

j=1;
discarded = [];
for i=1:length(cutoffs)
    discarded = [discarded sum(norm_range_LV_sorted < cutoffs(i))];
end


figure('DefaultAxesFontSize',30)
h = zeros(1,2);
h(1) = plot(cutoffs,100*discarded/max(discarded),'-','LineWidth',3,'MarkerSize',40,'Color','k')
hold

% 
x_line = [0.2 0.2];
y_line = [0 100];
% 
% Convert color code to 1-by-3 RGB array (0~1 each)
str = '#ffb142';
color = sscanf(str(2:end),'%2x%2x%2x',[1 3])/255;
 
h(2) = line(x_line,y_line,'Color','red','LineWidth',3,'LineStyle','--');

legend(h(2),'Cutoff chosen','Location','southeast');

title({"Percentage of LV phenotypes whose","normalised range is below a given threshold"},...
    'FontSize',50);
xlabel('Cutoff threshold')
ylabel('%')