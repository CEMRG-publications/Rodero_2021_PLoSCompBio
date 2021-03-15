F = [5,7,14,15,18,19];
M = [1,2,3,4,6,8,9,10,11,12,13,16,17,20];

%% Volume ALL

close all
h=figure;
hold on

set(h,'PaperOrientation','landscape');
set(h,'PaperPosition', [1 1 28 19]);
%% LV 
mean_vol = 143; % Petersen 2017
std_vol = 34;

patch([0.75 1.25 1.25 0.75], [mean_vol-std_vol mean_vol-std_vol,...
    mean_vol+std_vol mean_vol+std_vol], [0 1 0],...
    'FaceAlpha',0.3,'LineStyle','none')
%% RV
mean_vol = 154; % Petersen 2017
std_vol = 40;

patch([1+0.75 1+1.25 1+1.25 1+0.75], [mean_vol-std_vol mean_vol-std_vol,...
    mean_vol+std_vol mean_vol+std_vol], [0 1 0],...
    'FaceAlpha',0.3,'LineStyle','none')
%% LA
mean_vol = 74; % Petersen 2017
std_vol = 22;

patch([2+0.75 2+1.25 2+1.25 2+0.75], [mean_vol-std_vol mean_vol-std_vol,...
    mean_vol+std_vol mean_vol+std_vol], [0 1 0],...
    'FaceAlpha',0.3,'LineStyle','none')
%% RA
mean_vol = 80; % Petersen 2017
std_vol = 25;

patch([3+0.75 3+1.25 3+1.25 3+0.75], [mean_vol-std_vol mean_vol-std_vol,...
    mean_vol+std_vol mean_vol+std_vol], [0 1 0],...
    'FaceAlpha',0.3,'LineStyle','none')
%% Boxplots

LV_vol = [120.057 128.613 152.305 145.15 118.406 86.1567 113.198 142.959 154.357 154.98 142.291 101.844 139.445 177.252 103.964 140.259 118.204 92.8903 115.664 109.466]';
RV_vol = [147.499 178.024 191.621 189.982 101.945 88.1846 140.046 178.696 142.801 193.332 168.791 121.313 173.173 201.839 122.588 145.593 157.433 119.224 116.294 148.11]';
LA_vol = [68.5369 49.605 63.8641 63.427 47.4743 65.5117 67.8205 88.4302 96.2418 85.2434 96.6128 51.9832 74.0667 87.8193 65.5355 69.0059 68.9715 43.3464 72.399 54.6153]';
RA_vol = [112.379 91.5234 80.032 97.1162 50.4545 63.2971 110.446 108.549 89.184 111.711 95.6747 65.0821 85.029 98.6811 70.7274 80.2781 74.7992 52.8519 80.433 71.5133]';

bp=boxplot([LV_vol,RV_vol,LA_vol,RA_vol],'Labels',...
    {'LV','RV','LA','RA'});

ylabel('Volume (mL)','FontSize',44);

title('Chamber volumes for all the meshes','FontSize',44)
set(gca,'FontSize',70);
set(bp,'LineWidth',2);
print(gcf, '-dpdf', 'boxplots_vol_all.pdf');
hold off

% %% Volume MALES
% 
% close all
% h=figure;
% hold on
% 
% set(h,'PaperOrientation','landscape');
% set(h,'PaperPosition', [1 1 28 19]);
% %% LV 
% mean_vol = 166; % Petersen 2017
% std_vol = 32;
% 
% patch([0.75 1.25 1.25 0.75], [mean_vol-std_vol mean_vol-std_vol,...
%     mean_vol+std_vol mean_vol+std_vol], [0 1 0],...
%     'FaceAlpha',0.3,'LineStyle','none')
% %% RV
% mean_vol = 182; % Petersen 2017
% std_vol = 36;
% 
% patch([1+0.75 1+1.25 1+1.25 1+0.75], [mean_vol-std_vol mean_vol-std_vol,...
%     mean_vol+std_vol mean_vol+std_vol], [0 1 0],...
%     'FaceAlpha',0.3,'LineStyle','none')
% %% LA
% mean_vol = 78; % Petersen 2017
% std_vol = 23;
% 
% patch([2+0.75 2+1.25 2+1.25 2+0.75], [mean_vol-std_vol mean_vol-std_vol,...
%     mean_vol+std_vol mean_vol+std_vol], [0 1 0],...
%     'FaceAlpha',0.3,'LineStyle','none')
% %% RA
% mean_vol = 93; % Petersen 2017
% std_vol = 27;
% 
% patch([3+0.75 3+1.25 3+1.25 3+0.75], [mean_vol-std_vol mean_vol-std_vol,...
%     mean_vol+std_vol mean_vol+std_vol], [0 1 0],...
%     'FaceAlpha',0.3,'LineStyle','none')
% %% Boxplots
% 
% bp=boxplot([LV_vol(M),RV_vol(M),LA_vol(M),RA_vol(M)],'Labels',...
%     {'Left ventricle','Right ventricle','Left atrium','Right atrium'});
% 
% ylabel('Volume (mL)','FontSize',24);
% 
% title('Chamber volumes for meshes of males','FontSize',24)
% set(gca,'FontSize',20);
% set(bp,'LineWidth', 2);
% print(gcf, '-dpdf', 'boxplots_vol_males.pdf');
% hold off
% 
% %% Volume FEMALES
% 
% close all
% h=figure;
% hold on
% 
% set(h,'PaperOrientation','landscape');
% set(h,'PaperPosition', [1 1 28 19]);
% %% LV 
% mean_vol = 124; % Petersen 2017
% std_vol = 21;
% 
% patch([0.75 1.25 1.25 0.75], [mean_vol-std_vol mean_vol-std_vol,...
%     mean_vol+std_vol mean_vol+std_vol], [0 1 0],...
%     'FaceAlpha',0.3,'LineStyle','none')
% %% RV
% mean_vol = 130; % Petersen 2017
% std_vol = 24;
% 
% patch([1+0.75 1+1.25 1+1.25 1+0.75], [mean_vol-std_vol mean_vol-std_vol,...
%     mean_vol+std_vol mean_vol+std_vol], [0 1 0],...
%     'FaceAlpha',0.3,'LineStyle','none')
% %% LA
% mean_vol = 70; % Petersen 2017
% std_vol = 21;
% 
% patch([2+0.75 2+1.25 2+1.25 2+0.75], [mean_vol-std_vol mean_vol-std_vol,...
%     mean_vol+std_vol mean_vol+std_vol], [0 1 0],...
%     'FaceAlpha',0.3,'LineStyle','none')
% %% RA
% mean_vol = 69; % Petersen 2017
% std_vol = 17;
% 
% patch([3+0.75 3+1.25 3+1.25 3+0.75], [mean_vol-std_vol mean_vol-std_vol,...
%     mean_vol+std_vol mean_vol+std_vol], [0 1 0],...
%     'FaceAlpha',0.3,'LineStyle','none')
% %% Boxplots
% 
% bp=boxplot([LV_vol(F),RV_vol(F),LA_vol(F),RA_vol(F)],'Labels',...
%     {'Left ventricle','Right ventricle','Left atrium','Right atrium'});
% 
% ylabel('Volume (mL)','FontSize',24);
% 
% title('Chamber volumes for meshes of females','FontSize',24)
% set(gca,'FontSize',20);
% set(bp,'LineWidth', 2);
% print(gcf, '-dpdf', 'boxplots_vol_females.pdf');
% hold off

% %% Mass LV ALL, MALES AND FEMALES
% 
% close all
% h=figure;
% hold on
% 
% set(h,'PaperOrientation','landscape');
% set(h,'PaperPosition', [1 1 28 19]);
% %% ALL
% mean_mass_all = 85; % Petersen 2017
% std_mass_all = 24;
% 
% patch([0.75 1.25 1.25 0.75], [mean_mass_all-std_mass_all mean_mass_all-std_mass_all,...
%     mean_mass_all+std_mass_all mean_mass_all+std_mass_all], [0 1 0],...
%     'FaceAlpha',0.3,'LineStyle','none')
% %% MALES
% 
% mean_mass_males = 103;
% std_mass_males = 21;
% 
% patch([1+0.75 1+1.25 1+1.25 1+0.75], [mean_mass_males-std_mass_males mean_mass_males-std_mass_males,...
%     mean_mass_males+std_mass_males mean_mass_males+std_mass_males], [0 1 0],...
%     'FaceAlpha',0.3,'LineStyle','none')
% %% FEMALES
% 
% mean_mass_females = 70;
% std_mass_females = 13;
% 
% patch([2+0.75 2+1.25 2+1.25 2+0.75], [mean_mass_females-std_mass_females mean_mass_females-std_mass_females,...
%     mean_mass_females+std_mass_females mean_mass_females+std_mass_females], [0 1 0],...
%     'FaceAlpha',0.3,'LineStyle','none')
% 
% %% Boxplots
% 
% bp=boxplot([LV_mass',LV_mass(M)',LV_mass(F)'],...
%     [ones(1,20),2*ones(1,length(M)),3*ones(1,length(F))],'Labels',...
%     {'Whole cohort','Males','Females'});
% 
% ylabel('Mass (g)','FontSize',24);
% 
% title('Left ventricular mass','FontSize',24)
% set(gca,'FontSize',20);
% set(bp,'LineWidth', 2);
% print(gcf, '-dpdf', 'boxplots_mass_lv.pdf');
% hold off
% 
% close all

%% Mass LV ALL COMPARED W/ LITERATURE
% 
% close all
% h=figure;
% hold on
% 
% set(h,'PaperOrientation','landscape');
% set(h,'PaperPosition', [1 1 28 19]);
% 
% 
% %     mean_mass_females+std_mass_females mean_mass_females+std_mass_females], [0 1 0],...
% %     'FaceAlpha',0.3,'LineStyle','none')
% % 
% % %% Boxplots
% % 
% % bp=boxplot([LV_mass',LV_mass(M)',LV_mass(F)'],...
% %     [ones(1,20),2*ones(1,length(M)),3*ones(1,length(F))],'Labels',...
% %     {'Whole cohort','Males','Females'});
% % 
% % ylabel('Mass (g)','FontSize',24);
% % 
% % title('Left ventricular mass','FontSize',24)
% % set(gca,'FontSize',20);
% % set(bp,'LineWidth', 2);
% % print(gcf, '-dpdf', 'boxplots_mass_lv.pdf');
% 
% mean_klein2016 = 158;
% std_klein2016 =56.8;
% 
% patch([0.75 1.25 1.25 0.75], [mean_klein2016-std_klein2016 mean_klein2016-std_klein2016 ...
%     mean_klein2016+std_klein2016 mean_klein2016+std_klein2016],[0 1 0],...
% 'FaceAlpha',0.3,'LineStyle','none')
% 
% bp=boxplot([LV'*1.055],...
%     [ones(1,20)],'Labels',...
%     {'Whole cohort'});
% 
% ylabel('Mass (g)','FontSize',24);
% 
% title('Left ventricular mass','FontSize',24)
% set(gca,'FontSize',20);
% set(bp,'LineWidth', 2);
% print(gcf, '-dpdf', 'boxplots_mass_lv.pdf');
% hold off
% 
% close all