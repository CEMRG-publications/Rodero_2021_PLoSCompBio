close all

figure('DefaultAxesFontSize',20)
grid on

subplot(2,2,1)

plot(1:18,corEDVWT(:,1),'x-k','LineWidth',3,'MarkerSize',15)
hold on

% Convert color code to 1-by-3 RGB array (0~1 each)
str = '#67001F';
color = sscanf(str(2:end),'%2x%2x%2x',[1 3])/255;

plot(1:18,corEDVWT(:,2),'+-','LineWidth',3,'MarkerSize',15,...
'Color',color)

str = '#E98B6F';
color = sscanf(str(2:end),'%2x%2x%2x',[1 3])/255;

plot(1:18,corEDVWT(:,3),'*-','LineWidth',3,'MarkerSize',15,...
'Color',color)

legend('EDV','mass','EDV/mass')

xlim([1,18])
xticks(1:18)
xticklabels(1:18)
ylim([-1,1])
xlabel('Mode')
ylabel('Correlation')
title('Correlation plot for the LV')

hold off


subplot(2,2,2)


plot(1:18,corEDVWT(:,4),'x-k','LineWidth',3,'MarkerSize',15)
hold on

% Convert color code to 1-by-3 RGB array (0~1 each)
str = '#67001F';
color = sscanf(str(2:end),'%2x%2x%2x',[1 3])/255;

plot(1:18,corEDVWT(:,5),'+-','LineWidth',3,'MarkerSize',15,...
'Color',color)

str = '#E98B6F';
color = sscanf(str(2:end),'%2x%2x%2x',[1 3])/255;

plot(1:18,corEDVWT(:,6),'*-','LineWidth',3,'MarkerSize',15,...
'Color',color)

legend('EDV','mass','EDV/mass')

xlim([1,18])
xticks(1:18)
xticklabels(1:18)
ylim([-1,1])
xlabel('Mode')
ylabel('Correlation')
title('Correlation plot for the RV')

hold off




subplot(2,2,3)


plot(1:18,corEDVWT(:,7),'x-k','LineWidth',3,'MarkerSize',15)
hold on

% Convert color code to 1-by-3 RGB array (0~1 each)
str = '#67001F';
color = sscanf(str(2:end),'%2x%2x%2x',[1 3])/255;

plot(1:18,corEDVWT(:,8),'+-','LineWidth',3,'MarkerSize',15,...
'Color',color)

str = '#E98B6F';
color = sscanf(str(2:end),'%2x%2x%2x',[1 3])/255;

plot(1:18,corEDVWT(:,9),'*-','LineWidth',3,'MarkerSize',15,...
'Color',color)

legend('EDV','mass','EDV/mass')

xlim([1,18])
xticks(1:18)
xticklabels(1:18)
ylim([-1,1])
xlabel('Mode')
ylabel('Correlation')
title('Correlation plot for the LA')

hold off



subplot(2,2,4)


plot(1:18,corEDVWT(:,10),'x-k','LineWidth',3,'MarkerSize',15)
hold on

% Convert color code to 1-by-3 RGB array (0~1 each)
str = '#67001F';
color = sscanf(str(2:end),'%2x%2x%2x',[1 3])/255;

plot(1:18,corEDVWT(:,11),'+-','LineWidth',3,'MarkerSize',15,...
'Color',color)

str = '#E98B6F';
color = sscanf(str(2:end),'%2x%2x%2x',[1 3])/255;

plot(1:18,corEDVWT(:,12),'*-','LineWidth',3,'MarkerSize',15,...
'Color',color)

legend('EDV','mass','EDV/mass')

xlim([1,18])
xticks(1:18)
xticklabels(1:18)
ylim([-1,1])
xlabel('Mode')
ylabel('Correlation')
title('Correlation plot for the RA')

hold off