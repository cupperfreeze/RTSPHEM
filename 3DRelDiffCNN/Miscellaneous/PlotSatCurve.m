% make plots

close all 

figure('units', 'pixels', 'position', [0, 0, 2800, 1100], 'visible', 'on'); 
load('DiffSample2.mat') 
load('SampleSat2Data.mat') 
a = reshape(DiffData, 5000, 6); 
b = reshape(porosities, 5000, 6); 
for i = 1:18 
subplot(3, 6, i) 
plot(1:-0.1:0.5, a(50 + i, 1:6), '-*', 'MarkerSize', 15) 
hold on 
sw = 1:-0.1:0.5; 
Drel = a(50+i, 1:6); 
g = fittype('x.^a'); 
f0 = fit(sw', Drel(:)./Drel(1), g); 
fitalpha = f0.a; 
D_Brugge = sw.^fitalpha .* Drel(1); 

sf = nlinfit([repmat(b(50+i), 6, 1), kron(sw', ones(1, 1))], Drel', @(s, x) x(:, 1).^s(1).*x(:, 2).^(-s(2)), [0, 0]) 

 %D_MQ = (b(50+i)).^(sf(1))./(sw.^sf(2));
plot(1:-0.1:0.5, D_Brugge, '-*', 'MarkerSize', 15) 
hold on 
 %plot(1:-0.1:0.5,D_MQ,'-*','MarkerSize',15)
xlabel('$s_w$', 'Interpreter', 'Latex'); 
ylabel('$D_{rel}$', 'Interpreter', 'Latex'); 
set(gca, 'FontSize', 15) 
ylim = get(gca, 'Ylim'); 
text(0.5, 0.8*ylim(2)+0.2*ylim(1), ['$\alpha$=', num2str(round(fitalpha, 3))], 'FontSize', 15, 'Interpreter', 'Latex') 
 %text(0.5,0.8*ylim(2)+0.1*ylim(1),['$s$=',num2str(round(fitS,3))],'FontSize',15,'Interpreter','Latex')
 %text(0.5,0.8*ylim(2)+0*ylim(1),['$t$=',num2str(round(log(fitT)/log(b(50+i)),3))],'FontSize',15,'Interpreter','Latex')

end 
 % % % 
 % % AbsDiff=(abs(a(:,1)-a(:,6)))./a(:,1);
 % % mean(AbsDiff)
 % % std(AbsDiff)
 % % max(AbsDiff)
 % % min(AbsDiff)
 % % figure
 % % histogram(AbsDiff)
 % 
figure 
subplot(1, 2, 2) 
histogram(a(:, 1), 25, 'facealpha', 0.5) 

hold on 
histogram(a(:, 6), 25, 'facealpha', 0.5) 
legend('$s_w=1$', '$s_w=0.5$', 'Interpreter', 'Latex') 
xlabel('$D_{rel}$', 'Interpreter', 'Latex'); 
ylabel('frequency'); 
set(gca, 'FontSize', 15) 


subplot(1, 2, 1) 
plot(porosities(1:5000), a(:, 1), '*') 
hold on 
plot(porosities(25001:30000), a(:, 6), '*') 
 % 
 %  hold on
 %  D_MQ1 = (b(:,1)*1).^(10/3)./(b(:,1).^2);
 %  [temp1,sort1]=sort(b(:,1)); 
 %  plot(temp1,D_MQ1(sort1),'Color','black','LineWidth',3)
 %  
 %  hold on
 %  D_MQ1 = (b(:,6)).^(10/3)./(b(:,1).^2);
 %  [temp1,sort1]=sort(b(:,6)); 
 %  plot(temp1,D_MQ1(sort1),'Color','green','LineWidth',3)


legend('$s_w=1$', '$s_w=0.5$', 'Interpreter', 'Latex') 
xlabel('$\phi_w$', 'Interpreter', 'Latex'); 
ylabel('$D_{rel}$', 'Interpreter', 'Latex'); 
set(gca, 'FontSize', 15) 

load('CNNDiffPlain05.mat') 
load('DiffSample2.mat')
figure 
subplot(2, 2, 1) 
plot(trainingData, predictedTrain, '*') 
hold on 
plot([0, 0.3], [0, 0.3], 'LineWidth', 2) 
xlabel('$D_{cmp}$', 'Interpreter', 'Latex'); 
ylabel('$D_{prd}$', 'Interpreter', 'Latex'); 
set(gca, 'FontSize', 15) 
title('CNN 1 training') 

subplot(2, 2, 2) 
plot(validationData, predicted, '*') 
hold on 
plot([0, 0.3], [0, 0.3], 'LineWidth', 2) 
xlabel('$D_{cmp}$', 'Interpreter', 'Latex'); 
ylabel('$D_{prd}$', 'Interpreter', 'Latex'); 
set(gca, 'FontSize', 15) 
title('CNN 1 validation') 

 
DiffData1 = DiffData(1:4500); 
DiffData2 = DiffData(25000+(1:4500)); 
DiffData1 = DiffData1(DiffData2 < 0.3 & DiffData2 > 0.01); 

subplot(2, 2, 3) 
plot(trainingData, predictedTrain, '*') 
hold on 
plot([0, 0.17], [0, 0.17], 'LineWidth', 2) 
 % hold on
 % plot(trainingData,DiffData.* 0.5.^1.285,'.');
xlabel('$D_{cmp}$', 'Interpreter', 'Latex'); 
ylabel('$D_{prd}$', 'Interpreter', 'Latex'); 
set(gca, 'FontSize', 15) 
title('CNN 2 training') 

subplot(2, 2, 4) 
DiffData1 = DiffData(4501:5000); 
DiffData2 = DiffData(25000+(4501:5000)); 
DiffData1 = DiffData1(DiffData2 < 0.3 & DiffData2 > 0.01); 
plot(validationData, predicted, '*') 
hold on 
plot([0, 0.17], [0, 0.17], 'LineWidth', 2) 
 % hold on
 % plot(validationData,DiffData.* 0.5.^1.285,'.');
xlabel('$D_{cmp}$', 'Interpreter', 'Latex'); 
ylabel('$D_{prd}$', 'Interpreter', 'Latex'); 
set(gca, 'FontSize', 15) 
title('CNN 2 validation') 
 % 

sf = nlinfit([repmat(b(1:4500, 1), 1, 1), kron(sw(6), ones(4500, 1))], a(1:4500, 6), @(s, x) x(:, 1).^s(1).*x(:, 2).^(-s(2)), [0, 0]) 
predicted = repmat(b(4501:5000, 1), 1, 1).^sf(1) .* kron(sw(6), ones(500, 1)).^(-sf(2)); 
validationData = a(4501:5000, 6); 
R_squred = 1 - sum((predicted-(validationData)).^2) / (sum(((validationData)-mean((validationData))).^2)) 
