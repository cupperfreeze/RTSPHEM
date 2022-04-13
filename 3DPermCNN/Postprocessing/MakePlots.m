% % Plot sample data correlations and more
 
% Scalings arising from voxel size 2.25 \mu m 
scaleArea = (2.25*10^(-3)).^2;
scalePerm = ((225*10^(-6)).^2)*10^(12);
scaleVolume = (225*10^(-3))^3;

p=genpath('CNN');
addpath(p);
p=genpath('ComputedData');
addpath(p);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot correlation between permeability and pore space characteristics

load('SampleDataBENTHEIMER')
load('PermDataBENTHEIMEROrd1')
figure
subplot(1,4,1)
semilogy(areas(1:10000)*scaleArea, (PermData*scalePerm),'*')
xlabel('$A_\mathrm{cmp}$ [mm$^2$]','Interpreter','LaTeX');
ylabel('$k_\mathrm{cmp}$ [D]','Interpreter','LaTeX');
set(gca,'Fontsize',15)

subplot(1,4,2)
semilogy(porosities(1:10000),(PermData*scalePerm),'*')
xlabel('$\phi_\mathrm{cmp}$','Interpreter','LaTeX');
ylabel('$k_\mathrm{cmp}$ [D]','Interpreter','LaTeX');
set(gca,'Fontsize',15)

subplot(1,4,3)
semilogy(areas(1:10000)./(1-porosities(1:10000))*scaleArea/scaleVolume,(PermData*scalePerm),'*')
xlabel('$A_\mathrm{spec}$ [mm$^{-1}$]','Interpreter','LaTeX');
ylabel('$k_\mathrm{cmp}$ [D]','Interpreter','LaTeX');
set(gca,'Fontsize',15)

subplot(1,4,4)
semilogy(log10(MinCut(1:10000)),(PermData*scalePerm),'*')
xlabel('$\log_{10}(f_\mathrm{max})$','Interpreter','LaTeX');
ylabel('$k_\mathrm{cmp}$ [D]','Interpreter','LaTeX');
set(gca,'Fontsize',15)

Trafo=@(x) x.^1.407*10^(-8.183);
predicted = Trafo(MinCut(1:10^4));
R_squredLog = 1 - sum((log10(predicted)-log10(PermData)).^2) / (sum((log10(PermData)-mean(log10(PermData))).^2))
R_squred = 1 - sum(((predicted)-(PermData)).^2) / (sum(((PermData)-mean((PermData))).^2))


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot Traning success of PhyCNN

% load('CNNPhys2.mat')
% figure
% subplot(2,2,1)
% loglog((validationData*scalePerm),(predicted*scalePerm),'*')
% hold on
% loglog((validationData*scalePerm),(validationData*scalePerm),'LineWidth',3)
% xlabel('$k_\mathrm{cmp}$ [D]','Interpreter','LaTeX');
% ylabel('$k_\mathrm{prd}$ [D]','Interpreter','LaTeX');
% set(gca,'Fontsize',15);
% title('a)')
% 
% subplot(2,2,2)
% loglog((trainingData*scalePerm),(predictedTrain*scalePerm),'*')
% hold on
% loglog((trainingData*scalePerm),(trainingData*scalePerm),'LineWidth',3)
% xlabel('$k_\mathrm{cmp}$ [D]','Interpreter','LaTeX');
% ylabel('$k_\mathrm{prd}$ [D]','Interpreter','LaTeX');
% set(gca,'Fontsize',15);
% title('b)')
%  
% load('CompareBEREA')
% subplot(2,2,3)
% loglog((validationData*scalePerm),(predicted*scalePerm),'*')
% hold on
% loglog((validationData*scalePerm),(validationData*scalePerm),'LineWidth',3)
% xlabel('$k_\mathrm{cmp}$ [D]','Interpreter','LaTeX');
% ylabel('$k_\mathrm{prd}$ [D]','Interpreter','LaTeX');
% set(gca,'Fontsize',15);
% title('c)')
% 
% load('CompareCASTLE')
% subplot(2,2,4)
% loglog((validationData*scalePerm),(predicted*scalePerm),'*')
% hold on
% loglog((validationData*scalePerm),(validationData*scalePerm),'LineWidth',3)
% xlabel('$k_\mathrm{cmp}$ [D]','Interpreter','LaTeX');
% ylabel('$k_\mathrm{prd}$ [D]','Interpreter','LaTeX');
% set(gca,'Fontsize',15);
% title('d)')

% disp('Correlation Validation')
% R_squred = 1 - sum((predicted-(validationData)).^2) / (sum(((validationData)-mean((validationData))).^2))
% R_squredLog = 1 - sum((log10(predicted)-log10(validationData)).^2) / (sum((log10(validationData)-mean(log10(validationData))).^2))
% 
% disp('Correlation Training')
% R_squred = 1 - sum((predictedTrain-(trainingData)).^2) / (sum(((trainingData)-mean((trainingData))).^2))
% R_squredLog = 1 - sum((log10(predictedTrain)-log10(trainingData)).^2) / (sum((log10(trainingData)-mean(log10(trainingData))).^2))


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot distribution in porosity and permeability in Bentheimer data set
% 
% load('PermDataBENTHEIMEROrd1.mat');
% load('SampleDataBENTHEIMER.mat');
% subplot(1,2,1)
% histogram(porosities(1:10000),100)
% xlabel('$\phi_\mathrm{cmp}$','Interpreter','LaTeX');
% ylabel('frequency');
% set(gca,'Fontsize',15)
% Mean = mean(porosities(1:10000))
% STD = std(porosities(1:10000))
% 
% 
% subplot(1,2,2)
%  [~,edges] = histcounts(log10(PermData*scalePerm));
%  histogram(PermData*scalePerm,10.^edges);
%  set(gca, 'xscale','log');
% xlabel('$k_\mathrm{cmp}$ [D]','Interpreter','LaTeX');
% ylabel('frequency');
% set(gca,'Fontsize',15)
% Mean=mean(PermData*scalePerm)
% STD=std(PermData*scalePerm)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot computed and predicted permeability for artificially distorted
% samples
% 
% load('CompareVaryPoro.mat')
% MarkerSize=7;
% fig=figure
% set(fig,'defaultAxesColorOrder',[0,0,0;0,0.5,0.5])
% yyaxis left
% semilogy([flipud(porosities(1:17));porosities(18:end) ],[flipud(predicted(1:17));predicted(18:end) ].*scalePerm,'d','MarkerSize',MarkerSize,'Color',[195,0,0]/255)
% hold on
% semilogy(0.299856, 1.509882567916066e-04.*scalePerm,'d','MarkerSize',2*MarkerSize,'Color',[195,0,0]/255)
% hold on
% semilogy([flipud(porosities(1:17));porosities(18:end) ],[flipud(validationData(1:17));validationData(18:end)].*scalePerm,'d','MarkerSize',MarkerSize,'Color','b')
% hold on
% semilogy([0.299856 ],0.000147435.*scalePerm,'d','MarkerSize',2*MarkerSize,'Color','b')
% hold on
% rectangle('Position',[0,0.05,0.8,50],'FaceColor',[0.5,0.5,0.5,0.2]','EdgeColor',[0.5,0.5,0.5,0.2]')
% ylabel('$k$ [D]','Interpreter','LaTeX')
% ylim([0.001, 250])
% 
% 
% yyaxis right
% plot([flipud(porosities(1:17));porosities(18:end) ],abs([flipud(predicted(1:17));predicted(18:end) ]-[flipud(validationData(1:17));validationData(18:end)])./[flipud(validationData(1:17));validationData(18:end)],'s','MarkerSize',MarkerSize)
% hold on
% plot(0.299856,abs(1.509882567916066e-04 - 0.000147435)/0.000147435,'s','MarkerSize',2*MarkerSize)
% ylim([-0.1, 3])
% xlim([0, 0.75])
% xlabel('$\phi_\mathrm{cmp}$','Interpreter','LaTeX')
% ylabel('rel. deviation factor')
% legend('$k_\mathrm{prd}$','','$k_\mathrm{cmp}$','', '$|k_\mathrm{prd}-k_\mathrm{cmp}|/ k_\mathrm{cmp}$','Interpreter','LaTeX','Location','east')
% set(gca,'Fontsize',15);
% % 
% PoroLowerEnd = interp1([flipud(validationData(1:17));validationData(18:end)].*scalePerm,[flipud(porosities(1:17));porosities(18:end) ],0.05)
% PoroUpperEnd = interp1([flipud(validationData(1:17));validationData(18:end)].*scalePerm,[flipud(porosities(1:17));porosities(18:end) ],50)

