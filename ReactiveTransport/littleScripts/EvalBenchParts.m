%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%BenchPart1
%
% hydrogenC = hydrogenTransport.U.getdata(numel(timeSteps)-1);
%
% figure
% TrianglesVerticalIdx = find(abs(gridHyPHM.baryT(:,1)-0.05)<0.1/256/1);
% TrianglesVerticalY = gridHyPHM.baryT(TrianglesVerticalIdx,2);
% ConcentrationsVertical = hydrogenC(TrianglesVerticalIdx);
% plot(TrianglesVerticalY,-log10(ConcentrationsVertical*1000),'.');
% title('pH value vertical cut')
% axis([0 0.05 1.7 5])
% xlabel('y-Axis [cm]');
% ylabel('pH');
%
% figure
% TrianglesHorizontalIdx = find(abs(gridHyPHM.baryT(:,2)-0.025)<0.1/256/1);
% TrianglesHorizontalX = gridHyPHM.baryT(TrianglesHorizontalIdx,1);
% ConcentrationsHorizontal = hydrogenC(TrianglesHorizontalIdx);
% plot(TrianglesHorizontalX,-log10(ConcentrationsHorizontal*1000),'.');
% title('pH value horizontal cut')
% axis([0 0.1 1.7 5])
% xlabel('x-Axis [cm]');
% ylabel('pH');
%
% figure
% TrianglesDiagIdx = find(abs(gridHyPHM.baryT(:,2)-gridHyPHM.baryT(:,1)+0.025)<0.1/256/1.1);
% TrianglesDiagX = gridHyPHM.baryT(TrianglesDiagIdx,1);
% ConcentrationsDiag = hydrogenC(TrianglesDiagIdx);
% plot(TrianglesDiagX,-log10(ConcentrationsDiag*1000),'.');
% title('pH value 45 degree cut')
% axis([0.025 0.075 1.7 5])
% xlabel('x-Axis [cm]');
% ylabel('pH');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%BenchPart2
%
close all
figure
plot(timeSteps, Rate)
title('Average dissolution rate')
axis([-100, 8000, 0, 0.7 * 10^(-7)])
xlabel('t [s]');
ylabel('dissolution rate [mol/cm/s]');

figure
plot(timeSteps, surfaceArea{1})
title('Area')
axis([-100, 8000, 0, 0.07])
xlabel('t [s]');
ylabel('area in [cm²]');

figure
plot(timeSteps, Volume./Volume(2))
title('Volume')
axis([-100, 8000, 0, 1.1])
xlabel('t [s]');
ylabel('volume V/V0');

figure
[a, b] = meshgrid(0:lengthYAxis/numPartitionsMicroscale:lengthXAxis, 0:lengthYAxis/numPartitionsMicroscale:lengthYAxis);
contour(a, b, reshape(levelSet{1}(:, 24), 2 * numPartitionsMicroscale + 1, numPartitionsMicroscale + 1)', [-0.02, -0.01, 0, 1])
axis equal
