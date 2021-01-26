% Script to evaluate CNN prediction quality on SPE10 geometries

calculated = zeros(440, 15); % Stokes calculated permeabilities
predicted = zeros(2, 440, 15); % CNN predicted permeabilities

% import data and networks
load('SPE10.mat', 'permeabilityTensors', 'saveXi');
load('CNNPermGeneral.mat');
net1 = net;
load('CNNPermRectangle.mat');
net2 = net;
[x, y] = meshgrid(1:64);

% for 15 unit cells, extract geometry and calculated data in appropriate
% format
for i = 1:15
    container = zeros(64, 64, 1, 440);

    for timeIdx = 1:440

        currentImage = reshape(double(saveXi{i, timeIdx}(:, 1)), 101, 101);
        currentImage = interp2(currentImage, x*101/64, y*101/64, 'nearest');
        container(:, :, 1, timeIdx) = int8(currentImage > 1.5);
        calculated(timeIdx, i) = permeabilityTensors{30*i}(4, timeIdx) * 10^8;

    end
    % perform network predictions
    predicted(1, :, i) = exp(1/0.1*(net1.predict(turn(container, pi / 2)) - 0.9));
    predicted(2, :, i) = exp(1/0.2*(net2.predict(turn(container, pi / 2)) - 0.9));
end
%
% for i=1:15
%     subplot(3,5,i);
%     plot(1:440,log(predicted(1,:,i)));
%     hold on
%      plot(1:440,log(predicted(2,:,i)));
%     hold on
%      plot(1:440,log(predicted(3,:,i)));
%     hold on
%     plot(1:440,log(calculated(:,i)));
%     legend('predictedG','predictedSR','calculated','Location','northwest')
%     xlabel('timeStep')
%     ylabel('log(K_{1,1})')
%     set(gca, 'fontsize', 10)
% end

% plot results
figure
idx = [6, 13, 15];
for i = 1:3
    subplot(1, 3, i);
    plot(1:440, log(predicted(1, :, idx(i))), '.');
    hold on
    plot(1:440, log(predicted(2, :, idx(i))), '.');
    hold on
    plot(1:440, log(calculated(:, idx(i))), '.');
    legend('CNN', 'CNN special', 'calculated', 'Location', 'northwest')
    xlabel('timeStep')
    ylabel('log(K_{2,2})')
    set(gca, 'fontsize', 14)
end

% subplot(1,3,1);
% imshow(128-128*TrainingIm(:,:,1,141)); colormap gray;
% subplot(1,3,2);
% imshow(128-128*TrainingIm(:,:,1,297)); colormap gray;
% subplot(1,3,3);
% imshow(128-128*TrainingIm(:,:,1,35)); colormap gray;