% Time simulations using TauFactor
% For this script to run, please add path to TauFactor app

InidzesOfData = 1:5000;
imSize = [100, 100, 100];
Results = zeros(numel(InidzesOfData), 1);
% read in Data

disp('read in data');
tic
parfor i = InidzesOfData
    name = strcat('out_connected', num2str(i - 1), ['_', num2str(imSize(1)), 'x', num2str(imSize(2)), 'x', num2str(imSize(3))]);
    im = rwd2mat(name);
    out = TauFactor('InLine', 1, 0, 0, im, [0, 0, 0; 0, 0, 0; 1, 0, 0], [1, 1, 1]);
    Results(i) = out.Tau_W1.VolFrac / out.Tau_W1.Tau;
    if mod(i, 200) == 0
        disp([num2str(i / numel(InidzesOfData) * 100), ' Procent done'])
    end
end
disp('read in data done');
Time = toc;

save('EvalTaufactorEarly.mat', 'Results', 'Time', '-v7.3')
