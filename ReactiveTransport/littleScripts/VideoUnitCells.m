% Make video clips of the geometry evolution within the unit cells over
% simulation time of SPE

try
    load('SPE10.mat', 'saveXi', 'saveDist', 'timeIterationStep', 'microscaleGrid')
catch
    error('Data could not be loaded');
end


for j = [15] % Choose which unit cells shall be processed
    figure
    axis equal
    % axis tight manual
    % set(gca,'nextplot','replacechildren');


    cellNum = j;

    for i = 1:timeIterationStep
        [~, ~] = evaluateInterface(microscaleGrid, saveXi{j, i}, saveDist{j, i}, false)
        F(i) = getframe(gcf);
        clf;
        %     temp = F(i).cdata;
        %     F(i).cdata =temp(:,:,1);

    end


    % create the video writer with 1 fps
    writerObj = VideoWriter(char(strcat('CellNew', string(cellNum), '.avi')));
    writerObj.FrameRate = 20;
    % set the seconds per image
    % open the video writer
    open(writerObj);
    % write the frames to the video
    for i = 1:length(F)
        % convert the image to a frame
        frame = F(i);
        writeVideo(writerObj, frame);
    end
    % close the writer object
    close(writerObj);

    clear F
end
temp = 0;
% for i=1:numel(saveCellAdap)
%     temp = temp + numel(saveCellAdap{i});
% end
% saveAdap=temp/(numberOfSlices*timeIterationStep);
