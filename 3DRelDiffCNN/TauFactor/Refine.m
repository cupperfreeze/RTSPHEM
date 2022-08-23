% Compute convergence behavior of TauFactor with respect to mesh refinement
% For this script to run, please add path to TauFactor app

load('ForTauFactor.mat') %load example images

for i = 0:5
    for j = 0:2
        refinementLevel = j;
        A = im{i+1};
        B = zeros(size(A).*2^refinementLevel, 'uint8');
        for a = 1:2^refinementLevel
            for b = 1:2^refinementLevel
                for c = 1:2^refinementLevel
                    B(a:2^refinementLevel:end, b:2^refinementLevel:end, c:2^refinementLevel:end) = A;
                end
            end
        end
        Results = TauFactor('InLine', 1, 0, 0, B, [0, 0, 0; 0, 0, 0; 1, 0, 0], [100, 100, 100])
        save(3*i+j+1) = Results.Tau_W1.VolFrac / Results.Tau_W1.Tau;
    end
end

%tau: [0.129494528961026,0.131679729480527,0.132568405503143,0.196946486920997,0.199667422033846,0.200828221373451,0.0875835622973911,0.0895040496707393,0.0904736684579651,0.0579276781024580,0.0595021716776983,0.0601428098704709,0.0696462563237065,0.0725997056630481,0.0738356923529083,0.0244915837521070,0.0263903156567991,0.0272009686732979]
%sim: 0.1326  0.1330  0.1331
