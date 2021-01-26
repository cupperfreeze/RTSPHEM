function [isBalli, radiusEstimate] = isBall(LevelSet, tol)
%Check wheater a geometry corresponding to a LevelSet resembles a circle (Sensitivity tol)
%If yes, also return estimated radius

h = 1 / (sqrt(size(LevelSet, 1)) - 1);
H = h / 3;
[Xq, Yq] = meshgrid(-0.5:H:0.5, -0.5:H:0.5);
[x, y] = meshgrid(-0.5:h:0.5, -0.5:h:0.5);

%surf(x, y, reshape(LevelSet,round(1/h)+1,round(1/h)+1))
Vq = interp2(x, y, reshape(LevelSet, round(1 / h) + 1, round(1 / h) + 1), Xq, Yq);
inSolid = logical(Vq > 0);
Xq = Xq(inSolid);
Yq = Yq(inSolid);
Points = [Xq(:), Yq(:)]; %Coordinates of points within solid
baryCenter = [mean(Points(:, 1)), mean(Points(:, 2))];
Points(:, 1) = Points(:, 1) - baryCenter(1, 1); %Normalize
Points(:, 2) = Points(:, 2) - baryCenter(1, 2);

%scatter(Points(:,1),Points(:,2))

MainAxis = Points' * Points; %Correlation matrix


if norm((MainAxis-eye(2)*MainAxis(1, 1))) / norm(MainAxis) < tol
    radiusEstimate = (MainAxis(1, 1) * H^2 * 4 / pi).^(1 / 4);
    inBall = sum((Points(:, 1).^2+Points(:, 2).^2) <= radiusEstimate.^2);
    isBalli = logical((1-inBall/size(Points, 1)) < tol);
else
    isBalli = false;
    radiusEstimate = nan;
end
end