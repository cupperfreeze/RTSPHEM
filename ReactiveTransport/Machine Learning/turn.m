% turn a stack of n images [height x width x 1 x n] by angle using
% a periodic boundary (torus)

function out = turn(images, angles);

height = size(images, 1);
[x, y] = meshgrid(1:height);

num = size(images, 4);
out = images;
% Coordinate mapping
Xcoord = pixToCoord(x, height) .* cos(angles) + pixToCoord(y, height) .* sin(angles);
Ycoord = -pixToCoord(x, height) .* sin(angles) + pixToCoord(y, height) .* cos(angles);
indizes = sub2ind(size(x), CoordToPix(Xcoord, height), CoordToPix(Ycoord, height));

% image-wise rotation
for i = 1:num
    local = images(:, :, 1, i);
    temp = local(indizes);
    temp = reshape(temp, size(x))';
    out(:, :, 1, i) = temp;
end
% vectorizes version, slower!
% indizes = repmat(indizes,num,1)+ kron((0:num-1)',height^2*ones(height^2,1));
% out2 = permute(reshape(images(indizes),height,height,1,num),[2,1,3,4]);
% norm(single(out(:)-out2(:)))

%Convert pixel index to coordinate
    function out = pixToCoord(pix, height)
        out = (pix - 1) / (height - 1) - 0.5;
    end

%Convert coordinate to pixel index
    function out = CoordToPix(Coord, height)
        idx = find(Coord > 0.5);
        Coord(idx) = Coord(idx) - 1;
        idx = find(Coord < -0.5);
        Coord(idx) = Coord(idx) + 1;
        out = 1 + (Coord + 0.5) * (height - 1);
        out = round(out(:));
    end
end
