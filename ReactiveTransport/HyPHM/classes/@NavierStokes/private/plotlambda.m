j = 5;

[X, Y] = meshgrid(1:.2:4, 0:.2:4);

for k = 1:size(X, 1)
    for ell = 1:size(X, 2)
        fct = phi(verts, [X(k, ell); Y(k, ell)]);
        Z(k, ell) = fct(j);
    end
end

hold on

surf(X, Y, Z)
line(verts(1:2, 1), verts(1:2, 2), [0, 0], 'Color', 'k')
line(verts(2:3, 1), verts(2:3, 2), [0, 0], 'Color', 'k')
line(verts([1, 3], 1), verts([1, 3], 2), [0, 0], 'Color', 'k')

line(verts(1:2, 1), verts(1:2, 2), [1, 1], 'Color', 'k')
line(verts(2:3, 1), verts(2:3, 2), [1, 1], 'Color', 'k')
line(verts([1, 3], 1), verts([1, 3], 2), [1, 1], 'Color', 'k')

for k = [0, 1]
    text(verts(1, 1), verts(1, 2), k, '1')
    text(verts(2, 1), verts(2, 2), k, '2')
    text(verts(3, 1), verts(3, 2), k, '3')
end
