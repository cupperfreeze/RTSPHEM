%> @file SurfaceGrid/visualize.m Visualization of instances of the class SurfaceGrid.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function visualize(sg)

assert(isa(sg, 'SurfaceGrid'))

sg.grid.visualize

line(sg.coordV([1:end, 1], 1), sg.coordV([1:end, 1], 2), 'LineWidth', 2)

end % visualize
