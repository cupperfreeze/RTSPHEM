%> @file AbstactGrid.print.m Print grid data information to stdout.

function print(obj)

printline(2, 'Some grid properties')
printline(3, 'Number of triangles:                        %d', obj.numT)
printline(3, 'Number of edges:                            %d', obj.numE)
printline(3, '   ...thereof inner edges:               %d', sum(obj.idE == 0))
printline(3, '              boundary edges:            %d', sum(obj.idE ~= 0))
printline(3, 'Number of vertices:                         %d', obj.numV)

end % print
