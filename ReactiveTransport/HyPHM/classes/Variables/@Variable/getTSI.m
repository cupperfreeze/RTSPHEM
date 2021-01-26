%> @file Variable/getTSI.m Get data as instance of TriScatteredInterp.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%> @brief Get data as instance of TriScatteredInterp.
%>
%> @b Dependencies
%>   Variable.getdata
%>
%> @b Output
%>   varargout   [ TriScatteredInterp ] or [ TriScatteredInterp, TriScatteredInterp ]
%>   NOTE: RT0 is returned in Cartesian Coordinates!
%>
%> @b Howto
%>
%> @code
%>   tri = scalunk.getTSI(1);
%>   [QX, QY] = meshgrid(0:.1:10, 0:.1:10);
%>   mesh(QX, QY, tri(QX,QY))
%>
%>   [tri1, tri2] = vectunk.getTSI(1);
%>   [QX, QY] = meshgrid(0:.1:10, 0:.1:10);
%>   quiver(QX, QY, tri1(QX,QY), tri2(QX,QY))
%> @endcode


function varargout = getTSI(this, stepno)

msg = 'HyPHM: Wrong number of output arguments (1 for scalars, 2 for vectors).';

switch this.type
    case 'P0'
        assert(nargout == 0 || nargout == 1, msg)
        varargout{1} = TriScatteredInterp(this.grid.baryT(:, 1), ...
            this.grid.baryT(:, 2), ...
            this.getdata(stepno), 'nearest'); %changed by SG
    case 'P1'
        assert(nargout == 0 || nargout == 1, msg)
        varargout{1} = TriScatteredInterp(this.grid.coordV(:, 1), ...
            this.grid.coordV(:, 2), ...
            this.getdata(stepno));

    case 'RT0'
        assert(nargout == 2, msg)
        g = this.grid;
        dataRT0 = this.getdata(stepno);
        dataP2P2 = zeros(g.numT, 2);
        for kT = 1:g.numT
            dataP2P2(kT, 1:2) = RT0.getCartesianCoords(g, dataRT0(g.E0T(kT, :))', kT, g.baryT(kT, :)');
        end
        varargout{1} = TriScatteredInterp(g.baryT(:, 1), g.baryT(:, 2), dataP2P2(:, 1));
        varargout{2} = TriScatteredInterp(g.baryT(:, 1), g.baryT(:, 2), dataP2P2(:, 2));

    case 'P2P2'
        assert(nargout == 2, msg)
        data = this.getdata(stepno);
        u = data(:, 1);
        v = data(:, 2);
        varargout{1} = TriScatteredInterp([this.grid.coordV(:, 1); this.grid.baryE(:, 1)], ...
            [this.grid.coordV(:, 2); this.grid.baryE(:, 2)], ...
            u);
        varargout{2} = TriScatteredInterp([this.grid.coordV(:, 1); this.grid.baryE(:, 1)], ...
            [this.grid.coordV(:, 2); this.grid.baryE(:, 2)], ...
            v);

    otherwise
        error('HyPHM: Undefined type.')
end

end