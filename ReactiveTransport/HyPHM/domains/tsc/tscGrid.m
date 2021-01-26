%> @file tscGrid.m Generation of (semi-)periodic grids for the two-scale
%>convergence problem, ie either @f$\Omega@f$, @f$\Omega_\varepsilon@f$ or
%> @f$Y@f$.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%>
%> @params dom String specifying the domain: @f$\Omega@f$: <tt>'OM'</tt>, @f$\Omega_\varepsilon@f$: <tt>'OMe'</tt> or @f$Y@f$: <tt>'Y'</tt>.
%> @params hmax Grid fineness.
%> @params varargin Scaling parameter @f$\varepsilon@f$ (has to be a divisor of one).
%>

function g = tscGrid(dom, hmax, varargin)

assert(nargin == 2 || nargin == 3, 'HyPHM: Check number of arguments.')
assert(hmax < 1 || hmax > 0)

writehmax(hmax) % write file which is included by .geo-file

switch dom
    case {'OM', 'OMe'}
        assert(nargin == 3, 'HyPHM: Third agument has to be epsilon, required for vertical length.');
        eps = varargin{1};
        if eps == 1 / 6 || eps == 1 / 7, error('HyPHM: eps = 1/6, 1/7 is illshaped grid (bug?!).'), end
        assert(1/eps == round(1 / eps), 'HyPHM: eps has to be a divisor of one.')

        writeeps(eps) % write file which is included by .geo-file
        filename = ['tsc_', dom, '.geo']; % filename is tsc_OM.geo or tsc_OMe.geo
        g = Grid(filename, 1);
    case 'Y'
        writeeps(1) % write file which is included by .geo-file with eps == 1 for unity cell Y
        g = Grid('tsc_Y.geo', 1);
    otherwise
        error('HyPHM: Wrong arguments, read the documentation please.')
end

end

function writeeps(eps) % write ./domains/tsc/tsc_epsilon.geo

file = fopen('./domains/tsc/tsc_epsilon.geo', 'wt');
fprintf(file, 'eps = %.5f;\n', eps); % see Remark 1 below
fclose(file);

end

function writehmax(hmax) % write ./domains/tsc/tsc_meshwidth.geo

file = fopen('./domains/tsc/tsc_meshwidth.geo', 'wt');
fprintf(file, 'h = %f;\n', hmax);
fclose(file);

end

%% Remark 1
% There's a bug in gmsh causing a crash for eps=0.333333... which costed me
% roughly about 2 days of my life.  In any case, eps = 0.33333 (<- 5 post
% digits) works.
