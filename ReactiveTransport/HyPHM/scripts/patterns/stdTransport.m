%> @file stdTransport.m
%> @brief Script for one species Transport, see @ref stdTransportTutorial.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%> @page stdTransportTutorial Tutorial: One-species transport
%>
%> @section sec1 Discretization in Space
%>
%> Let the considered domain @f$\Omega@f$ be @f$[0,2]\times[0,1]@f$.  For
%> rectangular bounded domains, there's a pre-defined Grid generation
%> function.  With the mesh size @f$h = 0.08@f$, an instance of Grid is
%> created by
%> @code
%>   g = domainRectangle(0, 2,  0, 1,  0.08);
%> @endcode
%> Display the main properties by @c g.print.  The Grid may be visualized
%> by
%> @code
%>   g.visualize()
%>   g.visualize('idE')          % with edge ids
%>   g.visualize('numT', 'numE') % with triangle and edge numbers
%> @endcode
%> You find more visualization options in Grid.visualize.
%>
%> @section sec2 Discretization in Time
%>
%> Say we want to simulate in the time intervall @f$J = [0,0.5]@f$.  A
%> decomposition of @f$J@f$
%> @f[ t_0 = 0,\,0.1,\,0.2,\,\ldots,\,1.5=t_{\#\mathrm{steps}}@f]
%> is obtained via
%> @code
%>   st = Stepper(0 : .1 : 1.5);
%> @endcode
%> Similar to Grid, there's also a Stepper.print and Stepper.visualize
%> method.
%>
%> @section sec3 Definition and Initialization of Unknowns
%> The Unknowns live in discrete function spaces, see Variable.type.  The
%> number of unknowns as well as their types are defined by the respective
%> problem.  Transport requires a flux and a scalar unknown in @f$\vec{RT}_0(T)@f$
%> and @f$P_0(T)@f$, respectively:
%> @code
%>   conc  = Variable(g, st, 'concentration', 'P0');
%>   mflux = Variable(g, st, 'mass flux',     'RT0');
%> @endcode
%> Initialize it by
%> @code
%>   conc.setdata(0,  @(t, x) norm(x-[.5; .5]) < .2);
%>   mflux.setdata(0, zeros(g.numE, 1));
%> @endcode
%> where the first line defines a circle of concentration by one with center
%> @f$[.5, .5]@f$ and radius @f$.2@f$.
%>
%> more to come...
%>
%> @section secend The Complete Script
%> All together, the script may read
%> @code
%> g = domainRectangle(0, 2,  0, 1,  0.08);
%>
%> st = Stepper(0:.1:1.5);
%>
%> conc  = Variable(g, st, 'concentration', 'P0');
%> mflux = Variable(g, st, 'mass flux', 'RT0');
%> conc.setdata(0,  @(t, x) norm(x-[.5; .5]) < .2);
%> mflux.setdata(0, zeros(g.numE, 1));
%>
%> d = Transport(g, st);
%> d.id2F = {1,3,4}; d.id2N = {2};
%> d.Q = mflux;
%> d.U = conc;
%> d.A  = Variable(g, st, 'A',  'P0');    d.A.setdata(@(t,x) 1.0);
%> d.C  = Variable(g, st, 'C',  'RT0');   d.C.setdata(@(t,X,Y) {X==X; (1-X)/3}); % equivalent to 'd.C.setdata(@(t,x) [1; (1-x(1))/3]);' but much faster
%> d.D  = Variable(g, st, 'D',  'CCCC');  d.D.setdata(0.01*eye(2));
%> d.gF = Variable(g, st, 'gF', 'P0E');   d.gF.setdata(@(t,x) 0.0);
%>
%> while st.next
%>   d.computeLevel;
%> end
%>
%> conc.visualize
%> mflux.visualize
%> d.C.visualize
%> @endcode


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialization of Grids
% g = domainRectangle(0, 2,  0, 1,  1);
g = domainRectangle(0, 1, 0, 1, 0.05);
g.print
g.visualize()

%% Initialization of Timer
st = Stepper(0:.1:1);
st.print
st.status

%% Solution Vectors
conc = Variable(g, st, 'concentration', 'P0');
mflux = Variable(g, st, 'mass flux', 'RT0');
% conc.setdata(0,  @(t, x) norm(x-[.5; .5]) < .2); % circle with center [.5, .5] and radius .2.
conc.setdata(0, zeros(g.numT, 1)); %@(t, x) x(1)<=0.75 && x(1)>=0.25 && x(2)<=0.75 && x(2)>=0.25); % circle with center [.5, .5] and radius .2.
mflux.setdata(0, zeros(g.numE, 1));

%% Initialization of Problem Data
d = Transport(g, st, 'Transport problem');
d.Q = mflux;
d.U = conc;
d.balanceU = [];
d.isUpwind = 'full';
d.id2F = {1, 3};
d.id2N = {2};
d.id2D = {4};
d.A = Variable(g, st, '', 'P0');
d.A.setdata(@(t, x) 1.0);
d.C = Variable(g, st, 'C', 'RT0');
d.C.setdata(@(t, X, Y) {X, ==, X; 0 * X}); % equivalent to 'd.C.setdata(@(t,x) [1; (1-x(1))/3]);' but much faster
d.D = Variable(g, st, '', 'CCCC');
d.D.setdata(0.004*eye(2));
d.gF = Variable(g, st, '', 'P0E');
d.gF.setdata(@(t, x) 0.0); % no-flow boundary
d.uD.setdata(@(t, x) x(2) <= 0.75 && x(2) >= 0.25)
d.print

%% Assembly and Computation
printline(1, 'Assembly and Computation')

while st.next
    d.computeLevel('silent');
end % end stepper iteration

% %% Visualization
printline(1, 'Visualization')
conc.visualize
% mflux.visualize
% d.C.visualize
