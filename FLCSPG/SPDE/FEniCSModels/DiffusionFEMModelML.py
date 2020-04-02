from dolfin import *
from .FEMModel import *

__author__ = ["Benjamin, Bykowski", "Jean-Luc Bouchot"]
__copyright__ = "Copyright 2015-2020, Chair C for Mathematics (Analysis), RWTH Aachen; Seminar for Applied Mathematics, ETH Zurich; School of Mathematics and Statistics, Beijing Institute of Technology"
__credits__ = ["Jean-Luc Bouchot", "Benjamin, Bykowski", "Holger Rauhut", "Christoph Schwab"]
__license__ = "GPLv3"
__version__ = "0.1.0-dev"
__maintainer__ = "Jean-Luc Bouchot"
__email__ = "jlbmathit@gmail.com"
__status__ = "Development"
__lastmodified__ = "2020/04/02"

class DiffusionFEMModelML(FEMModel):
    def __init__(self, a, f, M_gen, mesh_size, onb_fname = None):
        self.a         = a
        self.f         = f
        self.M_gen     = M_gen

        self.mesh_size = mesh_size
        self.init_simple_mesh()
        if onb_fname is not None: 
            self.hilbert_basis = np.load(onb_fname) # This is your responsibility to make sure the file is correct.

    def solve(self, z):
        # Make FEniCS output only the most important messages
        CRITICAL  = 50 #, // errors that may lead to data corruption and suchlike
        ERROR     = 40 #, // things that go boom
        WARNING   = 30 #, // things that may go boom later
        INFO      = 20 #, // information of general interest
        PROGRESS  = 16 #, // what's happening (broadly)
        TRACE     = 13 #, // what's happening (in detail)
        DBG       = 10#   // sundry
        set_log_level(ERROR)

        # Create mesh if there is none
        if not hasattr(self, 'mesh'):
            self.init_simple_mesh()

        # Make sure we have a basis of functions for our Hilbert space expansion
        if not hasattr(self, 'hilbert_basis') or self.hilbert_basis is None: 
            self.generate_Hilbert_basis()

        # Create approximation space
        V = FunctionSpace(self.mesh, 'Lagrange', 1)

        # Define boundary conditions
        bc = DirichletBC(V, Constant(0.0), lambda x, on_boundary: on_boundary)

        # Define variational problem
        w = TrialFunction(V)
        v = TestFunction(V)

        params = self.split_params([self.a, self.f], z)

        x = SpatialCoordinate(self.mesh)
        A = self.a(x, Constant(params[0])) * inner(nabla_grad(w), nabla_grad(v)) * dx
        L = self.f(x, Constant(params[1])) * v * dx

        # print("Current parameters {} ".format(params))

        # Create goal-functional for error estimation
        u      = Function(V)
        self.M = self.M_gen(self, u, dx)

        # Create solver
        problem     = LinearVariationalProblem(A, L, u, bc)
        self.solver = LinearVariationalSolver(problem) #, solver_parameters={'linear_solver': 'iterative'})
        self.solver.parameters["linear_solver"] ="iterative"
        # y[k] = assemble(myAverage(mesh, u, dx))

        # Compute solution
        self.solver.solve()
        # print("Hello there! This is u: {}".format(u.vector().get_local()))

        return u


        

