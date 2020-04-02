from dolfin import *
import numpy as np
from ..     import SPDEModel

__author__ = ["Benjamin, Bykowski", "Jean-Luc Bouchot"]
__copyright__ = "Copyright 2015-2020, Chair C for Mathematics (Analysis), RWTH Aachen; Seminar for Applied Mathematics, ETH Zurich; School of Mathematics and Statistics, Beijing Institute of Technology"
__credits__ = ["Jean-Luc Bouchot", "Benjamin, Bykowski", "Holger Rauhut", "Christoph Schwab"]
__license__ = "GPLv3"
__version__ = "0.1.0-dev"
__maintainer__ = "Jean-Luc Bouchot"
__email__ = "jlbmathit@gmail.com"
__status__ = "Development"
__lastmodified__ = "2020/04/02"

class FEMModel(SPDEModel):
    def init_simple_mesh(self):
        # Create mesh and define function space
        if type(self.mesh_size) is tuple:
            if len(self.mesh_size) == 1:
                self.mesh = UnitIntervalMesh(*self.mesh_size)
            elif len(self.mesh_size) == 2:
                self.mesh = UnitSquareMesh(*self.mesh_size)
            elif len(self.mesh_size) == 3:
                self.mesh = UnitCubeMesh(*self.mesh_size)
            else:
                assert False, "Only one to three dimensional problems supported"
        else:
            self.mesh = UnitIntervalMesh(self.mesh_size)

    def generate_Hilbert_basis(self, doGS = False): 
        # This is basically a Gram Schmidt orthonormalization of the Finite Element basis


        V = FunctionSpace(self.mesh, 'Lagrange', 1)
        if (doGS): # Do Gram Schmidt!
            print("-->Computing orthonormal basis ... ")
            nb_dim = V.dim()
            running_vector = Function(V) # The one that is being orthonomalize
            projection_vector = Function(V) # Vector that will contain the "subspace" whose contribution we are getting rid
            ONBasis = np.zeros((nb_dim,nb_dim)) # Keep the coordinates on the FE basis of the ON Basis
            coord_vector = np.zeros(nb_dim) # Generate the canonical basis on the fly
        
            # Deal with the first basis vector
            coord_vector[0] = 1
            running_vector.vector().set_local(coord_vector)
            cur_energy = np.sqrt(assemble(inner(running_vector, running_vector)*dx))
            ONBasis[:,0] = coord_vector/cur_energy
            
            # Now on to the remaining vectors
            for i in range(nb_dim-1): 
                coord_vector = np.zeros(nb_dim) # Deletes the previous contribution
                coord_vector[i+1] = 1 # Creates the new canonical vector (with respect to the FE basis) 
                for j in range(i+1): 
                    cur_proj_coordinates = ONBasis[:,j]
                    projection_vector.vector().set_local(cur_proj_coordinates)
                    running_vector.vector().set_local(coord_vector) # That's the current function we'll be working with! 
                    
                    coord_vector = coord_vector - assemble(inner(projection_vector,running_vector)*dx)*cur_proj_coordinates
            
                running_vector.vector().set_local(coord_vector)
                cur_energy = np.sqrt(assemble(inner(running_vector,running_vector)*dx))
                ONBasis[:,i+1] = coord_vector/cur_energy
                # And the loop is done!

            print("-->DONE: Orthonormal basis computed!")
            np.save("ONB_big.npy", ONBasis)
            np.save("ONB_big_fctSpace.npy", V)
            print("-->SAVED: ONB in file ONB.npy")

        else: 
            nb_dim = V.dim()
            ONBasis = np.eye(nb_dim) # Keep the coordinates on the FE basis of the ON Basis

        self.hilbert_basis = ONBasis


    def refine_mesh(self, ratio=2): # Note, this can also be used to coarsen the mesh
        self.mesh_size = tuple(int(one_direction*ratio) for one_direction in self.mesh_size)
        self.init_simple_mesh()

    # @staticmethod
    def split_params(self, coeff, z):
        cur_split = 0
        params    = []

        for c in coeff:
            z_c = np.array(0)

            if c.num_params > 0:
                z_c = z[cur_split:(cur_split+c.num_params)]
                cur_split += c.num_params

            params.append(z_c)
        return params

    def sample(self, z):
        u = self.solve(z)
        # Return correlations with basis functions
        # return u.vector().get_coordinates() # XXX Change this to return the inner product with the basis /!\/!\
        # return self.get_coordinates_on_cell(u) 
        return self.get_coordinates(u) 
        # print(u.vector().get_local())
        # return u.vector().get_local()

    def get_coordinates_on_cell(self,a_sol):
        V = FunctionSpace(self.mesh, 'DG',0)
        projection = project(a_sol,V)
        return projection.vector().get_local()


    def get_coordinates(self,a_sol): # Surely we can do this better!
        # print("This is the current solution: {}".format(a_sol.vector().get_local())) 
        # return a_sol.vector().get_local()
        # local_data = self.find_coef_basis(a_sol,self.hilbert_basis,self.mesh,FunctionSpace(self.mesh, 'Lagrange', 1))
        # print("The computed coefficients are {} while the matrix based coefficients are {}".format(local_data,np.linalg.solve(self.hilbert_basis,a_sol.vector().get_local())))
        # print("ORIGINAL DATA = {} vs \n RECONSTRUCTED FROM INNER PRODUCTS = {}".format(np.linalg.norm(a_sol.vector().get_local()), np.linalg.norm( np.dot(self.hilbert_basis, local_data))))
        # print("ORIGINAL DATA = {} vs \n RECONSTRUCTED FROM INNER PRODUCTS = {}".format(a_sol.vector().get_local(), np.dot(self.hilbert_basis, local_data)))
        return self.find_coef_basis(a_sol,self.hilbert_basis,self.mesh,FunctionSpace(self.mesh, 'Lagrange', 1))



    def find_coef_basis(self, a_sol, some_functions, mesh, fctSpace):
        x = SpatialCoordinate(mesh)
        nb_dim,nb_fct = some_functions.shape # This only works if you are dealing with matrices of coordinate on the basis given by the fctSpace. We should be able to extend this to a list of generic vectors.
        # Note the the nmber of functions needs not be equal to the dimension. We could also use this for redundant basis of functions
        
        correlations = np.zeros(nb_fct) # Keeps the correlation with the basis functions somewhere. 
    
        cur_fct = Function(fctSpace)
        for one_fct in range(nb_fct): 
            cur_fct.vector().set_local(some_functions[:,one_fct])
            correlations[one_fct] = assemble(inner(a_sol,cur_fct)*dx)
    
        return correlations


    def __getstate__(self):
        odict = self.__dict__.copy()

        # Can't pickle mesh, solver and M
        del odict['mesh']
        del odict['solver']
        del odict['M']

        return odict
