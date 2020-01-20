import heterocl as hcl
import numpy as np

import math

class grid:
  def __init__(self, max, min, pts_each_dim, pDim):
      self.max = max
      self.min = min
      self.pts_each_dim = pts_each_dim
      self.pDim = pDim
      self.dx = (max - min)/(pts_each_dim - 1.0)


# NOTE: No information about structure of grid, dynamics of system passed in for now
# Hardcoded these information

# Global variables

def HJ_PDE_solver(V_new, V_init, thetas):
    # Calculate spatial derivative based on index and dimension number
    def spatial_derivative(i,j,k,dim):
        left = i*j - k
        right = i*j + k
        return left, right
    # Calculate Hamiltonian for every grid point in V_init
    with hcl.Stage("Hamiltonian"):
        with hcl.for_(1, V_init.shape[0], name="i") as i:
            with hcl.for_(1, V_init.shape[1], name="j") as j:
                with hcl.for_(1, V_init.shape[2], name="k") as k:
                    # Calculate dV_dx
                    dV_dx_L, dV_dx_R = spatial_derivative(i,j,k,0)
                    dV_dy_L, dV_dy_R = spatial_derivative(i,j,k,1)
                    dV_dtheta_L, dV_dtheta_R = spatial_derivative(i, j, k, 2)

                    # Calculate average gradient
                    dV_dx_C = (dV_dx_L + dV_dx_R)/2
                    dV_dy_C = (dV_dy_L + dV_dy_R) / 2
                    dV_dtheta_C = (dV_dtheta_L + dV_dtheta_R) / 2

                    # Get optimal control
                    uOpt = 1

                    # Velocity
                    v = 1

                    # Assume that mode is min
                    with hcl.if_(dV_dtheta_C > 0):
                        uOpt = -uOpt

                    # Calculate dynamics function
                    #V_new[i,j,k] = 1 * cos(thetas[k]) * dV_dx_C +1 * sin(thetas[k]) * dV_dy_C +uOpt * dV_theta_C
                    #angle = hcl.scalar(thetas[k], "angle")
                    V_new[i,j,k] = v * hcl.cos(thetas[k]) * dV_dx_C + v * hcl.sin(thetas[k]) * dV_dy_C +  dV_dtheta_C * uOpt

def main():
    hcl.init()
    V_f = hcl.placeholder((50, 50, 50), name="V_f", dtype = hcl.Float())
    V_init = hcl.placeholder((50, 50, 50), name="V_init", dtype=hcl.Float())
    thetas = hcl.placeholder((50,), name="thetas", dtype=hcl.Float())

    # Create schedule
    s = hcl.create_schedule([V_f, V_init, thetas], HJ_PDE_solver)

    # Inspect IR
    print(hcl.lower(s))

if __name__ == '__main__':
  main()