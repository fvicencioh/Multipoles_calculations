import bempp.api
import numpy as np
import os
from bempp.api.operators.boundary import sparse, laplace, modified_helmholtz
import inspect
from scipy.sparse.linalg import gmres


from .direct.direct import compute_induced_dipole, coulomb_energy_multipole
from .util.getData import *


def generate_nanoshaper_grid(filename):
    
    """
    filename: Name of mesh files, whitout .face or .vert extension
    
    Returns:
    
    grid: Bempp Grid
        
    """
    
    faces = np.loadtxt(filename+".face", dtype= int) - 1
    verts = np.loadtxt(filename+".vert", dtype= float)
    
    grid = bempp.api.Grid(verts.transpose(), faces.transpose())
    
    return grid

def getLHS(dirichl_space, neumann_space, ep_in, ep_ex, kappa, assembler="dense"):
    
    identity = sparse.identity(dirichl_space, dirichl_space, dirichl_space);
    VL       = laplace.single_layer(neumann_space, dirichl_space, dirichl_space, assembler=assembler);
    KL       = laplace.double_layer(dirichl_space, dirichl_space, dirichl_space, assembler=assembler);
    VY       = modified_helmholtz.single_layer(neumann_space, dirichl_space, dirichl_space, kappa, assembler=assembler);
    KY       = modified_helmholtz.double_layer(dirichl_space, dirichl_space, dirichl_space, kappa, assembler=assembler);
    
    blocked = bempp.api.BlockedOperator(2, 2);
    blocked[0, 0] = 0.5*identity + KL;
    blocked[0, 1] = -VL;
    blocked[1, 0] = 0.5*identity - KY;
    blocked[1, 1] = (ep_in/ep_ex)*VY;
    LHS = blocked.strong_form();
    
    return LHS

def getRHS(x_q, q, p, Q, ep_in, ep_ex, dirichl_space, neumann_space):
    
    @bempp.api.real_callable
    def multipolar_quadrupoles_charges_fun(x, n, i, result):
        T2 = np.zeros((len(x_q),3,3))
        phi = 0
        dist = x - x_q
        norm = np.sqrt(np.sum((dist*dist), axis = 1))
        T0 = 1/norm[:]
        T1 = np.transpose(dist.transpose()/norm**3)
        T2[:,:,:] = np.ones((len(x_q),3,3))[:]*dist.reshape((len(x_q),1,3))*np.transpose(np.ones((len(x_q),3,3))*dist.reshape((len(x_q),1,3)), (0,2,1))/norm.reshape((len(x_q),1,1))**5
        phi = np.sum(q[:]*T0[:]) + np.sum(T1[:]*p[:]) + 0.5*np.sum(np.sum(T2[:]*Q[:],axis=1))
        result[0] = (phi/(4*np.pi*ep_in))
        
    charged_grid_fun = bempp.api.GridFunction(dirichl_space, fun = multipolar_quadrupoles_charges_fun)
    RHS = np.concatenate([charged_grid_fun.coefficients, np.zeros(neumann_space.global_dof_count)])
    
    return RHS

def iteration_counter(x):
    global array_it, array_frame, it_count
    it_count += 1
    frame = inspect.currentframe().f_back
    array_it = np.append(array_it, it_count)
    array_frame = np.append(array_frame, frame.f_locals["resid"])
    print(F"Gmres iteration {it_count}, residual: {x}")

def solvation_energy_solvent(q, p, Q, phi, dphi, ddphi):
    
    cal2J = 4.184
    qe = 1.60217646e-19
    Na = 6.0221415e+23
    ep_vacc = 8.854187818e-12
    C0 = qe**2*Na*1e-3*1e10/(cal2J*ep_vacc)
    
    q_aux = 0
    p_aux = 0
    Q_aux = 0
    
    for i in range(len(q)):
        q_aux += q[i]*phi[i]
        
        for j in range(3):
            p_aux += p[i,j]*dphi[i,j]
            
            for k in range(3):
                Q_aux += Q[i,j,k]*ddphi[i,j,k]/6.
                
    solvent_energy = 0.5 * C0 * (q_aux + p_aux + Q_aux)
    
    return solvent_energy

def solvent_potential_first_derivate(xq, h, neumann_space, dirichl_space, solution_neumann, solution_dirichl):
    
    """
    Compute the first derivate of potential due to solvent
    in the position of the points
    Inputs:
    -------
        xq: Array size (Nx3) whit positions to calculate the derivate.
        h: Float number, distance for the central difference.

    Return:

        dpdr: Derivate of the potential in the positions of points.
    """

    dpdr = np.zeros([len(xq), 3])
    dist = np.array(([h,0,0],[0,h,0],[0,0,h]))
    # x axis derivate
    dx = xq[:] + dist[0]
    dx = np.concatenate((dx, xq[:] - dist[0]))
    slpo = bempp.api.operators.potential.laplace.single_layer(neumann_space, dx.transpose())
    dlpo = bempp.api.operators.potential.laplace.double_layer(dirichl_space, dx.transpose())
    phi = slpo.evaluate(solution_neumann) - dlpo.evaluate(solution_dirichl)
    dpdx = 0.5*(phi[0,:len(xq)] - phi[0,len(xq):])/h
    dpdr[:,0] = dpdx

    #y axis derivate
    dy = xq[:] + dist[1]
    dy = np.concatenate((dy, xq[:] - dist[1]))
    slpo = bempp.api.operators.potential.laplace.single_layer(neumann_space, dy.transpose())
    dlpo = bempp.api.operators.potential.laplace.double_layer(dirichl_space, dy.transpose())
    phi = slpo.evaluate(solution_neumann) - dlpo.evaluate(solution_dirichl)
    dpdy = 0.5*(phi[0,:len(xq)] - phi[0,len(xq):])/h
    dpdr[:,1] = dpdy

    #z axis derivate
    dz = xq[:] + dist[2]
    dz = np.concatenate((dz, xq[:] - dist[2]))
    slpo = bempp.api.operators.potential.laplace.single_layer(neumann_space, dz.transpose())
    dlpo = bempp.api.operators.potential.laplace.double_layer(dirichl_space, dz.transpose())
    phi = slpo.evaluate(solution_neumann) - dlpo.evaluate(solution_dirichl)
    dpdz = 0.5*(phi[0,:len(xq)] - phi[0,len(xq):])/h
    dpdr[:,2] = dpdz

    return dpdr

def solvent_potential_second_derivate(x_q, h, neumann_space, dirichl_space, solution_neumann, solution_dirichl):
    
    """
    Compute the second derivate of potential due to solvent
    in the position of the points

    xq: Array size (Nx3) whit positions to calculate the derivate.
    h: Float number, distance for the central difference.

    Return:

    ddphi: Second derivate of the potential in the positions of points.
    """
    ddphi = np.zeros((len(x_q),3,3))
    dist = np.array(([h,0,0],[0,h,0],[0,0,h]))
    for i in range(3):
        for j in np.where(np.array([0, 1, 2]) >= i)[0]:
            if i==j:
                dp = np.concatenate((x_q[:] + dist[i], x_q[:], x_q[:] - dist[i]))
                slpo = bempp.api.operators.potential.laplace.single_layer(neumann_space, dp.transpose())
                dlpo = bempp.api.operators.potential.laplace.double_layer(dirichl_space, dp.transpose())
                phi = slpo.evaluate(solution_neumann) - dlpo.evaluate(solution_dirichl)
                ddphi[:,i,j] = (phi[0,:len(x_q)] - 2*phi[0,len(x_q):2*len(x_q)] + phi[0, 2*len(x_q):])/(h**2)
      
            else:
                dp = np.concatenate((x_q[:] + dist[i] + dist[j], x_q[:] - dist[i] - dist[j], x_q[:] + dist[i] - dist[j], x_q[:] - dist[i] + dist[j]))
                slpo = bempp.api.operators.potential.laplace.single_layer(neumann_space, dp.transpose())
                dlpo = bempp.api.operators.potential.laplace.double_layer(dirichl_space, dp.transpose())
                phi = slpo.evaluate(solution_neumann) - dlpo.evaluate(solution_dirichl)
                ddphi[:,i,j] = (phi[0,:len(x_q)] + phi[0,len(x_q):2*len(x_q)] - phi[0, 2*len(x_q):3*len(x_q)] - phi[0, 3*len(x_q):])/(4*h**2)
                ddphi[:,j,i] = (phi[0,:len(x_q)] + phi[0,len(x_q):2*len(x_q)] - phi[0, 2*len(x_q):3*len(x_q)] - phi[0, 3*len(x_q):])/(4*h**2)
  
    return ddphi

def get_induced_dipole(x_q, q, p, p_pol, Q, alpha, ep_in,thole, polar_group, connections_12, connections_13, pointer_connections_12, pointer_connections_13, dphi_solvent):
    
    Nq = len(q)
    
    xq_temp = np.zeros(Nq)
    xq_temp[:] = x_q[:,0]
    
    yq_temp = np.zeros(Nq)
    yq_temp[:] = x_q[:,1]
    
    zq_temp = np.zeros(Nq)
    zq_temp[:] = x_q[:,2]
    
    dx_temp = np.zeros(Nq)
    dx_temp[:] = p[:,0]
    
    dy_temp = np.zeros(Nq)
    dy_temp[:] = p[:,1]
    
    dz_temp = np.zeros(Nq)
    dz_temp[:] = p[:,2]
    
    Qxx_temp = np.zeros(Nq)
    Qxx_temp[:] = Q[:,0,0]
    
    Qxy_temp = np.zeros(Nq)
    Qxy_temp[:] = Q[:,0,1]
    
    Qxz_temp = np.zeros(Nq)
    Qxz_temp[:] = Q[:,0,2]
    
    Qyx_temp = np.zeros(Nq)
    Qyx_temp[:] = Q[:,1,0]
    
    Qyy_temp = np.zeros(Nq)
    Qyy_temp[:] = Q[:,1,1]
    
    Qyz_temp = np.zeros(Nq)
    Qyz_temp[:] = Q[:,1,2]
    
    Qzx_temp = np.zeros(Nq)
    Qzx_temp[:] = Q[:,2,0]
    
    Qzy_temp = np.zeros(Nq)
    Qzy_temp[:] = Q[:,2,1]
    
    Qzz_temp = np.zeros(Nq)
    Qzz_temp[:] = Q[:,2,2]
    
    alphaxx_temp = np.zeros(Nq)
    alphaxx_temp[:] = alpha[:,0,0]
    
    alphaxy_temp = np.zeros(Nq)
    alphaxy_temp[:] = alpha[:,0,1]
    
    alphaxz_temp = np.zeros(Nq)
    alphaxz_temp[:] = alpha[:,0,2]
    
    alphayx_temp = np.zeros(Nq)
    alphayx_temp[:] = alpha[:,1,0]
    
    alphayy_temp = np.zeros(Nq)
    alphayy_temp[:] = alpha[:,1,1]
    
    alphayz_temp = np.zeros(Nq)
    alphayz_temp[:] = alpha[:,1,2]
    
    alphazx_temp = np.zeros(Nq)
    alphazx_temp[:] = alpha[:,2,0]
    
    alphazy_temp = np.zeros(Nq)
    alphazy_temp[:] = alpha[:,2,1]
    
    alphazz_temp = np.zeros(Nq)
    alphazz_temp[:] = alpha[:,2,2]
    
    px_pol = np.zeros(Nq)
    px_pol[:] = p_pol[:,0]
        
    py_pol = np.zeros(Nq)
    py_pol[:] = p_pol[:,1]
        
    pz_pol = np.zeros(Nq)
    pz_pol[:] = p_pol[:,2]
        
    dphix_temp = np.zeros(Nq)
    dphix_temp[:] = dphi_solvent[:,0]
        
    dphiy_temp = np.zeros(Nq)
    dphiy_temp[:] = dphi_solvent[:,1]
        
    dphiz_temp = np.zeros(Nq)
    dphiz_temp[:] = dphi_solvent[:,2]
    
    compute_induced_dipole(xq_temp, yq_temp, zq_temp, q, \
                               dx_temp, dy_temp, dz_temp, \
                               px_pol, py_pol, pz_pol, \
                               Qxx_temp, Qxy_temp, Qxz_temp, \
                               Qyx_temp, Qyy_temp, Qyz_temp, \
                               Qzx_temp, Qzy_temp, Qzz_temp, \
                               alphaxx_temp, alphaxy_temp, alphaxz_temp, \
                               alphayx_temp, alphayy_temp, alphayz_temp, \
                               alphazx_temp, alphazy_temp, alphazz_temp, \
                               thole, np.int32(polar_group), \
                               np.int32(connections_12), np.int32(pointer_connections_12), \
                               np.int32(connections_13), np.int32(pointer_connections_13), \
                               dphix_temp, dphiy_temp, dphiz_temp, ep_in)
    
    p_pol[:,0] = px_pol
    p_pol[:,1] = py_pol
    p_pol[:,2] = pz_pol
    
    return p_pol

def get_coulomb_energy(x_q, q, p, p_pol, Q, alpha, ep_in, thole, polar_group, connections_12, connections_13, pointer_connections_12, pointer_connections_13, p12scale, p13scale):
    
    
    Nq = len(q)
    
    xq_temp = np.zeros(Nq)
    xq_temp[:] = x_q[:,0]
    
    yq_temp = np.zeros(Nq)
    yq_temp[:] = x_q[:,1]
    
    zq_temp = np.zeros(Nq)
    zq_temp[:] = x_q[:,2]
    
    dx_temp = np.zeros(Nq)
    dx_temp[:] = p[:,0]
    
    dy_temp = np.zeros(Nq)
    dy_temp[:] = p[:,1]
    
    dz_temp = np.zeros(Nq)
    dz_temp[:] = p[:,2]
    
    Qxx_temp = np.zeros(Nq)
    Qxx_temp[:] = Q[:,0,0]
    
    Qxy_temp = np.zeros(Nq)
    Qxy_temp[:] = Q[:,0,1]
    
    Qxz_temp = np.zeros(Nq)
    Qxz_temp[:] = Q[:,0,2]
    
    Qyx_temp = np.zeros(Nq)
    Qyx_temp[:] = Q[:,1,0]
    
    Qyy_temp = np.zeros(Nq)
    Qyy_temp[:] = Q[:,1,1]
    
    Qyz_temp = np.zeros(Nq)
    Qyz_temp[:] = Q[:,1,2]
    
    Qzx_temp = np.zeros(Nq)
    Qzx_temp[:] = Q[:,2,0]
    
    Qzy_temp = np.zeros(Nq)
    Qzy_temp[:] = Q[:,2,1]
    
    Qzz_temp = np.zeros(Nq)
    Qzz_temp[:] = Q[:,2,2]
    
    alphaxx_temp = np.zeros(Nq)
    alphaxx_temp[:] = alpha[:,0,0]
    
    px_pol = np.zeros(Nq)
    px_pol[:] = p_pol[:,0]
        
    py_pol = np.zeros(Nq)
    py_pol[:] = p_pol[:,1]
        
    pz_pol = np.zeros(Nq)
    pz_pol[:] = p_pol[:,2]
    
    point_energy = np.zeros(Nq)
    
    cal2J = 4.184
    ep_vacc = 8.854187818e-12
    qe = 1.60217646e-19
    Na = 6.0221415e+23
    C0 = qe**2*Na*1e-3*1e10/(cal2J*ep_vacc)
    
    coulomb_energy_multipole(xq_temp, yq_temp, zq_temp, q, \
                             dx_temp, dy_temp, dz_temp, \
                             px_pol, py_pol, pz_pol, \
                             Qxx_temp, Qxy_temp, Qxz_temp, \
                             Qyx_temp, Qyy_temp, Qyz_temp, \
                             Qzx_temp, Qzy_temp, Qzz_temp, \
                             alphaxx_temp, thole, np.int32(polar_group), \
                             np.int32(connections_12), np.int32(pointer_connections_12), \
                             np.int32(connections_13), np.int32(pointer_connections_13), \
                             p12scale, p13scale, point_energy)
    
    coulomb_energy = np.sum(point_energy) * 0.5*C0/(4*np.pi*ep_in)
    
    return coulomb_energy

def multipole_calculations(filename, grid_filename, ep_in, ep_ex, k, h, maxiter=100, gmres_maxiter=500, mu="None", tol=1e-2, gmres_tol=1e-5):
    
    global array_it, array_frame, it_count
    
    x_q, q, d, Q, alpha, mass, polar_group, thole, \
           connections_12, connections_13, \
           pointer_connections_12, pointer_connections_13, \
           p12scale, p13scale, N = read_tinker(filename,float)
    
    Nq = len(q)
    
    if mu=="None":
        mu = np.zeros((Nq,3))
    
    grid = generate_nanoshaper_grid(grid_filename)
    
    dirichl_space = bempp.api.function_space(grid, "DP", 0)
    neumann_space = bempp.api.function_space(grid, "DP", 0)
    
    x_b = np.zeros((2*dirichl_space.global_dof_count))
    
    assembler="dense"
    
    if grid.number_of_elements>50000:
        assembler="fmm"
    
    lhs = getLHS(dirichl_space, neumann_space, ep_in, ep_ex, k, assembler=assembler)
    
    for iter_number in range(maxiter):
        
        print(F"------- Iteration {iter_number +1} -------")
        
        array_it, array_frame, it_count = np.array([]), np.array([]), 0
        
        p = d + mu #d is the permanent dipole, mu is the polarizable dipole
        
        rhs = getRHS(x_q, q, p, Q, ep_in, ep_ex, dirichl_space, neumann_space)
        
        x, info = gmres(lhs, rhs, x0=x_b, callback=iteration_counter, \
                        callback_type="pr_norm", tol=gmres_tol, maxiter=gmres_maxiter, restart = 1000)
        x_b = x.copy()
        
        solution_dirichl = bempp.api.GridFunction(dirichl_space, coefficients=x[:dirichl_space.global_dof_count])
        solution_neumann = bempp.api.GridFunction(neumann_space, coefficients=x[dirichl_space.global_dof_count:])
        
        #Calculation dphi due to solvent.
        
        dphi_solvent = solvent_potential_first_derivate(x_q, h, neumann_space, dirichl_space, solution_neumann, solution_dirichl)
        
        #Calculation of induced dipole.
        
        mu_b = mu.copy()
        
        mu = get_induced_dipole(x_q, q, d, mu, Q, alpha, ep_in, thole, polar_group, \
                                connections_12, connections_13, pointer_connections_12, pointer_connections_13, dphi_solvent)
        
        dipole_diff = np.max(np.sqrt(np.sum((np.linalg.norm(mu_b-mu,axis=1))**2)/len(mu)))
        if dipole_diff<tol:
            print(F"The induced dipole in dissolved state has converged in {iter_number+1} iterations")
            break
            
        print(F"Induced dipole residual: {dipole_diff}")
        
    #Calculation of phi and ddphi due to solvent once induced dipole has converged
    
    slpo = bempp.api.operators.potential.laplace.single_layer(neumann_space, x_q.transpose())
    dlpo = bempp.api.operators.potential.laplace.double_layer(dirichl_space, x_q.transpose())
    
    phi_solvent = slpo.evaluate(solution_neumann) - dlpo.evaluate(solution_dirichl)
    
    ddphi_solvent = solvent_potential_second_derivate(x_q, h, neumann_space, dirichl_space, solution_neumann, solution_dirichl)
    
    G_diss_solv = solvation_energy_solvent(q, d, Q, phi_solvent[0], dphi_solvent, ddphi_solvent)
    G_diss_mult = get_coulomb_energy(x_q, q, d, mu, Q, alpha, ep_in, thole, polar_group, connections_12, \
                                     connections_13, pointer_connections_12, pointer_connections_13, p12scale, p13scale)
    
    #Calcution of induced dipole in vacum:
    
    p_pol_vacc = np.zeros((Nq,3))
    
    dipole_diff_vacc = 1.
    iteration = 0
    
    while dipole_diff_vacc>tol:
        
        iteration += 1 
        
        p_pol_prev = p_pol_vacc.copy()
        
        p_pol_vacc = get_induced_dipole(x_q, q, d, p_pol_vacc, Q, alpha, ep_in, thole, polar_group, \
                                        connections_12, connections_13, pointer_connections_12, pointer_connections_13, np.zeros((Nq,3)))
        
        dipole_diff_vacc = np.max(np.sqrt(np.sum((np.linalg.norm(p_pol_prev-p_pol_vacc,axis=1))**2)/len(p_pol_vacc)))
        
        print(F"Induced dipole residual in vacuum: {dipole_diff_vacc}")
        
    print(F"{iteration} iterations for vacuum induced dipole to converge")
    
    G_vacc = get_coulomb_energy(x_q, q, d, p_pol_vacc, Q, alpha, ep_in, thole, polar_group, connections_12, \
                                     connections_13, pointer_connections_12, pointer_connections_13, p12scale, p13scale)
    
    print(F"Solvent contribution: {G_diss_solv} [kcal/Mol]")
    print(F"Multipoles contribution: {G_diss_mult} [kcal/Mol]")
    print(F"Coulomb vacuum energy: {G_vacc}")
    
    print("These values consider the polarization energy")
    
    total_energy = G_diss_solv + G_diss_mult - G_vacc
    
    print(F"Total solvation energy: {total_energy} [kcal/Mol]")
    
    return total_energy
    
    