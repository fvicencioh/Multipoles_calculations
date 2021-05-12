from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
import numpy as np

def read_pdb(pdb_filename, ff="amoeba2013.xml"):
    
    try:
        pdb = PDBFile(pdb_filename+'.pdb')
        
    except:
        raise ValueError(F"No pdb file whit {pdb_filename}.pdb filename")
        
    modeller = Modeller(pdb.topology, pdb.positions)
    forcefield = ForceField(ff)
    
    modeller.deleteWater()
    modeller.addHydrogens(ForceField('amber14-all.xml')) #Use amber only to add hydrogens
    modeller.addExtraParticles(forcefield)
    
    system = forcefield.createSystem(modeller.topology, nonbondedMethod=PME, nonbondedCutoff=1*nanometer,
    vdwCutoff=1.2*nanometer, polarization='mutual', mutualInducedTargetEpsilon=0.0001)
    
    return system, modeller

def get_xyzr_file(pdb_filename, ff="amoeba2013.xml"):
    
    system, modeller = read_pdb(pdb_filename, ff)
    
    xyzr_file = open(pdb_filename+".xyzr","w")
    
    Nq = system.getNumParticles()
    
    force_parameters = system.getForces()
    
    for i in range(Nq):
        
        pos = modeller.getPositions()[i].value_in_unit(nanometer)
        radius = force_parameters[8].getParticleParameters(i)[1].value_in_unit(nanometer)
        
        line = str(pos[0])+" "+str(pos[1])+" "+str(pos[2])+" "+str(radius)+"\n"
        xyzr_file.write(line)
        
    xyzr_file.close()
    
def get_parameters(pdb_filename, ff="amoeba2013.xml"):
    
    system, modeller = read_pdb(pdb_filename, ff)
    
    Nq = system.getNumParticles()
    force_parameters = system.getForces()
    
    q = np.zeros(Nq)
    p = np.zeros((Nq, 3))
    Q = np.zeros((Nq, 3, 3))
    thole = np.zeros(Nq)
    
    for i in range(Nq):
        
        param = force_parameters[9].getMultipoleParameters(i)
        
        q[i] = param[0].value_in_unit(elementary_charge)
        p[i,:] = param[1].value_in_unit(nanometer*elementary_charge)
        Q[i,0,:] = param[2].value_in_unit(nanometer**2*elementary_charge)[:3]
        Q[i,1,:] = param[2].value_in_unit(nanometer**2*elementary_charge)[3:6]
        Q[i,2,:] = param[2].value_in_unit(nanometer**2*elementary_charge)[6:]
        
        thole[i] = param[7]
        
        if Q[i]!=Q[i].transpose(): #Easy Test
            print(F"Bad assignment in Quadrupole {i}")
        
    return q, p, Q, thole

def get_induced_dipole_openmm(pdb_file, ff="amoeba2013.xml", maxiter=100):
    
    system, modeller = read_pdb(pdb_filename, ff)
    integrator = LangevinMiddleIntegrator(300*kelvin, 1/picosecond, 0.004*picoseconds)
    simulation = Simulation(modeller.topology, system, integrator)
    
    simulation.context.setPositions(modeller.positions)
    simulation.minimizeEnergy(maxIterations=maxiter)
    
    last_step = 0
    for i in range(maxiter):
        try:
            simulation.step(i)
        except:
            simulation.step(last_step)
            break
        last_step = i
    
    multipoles = system.getForces()[9]
    context = simulation.context
    
    dipoles = multipoles.getInducedDipoles(context)
    
    return dipoles