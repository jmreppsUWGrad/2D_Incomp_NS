# -*- coding: utf-8 -*-
"""
2D Incompressible Navier-Stokes solver

Boundary condition options:
    -'wall' for bc_type will imply no-slip and dp of 0
    -'zero_grad' will impose 0 normal gradient of that variable
    
Features to include (across all classes):
    -CoolProp library for material properties (track down needed functions)
    -Fix biasing meshing tools (this script and GeomClasses)
        ->Figure out biasing wrt dx and dy array sizes and mesh griding those (GeomClasses)
    -File reader for settings
    -periodic boundary conditions (SolverClasses)

@author: Joseph
"""

##########################################################################
# ----------------------------------Libraries and classes
##########################################################################
import numpy
from matplotlib import pyplot, cm
from mpl_toolkits.mplot3d import Axes3D
from datetime import datetime
import os
import CoolProp.CoolProp as CP

#from GeomClasses import OneDimLine as OneDimLine
from GeomClasses import TwoDimPlanar as TwoDimPlanar
#import MatClasses as Mat
import SolverClasses as Solvers
import FileClasses

##########################################################################
# ------------------------------ Geometry, Domain and BCs Setup
#    Reference directions:
#    left-smallest x coordinate
#    right-largest x value
#    north-largest y coordinate
#    south-smallest y coordinate
##########################################################################
settings={} # Dictionary of problem settings
BCs={} # Dictionary of boundary conditions
# Geometry details
settings['Length']                  = 2.0
settings['Width']                   = 2.0
settings['Nodes_x']                 = 41
settings['Nodes_y']                 = 41
settings['Fluid']                   = 'Air'
#CP.PropsSI('L','T', 300, 'P', 101325, settings['Fluid'])
settings['mu']                      = 0.1#CP.PropsSI('V','T', 300, 'P', 101325, settings['Fluid'])
settings['rho']                     = 1#CP.PropsSI('D','T', 300, 'P', 101325, settings['Fluid'])
settings['Gravity_x']               = 0
settings['Gravity_y']               = 0
settings['Pressure_grad_x']         = 0
settings['Pressure_grad_y']         = 0
settings['All_pressure_terms']      = False

# Meshing details
"""
Biasing options:
    -'OneWayUp'   for linearly increasing element sizes with increasing x/y
    -'OneWayDown' for linearly decreasing element sizes with increasing x/y
    -'TwoWayEnd'  for linearly increasing sizes till middle, then decrease again
    -'TwoWayMid'  for linearly decreasing sizes till middle, then increase again
    -size         is the smallest element size based on above selection
"""
settings['bias_type_x']             = None
settings['bias_size_x']             = 0.005 # Smallest element size (IN PROGRESS)
settings['bias_type_y']             = None
settings['bias_size_y']             = 0.00005 # Smallest element size (IN PROGRESS)

# Boundary conditions
BCs['bc_type_left']                 = 'wall'
BCs['bc_left_u']                    = None
BCs['bc_left_v']                    = None
BCs['bc_left_p']                    = None
# numpy.linspace(400, 900, settings['Nodes_y'])
BCs['bc_type_right']                = 'wall'
BCs['bc_right_u']                   = None
BCs['bc_right_v']                   = None
BCs['bc_right_p']                   = None
# numpy.linspace(400, 900, settings['Nodes_y'])
BCs['bc_type_south']                = 'wall'
BCs['bc_south_u']                   = None
BCs['bc_south_v']                   = None
BCs['bc_south_p']                   = None
# numpy.linspace(400, 900, settings['Nodes_x'])
BCs['bc_type_north']                = 'top'
BCs['bc_north_u']                   = 1
BCs['bc_north_v']                   = 0
BCs['bc_north_p']                   = 0
# numpy.linspace(400, 900, settings['Nodes_x'])

# Initial conditions ????

# Time advancement
settings['CFL']                     = 0.8
settings['dt']                      = 0.001
settings['total_time_steps']        = 3
settings['total_time']              = 'None'
settings['Convergence']             = 0.01
settings['Max_iterations']          = 1000
settings['Time_Scheme']             = 'Euler'
settings['Number_Data_Output']      = 1
#settings['Output_directory']        = 'C:\Users\mepps\Documents\Research\2D_compressible_NS\Tests'
settings['Output_directory']        = 'C:\Users\Joseph\Documents\School\Research\2D_Incomp_NS\Tests'

print('######################################################')
print('#      2D Incompressible Navier-Stokes Solver        #')
print('#              Created by J. Mark Epps               #')
print('#          Part of Masters Thesis at UW 2018-2020    #')
print('######################################################\n')
print 'Initializing geometry package...'
#domain=OneDimLine(L,Nx)
domain=TwoDimPlanar(settings)
domain.mesh()
print '################################'

##########################################################################
# -------------------------------------Initialize solver and domain
##########################################################################

print 'Initializing solver package...'
solver=Solvers.TwoDimPlanarSolve(domain, settings, BCs)
print '################################'
print 'Initializing domain...'

domain.u[:,:]=0
domain.v[:,:]=0
domain.p[:,:]=0
#solver.Apply_BCs(u, v, p, True)

print '################################'
##########################################################################
# -------------------------------------File setups
##########################################################################
#print 'Initializing files...'
#os.chdir('Tests')
#datTime=str(datetime.date(datetime.now()))+'_'+'{:%H%M}'.format(datetime.time(datetime.now()))
#isBinFile=False
#
##output_file=FileClasses.FileOut('Output_'+datTime, isBinFile)
#input_file=FileClasses.FileOut('Input_'+datTime, isBinFile)
#
## Write headers to files
#input_file.header('INPUT')
##output_file.header('OUTPUT')
#
## Write input file with settings
#input_file.input_writer(settings, BCs, domain.rho, domain.rhou, domain.rhov, domain.rhoE)
#input_file.close()
#print '################################'

##########################################################################
# -------------------------------------Solve
##########################################################################
print 'Solving:'
t,nt=0,0
output_data_t,output_data_nt=0,0
if settings['total_time_steps']=='None':
    settings['total_time_steps']=settings['total_time']*10**9
    output_data_t=settings['total_time']/settings['Number_Data_Output']
elif settings['total_time']=='None':
    settings['total_time']=settings['total_time_steps']
    output_data_nt=int(settings['total_time_steps']/settings['Number_Data_Output'])

#while nt<settings['total_time_steps'] and t<settings['total_time']:
for nt in range(settings['total_time_steps']):
#    print 'Time step %i of %i'%(nt+1, settings['total_time_steps'])
    err, dt=solver.Advance_Soln()
    print 'Time step %i, Step size=%.6f, Time elapsed=%f;'%(nt+1,dt, t)
    t+=dt
#    nt+=1
    if err==1:
        print '#################### Solver aborted #######################'
        break

#output_file.close()

##########################################################################
# ------------------------------------Post-processing
##########################################################################
X,Y,u,v,p=domain.X,domain.Y,domain.u,domain.v,domain.p
# 2D plot
#fig=pyplot.figure(figsize=(7, 7))
#ax = fig.gca(projection='3d')
#ax.plot_surface(domain.X, domain.Y, T, rstride=1, cstride=1, cmap=cm.viridis,linewidth=0, antialiased=True)
##ax.set_xlim(0,0.001)
##ax.set_ylim(0.005,0.006)
#ax.set_zlim(300, BCs['bc_south_T'])
#ax.set_xlabel('$x$ (m)')
#ax.set_ylabel('$y$ (m)')
#ax.set_zlabel('T (K)');
#fig.savefig('Plot1.png',dpi=300)

# 1D Plot
#fig2=pyplot.figure(figsize=(7,7))
#pyplot.plot(domain.Y[:,1]*1000, domain.T[:,1],marker='x')
#pyplot.xlabel('$y$ (mm)')
#pyplot.ylabel('T (K)')
#pyplot.title('Temperature distribution at 2nd x')
#pyplot.xlim(5,6);
#fig2.savefig('Plot2.png',dpi=300)

# Velocity Quiver plot and pressure contour
pl=5
fig3=pyplot.figure(figsize=(7, 7))
pyplot.quiver(X[::pl, ::pl], Y[::pl, ::pl], \
              u[::pl, ::pl], v[::pl, ::pl]) 
pyplot.contourf(X, Y, p, alpha=0.5, cmap=cm.viridis)  
pyplot.colorbar()
pyplot.xlabel('$x$ (m)')
pyplot.ylabel('$y$ (m)')
pyplot.title('Velocity plot and Pressure contours');
#fig3.savefig(datTime+'_Vel_Press.png',dpi=300)

print('Solver has finished its run')