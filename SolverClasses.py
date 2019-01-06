# -*- coding: utf-8 -*-
"""
Created on Sat Sep 29 13:17:11 2018

@author: Joseph

Solver classes for Compressible N-S equations. Takes in given object (geometry),
time step and convergence information and alters the object's temperature, 
velocity, pressure, density. BCs are applied as appropriate, but must be 
defined and copied into the solver object.

Assumptions:
    -equal discretization spacings in either x or y
    -constant thermal conductivity for conduction
    -constant viscosity for shear stress

Features:
    -Conservative Fourrier number correction based on smallest discretization
    -


"""

import numpy
#import GeomClasses
#import MatClasses
import CoolProp.CoolProp as CP
import temporal_schemes

# 1D Solvers (CURRENTLY ONLY FOR CONDUCTION)
class OneDimSolve():
    def __init__(self, geom, timeSize, timeSteps, conv):
        self.Domain=geom # Geometry object
        self.dt=timeSize
        self.Nt=timeSteps
        self.conv=conv
        self.T=self.Domain.T
        self.dx=self.Domain.dx
        self.maxCount=1000
        self.Fo=1.0*self.Domain.mat_prop['k']*self.dt\
        /(self.Domain.mat_prop['rho']*self.Domain.mat_prop['Cp'])
        self.BCs={'BCx1': ('T',600,(0,-1)),\
                 'BCx2': ('T',300,(0,-1)),\
                 'BCy1': ('T',600,(0,-1)),\
                 'BCy2': ('T',300,(0,-1))\
                 }
    
    # Convergence checker
    def CheckConv(self, Tprev, Tnew):
        diff=numpy.sum(numpy.abs(Tnew[:]-Tprev[:]))/numpy.sum(numpy.abs(Tprev[:]))
        print(diff)
        if diff<=self.conv:
            return True
        else:
            return False
    # Solve
    def SolveExpTrans(self):
        Tc=numpy.empty_like(self.T)
        for i in range(self.Nt):
            Tc=self.T.copy()
            self.T[1:-1]=2*self.Fo/(self.dx[:-1]+self.dx[1:])*(Tc[:-2]/self.dx[:-1]+Tc[2:]/self.dx[1:])\
            +(1-2*self.Fo/(self.dx[:-1]+self.dx[1:])*(1/self.dx[:-1]+1/self.dx[1:]))*Tc[1:-1]
        
    def SolveSS(self):
        Tc=numpy.empty_like(self.T)
        count=0
        print 'Residuals:'
        while count<self.maxCount:
            Tc=self.T.copy()
            self.T[1:-1]=(self.dx[1:]*Tc[:-2]+self.dx[:-1]*Tc[2:])\
            /(self.dx[1:]+self.dx[:-1])
            if self.CheckConv(Tc,self.T):
                break

# 2D solver
class TwoDimPlanarSolve():
    def __init__(self, geom_obj, settings, BCs):
        self.Domain=geom_obj # Geometry object
        self.CFL=settings['CFL']
        self.time_scheme=settings['Time_Scheme']
#        self.Nt=settings['total_time_steps']
        self.conv=settings['Convergence']
        self.maxCount=settings['Max_iterations']
        self.gx=settings['Gravity_x']
        self.gy=settings['Gravity_y']
        self.dpx=settings['Pressure_grad_x']
        self.dpy=settings['Pressure_grad_y']
        self.higher_p=settings['All_pressure_terms']
        self.dx,self.dy=numpy.meshgrid(geom_obj.dx,geom_obj.dy)
        self.BCs=BCs
    
    # Time step check with dx, dy, T and CFL number
    def getdt(self):
        dx=numpy.sqrt(self.dx**2+self.dy**2)
        
        dx=numpy.zeros_like(self.dx)
        dx[1:-1,1:-1]=0.5*numpy.sqrt((self.dx[1:-1,1:-1]+self.dx[1:-1,:-2])**2+\
                  (self.dy[1:-1,1:-1]+self.dy[:-2,1:-1])**2)
        dx[0,0]      =0.5*numpy.sqrt((self.dx[0,0])**2+(self.dy[0,0])**2)
        dx[0,1:-1]   =0.5*numpy.sqrt((self.dx[0,1:-1]+self.dx[0,:-2])**2+\
                  (self.dy[0,1:-1])**2)
        dx[1:-1,0]   =0.5*numpy.sqrt((self.dx[1:-1,0])**2+\
                  (self.dy[1:-1,0]+self.dy[:-2,0])**2)
        dx[0,-1]     =0.5*numpy.sqrt((self.dx[0,-1])**2+(self.dy[0,-1])**2)
        dx[-1,0]     =0.5*numpy.sqrt((self.dx[-1,0])**2+(self.dy[-1,0])**2)
        dx[-1,1:-1]  =0.5*numpy.sqrt((self.dx[-1,1:-1]+self.dx[-1,:-2])**2+\
                  (self.dy[-1,1:-1])**2)
        dx[1:-1,-1]  =0.5*numpy.sqrt((self.dx[1:-1,-1])**2+(self.dy[1:-1,-1]+\
                  self.dy[:-2,-1])**2)
        dx[-1,-1]    =0.5*numpy.sqrt((self.dx[-1,-1])**2+(self.dy[-1,-1])**2)
#        print(dx)
        c=CP.PropsSI('A','T',300,'P',101325,self.Domain.fluid)
#        print(c)
        return numpy.amin(self.CFL*dx/(c))

    # Convergence checker (REMOVE? NO IMPLICIT CALCULATIONS DONE)
    def CheckConv(self, Tprev, Tnew):
        diff=numpy.sum(numpy.abs(Tnew[:]-Tprev[:]))/numpy.sum(numpy.abs(Tprev[:]))
        print(diff)
        if diff<=self.conv:
            return True
        else:
            return False
    # Solve
    """ To do:
        - flux terms calculator
        - source terms
        - RK time advancement (eventually)
        
    """
    # Spatial derivatives
    # Calculates for entire domain and accounts for periodicity
    def compute_derivative(self, u, v, dx, dy):
        ddx=numpy.empty_like(u)
        ddy=numpy.empty_like(v)
        
        ddx[:,1:-1]=(u[:,2:]-u[:,:-2])/(dx[:,1:-1]+dx[:,:-2])
        ddy[1:-1,:]=(v[2:,:]-v[:-2,:])/(dy[1:-1,:]+dy[:-2,:])
        
        if (self.BCs['bc_type_left']=='periodic') or (self.BCs['bc_type_right']=='periodic'):
            ddx[:,0] =(u[:,1]-u[:,-1])/(dx[:,0]+dx[:,-1])
            ddx[:,-1]=(u[:,0]-u[:,-2])/(dx[:,-1]+dx[:,0])
        else:
            # Forward/backward differences for boundaries
            ddx[:,0] =(u[:,1]-u[:,0])/(dx[:,0])
            ddx[:,-1]=(u[:,-1]-u[:,-2])/(dx[:,-1])
        if (self.BCs['bc_type_north']=='periodic') or (self.BCs['bc_type_south']=='periodic'):
            ddy[0,:] =(v[1,:]-v[-1,:])/(dy[0,:]+dy[-1,:])
            ddy[-1,:]=(v[0,:]-v[-2,:])/(dy[-1,:]+dy[0,:])
        else:
            # Forward/backward differences for boundaries
            ddy[0,:] =(v[1,:]-v[0,:])/(dy[0,:])
            ddy[-1,:]=(v[-1,:]-v[-2,:])/(dy[-1,:])
        
        return ddx, ddy
    
    # Bondary condition handler (not including periodic BCs)
    def Apply_BCs(self, u, v, p, isPress):
        # Start with wall BCs (implied 0 gradients and no slip)
        # Sort BCs by pressure call or not
        
        # Pressure BCs
        if isPress:
            # Left face
            if self.BCs['bc_type_left']=='wall':
                p[:,0]  =p[:,1]
                
            elif self.BCs['bc_type_left']=='inlet':
                p[:,0]  =self.BCs['bc_left_p']
                            
            elif self.BCs['bc_type_left']=='outlet':
                p[:,0]=self.BCs['bc_left_p']
            elif self.BCs['bc_type_left']!='periodic':
                if (type(self.BCs['bc_left_p']) is str)\
                    and (self.BCs['bc_left_p']=='zero_grad'):
                    p[:,0]  =p[:,1]
                else:
                    p[:,0]  =self.BCs['bc_left_p']
                
            # Right face
            if self.BCs['bc_type_right']=='wall':
                p[:,-1]  =p[:,-2]
                
            elif self.BCs['bc_type_right']=='inlet':
                p[:,-1]  =self.BCs['bc_right_p']
                
            elif self.BCs['bc_type_right']=='outlet':
                p[:,-1]=self.BCs['bc_right_p']
            elif self.BCs['bc_type_right']!='periodic':
                if (type(self.BCs['bc_right_p']) is str)\
                    and (self.BCs['bc_right_p']=='zero_grad'):
                    p[:,-1]  =p[:,-2]  
                else:
                    p[:,-1]  =self.BCs['bc_right_p']
                
            # South face
            if self.BCs['bc_type_south']=='wall':
                p[0,:]  =p[1,:]
                
            elif self.BCs['bc_type_south']=='inlet':
                p[0,:]  =self.BCs['bc_south_p']
                            
            elif self.BCs['bc_type_south']=='outlet':
                p[0,:]=self.BCs['bc_south_p']
            elif self.BCs['bc_type_south']!='periodic':
                if (type(self.BCs['bc_south_p']) is str)\
                    and (self.BCs['bc_south_p']=='zero_grad'):
                    p[0,:]  =p[1,:]  
                else:
                    p[0,:]  =self.BCs['bc_south_p']
                
            # North face
            if self.BCs['bc_type_north']=='wall':
                p[-1,:]  =p[-2,:]
                
            elif self.BCs['bc_type_north']=='inlet':
                p[-1,:]  =self.BCs['bc_north_p']
                
            elif self.BCs['bc_type_north']=='outlet':
                p[-1,:]=self.BCs['bc_north_p']
            elif self.BCs['bc_type_north']!='periodic':
                if (type(self.BCs['bc_north_p']) is str)\
                    and (self.BCs['bc_north_p']=='zero_grad'):
                    p[-1,:]  =p[-2,:]  
                else:
                    p[-1,:]  =self.BCs['bc_north_p']
            
        # Velocity field BCs
        else:
            # Left face
            if self.BCs['bc_type_left']=='wall':
                u[:,0]  =0
                v[:,0]  =0
            elif self.BCs['bc_type_left']=='inlet':
                u[:,0]  =self.BCs['bc_left_u']
                v[:,0]  =self.BCs['bc_left_v']
                
    #        elif self.BCs['bc_type_left']=='outlet':
    #            p[:,0]=self.BCs['bc_left_p']
            elif self.BCs['bc_type_left']!='periodic':
                if (type(self.BCs['bc_left_u']) is str)\
                    and (self.BCs['bc_left_u']=='zero_grad'):
                    u[:,0]  =u[:,1]
                else:
                    u[:,0]  =self.BCs['bc_left_u']
                if (type(self.BCs['bc_left_v']) is str)\
                    and (self.BCs['bc_left_v']=='zero_grad'):
                    v[:,0]  =v[:,1]
                else:
                    v[:,0]  =self.BCs['bc_left_v']    
                        
            # Right face
            if self.BCs['bc_type_right']=='wall':
                u[:,-1]  =0
                v[:,-1]  =0
            elif self.BCs['bc_type_right']=='inlet':
                u[:,-1]  =self.BCs['bc_right_u']
                v[:,-1]  =self.BCs['bc_right_v']
                
    #        elif self.BCs['bc_type_right']=='outlet':
    #            p[:,-1]=self.BCs['bc_right_p']
            elif self.BCs['bc_type_right']!='periodic':
                if (type(self.BCs['bc_right_u']) is str)\
                    and (self.BCs['bc_right_u']=='zero_grad'):
                    u[:,-1]  =u[:,-2]  
                else:
                    u[:,-1]  =self.BCs['bc_right_u']
                if (type(self.BCs['bc_right_v']) is str)\
                    and (self.BCs['bc_right_v']=='zero_grad'):
                    v[:,-1]  =v[:,-2]  
                else:
                    v[:,-1]  =self.BCs['bc_right_v']
                        
            # South face
            if self.BCs['bc_type_south']=='wall':
                u[0,:]  =0
                v[0,:]  =0
            elif self.BCs['bc_type_south']=='inlet':
                u[0,:]  =self.BCs['bc_south_u']
                v[0,:]  =self.BCs['bc_south_v']
                
    #        elif self.BCs['bc_type_south']=='outlet':
    #            p[0,:]=self.BCs['bc_south_p']
            elif self.BCs['bc_type_south']!='periodic':
                if (type(self.BCs['bc_south_u']) is str)\
                    and (self.BCs['bc_south_u']=='zero_grad'):
                    u[0,:]  =u[1,:]  
                else:
                    u[0,:]  =self.BCs['bc_south_u']
                if (type(self.BCs['bc_south_v']) is str)\
                    and (self.BCs['bc_south_v']=='zero_grad'):
                    v[0,:]  =v[1,:]  
                else:
                    v[0,:]  =self.BCs['bc_south_v']
                    
            # North face
            if self.BCs['bc_type_north']=='wall':
                u[-1,:]  =0
                v[-1,:]  =0
            elif self.BCs['bc_type_north']=='inlet':
                u[-1,:]  =self.BCs['bc_north_u']
                v[-1,:]  =self.BCs['bc_north_v']
                
    #        elif self.BCs['bc_type_north']=='outlet':
    #            p[-1,:]=self.BCs['bc_north_p']
            elif self.BCs['bc_type_north']!='periodic':
                if (type(self.BCs['bc_north_u']) is str)\
                    and (self.BCs['bc_north_u']=='zero_grad'):
                    u[-1,:]  =u[-2,:]  
                else:
                    u[-1,:]  =self.BCs['bc_north_u']
                if (type(self.BCs['bc_north_v']) is str)\
                    and (self.BCs['bc_north_v']=='zero_grad'):
                    v[-1,:]  =v[-2,:]  
                else:
                    v[-1,:]  =self.BCs['bc_north_v']
    
    # Calculate pressure field
    def compute_pressure(self,dt,u,v,dudx,dvdy,dvdx,dudy,d2udx,d2vdy,d2vdx,d2udy,dx,dy):
        rhs=numpy.zeros_like(u)
        p_c=self.Domain.p.copy()
        print(p_c)
        count=0
        
        rhs =-self.Domain.rho/dt*(dudx+dvdy)
        # Calculate higher/mixed derivatives
        if self.higher_p:
            rhs*=-1
            d3udx,d3vdy=self.compute_derivative(d2udx, d2vdy, self.dx, self.dy)
            d3udyyx,d3vdxxy=self.compute_derivative(d2udy, d2vdx, self.dx, self.dy)
            exp, dummy = self.compute_derivative(u*dvdy+v*dudy, numpy.zeros_like(v), self.dx, self.dy)
            rhs-=self.Domain.rho*(u*d2udx + v*d2vdy)
            rhs-=self.Domain.rho*exp
            rhs+=self.Domain.mu*(d3udx+d3vdy+d3udyyx+d3vdxxy)
        rhs-=self.Domain.rho*(dudx**2 + dvdy**2 + 2*dudy*dvdx)
        
        while (count<self.maxCount):
            # Central difference for bulk
            self.Domain.p[1:-1,1:-1]=((dx[1:-1,1:-1]**2*(p_c[2:,1:-1]+p_c[:-2,1:-1])\
                         +dy[1:-1,1:-1]**2*(p_c[1:-1,2:]+p_c[1:-1,:-2]))\
                         -(dy[1:-1,1:-1]*dx[1:-1,1:-1])**2*rhs[1:-1,1:-1])\
                        /2/(dy[1:-1,1:-1]**2+dx[1:-1,1:-1]**2)
            
            # Forward/backward difference at boundaries
            if self.BCs['bc_type_south']=='periodic' or self.BCs['bc_type_north']=='periodic':
                self.Domain.p[0,1:-1]=((dx[0,1:-1]**2*(p_c[1,1:-1]+p_c[-1,1:-1])\
                             +dy[0,1:-1]**2*(p_c[0,2:]+p_c[0,:-2]))\
                             -(dy[0,1:-1]*dx[0,1:-1])**2*rhs[0,1:-1])\
                            /2/(dy[0,1:-1]**2+dx[0,1:-1]**2)
                self.Domain.p[-1,1:-1]=((dx[-1,1:-1]**2*(p_c[0,1:-1]+p_c[-2,1:-1])\
                             +dy[-1,1:-1]**2*(p_c[-1,2:]+p_c[-1,:-2]))\
                             -(dy[-1,1:-1]*dx[-1,1:-1])**2*rhs[-1,1:-1])\
                            /2/(dy[-1,1:-1]**2+dx[-1,1:-1]**2)
            else:
                self.Domain.p[0,1:-1]=((dx[0,1:-1]**2*(p_c[2,1:-1]-2*p_c[1,1:-1])\
                             +dy[0,1:-1]**2*(p_c[0,2:]+p_c[0,:-2]))\
                             -(dy[0,1:-1]*dx[0,1:-1])**2*rhs[0,1:-1])\
                            /(2*dy[0,1:-1]**2-dx[0,1:-1]**2)
                self.Domain.p[-1,1:-1]=((dx[-1,1:-1]**2*(p_c[-3,1:-1]-2*p_c[-2,1:-1])\
                             +dy[-1,1:-1]**2*(p_c[-1,2:]+p_c[-1,:-2]))\
                             -(dy[-1,1:-1]*dx[-1,1:-1])**2*rhs[-1,1:-1])\
                            /(2*dy[-1,1:-1]**2-dx[-1,1:-1]**2)
            if self.BCs['bc_type_left']=='periodic' or self.BCs['bc_type_right']=='periodic':
                self.Domain.p[1:-1,0]=((dx[1:-1,0]**2*(p_c[2:,0]+p_c[:-2,0])\
                             +dy[1:-1,0]**2*(p_c[1:-1,1]+p_c[1:-1,-1]))\
                             -(dy[1:-1,0]*dx[1:-1,0])**2*rhs[1:-1,0])\
                            /2/(dy[1:-1,0]**2+dx[1:-1,0]**2)
                self.Domain.p[1:-1,-1]=((dx[1:-1,-1]**2*(p_c[2:,-1]+p_c[:-2,-1])\
                             +dy[1:-1,-1]**2*(p_c[1:-1,0]+p_c[1:-1,-2]))\
                             -(dy[1:-1,-1]*dx[1:-1,-1])**2*rhs[1:-1,-1])\
                            /2/(dy[1:-1,-1]**2+dx[1:-1,-1]**2)
            else:
                self.Domain.p[1:-1,0]=((dx[1:-1,0]**2*(p_c[2:,0]+p_c[:-2,0])\
                             +dy[1:-1,0]**2*(p_c[1:-1,2]-2*p_c[1:-1,1]))\
                             -(dy[1:-1,0]*dx[1:-1,0])**2*rhs[1:-1,0])\
                            /(-dy[1:-1,0]**2+2*dx[1:-1,0]**2)
                self.Domain.p[1:-1,-1]=((dx[1:-1,-1]**2*(p_c[2:,-1]+p_c[:-2,-1])\
                             +dy[1:-1,-1]**2*(p_c[1:-1,-3]-2*p_c[1:-1,-2]))\
                             -(dy[1:-1,-1]*dx[1:-1,-1])**2*rhs[1:-1,-1])\
                            /(-dy[1:-1,-1]**2+2*dx[1:-1,-1]**2)
            # Corner treatments
            self.Domain.p[0,0]=((dx[0,0]**2*(p_c[2,0]-2*p_c[1,0])\
                         +dy[0,0]**2*(p_c[0,2]-2*p_c[0,1]))\
                         -(dy[0,0]*dx[0,0])**2*rhs[0,0])\
                        /(2*dy[0,0]**2-dx[0,0]**2)
            self.Domain.p[-1,0]=((dx[-1,0]**2*(p_c[-3,0]-2*p_c[-2,0])\
                         +dy[-1,0]**2*(p_c[-1,2]-2*p_c[-1,1]))\
                         -(dy[-1,0]*dx[-1,0])**2*rhs[-1,0])\
                        /(2*dy[-1,0]**2-dx[-1,0]**2)
            self.Domain.p[0,-1]=((dx[0,-1]**2*(p_c[2,-1]-2*p_c[1,-1])\
                         +dy[0,-1]**2*(p_c[0,-3]-2*p_c[0,-2]))\
                         -(dy[0,-1]*dx[0,-1])**2*rhs[0,-1])\
                        /(2*dy[0,-1]**2-dx[0,-1]**2)
            self.Domain.p[-1,-1]=((dx[-1,-1]**2*(p_c[-3,-1]-2*p_c[-2,-1])\
                         +dy[-1,-1]**2*(p_c[-1,-3]-2*p_c[-1,-2]))\
                         -(dy[-1,-1]*dx[-1,-1])**2*rhs[-1,-1])\
                        /(2*dy[-1,-1]**2-dx[-1,-1]**2)
            
            self.Apply_BCs(u,v,self.Domain.p, True)
            if self.CheckConv(p_c, self.Domain.p):
                return 0
            p_c=self.Domain.p.copy()
            count+=1
        
        return -1
                    
    # Main compressible solver (1 time step)
    def Advance_Soln(self):
        u_0=self.Domain.u.copy()
        v_0=self.Domain.v.copy()
        u_c=u_0.copy()
        v_c=v_0.copy()
        mu=self.Domain.mu
        rho=self.Domain.rho
                
        if self.time_scheme=='Euler':
            rk_coeff = numpy.array([1,0])
            rk_substep_fraction = numpy.array([1,0])
            Nstep = 1
            dudt=[0]*Nstep
            dvdt=[0]*Nstep
            
        else:
            RK_info=temporal_schemes.runge_kutta(self.time_scheme)
            Nstep = RK_info.Nk
            if Nstep<0:
                return 1, -1 # Scheme not recognized; abort solver
            rk_coeff = RK_info.rk_coeff
            rk_substep_fraction = RK_info.rk_substep_fraction

            dudt=[0]*Nstep
            dvdt=[0]*Nstep
            
        dt=self.getdt()
        if (numpy.isnan(dt)) or (dt<=0):
            print '*********Diverging time step***********'
            return 1, dt
#        print 'Time step size: %f'%dt
        
        for step in range(Nstep):
            
            ###################################################################
            # Compute time deriviatives of conservatives (2nd order central schemes)
            ###################################################################
            # Calculate velocity spatial derivatives
            dudx,dvdy=self.compute_derivative(u_c,v_c,self.dx,self.dy)
            dvdx,dudy=self.compute_derivative(v_c, u_c, self.dx, self.dy)
            d2udx,d2vdy=self.compute_derivative(dudx, dvdy, self.dx,self.dy)
            d2vdx,d2udy=self.compute_derivative(dvdx, dudy, self.dx, self.dy)
            
            # Get pressure field
            err=self.compute_pressure(dt,u_c,v_c,dudx,dvdy,dvdx,dudy,d2udx,d2vdy,d2vdx,d2udy,self.dx, self.dy)
            if err==-1:
                print '********Pressure field problem**************'
                return 1, dt
            # x-momentum (pressure, flux, shear stress, gravity, pressure gradient)
            dudt[step], dummy =self.compute_derivative(self.Domain.p, numpy.zeros_like(v_c), self.dx, self.dy)
            dudt[step]       *=-1
            dudt[step]       -=u_c*dudx+v_c*dudy
            dudt[step]       +=mu/rho*(d2udx+d2udy)
            dudt[step]       +=self.gx
            dudt[step]       -=self.dpx/rho
    
            # y-momentum (pressure, flux, shear stress, gravity, pressure gradient)
            dummy, dvdt[step] =self.compute_derivative(numpy.zeros_like(u_c), self.Domain.p, self.dx, self.dy)
            dvdt[step]       *=-1
            dvdt[step]       -=u_c*dvdx+v_c*dvdy
            dvdt[step]       +=mu/rho*(d2vdx+d2vdy)
            dvdt[step]       +=self.gy
            dvdt[step]       -=self.dpy/rho
            
            # Compute intermediate conservative values for RK stepping
            u_c=u_0.copy()
            v_c=v_0.copy()
            
            if step < (Nstep - 1):
                for rk_index in range(step + 1):
                    
                    u_c+= dt*rk_coeff[step+1][rk_index]*dudt[rk_index]
                    v_c+= dt*rk_coeff[step+1][rk_index]*dvdt[rk_index]
                    
                ###################################################################
                # Apply boundary conditions
                ###################################################################
                self.Apply_BCs(u_c, v_c, self.Domain.p, False)
            
            ###################################################################
            # END OF TIME STEP CALCULATIONS
            ###################################################################
            
        ###################################################################
        # Compute new conservative values at new time step
        ###################################################################
        for step in range(Nstep):    
            self.Domain.u+= dt * rk_substep_fraction[step] * dudt[step]
            self.Domain.v+= dt * rk_substep_fraction[step] * dvdt[step]
            
        ###################################################################
        # Apply boundary conditions
        ###################################################################
        self.Apply_BCs(self.Domain.u, self.Domain.v, self.Domain.p, False)
        
        ###################################################################
        # Divergence check
        ###################################################################
        
        if (numpy.isnan(numpy.amax(self.Domain.u))) or \
            (numpy.isnan(numpy.amax(self.Domain.v))):
            print '**************Divergence detected****************'
            return 1, dt
        
        ###################################################################
        # Output data to file?????
        ###################################################################
        
        
        
        return 0, dt