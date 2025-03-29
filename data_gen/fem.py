import sys         
from ngsolve import *
#from ngsolve.webgui import Draw
#import netgen.gui
from netgen.geom2d import SplineGeometry
import numpy as np

# FEM Functions to generate far field data
def farf2d(u,mesh,farp,kappa,porder):
    ntheta=farp["n"]
    ctheta=farp["cent"]
    apptheta=farp["app"]
    nv=specialcf.normal(mesh.dim)
    theta=np.zeros(ntheta)
    uinf=np.zeros(ntheta,dtype=complex)
    fesa=H1(mesh, order=porder, complex=True, definedon=mesh.Materials("air"))
    #print(Integrate(CoefficientFunction(nv[0]*(x-x0[0])+nv[1]*(y-x0[1])),mesh,order=porder+1,definedon=mesh.Boundaries("scatterer")))
    for jp in range(0,ntheta):
        theta[jp]=(ctheta-apptheta/2)+apptheta*jp/(ntheta-1)
        # phi[jp]=2*np.pi*jp/(nphi-1)  ## Special for testing
        xhat=CoefficientFunction((np.cos(theta[jp]),np.sin(theta[jp])))
        Eout = exp(-1J*kappa*(x*xhat[0]+y*xhat[1]))
        func1=CoefficientFunction(-1j*kappa*(xhat*nv)*Eout * u)
        uinf1=Integrate(tuple(func1),mesh,order=porder+1,definedon=
                            mesh.Boundaries("scatterer"))
        vv=GridFunction(fesa)
        vv.vec[:]=0
        vv.Set(Eout,BND,definedon=mesh.Boundaries("scatterer"))
        fvv=CoefficientFunction(grad(vv)*grad(u)-kappa*kappa*vv*u)
        uinf2=Integrate(tuple(fvv),mesh,order=porder+1,definedon=mesh.Materials("air"))
        uinf[jp]=exp(1J*np.pi/4)/np.sqrt(8*np.pi*kappa)*(uinf1+uinf2)
    return(uinf,theta)

def helmsol(mesh,porder,ncoef,kappa,incp,farp):
    fes = H1(mesh, order=porder, complex=True)
    u = fes.TrialFunction()
    v = fes.TestFunction()
    a = BilinearForm(fes)
    a += SymbolicBFI(grad(u)*grad(v) - kappa**2*ncoef*u*v)
    a += SymbolicBFI(-1j*kappa*u*v,definedon=mesh.Boundaries("outerbnd"))
    print('Number of DoFs: ',fes.ndof)
    gfu = GridFunction(fes)
    #Draw(gfu,mesh,'us')
    with TaskManager():
        a.Assemble()
        Ainv=a.mat.Inverse()
    uinf=np.zeros((farp["n"],incp["n"]),dtype=complex)
    phi=np.zeros(incp["n"]);
    center=incp["cent"]
    app=incp["app"]
    for ip in range(0,incp["n"]):
        if ip%10==0:
            print("Done ip = ", ip,' of ',incp["n"])
        if incp["n"]==1:
            phi[0]=0.0
        else:
            phi[ip]=(center-app/2)+app*ip/(incp["n"]-1)
        d=[np.cos(phi[ip]),np.sin(phi[ip])]
        with TaskManager():
            b = LinearForm(fes)
            ui=exp(1J*kappa*(d[0]*x+d[1]*y))
            b += SymbolicLFI(kappa*kappa*(ncoef-1)*ui * v)
            b.Assemble()
            gfu.vec.data =  Ainv * b.vec
            Redraw()
        uinf[:,ip],theta=farf2d(gfu,mesh,farp,kappa,porder)
    return(uinf,phi,theta)
    
# Discretize the born operator
def discretize_born(k, xlim, phi, Ngrid, theta):
    vert_step = 2 * xlim / Ngrid  # Vertical discretization step size
    hor_step = 2 * xlim / Ngrid  # Horizontal discretization step size
    Cfac = vert_step * hor_step * np.exp(1j * np.pi / 4) * np.sqrt(k**3 / (np.pi * 8))  # Compute constant factor
    y1 = np.linspace(-xlim, xlim, Ngrid)  # Discretize horizontal domain
    y2 = np.linspace(-xlim, xlim, Ngrid)  # Discretize vertical domain
    Y1, Y2 = np.meshgrid(y1, y2)  # Produce matrices of coordinates on square domain
    grid_points = np.column_stack((Y1.ravel(), Y2.ravel()))  # Combine matrices into N^2 x 2 (matrix of coordinate pairs)
    xhat = np.array([np.cos(theta), np.sin(theta)]).T # Measurement angle values to evaluate far field at
    d = np.array([np.cos(phi), np.sin(phi)]) # Incident wave
    diff = xhat - d
    dot_products = np.dot(diff, grid_points.T)
    Exp = np.exp(1j * k * dot_products)
    A = Cfac * Exp # Discretize born operator
    return A
