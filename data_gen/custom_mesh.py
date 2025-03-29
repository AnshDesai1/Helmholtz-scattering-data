import sys         
from ngsolve import *
#from ngsolve.webgui import Draw
#import netgen.gui
from netgen.geom2d import SplineGeometry
import numpy as np

def add_ellipse(geo,R1,R2,xcen,ang,ind,ind_out):
    # Here the curve defines an ellipse.  note t isin [0,1]
    mat=np.array([[np.cos(ang),np.sin(ang)],[-np.sin(ang),np.cos(ang)]])
    print('ind',ind,'ang=',ang,'x0',xcen,'R1',R1,'R2',R2)
    Curve= lambda t: (
        xcen[0]+mat[0,0]*R1*cos(t*2*np.pi)+mat[0,1]*R2*sin(t*2*np.pi),
        xcen[1]+mat[1,0]*R1*cos(t*2*np.pi)+mat[1,1]*R2*sin(t*2*np.pi))
    if ind_out==1:
        geo.AddCurve(Curve,leftdomain=ind,rightdomain=ind_out,bc="scatterer")
    else:
        geo.AddCurve(Curve,leftdomain=ind,rightdomain=ind_out)
    return(geo)

def circle(hmax_s=0.3,hmax_a=0.3,pml_rad=2,R=1,pml_delta=.8,order=3,x0=(0,0),
           circle_data={'numcircle':3,
                             'R1':0.2,'xcen1':(0.3,0.2),'ind1':3,
                             'R2':0.2,'xcen2':(-0.5,-0.2),'ind2':4,
                             'R3':0.2,'xcen3':(0.4,-0.4),'ind3':5}):
    geo = SplineGeometry()
    geo.AddCircle( (0,0), pml_rad+pml_delta, leftdomain=2, bc="outerbnd")
    geo.AddCircle( (0,0), pml_rad, leftdomain=1, rightdomain=2, bc="innerbnd")
    numcircle = circle_data['numcircle']
    for i in range(numcircle):
        R = circle_data[f'R{i+1}']
        xcen = circle_data[f'xcen{i+1}']
        ind = circle_data[f'ind{i+1}']
        geo.AddCircle(xcen,R,leftdomain=ind,rightdomain=1, bc="scatterer")
        print(ind)
        geo.SetMaterial(ind,f"circle{i+1}")
        geo.SetDomainMaxH(ind,hmax_s)
        print(f"Circle {i+1} Done")
    geo.SetMaterial(1, "air")
    geo.SetMaterial(2, "pmlregion")
    geo.SetDomainMaxH(2,hmax_a)
    geo.SetDomainMaxH(1,hmax_a)
    mesh = Mesh(geo.GenerateMesh ())
    mesh.Curve(order)
    return(mesh)
    
def ellipses(hmax_a=0.1,hmax_s=.8,pml_delta=.5,pml_rad=1.8,order=3,
                 ellip_data={'numellip':3,
                             'R1a':0.5,'R1b':.2,'xcen1':(0.3,0.2),'ang1':np.pi/8,'ind1':3,
                             'R2a':0.5,'R2b':.2,'xcen2':(-0.5,-0.2),'ang2':np.pi/8,'ind2':4,
                             'R3a':0.5,'R3b':.2,'xcen3':(0.4,-0.4),'ang3':-np.pi/8,'ind3':5}):
    geo = SplineGeometry()
    geo.AddCircle( (0,0), pml_rad+pml_delta, leftdomain=2, bc="outerbnd")
    geo.AddCircle( (0,0), pml_rad, leftdomain=1, rightdomain=2)
    # first ellipse
    ang=ellip_data["ang1"]
    R1=ellip_data["R1a"]
    R2=ellip_data["R1b"]
    xcen=ellip_data["xcen1"]
    ind=ellip_data["ind1"]
    geo=add_ellipse(geo,R1,R2,xcen,ang,ind,1)
    geo.SetMaterial(ind,"ellip1")
    geo.SetDomainMaxH(ind,hmax_s)
    if ellip_data["numellip"] > 1:
    # second ellipse
        ang=ellip_data["ang2"]
        R1=ellip_data["R2a"]
        R2=ellip_data["R2b"]
        xcen=ellip_data["xcen2"]
        ind=ellip_data["ind2"]
        geo=add_ellipse(geo,R1,R2,xcen,ang,ind,1)
        geo.SetMaterial(ind,"ellip2")
        geo.SetDomainMaxH(ind,hmax_s)
    if ellip_data["numellip"] > 2:
    # third ellipse
        ang=ellip_data["ang3"]
        R1=ellip_data["R3a"]
        R2=ellip_data["R3b"]
        xcen=ellip_data["xcen3"]
        ind=ellip_data["ind3"]
        geo=add_ellipse(geo,R1,R2,xcen,ang,ind,1)
        geo.SetMaterial(ind,"ellip3")
        geo.SetDomainMaxH(5,hmax_s)
    geo.SetMaterial(1, "air")
    geo.SetMaterial(2, "pmlregion")
    geo.SetDomainMaxH(2,hmax_a)
    geo.SetDomainMaxH(1,hmax_a)
    mesh = Mesh(geo.GenerateMesh())
    mesh.Curve(order)
    return(mesh)

def ellipses_in_circle(R_circ=1.5,hmax_a=0.2,hmax_s=.1,pml_delta=.5, pml_rad=1.8,order=3,
                ellip_data={'numellip':3,
                            'R1a':0.5,'R1b':.2,'xcen1':(0.3,0.2),'ang1':np.pi/8,'ind1':3,
                            'R2a':0.5,'R2b':.2,'xcen2':(-0.5,-0.2),'ang2':np.pi/8,'ind2':4,
                            'R3a':0.5,'R3b':.2,'xcen3':(0.4,-0.4),'ang3':-np.pi/8,'ind3':5}):
    geo = SplineGeometry()
    geo.AddCircle( (0,0), pml_rad+pml_delta, leftdomain=2, bc="outerbnd")
    geo.AddCircle( (0,0), pml_rad, leftdomain=1, rightdomain=2)
    geo.AddCircle( (0,0), R_circ, leftdomain=6, rightdomain=1, bc="scatterer")
    # first ellipse
    ang=ellip_data["ang1"]
    R1=ellip_data["R1a"]
    R2=ellip_data["R1b"]
    xcen=ellip_data["xcen1"]
    ind=ellip_data["ind1"]
    geo=add_ellipse(geo,R1,R2,xcen,ang,ind,6)
    geo.SetMaterial(ind,"ellip1")
    geo.SetDomainMaxH(ind,hmax_s)
    if ellip_data["numellip"] > 1:
    # second ellipse
        ang=ellip_data["ang2"]
        R1=ellip_data["R2a"]
        R2=ellip_data["R2b"]
        xcen=ellip_data["xcen2"]
        ind=ellip_data["ind2"]
        geo=add_ellipse(geo,R1,R2,xcen,ang,ind,6)
        geo.SetMaterial(ind,"ellip2")
        geo.SetDomainMaxH(ind,hmax_s)
    if ellip_data["numellip"] > 2:
    # third ellipse
        ang=ellip_data["ang3"]
        R1=ellip_data["R3a"]
        R2=ellip_data["R3b"]
        xcen=ellip_data["xcen3"]
        ind=ellip_data["ind3"]
        geo=add_ellipse(geo,R1,R2,xcen,ang,ind,6)
        geo.SetMaterial(ind,"ellip3")
        geo.SetDomainMaxH(5,hmax_s)
    geo.SetMaterial(1, "air")
    geo.SetMaterial(2, "pmlregion")
    geo.SetMaterial(6,"big_circle")
    geo.SetDomainMaxH(6,hmax_s)
    geo.SetDomainMaxH(2,hmax_a)
    geo.SetDomainMaxH(1,hmax_a)
    mesh = Mesh(geo.GenerateMesh())
    mesh.Curve(order)
    return mesh

def random_scatterer(hmax=0.1,pml_rad=2.,pml_width=.8,xcen=(0,0), 
                     Rmax=1, Rmin=.2, 
                     npoints=10,
                     angs=2*np.pi*np.sort(np.random.rand(10)),
                     mags=np.random.rand(10),
                     order=3
                    ):
    geo = SplineGeometry()
    geo.AddCircle( (0,0), pml_rad+pml_width, leftdomain=2, bc="outerbnd")
    geo.AddCircle( (0,0), pml_rad, leftdomain=1, rightdomain=2)
    print('Number of points= ',npoints)
    Rads=(Rmax-Rmin)*mags+Rmin
    xav=0
    yav=0
    for j in range(0,npoints):
        xav=xav+Rads[j]*np.cos(angs[j])
        yav=yav+Rads[j]*np.sin(angs[j])
    # used to place origin inside the scatterer
    xav=xav/npoints
    yav=yav/npoints
    for j in range(0,npoints):
        x=Rads[j]*np.cos(angs[j])-xav
        y=Rads[j]*np.sin(angs[j])-yav
        Rads[j]=np.sqrt(x*x+y*y)
        angs[j]=np.arctan2(x,y)
    Ra=np.zeros((npoints,2))
    Ra[:,0]=angs
    Ra[:,1]=Rads
    Ras=np.sort(Ra,0)
    Rads=Ras[:,1]
    angs=Ras[:,0]
    print('Rads = ',Rads)
    print('angs = ',angs)
    for j in range(0,npoints-1):
        x1=Rads[j]*np.cos(angs[j])+xcen[0]
        y1=Rads[j]*np.sin(angs[j])+xcen[1]
        x2=Rads[j+1]*np.cos(angs[j+1])+xcen[0]
        y2=Rads[j+1]*np.sin(angs[j+1])+xcen[1]
        p1=geo.AppendPoint(x1,y1)
        p2=geo.AppendPoint(x2,y2)
        geo.Append (["line", p1, p2],leftdomain=3,rightdomain=1,bc="scatterer",
                    maxh=hmax/2)
    x1=Rads[npoints-1]*np.cos(angs[npoints-1])+xcen[0]
    y1=Rads[npoints-1]*np.sin(angs[npoints-1])+xcen[1]
    x2=Rads[0]*np.cos(angs[0])+xcen[0]
    y2=Rads[0]*np.sin(angs[0])+xcen[1]
    p1=geo.AppendPoint(x1,y1)
    p2=geo.AppendPoint(x2,y2)
    geo.Append (["line", p1, p2],leftdomain=3,rightdomain=1,bc="scatterer",
                    maxh=hmax/2)
    geo.SetMaterial(1, "air")
    geo.SetMaterial(2, "pmlregion")
    geo.SetMaterial(3, "D")
    mesh = Mesh(geo.GenerateMesh(maxh=hmax))
    mesh.Curve(order)
    return(mesh,angs,Rads)

def Lshape(hmax=0.1,pml_rad=2.,pml_width=.8,L2=np.sqrt(np.pi)/2,xcen=(0,0),
               ang=0, order=3):
    geo = SplineGeometry()
    geo.AddCircle( (0,0), pml_rad+pml_width, leftdomain=2, bc="outerbnd")
    geo.AddCircle( (0,0), pml_rad, leftdomain=1, rightdomain=2)
    #L2=np.sqrt(np.pi)/2 # for equal area
    mat=np.array([[np.cos(ang),np.sin(ang)],[-np.sin(ang),np.cos(ang)]])
    X1=np.matmul(mat,np.array([-L2,-L2]))+xcen
    X2=np.matmul(mat,np.array([-L2,L2]))+xcen
    X3=np.matmul(mat,np.array([L2,L2]))+xcen
    X4=np.matmul(mat,np.array([L2,0]))+xcen
    X5=np.matmul(mat,np.array([0,0]))+xcen
    X6=np.matmul(mat,np.array([0,-L2]))+xcen
    p1,p2,p3,p4,p5,p6 = [ geo.AppendPoint(x,y) for x,y in
                    [ X1, X2,X3, X4,X5,X6]]
    geo.Append (["line", p2, p1],leftdomain=3,rightdomain=1,bc="scatterer",
                    maxh=hmax/2)
    geo.Append (["line", p3, p2],leftdomain=3,rightdomain=1,bc="scatterer",
                    maxh=hmax/2)
    geo.Append (["line", p4, p3],leftdomain=3,rightdomain=1,bc="scatterer",
                    maxh=hmax/2)
    geo.Append (["line", p5, p4],leftdomain=3,rightdomain=1,bc="scatterer",
                    maxh=hmax/2)
    geo.Append (["line", p6, p5],leftdomain=3,rightdomain=1,bc="scatterer",
                    maxh=hmax/2)
    geo.Append (["line", p1, p6],leftdomain=3,rightdomain=1,bc="scatterer",
                    maxh=hmax/2)
    geo.SetMaterial(1, "air")
    geo.SetMaterial(2, "pmlregion")
    geo.SetMaterial(3, "D")
    mesh = Mesh(geo.GenerateMesh(maxh=hmax))
    mesh.Curve(order)
    return(mesh)

def Ushape(hmax=0.1,pml_rad=2.,pml_width=.8,Box=[-.5,.3,-.3,.3],t=0.1, xcen=(0.5,0.5),
               ang=np.pi/4, order=3):
    geo = SplineGeometry()
    geo.AddCircle( (0,0), pml_rad+pml_width, leftdomain=2, bc="outerbnd")
    geo.AddCircle( (0,0), pml_rad, leftdomain=1, rightdomain=2)
    #
    mat=np.array([[np.cos(ang),np.sin(ang)],[-np.sin(ang),np.cos(ang)]])
    Lmin=Box[0]
    Lmax=Box[1]
    Hmin=Box[2]
    Hmax=Box[3]
    print(Lmin,Lmax,Hmin,Hmax)
    X1=np.matmul(mat,np.array([Lmin,Hmin]))+xcen
    X2=np.matmul(mat,np.array([Lmax,Hmin]))+xcen
    X3=np.matmul(mat,np.array([Lmax,Hmin+t]))+xcen
    X4=np.matmul(mat,np.array([Lmin+t,Hmin+t]))+xcen
    X5=np.matmul(mat,np.array([Lmin+t,Hmax-t]))+xcen
    X6=np.matmul(mat,np.array([Lmax,Hmax-t]))+xcen
    X7=np.matmul(mat,np.array([Lmax,Hmax]))+xcen
    X8=np.matmul(mat,np.array([Lmin,Hmax]))+xcen
    p1,p2,p3,p4,p5,p6,p7,p8 = [ geo.AppendPoint(x,y) for x,y in
                    [ X1, X2,X3, X4,X5,X6,X7,X8]]
    geo.Append (["line", p1, p2],leftdomain=3,rightdomain=1,bc="scatterer",
                    maxh=hmax/2)
    geo.Append (["line", p2, p3],leftdomain=3,rightdomain=1,bc="scatterer",
                    maxh=hmax/2)
    geo.Append (["line", p3, p4],leftdomain=3,rightdomain=1,bc="scatterer",
                    maxh=hmax/2)
    geo.Append (["line", p4, p5],leftdomain=3,rightdomain=1,bc="scatterer",
                    maxh=hmax/2)
    geo.Append (["line", p5, p6],leftdomain=3,rightdomain=1,bc="scatterer",
                    maxh=hmax/2)
    geo.Append (["line", p6, p7],leftdomain=3,rightdomain=1,bc="scatterer",
                    maxh=hmax/2)
    geo.Append (["line", p7, p8],leftdomain=3,rightdomain=1,bc="scatterer",
                    maxh=hmax/2)
    geo.Append (["line", p8, p1],leftdomain=3,rightdomain=1,bc="scatterer",
                    maxh=hmax/2)
    geo.SetMaterial(1, "air")
    geo.SetMaterial(2, "pmlregion")
    geo.SetMaterial(3, "D")
    mesh = Mesh(geo.GenerateMesh(maxh=hmax))
    mesh.Curve(order)
    return(mesh)
    
def Venn(hmax=0.1,pml_rad=3.,pml_width=.8,x1=0.8,R1=1.1,order=3):
# intersecting circles centered at (-x1,R1) and (x1,R1) and radius R1
    geo = SplineGeometry()
    geo.AddCircle( (0,0), pml_rad+pml_width, leftdomain=2, bc="outerbnd")
    geo.AddCircle( (0,0), pml_rad, leftdomain=1, rightdomain=2)
    #
    p1=geo.AppendPoint(0,np.sqrt(R1**2-x1**2))
    print(x1+R1,-x1+R1)
    alpha=(x1+R1)/np.sqrt(R1**2-x1**2)
    p1a=geo.AppendPoint(x1+R1,np.sqrt(R1**2-x1**2)+alpha*x1)
    p1aR=geo.AppendPoint(-(x1+R1),np.sqrt(R1**2-x1**2)+alpha*x1)
    p2=geo.AppendPoint(x1+R1,0)
    p2R=geo.AppendPoint(-(x1+R1),0)
    p2a=geo.AppendPoint(x1+R1,-(np.sqrt(R1**2-x1**2)+alpha*x1))
    p2aR=geo.AppendPoint(-(x1+R1),-(np.sqrt(R1**2-x1**2)+alpha*x1))
    p3=geo.AppendPoint(0,-np.sqrt(R1**2-x1**2))
    beta=(-x1+R1)/np.sqrt(R1**2-x1**2)
    p3a=geo.AppendPoint(-x1+R1,np.sqrt(R1**2-x1**2)-beta*x1)
    p3aR=geo.AppendPoint(x1-R1,np.sqrt(R1**2-x1**2)-beta*x1)
    print(-x1+R1,np.sqrt(R1**2-x1**2)-beta*x1)
    p4=geo.AppendPoint(-x1+R1,0)
    p4a=geo.AppendPoint(-x1+R1,-(np.sqrt(R1**2-x1**2)-beta*x1))
    p4aR=geo.AppendPoint(x1-R1,-(np.sqrt(R1**2-x1**2)-beta*x1))
    p4R=geo.AppendPoint(x1-R1,0)
    geo.Append(["spline3",p2,p1a,p1],leftdomain=3,rightdomain=1,bc="scatterer")
    geo.Append(["spline3",p3,p2a,p2],leftdomain=3,rightdomain=1,bc="scatterer")
    geo.Append(["spline3",p1,p3a,p4],leftdomain=3,rightdomain=4)
    geo.Append(["spline3",p4,p4a,p3],leftdomain=3,rightdomain=4)
    geo.Append(["spline3",p3,p4aR,p4R],leftdomain=5,rightdomain=4)
    geo.Append(["spline3",p4R,p3aR,p1],leftdomain=5,rightdomain=4)
    geo.Append(["spline3",p1,p1aR,p2R],leftdomain=5,rightdomain=1,bc="scatterer")
    geo.Append(["spline3",p2R,p2aR,p3],leftdomain=5,rightdomain=1,bc="scatterer")
    geo.SetMaterial(1, "air")
    geo.SetMaterial(2, "pmlregion")
    geo.SetMaterial(3, "DL")
    geo.SetMaterial(4, "DLR")
    geo.SetMaterial(5, "DR")
    mesh=Mesh(geo.GenerateMesh(maxh=hmax))
    mesh.Curve(order)
    return(mesh)

def Shepp_Logan(hmax=0.02,pml_rad=1.,pml_width=.8, order=3, data=[
        ((0,0), 0, .69, .92, 3, 1),
        ((0,-0.0184), 0, .6624, .874, 4, 3),
        ((0.28,0), 18./360.*2.*np.pi, .11, .31, 5, 4),
        (((-0.22,0)), -18./360.*2.*np.pi, .16, .41, 6, 4),
        ((0,0.45), 0, .21, .25, 7, 4),
        ((0,0.1), 0, 0.046, 0.046, 8, 4),
        ((0,-0.1+.1), 0, 0.046, 0.046, 9, 4),
        ((-0.08,-0.605), 0, .046, .023, 10, 4),
        (((0.0,-0.605)), 0, .023, .023, 11, 4),
        ((0.06,-0.605), 0, .023, .046, 12, 4)
    ]):
    geo = SplineGeometry()
    geo.AddCircle( (0,0), pml_rad+pml_width, leftdomain=2, bc="outerbnd")
    geo.AddCircle( (0,0), pml_rad, leftdomain=1, rightdomain=2)
    for (xcen, ang, R1, R2, ind, ind_out) in data:
        add_ellipse(geo, R1, R2, xcen, ang, ind, ind_out)
    geo.SetMaterial(1, "air")
    geo.SetMaterial(2, "pmlregion")
    geo.SetMaterial(3, "a")
    geo.SetMaterial(4, "b")
    geo.SetMaterial(5, "c")
    geo.SetMaterial(6, "d")
    geo.SetMaterial(7, "e")
    geo.SetMaterial(8, "f")
    geo.SetMaterial(9, "g")
    geo.SetMaterial(10, "h")
    geo.SetMaterial(11, "i")
    geo.SetMaterial(12, "j")
    mesh=Mesh(geo.GenerateMesh(maxh=hmax))
    mesh.Curve(order)
    return(mesh)

def Bullseye(hmax=0.1,pml_rad=2.,pml_width=.8, order=3):
    geo = SplineGeometry()
    geo.AddCircle( (0,0), pml_rad+pml_width, leftdomain=2, bc="outerbnd")
    geo.AddCircle( (0,0), pml_rad, leftdomain=1, rightdomain=2)
    # a
    data = [((0,0), 0, 1., 1., 3, 1),
            ((0,0), 0, .8, .8, 4, 3), 
            ((0,-0.0184), 0, .6, .4, 5, 4)
           ]
    for (xcen, ang, R1, R2, ind, ind_out) in data:
        add_ellipse(geo,R1,R2,xcen,ang,ind,ind_out)
    geo.SetMaterial(1, "air")
    geo.SetMaterial(2, "pmlregion")
    geo.SetMaterial(3, "a")
    geo.SetMaterial(4, "b")
    geo.SetMaterial(5, "c")
    mesh=Mesh(geo.GenerateMesh(maxh=hmax))
    mesh.Curve(order)
    return(mesh)

def Squares(hmax=0.1,pml_rad=2.,pml_width=.8):
    geo = SplineGeometry()
    geo.AddCircle( (0,0), pml_rad+pml_width, leftdomain=2, bc="outerbnd")
    geo.AddCircle( (0,0), pml_rad, leftdomain=1, rightdomain=2)
    # points
    # need x1<x2<x3
    x1=-.4
    x2=.4
    x3=.6
    # need y1<y2<y3<y4
    y1=-.4
    y2=-.2
    y3=.1
    y4=.4
    p1=geo.AppendPoint(0,y2)
    p2=geo.AppendPoint(0,y3)
    p3=geo.AppendPoint(x2,y3)
    p4=geo.AppendPoint(x2,y2)
    # a
    geo.Append(["line",p1,p2],leftdomain=3,rightdomain=4)
    geo.Append(["line",p2,p3],leftdomain=3,rightdomain=4)
    geo.Append(["line",p3,p4],leftdomain=5,rightdomain=4)
    geo.Append(["line",p4,p1],leftdomain=3,rightdomain=4)
    #
    # b
    p5=geo.AppendPoint(x3,y3)
    p6=geo.AppendPoint(x3,y2)
    geo.Append(["line",p5,p3],leftdomain=5,rightdomain=1,bc="scatterer")
    geo.Append(["line",p6,p5],leftdomain=5,rightdomain=1,bc="scatterer")
    geo.Append(["line",p4,p6],leftdomain=5,rightdomain=1,bc="scatterer")
    #c
    p7=geo.AppendPoint(x2,y1)
    p8=geo.AppendPoint(x2,y4)
    p9=geo.AppendPoint(x1,y4)
    p10=geo.AppendPoint(x1,y1)
    geo.Append(["line",p3,p8],leftdomain=3,rightdomain=1,bc="scatterer")
    geo.Append(["line",p8,p9],leftdomain=3,rightdomain=1,bc="scatterer")
    geo.Append(["line",p9,p10],leftdomain=3,rightdomain=1,bc="scatterer")
    geo.Append(["line",p10,p7],leftdomain=3,rightdomain=1,bc="scatterer")
    geo.Append(["line",p7,p4],leftdomain=3,rightdomain=1,bc="scatterer")
    #
    geo.SetMaterial(1, "air")
    geo.SetMaterial(2, "pmlregion")
    geo.SetMaterial(3, "a")
    geo.SetMaterial(4, "b")
    geo.SetMaterial(5, "c")
    mesh=Mesh(geo.GenerateMesh(maxh=hmax))
    return(mesh)
