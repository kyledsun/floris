#--- Legacy python 2.7
from __future__ import division
from __future__ import print_function
# --- General
import unittest
import numpy as np
# --- Local
from .VortexCylinder import vc_tang_u, vc_longi_u, vc_root_u, vcs_tang_u, vcs_longi_u, vc_tang_u_doublet
from .VortexDoublet  import doublet_line_u
from .SelfSimilar    import ss_u
from .VortexCylinderSkewed import svc_tang_u, svc_longi_u, svc_root_u, svcs_tang_u, svcs_longi_u
from .Solver import Ct_const_cutoff, WakeVorticityFromCt, WakeVorticityFromGamma

# --------------------------------------------------------------------------------}
# --- Helper functions for geometry 
# --------------------------------------------------------------------------------{
def RotMat_AxisAngle(u,theta):
    """ Returns the rotation matrix for a rotation around an axis u, with an angle theta """
    R=np.zeros((3,3))
    ux,uy,uz=u
    c,s=np.cos(theta),np.sin(theta)
    R[0,0]=ux**2*(1-c)+c    ; R[0,1]=ux*uy*(1-c)-uz*s; R[0,2]=ux*uz*(1-c)+uy*s;
    R[1,0]=uy*ux*(1-c)+uz*s ; R[1,1]=uy**2*(1-c)+c   ; R[1,2]=uy*uz*(1-c)-ux*s
    R[2,0]=uz*ux*(1-c)-uy*s ; R[2,1]=uz*uy*(1-c)+ux*s; R[2,2]=uz**2*(1-c)+c;
    return R

def orth_vect(u):
    """ Given one vector, returns a 3x1 orthonormal vector to it"""
    u=(u/np.linalg.norm(u)).ravel()
    if abs(u[0])>=1/np.sqrt(3):
        v = np.array([[-u[1]],[u[0]] ,[0]]    )
    elif abs(u[1])>=1/np.sqrt(3):
        v = np.array([[0]    ,[-u[2]],[u[1]]] )
    elif abs(u[2])>=1/np.sqrt(3):
        v = np.array([[u[2]] ,[0]    ,[-u[0]]])
    else:
        raise Exception('Cannot find orthogonal vector to a zero vector: {} '.format(u))
    return v/np.linalg.norm(v)

def transform_T(T_a2b,Xb,Yb,Zb):
    Xa=T_a2b[0,0]*Xb+T_a2b[1,0]*Yb+T_a2b[2,0]*Zb
    Ya=T_a2b[0,1]*Xb+T_a2b[1,1]*Yb+T_a2b[2,1]*Zb
    Za=T_a2b[0,2]*Xb+T_a2b[1,2]*Yb+T_a2b[2,2]*Zb
    return Xa,Ya,Za
def transform(T_a2b,Xa,Ya,Za):
    Xb=T_a2b[0,0]*Xa+T_a2b[0,1]*Ya+T_a2b[0,2]*Za
    Yb=T_a2b[1,0]*Xa+T_a2b[1,1]*Ya+T_a2b[1,2]*Za
    Zb=T_a2b[2,0]*Xa+T_a2b[2,1]*Ya+T_a2b[2,2]*Za
    return Xb,Yb,Zb


# --------------------------------------------------------------------------------}
# --- Main Class 
# --------------------------------------------------------------------------------{
class WindTurbine:
    def __init__(self,R,r_hub=[0,0,0],e_shaft_yaw0=[0,0,1],e_vert=[0,1,0],U0=[0,0,10],Ct=None,Lambda=np.inf,name='',Ground=False,Model='VC'):
        """ 
        INPUTS
         - R: rotor radius
         - r_hub: position of the turbine hub in global coordinate system
         - e_shaft_yaw0: unit vector along the shaft (untitled for now), going downwind, when the turbine has zero yaw
         - e_vert: unit vertical vector, about which positive yawing is done
         - U0: Free stream velocity in global coordinates (can be changerd with `update_wind`)
         - Ct: Thrust coefficient (can be changed with `update_loading`)
         - Ground: Include ground effect in calculations
         - Model: one of ['VC','VCFF', 'VD', 'SS']
                   'VCFF': Vortex cylinder with far-field approximation (fastest)
                   'VC': Vortex cylinder
                   'SS': Self similar model of Troldborg et al.  (not good close to rotor)
                   'VD': Self similar model of Troldborg et al.  (not good close to rotor)
        """
        self.set_yaw0_coord(e_shaft_yaw0,e_vert)
        self.update_position(r_hub)
        self.update_wind(U0)
        self.name=name
        self.R=R
        self.r=None
        self.gamma_t=None
        self.gamma_t=None
        self.Gamma_r=None
        self.Lambda=Lambda
        self.Ground=Ground # Ground effect will be included in calculation of induced velocity
        self.chi=None
        self.Model=Model

    def set_yaw0_coord(self,e_shaft_yaw0,e_vert):
        self.e_shaft_g = np.asarray(e_shaft_yaw0).ravel().reshape(3,1)
        self.e_shaft_g0= self.e_shaft_g
        self.e_vert_g  = np.asarray(e_vert).ravel().reshape(3,1)
        v_horz0 = self.e_shaft_g -np.dot(self.e_shaft_g.T,self.e_vert_g)
        self.e_horz_g = v_horz0/np.linalg.norm(v_horz0)
        self.e_FAST_z = self.e_vert_g
        self.e_FAST_x = self.e_horz_g
        self.e_FAST_y = np.cross(self.e_FAST_z.T,self.e_FAST_x.T).T
        # Transformation matrix from global to FAST coordinate system
        self.T_F2g   = np.column_stack((self.e_FAST_x,self.e_FAST_y,self.e_FAST_z))
        # Transformation matrix from cylinder coordinate system to wind turbine
        # TODO: handle tilt
        e_c_x=np.cross(self.e_vert_g.T,self.e_shaft_g.T).T
        self.T_c2wt = np.column_stack((e_c_x,self.e_vert_g,self.e_shaft_g))
        self.update_yaw_pos(0)

    def update_yaw_pos(self,yaw):
        """ 
            yaw : radians, turbine yaw angle (with respect to the 0 yaw position defined at init)
        """
        self.yaw_pos = yaw
        self.T_wt2g  = RotMat_AxisAngle(self.e_vert_g,yaw)
        # Rotating the shaft vector so that its coordinate follow the new yaw position
        self.e_shaft_g=np.dot(self.T_wt2g , self.e_shaft_g0)

    def update_position(self,r_hub):
        self.r_hub=np.asarray(r_hub).ravel().reshape(3,1)

    def update_wind(self,U0_g):
        self.U0_g = np.asarray(U0_g).ravel().reshape(3,1)

    def update_loading(self,r=None,Ct=None,Gamma=None,Lambda=None,nCyl=1,gamma_t_Ct=None):
        """ 
        Computes relevant parameters when the turbine loading is updated, mainly, gamma_t, 
        the intensity of the tangential vorticity sheet.
        The ditributon will be determined based on the inputs, with one these three approaches:
           1. Ct(r) distribution
           2. Gamma(r) distribution
           3. gamma_t(Ct(r)) function

        INPUTS:
          r: radial coordinates at which Ct or Gamma are provided
          Ct: local thrust coefficient (Ct(r), array), or total thrust coefficient (CT, scalar)
          Gamma:  bound circulation (Gamma(r), array), or total rotor circulation (Gamma_tot, scalar)
          Lambda: tip speed ratio (assumed infinite if None)
          nCyl : number of cylindrical model used in the spanwise direction (default is 1)
                 The circulation (gamma_t) will be determined for each of the radial cylinder
          gamma_t_Ct: function that provides gamma_t as function of Ct (or gamma_t as function of CT)
        """
        U0=np.linalg.norm(self.U0_g)

        # --- Reinterpolating loading to number of cylinders if needed
        if nCyl is not None:
            if nCyl==1:
                vr0= np.array([0.995*self.R])
                if Ct is not None:
                    Ct =np.array([np.mean(Ct)])
                if Gamma is not None:
                    Gamma =np.array([np.mean(Gamma)])
            else:
                vr0= np.linspace(0.005,0.995,nCyl)*self.R
                if Ct is not None:
                    Ct = np.interp(vr0,r,Ct)
                else:
                    Gamma = np.interp(vr0,r,Gamma)
            r=vr0

        # Updating Lambda
        if Lambda is None:
            Lambda=self.Lambda
        if Lambda is None:
            raise Exception('Provide `Lambda` for update_loading. (Note: `Lambda=np.Inf` supported) ')
        Omega = Lambda*U0/self.R
        #print('U0',U0)
        #print('Ct',Ct)

        # Computing and storing gamma distribution and loading
        if gamma_t_Ct is not None:
            if Ct is None:
                raise Exception('Provide `Ct` along `gamma_t_Ct`')
            self.gamma_t = gamma_t_Ct(Ct)
            self.gamma_l=None # TODO
            self.Gamma_r=None # TODO
        elif Ct is not None:
            self.gamma_t,self.gamma_l,self.Gamma_r,misc=WakeVorticityFromCt(r,Ct,self.R,U0,Omega)
        elif Gamma is not None:
            self.gamma_t,self.gamma_l,self.Gamma_r,misc=WakeVorticityFromGamma(r,Gamma,self.R,U0,Omega)
        else:
            raise Exception('Unknown loading spec')
        #self.gamma_t=self.gamma_t*1.06
        #print('gamma_t    ',self.gamma_t)
        #print('gamma_l    ',self.gamma_l)
        #print('Gamma_r    ',self.Gamma_r)
        #print('Gamma_/2piR',-self.Gamma_r/(2*np.pi*self.R))
        #print(misc)
        self.Lambda=Lambda
        self.r=r
        self.Ct=Ct

    @property
    def yaw_wind(self):
        """ NOTE: this is wind angle not wind direction, measured with same convention as yaw:
            - around the axis e_vert
        """
        u_horz = self.U0_g - np.dot(self.U0_g.T,self.e_vert_g)*self.e_vert_g
        e_w    = u_horz/np.linalg.norm(u_horz)
        sign = np.sign  ( np.dot(np.cross(self.e_horz_g.T, self.U0_g.T),self.e_vert_g) )
        if sign==0:
            yaw_wind = np.arccos(np.dot(e_w.T,self.e_horz_g))
        else:
            yaw_wind = sign * np.arccos(np.dot(e_w.T,self.e_horz_g))
        return yaw_wind.ravel()[0]

    @property
    def yaw_error(self):
        v_horz = self.e_shaft_g -np.dot(self.e_shaft_g.T,self.e_vert_g)
        e_horz = v_horz/np.linalg.norm(v_horz)
        u_skew  = self.U0_g-np.dot(self.U0_g.T,e_horz)*e_horz-np.dot(self.U0_g.T,self.e_vert_g)*self.e_vert_g
        e_w=self.U0_g/np.linalg.norm(self.U0_g)
        yaw_sign=np.sign  ( np.dot(np.cross(e_w.T,self.e_shaft_g.T),self.e_vert_g) )
        yaw_error=yaw_sign * np.arccos(np.dot(e_w.T,self.e_shaft_g))
#         U0_wt = np.dot( self.T_wt2g.T , self.U0_g)
#         e_skew=u_skew/np.linalg.norm(u_skew)
#         print('e_skew',e_skew.T)
        return yaw_error.ravel()[0]

    #@property
    #def tilt(self):
    #    shaft_vert=np.dot(self.e_shaft.T,self.e_vert)
    #    return np.arcsin(shaft_vert)*180/np.pi

    #@property
    #def yaw(self):
    #    if self.tilt==0:
    #        u_skew=self.U0-np.dot(self.U0,self.e_shaft)
#   #          e_skew=self.

    #    else:
    #        raise Exception('Tilt and yaw not supported yet')
    #    return np.arcsin(shaft_vert)*180/np.pi

    def rotor_disk_points(self,nP=100):
        points=np.zeros((3,nP))
        theta=np.linspace(0,2*np.pi,nP) # dTheta=2 pi /(np-1)
        e_r = self.R*orth_vect(self.e_shaft_g)
        for i,t in enumerate(theta):
            T=RotMat_AxisAngle(self.e_shaft_g,t)
            points[:,i]= self.r_hub.ravel()+np.dot(T,e_r).ravel()
        return points


    def set_chi(self,chi):
        self.chi=chi
    
    def compute_u(self, Xg, Yg, Zg, only_ind=False, longi=False, tang=True, root=False, no_wake=False, ground=None, Model=None, R_far_field=6): 
        """ 
        INPUTS:
            Xg, Yg, Zg: Control points in global coordinates where the flow is to be computed. 
            only_ind: if true, only induction is returned (without the free stream)

            longi, tang, root: booleans specifying which component of vorticity is considered.
                               Default is `tang` only
            no_wake: boolean, if true: the induced velocity in the wake is set to 0. 
                     Typically set to true when combining with wake models.

            Model : string in ['VC','VCFF','SS','VD']
                   'VCFF': Vortex cylinder with far-field approximation (fastest)
                   'VC': Vortex cylinder
                   'SS': Self similar model of Troldborg et al.  (not good close to rotor)
                   'VD': Self similar model of Troldborg et al.  (not good close to rotor)
        """
        # --- Optional argument overriding self
        if ground is None:
            ground=self.Ground
        if Model is None:
            Model=self.Model

        # Control points in "Cylinder coordinate system" (rotation only)
        T_c2g=np.dot(self.T_wt2g,self.T_c2wt)
        Xc,Yc,Zc = transform_T(T_c2g, Xg,Yg,Zg)
        # Detecting whether our vertical convention match, and define chi
        e_vert_c = np.dot(T_c2g.T , self.e_vert_g)
        if self.chi is None:
            # TODO TODO chi needs induction effect!
            self.chi= np.sign(e_vert_c.ravel()[1])* (self.yaw_wind-self.yaw_pos)
        if self.gamma_t is None:
            raise Exception('Please set loading with `update_loading` before calling `compute_u`')

        uxc = np.zeros(Xg.shape)
        uyc = np.zeros(Xg.shape)
        uzc = np.zeros(Xg.shape)
        m=np.tan(self.chi)
        # Cylinder position in "Cylinder coordinate system) (rotation only)
        Xcyl, Ycyl, Zcyl = transform_T(T_c2g,np.array([self.r_hub[0]]), np.array([self.r_hub[1]]),  np.array([self.r_hub[2]]))
        # Translate control points such that origin is at rotor center. NOTE: not all routines use this
        Xc0,Yc0,Zc0=Xc-Xcyl[0],Yc-Ycyl[0],Zc-Zcyl[0]
        if ground:
            # Mirror control points are two time the hub height above the cylinder
            Yc0mirror=Yc0+2*Ycyl[0]
            Ylist=[Yc0,Yc0mirror]
            #print('>>> Ground effect',Ycyl[0])
        else:
            Ylist=[Yc0]

        # --- Root vortex influence
        if root and  (self.Gamma_r is not None) and self.Gamma_r!=0:
            for Y in Ylist:
                if np.abs(self.chi)>1e-7:
                    uxc0,uyc0,uzc0 = svc_root_u(Xc0,Y,Zc0,Gamma_r=self.Gamma_r,m=m,polar_out=False)
                else:
                    uxc0,uyc0,uzc0 =  vc_root_u(Xc0,Y,Zc0,Gamma_r=self.Gamma_r,polar_out=False)
                uxc += uxc0
                uyc += uyc0
                uzc += uzc0


        if len(self.gamma_t)==1:
            # --- Tangential and longi - ONE Cylinder only
            for iY,Y in enumerate(Ylist):
                if tang and (self.gamma_t!=0):
                    if np.abs(self.chi)>1e-7:
                        if Model =='VC':
                            uxc0,uyc0,uzc0 = svc_tang_u(Xc0,Y,Zc0,gamma_t=self.gamma_t,R=self.R,m=m,polar_out=False)
                        else:
                            raise NotImplementedError('Model '+Model + ', with yaw.')
                    else:
                        if Model =='VC':
                                uxc0,uyc0,uzc0 = vc_tang_u        (Xc0,Y,Zc0, gamma_t=self.gamma_t, R=self.R, polar_out=False)
                        elif Model =='VCFF':
                            uxc0,uyc0,uzc0 = vc_tang_u_doublet(Xc0,Y,Zc0, gamma_t=self.gamma_t, R=self.R, polar_out=False,r_bar_Cut=R_far_field)
                        elif Model =='VD':
                            uxc0,uyc0,uzc0 = doublet_line_u(Xc0, Y, Zc0, dmz_dz = self.gamma_t * self.R**2 * np.pi)
                        elif Model =='SS':
                            uzc0 = ss_u          (Xc0, Y, Zc0, gamma_t=self.gamma_t, R=self.R)
                            uxc0=uzc0*0
                            uyc0=uzc0*0
                        else:
                            raise NotImplementedError('Model'+Model)
                    uxc += uxc0
                    uyc += uyc0
                    uzc += uzc0
                if longi and (self.gamma_l is not None) and self.gamma_l!=0 :
                    if np.abs(self.chi)>1e-7:
                        if Model =='VC':
                            uxc0,uyc0,uzc0 = svc_longi_u(Xc0,Y,Zc0,gamma_l=self.gamma_l,R=self.R,m=m,polar_out=False)
                        else:
                            raise NotImplementedError('Model '+Model + ', longi component.')
                    else:
                        if Model =='VC':
                            uxc0,uyc0,uzc0 = vc_longi_u (Xc0,Y,Zc0,gamma_l=self.gamma_l,R=self.R    ,polar_out=False)
                        else:
                            raise NotImplementedError('Model'+Model + ', longi component.')
                    uxc += uxc0
                    uyc += uyc0
                    uzc += uzc0
        else:
            # --- Tangential and longi - MULTI Cylinders
            if Model =='VC':
                nr   = len(self.r)
                nWT = 1
                # Control points are directly translated by routine
                gamma_t = self.gamma_t.reshape((nWT,nr))
#                 print('r      ',self.r)
#                 print('gamma_t',gamma_t)
                if self.gamma_l is not None:
                    gamma_l = self.gamma_l.reshape((nWT,nr))
                vR      = self.r.reshape((nWT,nr))
                vm       = m* np.ones((nWT,nr))
                if tang:
                    if np.abs(self.chi)>1e-7:
                        uxc0,uyc0,uzc0 = svcs_tang_u(Xc,Yc,Zc,gamma_t=gamma_t,R=vR,m=vm,Xcyl=Xcyl,Ycyl=Ycyl,Zcyl=Zcyl,Ground=ground)
                    else:
                        uxc0,uyc0,uzc0 = vcs_tang_u (Xc,Yc,Zc,gamma_t=gamma_t,R=vR    ,Xcyl=Xcyl,Ycyl=Ycyl,Zcyl=Zcyl, Ground=ground)
                    uxc += uxc0
                    uyc += uyc0
                    uzc += uzc0
                if longi and (self.gamma_l is not None):
                    if np.abs(self.chi)>1e-7:
                        uxc0,uyc0,uzc0 = svcs_longi_u(Xc,Yc,Zc,gamma_l=gamma_l,R=vR,m=vm,Xcyl=Xcyl,Ycyl=Ycyl,Zcyl=Zcyl, Ground=ground)
                    else:
                        uxc0,uyc0,uzc0 = vcs_longi_u (Xc,Yc,Zc,gamma_l=gamma_l,R=vR      ,Xcyl=Xcyl,Ycyl=Ycyl,Zcyl=Zcyl, Ground=ground)
                    uxc += uxc0
                    uyc += uyc0
                    uzc += uzc0
            else:
                raise NotImplementedError('Model'+Model, 'with multiple cylinders')
        if no_wake:
#             uxc[:]=0
#             uyc[:]=0
#             uzc[:]=1
            # Zero wake induction
            bDownStream=Zc0>=-0.20*self.R
#             bDownStream=Zc0>=0
            Rc = np.sqrt(Xc0**2 + Yc0**2)
            bRotorTube = Rc<self.R*1.001 # we give a margin since VD and VC have fields very dissimilar at R+/-eps
            bSelZero = np.logical_and(bRotorTube,bDownStream)
            uxc[bSelZero]=0
            uyc[bSelZero]=0
            uzc[bSelZero]=0
            # Decay
#             rc = np.sqrt(Xc**2 + Yc**2 + Zc**2)
#             RadialDecay = np.ones(rc.shape)
#             bSelDecay = np.logical_and(rc>self.R, bDownStream)
#             RadialDecay[bSelDecay] = 1/(rc[bSelDecay]/self.R)**2
#             ZDecay = np.ones(uxc.shape)
            #bSelDecay = np.logical_and(Zc>self.R, bDownStream)
#             bSelDecay = bDownStream
            #ZDecay[bSelDecay] = 1/(Zc[bSelDecay]/self.R)**2
#             ZDecay[bSelDecay] = np.exp(-(Zc[bSelDecay]/self.R)**2)
#             uxc*=ZDecay
#             uyc*=ZDecay
#             uzc*=ZDecay

        # Back to global
        uxg,uyg,uzg = transform(T_c2g, uxc, uyc, uzc)
        # Add free stream if requested
        if not only_ind:
            uxg += self.U0_g[0]
            uyg += self.U0_g[1]
            uzg += self.U0_g[2]
        return uxg,uyg,uzg

    def tostring(self,short=True):
        s ='class WindTurbine({}), with attributes:\n'.format(self.name)
        s+=' - R      : {}\n'.format(self.R)
        s+=' - r_hub  : {}\n'.format(self.r_hub.T)
        s+=' - U0     : {}\n'.format(self.U0_g.T)
        s+=' - yaw_err: {} deg\n'.format(self.yaw_error*180/np.pi)
        if not short:
            s+=' - wd     : {} deg\n'.format(self.yaw_wind*180/np.pi)
            s+=' - e_shaft: {}\n'.format(np.round(self.e_shaft_g.T,4))
            s+=' - yaw_pos: {} deg\n'.format(self.yaw_pos  *180/np.pi)
            s+=' - yaw_err: {} deg\n'.format((-self.yaw_wind+self.yaw_pos)*180/np.pi)
            s+=' - T_F2g  : \n{}\n'.format(np.round(self.T_F2g ,4))
            s+=' - T_c2wt : \n{}\n'.format(np.round(self.T_c2wt,4))
            s+=' - T_wt2g : \n{}\n'.format(np.round(self.T_wt2g,4))
            s+=' - T_c2g  : \n{}\n'.format(np.round(np.dot(self.T_wt2g,self.T_c2wt),4))
        return s

    def __repr__(self):
        s=self.tostring(short=False)
        #s+=' - tilt   :  {}\n'.format(self.tilt)
        #s+=' - e_shaft: {}\n'.format(self.e_shaft)
        #s+=' - e_vert: {}\n'.format(self.e_vert)
        return s

if __name__ == "__main__":
    unittest.main()