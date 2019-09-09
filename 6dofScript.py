def compute6dof(mtotal, ftotal, H, txnmat):
    """Computes 6 dof dynamics on the given body
       mg - znGroupInfo object 
       mtotal - Total moments in body frame.
       ftotal - Total force in inertial frame.
       H - The timestep to be used."""

    K1 = ftotal[0,0]/mass 
    K2 = ftotal[0,1]/mass
    K3 = ftotal[0,2]/mass
    
    if not simple_euler_equations:
        K4 = (1.0/Ixx)*((Iyy-Izz)*q*r + mtotal[0,0])
        K5 = (1.0/Iyy)*((Izz-Ixx)*p*r + mtotal[0,1])
        K6 = (1.0/Izz)*((Ixx-Iyy)*p*q + mtotal[0,2])
    else:    
        K4 = (1.0/Ixx)*(mtotal[0,0])
        K5 = (1.0/Iyy)*(mtotal[0,1])
        K6 = (1.0/Izz)*(mtotal[0,2])

    K7 = (p + q*(sin(PHI))*tan(THETA) + r*(cos(PHI))*tan(THETA))
    K8 = (q*cos(PHI) - r*sin(PHI))
    K9 = (q*(sin(PHI))/cos(THETA) + r*(cos(PHI))/cos(THETA))

    # use the average of previous and current euler angle rates to update the
    # euler angles.
    K7dash = K7 + phidot
    K8dash = K8 + thetadot
    K9dash = K9 + psidot

    K7dash /= 2.0
    K8dash /= 2.0
    K9dash /= 2.0

    phidot = K7
    thetadot = K8
    psidot = K9
    
    K7 = K7dash
    K8 = K8dash
    K9 = K9dash

    # use the average of current and next velocities to update the 
    # positions
    K10 = (2*u + K1*H)/2.0
    K11 = (2*v + K2*H)/2.0
    K12 = (2*w + K3*H)/2.0
  
    u = u + K1*H
    v = v + K2*H
    w = w + K3*H
    
    # compute the velocities in fixed coordinates
    #uf = txnmat[0,0]*u + txnmat[1,0]*v + txnmat[2,0]*w
    #vf = txnmat[0,1]*u + txnmat[1,1]*v + txnmat[2,1]*w
    #wf = txnmat[0,2]*u + txnmat[1,2]*v + txnmat[2,2]*w
    uf = u
    vf = v
    wf = w

    p = p + K4*H
    q = q + K5*H
    r = r + K6*H

    dx = K10*H
    dy = K11*H
    dz = K12*H

    dPHI = K7*H
    dTHETA = K8*H
    dPSI = K9*H

