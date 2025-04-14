def getvm(al, vw, cmsq):
    """
    Parameters
    ----------
    al : TYPE
        The phase transition strength parameter \alpha.
    vw : TYPE
        The bubble wall velocity
    cmsq : TYPE
        The sound speed of the broken phase

    Returns
    -------
    Tuple
        The fluid velocity just behind the wall;
        The type of hydrodynamical mode
        0: deflagration
        1: hybrid
        2: detonation
   
    """
    if vw**2 < cmsq:
        return (vw, 0)
    cc = 1. - 3.*al + vw**2*(1./cmsq + 3.*al)
    disc = -4.*vw**2/cmsq + cc**2
    if (disc < 0.)|(cc < 0.):
        return (np.sqrt(cmsq), 1)
    return ((cc + np.sqrt(disc))/2.*cmsq/vw, 2)

def vJ_pseudo(alpha,cmsq):
    r1 = 1 + np.sqrt(3*alpha*(1 - cmsq + 3*cmsq*alpha))
    r2 = 3*np.sqrt(cmsq)*alpha + 1/np.sqrt(cmsq)
    r3 = 1 - np.sqrt(3*alpha*(1 - cmsq + 3*cmsq*alpha))
    R1 = r1/r2
    R2 = r3/r2
    return R1,R2

def get_K_Wow(vw, v0, cssq,r=1,profile=False):
    if v0==0:
        return 0,1
    n = 8*1024 # change accuracy here
    vs = np.linspace(v0, 0, n)
    sol = odeint(dfdv, [vw,1.,1.], vs, args=(cssq,))
    xis, wows, ToTs = (sol[:,0], sol[:,1], sol[:,2])
    if mu(vw,v0)*vw <= cssq:
        ll = max(int(sum(np.heaviside(cssq - (mu(xis,vs)*xis),0.))), 1)
        vs = vs[:ll]
        xis = xis[:ll]
        wows = wows[:ll]/wows[ll-1]*getwow(xis[-1], mu(xis[-1],vs[-1])) 
        # ToTs = ToTs[:ll]/ToTs[ll-1]*getwow(xis[-1], mu(xis[-1],vs[-1]))/r
    Kint = simpson(wows*(xis*vs)**2/(1.-vs**2), x=xis)
    if profile==False:
        return (Kint*4./vw**3, wows[0])#, xis, wows)
    else:
        return (xis,vs,wows,ToTs)

def alN(al, wow, cmsq, cpsq):
    da = (1./cmsq - 1./cpsq)/(1./cpsq + 1.)/3.
    return (al+da)*wow -da

def get_alN_wow(vp,vm,vw,cmsq,cpsq):
    Ksh,wow = get_K_Wow(vw,mu(vw,vp),cpsq)
    al = (vp/vm-1.)*(vp*vm/cmsq - 1.)/(1-vp**2)/3.
    return (alN(al,wow,cmsq,cpsq), wow) 

def solv(vw, v0, cssq, n):
    if v0==0:
        return 0,1
    vs = np.linspace(v0, 0, n)
    sol = odeint(dfdv, [vw,1.,1.], vs, args=(cssq,))
    xis, wows, ToTs = (sol[:,0], sol[:,1], sol[:,2])
    if mu(vw,v0)*vw <= cssq:
        ll = max(int(sum(np.heaviside(cssq - (mu(xis,vs)*xis),0.))), 1)#find the shock front
        vs = vs[:ll]
        xis = xis[:ll]
        wows = wows[:ll]/wows[ll-1]*getwow(xis[-1], mu(xis[-1],vs[-1]))#matching condition at shock front to give the w/w_n
        ToTs = ToTs[:ll]/ToTs[ll-1]*(getwow(xis[-1], mu(xis[-1],vs[-1])))**(0.25)#for the shock front, the relation of enthalpy and temperature
    return (xis,vs,wows,ToTs) 

def lbdDT(alpha, wpro, cmsq):
    r = wpro - (1. + 3.*cmsq*alpha)
    r *= 1./(1. + cmsq)
    return r 

def lbdDF(alpha, wpro, cpsq, cmsq):
    a = 1. + 1./cmsq
    b = (1. - cpsq/cmsq)/(1. + cpsq)
    d = 1. + cpsq
    r = (wpro + d*(3.*alpha - b)/a)/(1. + d*(3.*alpha - b)/a) -1
    r *= 1./d + (3.*alpha - b)/a
    return r

def profile(vw, al, cpsq, cmsq, r=1, n = 8*1024): # change accuracy here):
    vm, mode = getvm(al, vw, cmsq)
    m = 1024#512 # the precision for the zero v part of the full profile
    # print(vm)
    # print(mode)
    if mode<2:
        almax,wow = get_alN_wow(0,vm,vw,cmsq,cpsq)
        # print('a1 =',almax)
        if almax<al:
            print ("alpha too large for shock")
            return 0
        vp = min(cpsq/vw,vw)
        almin,wow = get_alN_wow(vp,vm,vw,cmsq,cpsq)
        # print('a2 =',almin)
        if almin>al:
            print ("alpha too small for shock")
            return 0
        iv = [[vp,almin],[0,almax]]
        while (abs(iv[1][0]-iv[0][0])>1e-7):
            vpm = (iv[1][0]+iv[0][0])/2.
            #print(vpm)
            alm = get_alN_wow(vpm,vm,vw,cmsq,cpsq)[0]
            #print(alm)
            if alm>al:
                iv = [iv[0],[vpm,alm]]
            else:
                iv = [[vpm,alm],iv[1]]
        vp = (iv[1][0]+iv[0][0])/2.
        xis,vs,wows,ToTs = solv(vw, mu(vw,vp),cpsq, n)
        wm = wows[0]*getwow(vp, vm)
        # print(wm)
        # print(vp)
        xi0 = np.linspace(0.01,vw-1e-3,m)
        v0 = np.zeros_like(xi0)
        w0 = wm*np.ones_like(xi0)
        xi1 = np.linspace(xis[-1]+1e-3, 0.99, m)
        v1 = np.zeros_like(xi1)
        w1 = np.ones_like(xi1)
        # print(xi1)
   
    if mode == 2:
        # print('detonation')
        # detionation
        xirf,vrf,wowrf,ToTrf = solv(vw, mu(vw, vm), cmsq, n)
        wowrf *= getwow(vw, vm)
        # print(getwow(vp, vm))
        ToTrf *= (wowrf[0]*((1. + cpsq)/(1. + cmsq))/r)**(0.25)
        xi0 = np.linspace(0.01, np.sqrt(cmsq)-1e-3, m)
        v0 = np.zeros_like(xi0)
        w0 = wowrf[-1]*np.ones_like(xi0)
        xi1 = np.linspace(vw+1e-3, 0.99, m)
        v1 = v0
        w1 = np.ones_like(xi1)
        xi = np.concatenate((xi0, xirf[::-1], xi1))
        vpro = np.concatenate((v0, vrf[::-1], v1))
        wpro = np.concatenate((w0, wowrf[::-1], w1))

        # lbdpro = np.concatenate((w0, wowrf[::-1])) - (1 + 3*cmsq*al)
        # lbdpro *= 1/(1 + cmsq)
        # lbdpro = np.concatenate((lbdpro, np.zeros_like(xi1)))
        wp = np.concatenate((w0, wowrf[::-1]))
        lbdpro = lbdDT(al, wp, cmsq)
        lbdpro = np.concatenate((lbdpro, np.zeros_like(xi1)))

    elif mode ==1:
        # print('hybrid')/
        # hybrid
        # xirf,vrf,wowrf,ToTrf = (0,0,0,0)
        xirf,vrf,wowrf,ToTrf = solv(vw, mu(vw, np.sqrt(cmsq)), cmsq, n)
        wowrf *= wows[0]*getwow(vp, np.sqrt(cmsq))
        ToTrf *= (wowrf[0]*((1. + cpsq)/(1. + cmsq))/r)**(0.25)
        xihs = np.concatenate((xirf[::-1],xis[1:]))
        vhs = np.concatenate((vrf[::-1],vs[1:]))
        whs = np.concatenate((wowrf[::-1],wows[1:]))
        xi0 = np.linspace(0.01, np.sqrt(cmsq)-1e-3, m)
        v0 = np.zeros_like(xi0)
        w0 = wowrf[-1]*np.ones_like(xi0)
        xi1 = np.linspace(xis[-1]+1e-3, 0.99, m)
        v1 = np.zeros_like(xi1)
        w1 = np.ones_like(xi1)
        xi = np.concatenate((xi0, xihs, xi1))
        vpro = np.concatenate((v0, vhs, v1))
        wpro = np.concatenate((w0, whs, w1))

        # lbdpro = np.zeros_like((1000,))
        wpdt = np.concatenate((w0, wowrf[::-1]))
        lbddt = lbdDT(al, wpdt, cmsq)
        wpdf = np.concatenate((wows[1:], w1))
        lbddf = lbdDF(al, wpdf, cpsq, cmsq)
        lbdpro = np.concatenate((lbddt, lbddf))

    elif mode == 0:
        # deflagration 
        xi = np.concatenate((xi0, xis, xi1))
        vpro = np.concatenate((v0, vs, v1))
        wpro = np.concatenate((w0, wows, w1))
        
        wp = np.concatenate((wows, w1))
        lbdpro = lbdDF(al, wp, cpsq, cmsq)
        lpm = lbdDT(al, wm, cmsq)
        lbdpro = np.concatenate((lpm*np.ones_like(xi0), lbdpro))
        # lbdpro = lbdDF(al, wpro, cpsq, cmsq)

    def smooth(x, y, n):
        new_x = np.linspace(x[0], x[-1], n)
        return new_x, np.interp(new_x, x, y)
    
    xi_new, vip_new = smooth(xi, vpro, len(xi))
    xi_new, wpro_new = smooth(xi, wpro, len(xi))
    xi_new, lamip_new = smooth(xi, lbdpro, len(xi))

    return xi_new, vip_new, wpro_new, lamip_new
    
    return xi, vpro, wpro, lbdpro




def EkinoRs(qR,T,vw,alpha,cpsq=1/3.,cmsq=1/3.,n_prof=4000, ltdist='exp',bk='cpu'):
    """
    The velocity power spectrum of Sound Shell Model
    
    parameter
    ---------
    qR: 1d array
        the frequency band of velocity spectrum.
    T: 1d array
        the time scaled by beta, accoording to the bubble life time distribution
        for the exponential nucleation rate, the bubble life time distribution 
        exponentially suppressed, the uper limit of T shoulb be large enough, 
        basically [0, 20] is enough.
    vw: float
        bubble wall velocity
    alpha: float
        strength parameter
    cpsq: float
        square of sound velocity of symmetric phase
    cmsq: float
        square of sound velocity of symmetric phase
    ltdist:
        bubble lifetime distribution
        
    return
    ------
    Pv: 1d array
        the velocity power specturm 
    A: 2d array
    realA: 2d array
    imgA: 2d array
    realAsq: 2d array
    imgAsq: 2d array
    """

    
    qR0 = qR
    qR = qR.reshape((qR.shape[0],1)) # qR.shape = (n,)--->(n,1), T.shape = (m,)
    z0 = qR*T/(vw*(np.pi*8)**(1/3.)) # 2-d array of z, z.shape(n,m)
    #z = z0.reshape((qR.shape[0],T.shape[0])) # get a 3d array, (n,m,1)
    xi, vip, wpro, lamip = profile(vw, alpha, cpsq, cmsq, n=n_prof)
    
    fz = fz_gen( vip,  z0,  xi)
    fz *= 4*np.pi/z0
    lz = lz_gen( vip,  z0,  lamip, xi)
    lz *= 4*np.pi/z0

    realA = np.zeros_like(fz)
    for i in range(fz.shape[0]):
        rA = 0.5*np.gradient(fz[i,::],z0[i,::])
        realA[i,...] = rA
    
    realA = np.array(realA)

    realA = np.array(realA)
    imgA = 0.5*lz*np.sqrt(cmsq)
    realAsq =  realA**2
    imgAsq = imgA**2   
    A = realAsq + imgAsq # (n,m) array
    if ltdist == 'exp':
        mublife = np.exp(-T)
    elif ltdist == 'sim':
        mublife = 0.5*T**2*np.exp(-T**3/6.)
    Pvintegrand = T**6*mublife*np.abs(A) # integrand of power spectrum
    Pv = np.trapezoid(Pvintegrand,T)
    Pv *= qR0**2/(128.*np.pi**4*vw**6)

    
    return Pv#, A, realA, imgA, realAsq, imgAsq







def Deltamn_flat(delt, pmn):
    r = (1 - np.cos(pmn*delt))/(2.*pmn**2)
    return r

def Deltamn(Pht, Htau, HRs):
    """
    The function controls the shape of the spectrum in the low frequency range

    Pht = pmn * Rs:
        R_* normalised momentum 

    Htau = dtau * Hs:
        H_* normalised source lifetime

    HRs = Hs * Rs
        H_* normalised mean bubble seperation
    """
    Si0, Ci0 = sici(Pht*(Htau + 1.)/HRs) # Htau + 1 = tau_fin * Hs (see under Fig 3)
    Si, Ci = sici(Pht/HRs)
    DSi = Si0 - Si
    DCi = Ci0 - Ci
    Dmn = 0.25*(DCi**2 + DSi**2)
    return Dmn

def Delt0(Pht, Htau, HRs, cs = 1/np.sqrt(3)):
    part1 = np.log(Htau + 1.)**2
    Pht = 2.*Pht*cs
    Si0, Ci0 = sici(Pht*(Htau + 1.)/HRs)
    Si, Ci = sici(Pht/HRs)
    DSi = Si0 - Si
    DCi = Ci0 - Ci
    r = DSi**2 + DCi**2 + part1
    return 0.5*r 

def Delta(Htau, HRs, K3, P3, z3, Pt, cs=1/np.sqrt(3)):
    """
    The full delta function
    """

    Dlt = 0
    for m in [1,-1]:
        for n in [1,-1]:
            Pht = (P3 + m*Pt)*cs + n*K3 
            #Dlt += Deltamn(Pht, Htau, HRs)
            Dlt += np.array(Deltamn_cy(Pht, Htau, HRs))
    return Dlt


# this function is unused! zetakin is just EkinoRs (I think)
def tildeDlt0(Htau, HRs, K, zetakin):
    """
    K and zetakin should be 1d array
    """
    integrand1 = zetakin**2/K**2
    integrand2 = Delt0(K, Htau, HRs)
    integrand = integrand1*integrand2
    num = trapezoid(integrand, K)
    den = trapezoid(integrand1, K)
    return num/den

# def tildeDlt():
def tildeDlt(K, P, z, vw, alpha, Htau, HRs, cs=1/np.sqrt(3.), n_prof=5000, n_ikR=500,
             log_kR_min=9e-6, log_kR_max=650, nsize=100, T=np.logspace(-2, np.log10(20), 1000)):
    """
    To calculate the integral part for the normalized UETC of shear stress 
    We neglect the prefactor part

    Parameters:
        K: 1d array (m,)
            kR_*: the scaled frequency
        P: 1d array (n,)
            pR_*: the sclaled integration variable
        z: 1d array (i,)
            the integration variable
    """
    tot_t = 0

    tmp_kR = np.logspace(np.log10(log_kR_min), np.log10(log_kR_max), n_ikR)
    
    tmp_Tt = T
    tmp_Pv = EkinoRs(tmp_kR, tmp_Tt, vw, alpha,n_prof=n_prof)
    EkinoRs_f = CubicSpline(tmp_kR, tmp_Pv)
    
    integrand1 = EkinoRs_f(P)
    
    EstoRst = np.max(integrand1)
    integrand1 *= 1/EstoRst
    integrand1 *= P**2 # we got the 1d array with shape (n,)
    
    r = np.zeros_like(K)
    integrand = np.zeros_like(P)
    for Ki in range(len(K)):
        for Pi in range(len(P)):
            Pt = np.sqrt(np.abs(K[Ki]**2 + P[Pi]**2 - 2. * P[Pi] * K[Ki] * z) )
            Ptt = Pt
            
            zetaPtt = EkinoRs_f(Ptt)/EstoRst
            
            
            integrand2 = (1. - z**2)**2/Pt**4 # 3d array
            integrand2 *= zetaPtt
            
            integrand2 *= Delta(Htau, HRs, K[Ki], P[Pi], z, Pt)
            
            tot_t += t1 - t0
            integral2 = trapezoid(integrand2, z) 
            #integral2 = simpson(integrand2, x=z) 
            
            integrand[Pi] = integral2*integrand1[Pi]
        r[Ki] = trapezoid(integrand, P)
        #r[Ki] = simpson(integrand, x=P)
    
    return r, integrand#r*prefactor, prefactor, r, EstoRst


def Delta_intp(Htau, HRs, K3, P3, z3, Pt, Dmn_f, cs=1/np.sqrt(3)):
    """
    The full delta function
    """

    Dlt = 0
    for m in [1,-1]:
        for n in [1,-1]:
            Pht = (P3 + m*Pt)*cs + n*K3 
            Dlt += Dmn_f(np.abs(Pht))
    return Dlt

# eq 95: tildeDlt = zetaGW / zetaPi
# tildeDlt_intp[0] = eq 51 without prefactor (integral part only) and divided by Ekin_max**2
# tildeDlt_intp[1] = Ekin_max
def tildeDlt_intp(K, P, z, vw, alpha, Htau, HRs, cs=1/np.sqrt(3.), n_prof=5000, n_ikR=500,
             log_kR_min=1e-6, log_kR_max=650, nsize=100, T=np.logspace(-2, np.log10(20), 1000)):
    """
    To calculate the integral part for the normalized UETC of shear stress 
    We neglect the prefactor part

    Parameters:
        K: 1d array (m,)
            kR_*: the scaled frequency
        P: 1d array (n,)
            pR_*: the sclaled integration variable
        z: 1d array (i,)
            the integration variable
    """
    # creating interpolation function for Dmn
    Pht_max = ((np.sqrt(K.max() + P.max())**2 + P.max() ) * (1/np.sqrt(3)) + K.max()) * 1.05
    nPhs = 5000000
    
    Phs = np.logspace(-11, np.log10(Pht_max), nPhs)
    Phs = np.concatenate( ( np.flip(-Phs), Phs) )
    
    Si0, Ci0 = sici(Phs*(Htau + 1.)/HRs)
    Si, Ci = sici(Phs/HRs)
    
    Dmn = 0.25*((Si0 - Si)**2 + (Ci0 - Ci)**2)
    del Si0, Ci0, Si, Ci
    #Dmn_f = interp1d(Phs, Dmn)


    tmp_kR = np.logspace(np.log10(log_kR_min), np.log10(log_kR_max), n_ikR)
    
    tmp_Tt = T
    tmp_Pv = EkinoRs(tmp_kR, tmp_Tt, vw, alpha,n_prof=n_prof) ## update this to include option for nucleation type (default exp)
    #EkinoRs_f = CubicSpline(tmp_kR, tmp_Pv)
    EkinoRs_f = interp1d(tmp_kR, tmp_Pv)
    
    integrand1 = EkinoRs_f(P)
    
    EstoRst = np.max(integrand1)
    integrand1 *= 1/EstoRst
    integrand1 *= P**2 # we got the 1d array with shape (n,)
    
    r = np.zeros_like(K)
    integrand = np.zeros_like(P)
    #tot_t = 0
    for Ki in range(len(K)):
        for Pi in range(len(P)):
            Pt = np.sqrt(np.abs(K[Ki]**2 + P[Pi]**2 - 2. * P[Pi] * K[Ki] * z) ) # Ptilde
            Ptt = Pt
            #zetaPtt = EkinoRs_f(Ptt)/EstoRst
            zetaPtt = np.array(c_interp(Ptt, tmp_kR, tmp_Pv))/EstoRst # Ekin/Ekin_max as a function of Ptilde
            
            integrand2 = zetaPtt * (1. - z**2)**2/Pt**4 # 3d array
            
            #integrand2 *= Delta_intp(Htau, HRs, K[Ki], P[Pi], z, Pt, Dmn_f)
            integrand2 *= np.array(Delta_interp_cy(Htau, HRs, K[Ki], P[Pi], z, Pt, Phs, Dmn))
            
            integral2 = trapezoid(integrand2, z) 

            integrand[Pi] = integral2*integrand1[Pi]
        r[Ki] = trapezoid(integrand, P)
        #r[Ki] = simpson(integrand, x=P)
    
    OK_over_K = EstoRst 
    return r, OK_over_K