import numpy as np
from scipy.misc import factorial

def uq_get1DNodesAndWeights(nclp,pc_type):
    if pc_type == 'LEGENDRE':
        x,w=np.polynomial.legendre.leggauss(nclp)
        w=0.5*w
    elif pc_type == 'HERMITE':
        x,w=np.polynomial.hermite.hermgauss(nclp)
        x=x*np.sqrt(2)
        w=w/np.sqrt(np.pi)
    elif pc_type == 'LAGUERRE':
        x,w=np.polynomial.laguerre.laggauss(nclp)
    else:
        print 'Unknown Quadrature rule!'
        return
    return x,w

def uq_computeNPCTerms(ndim,nord):
    if nord ==-1:
        return 0
    enume=1
    denom=1
    minNO=min(nord,ndim)

    for k in range(minNO):
        enume=enume*(nord+ndim-k)
        denom=denom*(k+1)

    nPCTerms=enume/denom
    return nPCTerms
    
def uq_psi(x,nord,pc_type):
    nclp=len(x)
    psi=np.zeros((nclp,nord+1))

    if pc_type == 'LEGENDRE':
        psi[:,0]=1
        if (nord>0):
            psi[:,1]=x
            for i in range(2,nord+1):
                ca=(2.0*(i-1)+1.)*x
                cb=1.0*(i-1)
                psi[:,i]=(ca * psi[:, i-1] - cb* psi[:, i-2])/(1.0*i);
    elif pc_type == 'HERMITE':
        psi[:,0]=1
        if (nord>0):
            psi[:,1]=x
            for i in range(2,nord+1):
                psi[:,i]=x*psi[:,i-1]-1.0*(i-1)*psi[:,i-2]
##    elif pc_type == 'LAGUERRE':
##        psi[:,0]=1
##        if (nord>0):
##            psi[:,1]=1-x
##            for i in range(2,nord+1):          
    else:
        print ' Unknown polynomial type!'
        return
    
    return psi        

def uq_initMultiIndex(nord,ndim):
    if (ndim ==0): return 1
    npc=uq_computeNPCTerms(nord,ndim)
    iup=0
    isum=0
    ic=np.zeros(ndim,dtype=np.int)
    ict=np.zeros(ndim,dtype=np.int)

    mi=np.zeros((npc,ndim),dtype=int)

    iup=0
    for idim in range(ndim):
        mi[iup,idim]=0

    if nord > 0:
        for idim in range(ndim): ##check
            iup = iup+1
            mi[iup,idim]=1
            ic[idim]=1
    if nord >1:
        for iord in range(2,nord+1):
            lessiord=iup
            for idim in range(ndim):
                isum=0
                for ii in range(idim,ndim):
                    isum=isum+ic[ii]
                ict[idim]=isum
            for idim in range(ndim):
                ic[idim]=ict[idim]
            
            for idimm in range(ndim):
                a=lessiord-ic[idimm]+1
                b=lessiord+1
                for ii in range(a,b):
                    iup=iup+1
                    for idim in range(ndim):
                        mi[iup,idim]=mi[ii,idim]
                    mi[iup,idimm]=mi[iup,idimm]+1
            
    
    return mi,npc

def uq_evalBasisNormsSquared(nord,ndim,pc_type):
    multiIndex,npc=uq_initMultiIndex(nord,ndim)
    psiMultiDsq=np.zeros((npc,1))

    def funleg(n):
        return 1.0/(2.0*n+1.0)

    def funhermite(n):
        return factorial(n)

    def funlag(n):
        return np.ones((npc,1))

    for k in range(npc):
        pprod=1.0
        for j in range(ndim):
            if pc_type == 'LEGENDRE':
                a=funleg(multiIndex[k,j])
            elif pc_type == 'HERMITE':
                a=funhermite(multiIndex[k,j])
            else:
                print ' Unknown quadrature'
                return
            pprod=pprod * a
            
        psiMultiDsq[k]=pprod
        
    return psiMultiDsq

##def uq_settripleprod(pcdata): # TODO: Need to implement
##
##    nPCTerms=pcdata['nPCTerms']
##    nclp=pcdata['nclp']
##    w=pcdata['w']
##    x=pcdata['x']
##    psi=pcdata['psi']
##    
##    isparse={}#np.zeros((nPCTerms),dtype=np.int)
##    jsparse={}#np.zeros((nPCTerms),dtype=np.int)
##    csparse={}#np.zeros((nPCTerms))
##    nsparse=[]#np.zeros(nPCTerms,dtype=np.int)
##
##    tol=1e-4
##    for k in range(nPCTerms):
##        isparse[k]=[]
##        jsparse[k]=[]
##        csparse[k]=[]
##        for j in range(nPCTerms):
##            for i in range(nPCTerms):
##                wksum=0.0
##                for iqp in range(nclp):
##                    wksum += psi[iqp,i]* psi[iqp,j]* psi[iqp,k]*w[iqp]
##
##                if (wksum > tol):
##                    isparse[k].append(i)
##                    jsparse[k].append(j)
##                    csparse[k].append(wksum)
##        nsparse.append(len(isparse[k]))
##    pcdata['nsparse']=nsparse
##    pcdata['isparse']=isparse
##    pcdata['jsparse']=jsparse
##    pcdata['csparse']=csparse
##    return pcdata

def uq_PCBasis(pcdata,xi):
    if len(xi) <> pcdata['ndim']:
        print 'Input vector must be {} dimensional'.format(pcdata['ndim'])
        return 0
    psi=uq_psi(xi,pcdata['nord'],pcdata['pc_type'])
    ndim=pcdata['ndim']
    multiIndex=pcdata['multiIndex']
    nPCTerms=pcdata['nPCTerms']
    Psi=np.ones((1,nPCTerms))
    for i in range(ndim):
        Psi=Psi * psi[i,multiIndex[:nPCTerms,i]]
    return Psi
    
def uq_pcset(nord,ndim,pc_type):
    
    pcdata=dict()
    pcdata['ndim']=ndim
    pcdata['nord']=nord
    pcdata['nclp']=2*nord+1
    pcdata['pc_type']=pc_type

    x,w = uq_get1DNodesAndWeights(pcdata['nclp'],pc_type)

    pcdata['x']=x
    pcdata['w']=w

    pcdata['psi']=uq_psi(pcdata['x'], nord, pc_type)

    pcdata['multiIndex'],pcdata['nPCTerms'] = uq_initMultiIndex(nord,ndim)

    pcdata['psiMultiDsq']= uq_evalBasisNormsSquared(nord,ndim,pc_type)

    pcdata=uq_settripleprod(pcdata)
    
    return pcdata

def uq_evalpce(pcdata,pce,xi):
    Psi=uq_PCBasis(pcdata,xi)
    Upc=np.dot(pce,Psi)
    return Upc

    
def uq_sample(pcdata,pce,Nsamp):
    U=np.zeros(Nsamp)
    ndim=pcdata['ndim']
    pc_type=pcdata['pc_type']

    for s in range(Nsamp):
        if pc_type == 'LEGENDRE':
            xi=2.0*np.random.rand(ndim)-1.0
        elif pc_type == 'HERMITE':
            xi=np.random.randn(ndim)
        else:
            print ' Unknown quadrature'
            return
        Psi=uq_PCBasis(pcdata,xi)
        U[s]=np.dot(pce.T,Psi.T)
    return U

def uq_getNISP(pcdata,quadrature):
    ndim=pcdata['ndim']
    nord=pcdata['nord']
    pc_type=pcdata['pc_type']
    nPCTerms=pcdata['nPCTerms']
    psiMultiDsq=pcdata['psiMultiDsq']

    X=quadrature['nodes']
    W=quadrature['weights']
    nquad=quadrature['nquad']

    K= np.zeros((nPCTerms,nquad))
    for j in range(nquad):
        Psi_Xj=uq_PCBasis(pcdata,X[j,:])
        K[:,j]=W[j]*Psi_Xj / psiMultiDsq.T
    
    return K

def uq_meanvar(pcdata,pce):
    psiMultiDsq=pcdata['psiMultiDsq']
    mu=pce[:,1]
    A=pce[:,1:]**2
    b=psiMultiDsq[1:]
    sigma2=np.dot(A,b)
    return mu,sigma2


def uq_quadtable(nclp,ndim):
    nquad=nclp**ndim
    ind=np.zeros((nquad,ndim),dtype=int)
    ind[0,:]=0

    for m in range(1,nquad):
        for i in range(ndim):
            ind[m,i]=ind[m-1,i]

        nflag=0

        ind[m,0]=ind[m,0]+1

        if (ind[m,0] >nclp-1):
            nflag=1
            ind[m,0]=0

            for i in range(1,ndim):
                if(nflag >0):
                    nflag=0
                    ind[m,i]=ind[m,i]+1

                    if(ind[m,i]>nclp-1):
                        ind[m,i]=0
                        nflag=1

    return ind

def uq_quadrature(ndim,num1dnodes,pc_type):
    ind=uq_quadtable(num1dnodes,ndim)
    nquad=ind.shape[0]

    x,w = uq_get1DNodesAndWeights(num1dnodes,pc_type)
    X=np.zeros((nquad,ndim))
    for i in range(nquad):
        X[i,:]=x[ind[i,:]]

    W=np.zeros((nquad,1))
    for i in range(nquad):
        prodw=1
        for j in range(ndim):
            prodw=prodw*w[ind[i,j]]
        W[i]=prodw

    quadrature={}
    quadrature['nodes']=X
    quadrature['weights']=W
    quadrature['nquad']=nquad
    return quadrature


    
##def uq_initMultiIndex(nord,ndim):
##    maxP=200
##    multiIndex=np.zeros((maxP,ndim))
##    nnp=np.zeros(maxP)
##    nn=np.zeros(maxP)
##    ic=np.zeros(maxP)
##    ict=np.zeros(maxP)
##    
##    nPCTerms=0
##    nnp[nPCTerms]=0
##    for i in range(ndim):
##        multiIndex[nPCTerms,i]=0
##
##    nn[0]=1
##
##    if (nord>0):
##
##        for j in range(ndim):
##            nPCTerms =nPCTerms+ 1
##            nnp[nPCTerms]=1
##            for i in range(j-1):
##                multiIndex[nPCTerms,i]=0
##
##            multiIndex[nPCTerms,j]=1
##
##            for i in range(j+1,ndim):
##                multiIndex[nPCTerms,i]=0
##        
##        for i in range(ndim):
##            ic[i]=1
##
##        nn[1]=ndim
##
##        if (nord >1):
##            for k in range(1,nord):
##                IP=nPCTerms
##                for j in range(ndim):
##                    isum=0
##                    for l in range(ndim-1,j,-1):
##                        isum=isum+ic[l]
##                    ic[j]=isum
##
##                for j in range(ndim):
##                    ic[j]=ict[j]
##
##                print IP
##                for jj in range(ndim):
##                    for j in range(IP-ic[jj]+1,IP):
##                        nPCTerms=nPCterms+1
##                        nnp[nPCTerms]=k
##                        for i in range(ndim):
##                            multiIndex[nPCTerms,i]=multiIndex[j,i]
##                        multiIndex[nPCTerms,jj]=multiIndex[nPCTerms,jj]+1
##                nnp[k+1]=nPCTerms-IP
##    
##    multiIndex=multiIndex[:nPCTerms,:]
##    return multiIndex,nPCTerms    


def uq_settripleprod(pcdata):
    nord=pcdata['nord']
    nclp=pcdata['nclp']
    ndim=pcdata['ndim']
    multiIndex=pcdata['multiIndex']
    nPCTerms=pcdata['nPCTerms']
    w=pcdata['w']
    x=pcdata['x']
    psi=pcdata['psi']

    apow=np.zeros((nord+1,nord+1,nord+1))
    tol=1e-4
    for k in range(nord+1):
        for j in range(nord+1):
            for i in range(nord+1):
                sum=0
                for m in range(nclp):
                    sum= sum+psi[m,i]*psi[m,j]*psi[m,k]*w[m]
                apow[i,j,k]=sum
    apow[apow<tol]=0.0

    nsize=1000
    nsparse=np.zeros(nPCTerms,dtype=np.int)
    isparse=np.zeros((nPCTerms,nsize),dtype=np.int)
    jsparse=np.zeros((nPCTerms,nsize),dtype=np.int)
    csparse=np.zeros((nPCTerms,nsize))

    
    for k in range(nPCTerms):
        for j in range(nPCTerms):
            for i in range(nPCTerms):
                aprod=1
                for m in range(ndim):
                    l1=multiIndex[k,m]
                    l2=multiIndex[j,m]
                    l3=multiIndex[i,m]
                    aprod=aprod*apow[l1,l2,l3]
                if (aprod > tol):
                    nsparse[k]=nsparse[k]+1
                    isparse[k,nsparse[k]-1]=i
                    jsparse[k,nsparse[k]-1]=j
                    csparse[k,nsparse[k]-1]=aprod

    nsize=nsparse.max()
    isparse=isparse[:,:nsize]
    jsparse=jsparse[:,:nsize]
    csparse=csparse[:,:nsize]
    pcdata['nsparse']=nsparse
    pcdata['isparse']=isparse
    pcdata['jsparse']=jsparse
    pcdata['csparse']=csparse
        
    return pcdata


def uq_product(pcdata, pce1, pce2):
    nsparse=pcdata['nsparse']
    isparse=pcdata['isparse']
    jsparse=pcdata['jsparse']
    csparse=pcdata['csparse']
    psiMultiDsq=pcdata['psiMultiDsq']

    pce=np.zeros(pce1.shape)

    # perform sparse galerkin product
    for k in range(pcdata['nPCTerms']):
        sumpc=0.0
        for l in range(nsparse[k]):
            sumpc=sumpc+pce1[isparse[k,l]]*pce2[jsparse[k,l]]*csparse[k,l]
        pce[k]=sumpc/psiMultiDsq[k]
    return pce

def uq_div(pcdata, pce1, pce2):
    nsparse=pcdata['nsparse']
    isparse=pcdata['isparse']
    jsparse=pcdata['jsparse']
    csparse=pcdata['csparse']
    psiMultiDsq=pcdata['psiMultiDsq']
    nPCTerms=pcdata['nPCTerms']

    A=np.zeros((nPCTerms,nPCTerms))
    for k in range(nPCTerms):
        for l in range(nsparse[k]):
            m1=isparse[k,l]
            m2=jsparse[k,l]
            c=csparse[k,l]/psiMultiDsq[k]
            A[k,m1]=A[k,m1]+pce2[m2]*c
    pce = np.linalg.solve(A,pce1)
    return pce


def uq_log(pcdata, f):
    x=np.zeros(f.shape)
    u=x.copy()
    x[0]=f[0]
    x0=x
    u[0]=np.log(x[0])

    nsteps=2**8
    dx=(f-x)/nsteps

    for i in range(nsteps):
        dxox=uq_div(pcdata,dx,x)
        u=u+dxox
        x=x0+i*dx
    return u

def uq_exp(pcdata, f):
    x=np.zeros(f.shape)
    u=x.copy()
    x[0]=f[0]
    u[0]=np.exp(x[0])
    
    nsteps=2**8
    dx=(f-x)/nsteps

    for i in range(nsteps):
        udx=uq_product(pcdata,u,dx)
        u=u+udx
    return u

    
    
##### tests

def test_uq_exp():
    pc = uq_pcset(2,1,'LEGENDRE')
    u=np.array([1.,2.,3.])
    v=np.array([4.,5.,6.])
    logu=np.array([15.0466,30.5587,38.0566])
    logv=np.array([12245.2636365558,27738.8213183231,31940.6081648767])
    pclogu=uq_exp(pc, u)
    pclogv=uq_exp(pc, v)
    np.testing.assert_allclose(pclogu,logu,atol=0.01)
    np.testing.assert_allclose(pclogv,logv,atol=0.01)
    
def test_uq_log():
    pc = uq_pcset(2,1,'LEGENDRE')
    u=np.array([1.,2.,3.])
    v=np.array([4.,5.,6.])
    logu=np.array([-2.3007,1.9222,4.2905])
    logv=np.array([1.0108,0.9808,1.2261])
    pclogu=uq_log(pc, u)
    pclogv=uq_log(pc, v)
    np.testing.assert_allclose(pclogu,logu,atol=0.1)
    np.testing.assert_allclose(pclogv,logv,atol=0.01)


def test_uq_product():
    pc = uq_pcset(2,1,'LEGENDRE')
    u=np.array([1.,2.,3.])
    v=np.array([4.,5.,6.])
    expectedres=np.array([10.9333,23.8000,29.8095])
    exp=uq_product(pc, u, v)
    np.testing.assert_almost_equal(expectedres,exp,decimal=4)
    

def test_uq_div():
    pc = uq_pcset(2,1,'LEGENDRE')
    u=np.array([1.,2.,3.])
    v=np.array([4.,5.,6.])
    expectedres=np.array([ 0.0882,0.1326,0.3550])
    exp=uq_div(pc, u, v)
    np.testing.assert_almost_equal(expectedres,exp,decimal=4)    
    


def test_uq_quadtable():
    ind1d=array([[ 0.],
       [ 1.],
       [ 2.],
       [ 3.],
       [ 4.]])
    ind=uq_quadtable(5,1)
    assert(np.array_equal(ind1d,ind) == True)
    ind=uq_quadtable(5,2)
    assert(ind.shape ==(25,2))
    assert(np.array_equal(ind[-1,:],np.array([5,5])) == True)

def test_uq_quadrature():
    dd=uq_quadrature(1,6,'LEGENDRE')
    wht=np.array([[ 0.08566225],
      [ 0.18038079],
      [ 0.23395697],
      [ 0.23395697],
      [ 0.18038079],
      [ 0.08566225]])
    node= np.array([[-0.93246951],
       [-0.66120939],
       [-0.23861919],
       [ 0.23861919],
       [ 0.66120939],
       [ 0.93246951]])
    
    
    assert(np.array_equal(dd['nodes'],node) == True)
    assert(np.array_equal(dd['weights'],wht) == True)
    assert(dd['nquad']==6)    

def test_pcset():
    a = uq_pcset(2,1,'LEGENDRE')
    psi2 = np.array([
       [ 1.        , -0.90617985,  0.73174287],
       [ 1.        , -0.53846931, -0.0650762 ],
       [ 1.        ,  0.        , -0.5       ],
       [ 1.        ,  0.53846931, -0.0650762 ],
       [ 1.        ,  0.90617985,  0.73174287]])

    psinorm = np.array([
       [ 1.        ],
       [ 0.33333333],
       [ 0.2       ]])

    assert(np.array_equal(psi2,a['psi']) == True)
    assert(a['w'].sum() == 1.0)
    assert(np.array_equal(psinorm,a['psiMultiDsq']) == True)


def test_simple():
    nsamples=5000
    f=np.random.standard_normal(nsamples)
    Y=f**2

    polyOrder=2
    nvar=1
    ndim=nvar
    pc_type='HERMITE'
    pc=uq_pcset(polyOrder,nvar,pc_type)

    quadrature=uq_quadrature(1,3,pc_type)
    
    K=uq_getNISP(pc,quadrature)

    Ygrid=quadrature['nodes']**2

    c=np.dot(K,Ygrid)

    U=uq_sample(pc,c,nsamples)

    print 'Mean: ',Y.mean(),U.mean()
    print 'St dev: ',Y.std(),U.std()

    print 'Mean ',c[0]
    print 'stddev ',np.sqrt(c.sum()) 


