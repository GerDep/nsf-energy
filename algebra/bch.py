"""Canonical normal-polynomial algebra, independent of occupation action."""
from collections import Counter,defaultdict
from itertools import product
import math

def add(*ps):
    out=defaultdict(complex)
    for p in ps:
        for k,z in p.items():out[k]+=z
    return {k:z for k,z in out.items() if abs(z)>1e-15}

def scale(p,z):return {k:z*v for k,v in p.items()}

def mul(p,q):
    out=defaultdict(complex)
    for (cre,ann),a in p.items():
      for (cre2,ann2),b in q.items():
        A=Counter(ann);C=Counter(cre2);common=sorted(A.keys()&C.keys())
        for ks in product(*(range(min(A[i],C[i])+1) for i in common)):
            z=a*b;r=list(cre+cre2);s=list(ann+ann2)
            for i,k in zip(common,ks):
                z*=math.comb(A[i],k)*math.comb(C[i],k)*math.factorial(k)
                for j in range(k):r.remove(i);s.remove(i)
            out[(tuple(sorted(r)),tuple(sorted(s)))]+=z
    return add(out)

def comm(p,q):return add(mul(p,q),scale(mul(q,p),-1))

def norm(p):return math.sqrt(sum(abs(z)**2 for z in p.values()))

def a(i):return {((),(i,)):1.}

def outgoing4(T1,T2,T3,modes):
    result={};iterations={};alternate_errors=[]
    for i in range(modes):
        ai=a(i);a2=scale(comm(ai,T1),1j)
        a3=add(scale(comm(ai,T2),1j),scale(comm(comm(ai,T1),T1),-.5))
        I4=add(scale(comm(comm(ai,T1),T2),-.5),scale(comm(comm(ai,T2),T1),-.5),
               scale(comm(comm(comm(ai,T1),T1),T1),-1j/6))
        I4_alt=add(scale(add(comm(a2,T2),comm(a3,T1)),.5j),
                   scale(comm(comm(a2,T1),T1),1/12))
        alternate_errors.append(norm(add(I4,scale(I4_alt,-1))))
        iterations[i]=I4;result[i]=add(scale(comm(ai,T3),1j),I4)
    return result,iterations,max(alternate_errors)

def lift_residual(R,degree=4):
    """Canonical CCR=1 lift: -i/(r+1) c_i† R_i. CCR=1/2 uses -2i."""
    out=defaultdict(complex)
    for i,p in R.items():
        for (cre,ann),z in p.items():
            if len(cre)+len(ann)!=degree:continue
            out[(tuple(sorted(cre+(i,))),ann)]+=-1j*z/(len(cre)+1)
    return add(out)

def reconstruct_quintic(independent_out4,T1,T2,modes,tolerance=1e-11):
    """Requires canonical delta c4 from an INDEPENDENT NSF source.

Supplying its BCH iteration instead makes a zero residual by construction and
does not establish that the physical five-leg contact vanishes.
"""
    if independent_out4 is None:raise ValueError('Independent fourth-order outgoing coefficient is missing')
    _,iterations,_=outgoing4(T1,T2,{},modes)
    residual={i:{k:z for k,z in add(independent_out4[i],scale(iterations[i],-1)).items()
                    if len(k[0])+len(k[1])==4} for i in range(modes)}
    T3=lift_residual(residual)
    error=max(norm(add(scale(comm(a(i),T3),1j),scale(residual[i],-1))) for i in range(modes))
    adjoint={(ann,cre):z.conjugate() for (cre,ann),z in T3.items()}
    hermitian_error=norm(add(T3,scale(adjoint,-1)))
    vacuum_error=norm({k:z for k,z in T3.items() if not k[1]})
    scale0=max(1.,norm(T3),max(norm(p) for p in residual.values()))
    if max(error,hermitian_error,vacuum_error)>tolerance*scale0:
        raise ValueError('Fourth-order input fails the inverse-commutator, adjoint or vacuum checks')
    return T3,{'commutator_residual':error,'hermiticity_residual':hermitian_error,'vacuum_residual':vacuum_error}
