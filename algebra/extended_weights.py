"""Extended NSF P/W coefficients from independently supplied amplitudes.

Does not reconstruct missing physical kernels. Arrays have axes
(momentum_nodes, final_helicities, four_initial_helicities). Quadrature weights
must already contain the final-state measure, identical-particle factor,
sampling-density correction and averaging over momentum nodes. No extra
(2*pi)^3 or Bose factor is inserted here. Initial helicities are averaged.
"""
from itertools import product
from math import factorial
import numpy as np

class MissingPhysicalInput(ValueError):
    pass

def compositions(n):
    if n==0:
        yield ()
    else:
        for first in range(1,n+1):
            for tail in compositions(n-first):
                yield (first,)+tail

def exponential_coefficient(generators,order):
    """S_n in exp(i sum epsilon^j T_j); ordered noncommuting matrices.

All T_1,...,T_n must be explicitly supplied. Explicit zero matrices are
allowed only when independently justified by the caller. Finite matrices are
an algebraic interface, not a discretization prescription for NSF loops.
"""
    absent=[j for j in range(1,order+1) if generators.get(j) is None]
    if absent:
        raise MissingPhysicalInput('Missing generators: '+str(absent))
    shape=np.asarray(generators[1]).shape
    if len(shape)!=2 or shape[0]!=shape[1]:
        raise ValueError('Generators must be square matrices.')
    matrices={j:np.asarray(generators[j],complex) for j in range(1,order+1)}
    if any(v.shape!=shape for v in matrices.values()):
        raise ValueError('Use a common Hilbert-space matrix representation.')
    out=np.zeros(shape,complex)
    for word in compositions(order):
        term=np.eye(shape[0],dtype=complex)
        for j in word:term=term@matrices[j]
        out+=(1j**len(word)/factorial(len(word)))*term
    return out

def checked_inputs(amplitudes,energy,weights,Ein,source_record):
    if not source_record:
        raise MissingPhysicalInput('Describe the independent amplitude sources and prescriptions.')
    missing=[name for name,value in amplitudes.items() if value is None]
    if missing:
        raise MissingPhysicalInput('Missing amplitudes: '+', '.join(missing))
    arrays={name:np.asarray(value,complex) for name,value in amplitudes.items()}
    shape=next(iter(arrays.values())).shape
    if len(shape)!=3 or shape[-1]!=4 or any(v.shape!=shape for v in arrays.values()):
        raise ValueError('Common shape (nodes, final helicities, 4 initial helicities) required.')
    E=np.asarray(energy,float);w=np.asarray(weights,float)
    if E.shape!=(shape[0],) or w.shape!=E.shape or Ein<=0:
        raise ValueError('Energy/weight arrays must match nodes; Ein must be positive.')
    if any(not np.isfinite(a).all() for a in list(arrays.values())+[E,w]):
        raise ValueError('Nonfinite samples are not discarded or replaced by zero.')
    if np.any(E<=0) or np.any(w<0):
        raise ValueError('Positive physical energies and nonnegative integration weights required.')
    return arrays,E,w

def square(a):return np.sum(abs(a)**2,axis=(1,2))/4
def interference(a,b):return 2*np.real(np.sum(a.conj()*b,axis=(1,2)))/4

def moments(density,E,w,Ein):
    """Density may be signed for an interference coefficient; never clip it."""
    return {'P':float(w@density),'W':float(w@(E*density)),
            'D':float(w@((E/Ein-1)*density))}

def channel22(a2,a4,energy,weights,*,Ein=1.,source_record):
    """P2=eps^4 P[4]+eps^6 P[6]+O(eps^8), similarly W2/D2.

    a2=<f|S2|i> connected; a4=<f|S4|i> connected, including the required
    internal momentum integrations. An external-sum tree pool is insufficient
    to evaluate a physical loop in S4.
    """
    a,E,w=checked_inputs({'a2_S2':a2,'a4_S4':a4},energy,weights,Ein,source_record)
    return {'channel':22,'Ein':Ein,'orders':{
        4:moments(square(a['a2_S2']),E,w,Ein),
        6:moments(interference(a['a2_S2'],a['a4_S4']),E,w,Ein)},
        'source_record':source_record,'remainder_probability_order':8}

def channel23(iteration,mixed12,mixed21,contact,energy,weights,*,Ein=1.,source_record,a5=None):
    """Complete leading P3/W3 at eps^6; optional interference at eps^8.

    Every leading block is required explicitly. In the compositional BCH
    model the caller can supply certified zero mixed blocks, documenting that
    scope in source_record. The unknown quintic contact is never inferred zero.
    a5, when supplied, is the independent complete <f|S5|i>.
    """
    a,E,w=checked_inputs({'iteration':iteration,'mixed12':mixed12,'mixed21':mixed21,
                          'contact_T3':contact},energy,weights,Ein,source_record)
    A=a['iteration'];B=a['mixed12']+a['mixed21'];C=a['contact_T3'];full=A+B+C
    density=square(full)
    terms={'AA':square(A),'BB':square(B),'CC':square(C),
           'AB':interference(A,B),'AC':interference(A,C),'BC':interference(B,C)}
    if not np.allclose(sum(terms.values()),density,rtol=1e-11,atol=1e-24):
        raise ArithmeticError('Coherent square does not match six interference blocks.')
    result={'channel':23,'Ein':Ein,'orders':{6:moments(density,E,w,Ein)},
            'six_blocks':{k:moments(v,E,w,Ein) for k,v in terms.items()},
            'source_record':source_record,'remainder_probability_order':8}
    if a5 is not None:
        arr,E,w=checked_inputs({'a3_S3':full,'a5_S5':a5},energy,weights,Ein,source_record)
        result['orders'][8]=moments(interference(full,arr['a5_S5']),E,w,Ein)
        result['remainder_probability_order']=10
    return result

def combine_through_six(c22,c23,epsilon=1.):
    """Common probability order six, retaining the final ratio unexpanded.

    Including the optional order-eight 23 alone would not complete order eight
    for the sum: the 22 order-eight coefficient also requires S6.
    """
    if c22['Ein']!=c23['Ein']:
        raise ValueError('Common incoming energy required.')
    if set(c22['orders'])!={4,6} or 6 not in c23['orders']:
        raise MissingPhysicalInput('Need 22 orders 4,6 and complete 23 order 6.')
    v2={k:sum(epsilon**n*row[k] for n,row in c22['orders'].items()) for k in ['P','W','D']}
    v3={k:epsilon**6*c23['orders'][6][k] for k in ['P','W','D']}
    total={k:v2[k]+v3[k] for k in ['P','W','D']}
    positive=total['P']>0 and v2['P']>=0 and v3['P']>=0
    r=total['D']/total['P'] if total['P']>0 else None
    return {'22':v2,'23':v3,'total':total,'r_E':r,
            'mean_E_over_Ein':None if r is None else 1+r,
            'positive_integrated_channel_weights':positive,
            'ratio_expanded':False,'probability_order':6,
            'warning':None if positive else 'Signed truncation has nonpositive channel/total weight; no clipping performed.'}
