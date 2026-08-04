import numpy as np
N=16
th=np.pi/N
z=np.exp(1j*np.arange(N+1)*th)
e=z[1:]-z[:-1]
R=.0065; C=26.1
j=np.arange(1,N)
jj,kk=np.meshgrid(j,j,indexing='ij')
Dinv=np.minimum(jj,kk)*(N-np.maximum(jj,kk))/N
surv=[]; margins=[]
for bits in range(1<<(N-1)):
    c=np.ones(N)
    for k in range(1,N):
        if (bits>>(k-1))&1: c[k]=-1
    b=abs(np.dot(c,e))
    a=c[:-1]-c[1:]
    B=np.vstack([-a*np.sin(j*th),a*np.cos(j*th)])
    sigma=np.sqrt(np.max(np.linalg.eigvalsh(B@Dinv@B.T)))
    mar=b-R*sigma-C*R*R
    code=''.join('+' if x>0 else '-' for x in c)
    if mar<=1e-13: surv.append((code,b,sigma,mar))
    else: margins.append((mar,code,b,sigma))
print('survivors',len(surv))
for x in sorted(surv): print(x)
print('worst excluded',min(margins))
REP='+--+-++-+--+-++-'
def neg(s): return ''.join('+' if ch=='-' else '-' for ch in s)
def orbit(code):
 full=code+neg(code); out=set()
 for base in (full,full[::-1]):
  for r in range(2*N):
   h=(base[r:]+base[:r])[:N]
   if h[0]=='-': h=neg(h)
   out.add(h)
 return out
print('orbit size',len(orbit(REP)),'equal',set(c for c,*_ in surv)==orbit(REP))
