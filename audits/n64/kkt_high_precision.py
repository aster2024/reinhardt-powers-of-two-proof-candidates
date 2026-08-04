import mpmath as mp, importlib.util
mp.mp.dps=100
N=64
code='-++++++-----+--+-+++--+---+--+++---++-+++-++---+-++-+++++------+'
c=[mp.mpf(1) if z=='+' else mp.mpf(-1) for z in code]
a=[mp.mpf(0)]+[c[j-1]-c[j] for j in range(1,N)]
theta=mp.pi/N
# initial from certified axial construction
S_INTS=[113553801169565658142384194265, -9880196998990015085044055410, -133314195167545688312472305084, -256748193336101361539900554759, -380182191504657034767328804433, -503616189673212707994757054108, -627050187841768381222185303782, -527080630656387284719934645129, -427111073471006188217683986476, -327141516285625091715433327823, -227171959100243995213182669170, -127202401914862898710932010517, -224518637743365449610080769943, -131254926381658612514512915081, -37991215019951775418945060220, -120535648883632559553294480236, -35302307098282933183053776754, -109412180028470470567041786873, -183522052958658007951029796992, -257631925888845545335017807111, -199919646362187568673185756696, -142207366835529592011353706281, -196345090812107605367512622967, -149035276097607335865970026399, -101725461383107066364427429832, -54415646668606796862884833264, -75982724031176814146623311785, -39897017973167708772317662954, -3811311915158603398012014123, -2540874610105735598674676082, -1270437305052867799337338041]
t0=mp.mpf('0.9999999999878733515')
s=[mp.mpf('0')]*(N+1)
for j,v in enumerate(S_INTS,1): s[j]=t0*mp.mpf(v)/mp.mpf(10)**40
s[32]=0
for j in range(33,N): s[j]=-s[N-j]
y1=y2=mp.mpf('0')

def eval_system(s,y1,y2):
    phi=[j*theta+s[j] for j in range(N+1)]
    alpha=[phi[j+1]-phi[j] for j in range(N)]
    grad=[None]*N
    for j in range(1,N):
        grad[j]=mp.cos(alpha[j]/2)-mp.cos(alpha[j-1]/2) + a[j]*(-y1*mp.sin(phi[j])+y2*mp.cos(phi[j]))
    G=sum(c[j]*(mp.e**(1j*phi[j+1])-mp.e**(1j*phi[j])) for j in range(N))
    r=mp.matrix(N+1,1)
    for j in range(1,N): r[j-1]=grad[j]
    r[N-1]=mp.re(G); r[N]=mp.im(G)
    return r,phi,alpha,G

for it in range(12):
    r,phi,alpha,G=eval_system(s,y1,y2)
    maxr=max(abs(r[k]) for k in range(N+1))
    print('iter',it,'maxres',mp.nstr(maxr,12))
    if maxr < mp.mpf('1e-85'): break
    J=mp.matrix(N+1,N+1)
    # Hessian of f=-F: path matrix with edge weights .5 sin(alpha/2)
    w=[mp.mpf('0.5')*mp.sin(alpha[j]/2) for j in range(N)]
    for j in range(1,N):
        row=j-1
        # diagonal from adjacent gaps + multiplier curvature
        J[row,j-1]=w[j-1]+w[j] + a[j]*(-y1*mp.cos(phi[j])-y2*mp.sin(phi[j]))
        if j>1: J[row,j-2]=-w[j-1]
        if j<N-1: J[row,j]=-w[j]
        J[row,N-1]=a[j]*(-mp.sin(phi[j]))
        J[row,N]=a[j]*mp.cos(phi[j])
    # closure Jacobian rows
    for j in range(1,N):
        J[N-1,j-1]=a[j]*(-mp.sin(phi[j]))
        J[N,j-1]=a[j]*mp.cos(phi[j])
    dz=mp.lu_solve(J,-r)
    for j in range(1,N): s[j]+=dz[j-1]
    y1+=dz[N-1]; y2+=dz[N]

r,phi,alpha,G=eval_system(s,y1,y2)
F=sum(2*mp.sin(x/2) for x in alpha)
U=2*N*mp.sin(mp.pi/(2*N))
print('max residual',mp.nstr(max(abs(r[k]) for k in range(N+1)),30))
print('closure',mp.nstr(G,30))
print('perimeter',mp.nstr(F,100))
print('upper deficit',mp.nstr(U-F,100))
print('y norm',mp.nstr(mp.sqrt(y1*y1+y2*y2),40))
print('min/max gaps',mp.nstr(min(alpha),30),mp.nstr(max(alpha),30))
print('axial symmetry max',mp.nstr(max(abs(s[N-j]+s[j]) for j in range(N+1)),30))
