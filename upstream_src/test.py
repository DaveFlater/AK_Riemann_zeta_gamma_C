import numpy as np
from Riemann_zeta_python import Riemann_zeta, Riemann_zeta_prime
from ln_gamma import ln_gamma
import matplotlib.pyplot as plt


print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
print('Test 1 : compute zeta(s) for s=2,4,6,8,10')
f=np.zeros(10,dtype=np.complex128)
for k in range(1,6):
    s=2*k
    f[k-1]=Riemann_zeta(s)
    print(f[k-1])
g=np.zeros(10,dtype=np.complex128)
g[0]=np.pi**2/6
g[1]=np.pi**4/90
g[2]=np.pi**6/945
g[3]=np.pi**8/9450
g[4]=np.pi**10/93555
print('maximum error=',np.max(np.abs(f[:4]-g[:4])))

print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
print('Test 2 : compute zeta(s) for s=0,-1,-3,-5,-7,-9')
f=np.zeros(10,dtype=np.complex128)
s=0
f[0]=Riemann_zeta(s)
print(f[0])
for k in range(1,6):
    s=1-2*k
    f[k]=Riemann_zeta(s)
    print(f[k])
g=np.zeros(10,dtype=np.complex128)
g[0]=-1.0/2
g[1]=-1.0/12
g[2]=1.0/120
g[3]=-1.0/252
g[4]=1.0/240
g[5]=-1.0/132
print('maximum error=',np.max(np.abs(f[:6]-g[:6])))


# the zeros of zeta(s) used in the next two tests
# were downloaded from https://www-users.cse.umn.edu/~odlyzko/zeta_tables/index.html
print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
print('Test 3 : compute zeta(1/2+i*t_k) for the first five non-trivial zeros 1/2+i*t_k of zeta(s)')
t=np.zeros(5,dtype=np.float64)
t[0]=14.13472514173469379045725198356247027
t[1]=21.02203963877155499262847959389690277
t[2]=25.01085758014568876321379099256282181
t[3]=30.42487612585951321031189753058409132
t[4]=32.93506158773918969066236896407490348
print(Riemann_zeta(0.5+1j*t))

print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
print('Test 4 : compute zeta(1/2+i*t_k) for the first five non-trivial zeros 1/2+i*t_k of zeta(s) at height t=10^6')
print('note that these zeros are only correct to 10^{-9}, so the values of zeta(1/2+i*t_k) should be of the same order of magnitude')
t[0]=1000000.584097696
t[1]=1000000.828343490
t[2]=1000001.435265267
t[3]=1000001.905648406
t[4]=1000002.877617740
print(Riemann_zeta(0.5+1j*t))

    
# the zeros of zeta(s) used in the next test
# were downloaded from https://www.lmfdb.org/zeros/zeta/
print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
print('Test 5 : compute zeta(1/2+i*t_k) for the first five non-trivial zeros 1/2+i*t_k of zeta(s) at height t=10^9')
t[0]=1000000000.115650890
t[1]=1000000000.434026895
t[2]=1000000000.530342856
t[3]=1000000001.019020177
t[4]=1000000001.293739492
print(Riemann_zeta(0.5+1j*t))



print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
print('Test 6 : compare the vaules of zeta(s_j) with the values pre-computed using Fortran quadruple precision implementation')
s=np.zeros(10,dtype=np.complex128)
zeta_exact=np.zeros(10,dtype=np.complex128)
s[0]=-15.5+1.0j
s[1]=0.5+20.0j
s[2]=1.5-150.0j
s[3]=0+500.0j
s[4]=10+3000.0j
s[5]=1-100000.0j
s[6]=-0.5+5000000.0j
s[7]=0.5+20000000.0j
s[8]=0.75+1.0e+10j
s[9]=0.4-1.0e+12j
zeta_exact[0]=1.607209928015388936363-0.314611728846899873791j
zeta_exact[1]=0.429913860437843372157-1.064291443080589112727j
zeta_exact[2]=0.643164286756546270307+0.138652708369380808021j
zeta_exact[3]=-9.82537217473325310405-2.187513665314884630833j
zeta_exact[4]=1.000919537804919266834+2.874631032437997010785e-4j
zeta_exact[5]=1.618122122846936796567-1.070441041470623686626j
zeta_exact[6]=1537232.276423317013525-590208.9908762055608925j
zeta_exact[7]=0.281689403857957175927+2.350684988978597606149j
zeta_exact[8]=0.315056654699321886152+0.561446548967376614715j
zeta_exact[9]=14.31909665138080782830+17.85833478123989396854j
print(Riemann_zeta(s)-zeta_exact)

print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
print("Test 7: computing the first non-trivial zero with Newton's method")
s = 0.5 + 14.1j # initial approximation
for i in range(5):
    f, f1 = Riemann_zeta_prime(s)
    s = s - f / f1
    print(f'  s = {s}')
s1 = 0.5 + 14.13472514173469379j    
print(f'  error = {s - s1}')

print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
print("Test 8: computing the non-trivial zero close to 1/2+401.8i with Newton's method")
s = 0.5 + 401.8j # initial approximation
for i in range(5):
    f, f1 = Riemann_zeta_prime(s)
    s = s - f / f1
    print(f'  s = {s}')
s1 = 0.5 + 401.8392286005332165399113j    
print(f'  error = {s - s1}')


print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
print('Test 9 : plot the Riemann-Siegel function Z(t) for 0<t<100, 1000<t<1100 and 10000<t<10100')
t=np.linspace(0,100,10000)
g=ln_gamma(0.25+0.5j*t)
f0=Riemann_zeta(0.5+1j*t)*np.exp(1j*g.imag-0.5j*t*np.log(np.pi))

t1=t+1000
g1=ln_gamma(0.25+0.5j*t1)
f1=Riemann_zeta(0.5+1j*t1)*np.exp(1j*g1.imag-0.5j*t1*np.log(np.pi))

t2=t+10000
g2=ln_gamma(0.25+0.5j*t2)
f2=Riemann_zeta(0.5+1j*t2)*np.exp(1j*g2.imag-0.5j*t2*np.log(np.pi))

fig,axes=plt.subplots(3,1,figsize=(8,12))  
axes[0].plot(t,f0.real)
axes[0].set_title('Z(t) for 0<t<100')
axes[0].grid(True)
axes[1].plot(t1,f1.real)
axes[1].set_title('Z(t) for 1000<t<1100')
axes[1].grid(True)
axes[2].plot(t2,f2.real)
axes[2].set_title('Z(t) for 10000<t<10100')
axes[2].grid(True)
plt.tight_layout()
plt.show()

