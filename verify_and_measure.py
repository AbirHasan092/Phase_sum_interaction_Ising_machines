"""Independent evaluation of main.m; writes results without requiring MATLAB."""
from pathlib import Path
import sys, csv
import numpy as np
root=Path(__file__).resolve().parent
d=np.loadtxt(sys.argv[1],delimiter=',',skiprows=1)
out=root/'results_csv_gain'; out.mkdir(exist_ok=True)
u=(d[:,0]-d[0,0])/(d[-1,0]-d[0,0]); H=(len(u)-1)//4
h=np.arange(-H,H+1); E=np.exp(-2j*np.pi*u[:,None]*h)
z=np.trapezoid(d[:,1,None]*E,u,axis=0); q=np.trapezoid(d[:,3,None]*E,u,axis=0)
lz=20*np.log10(np.maximum(abs(z),1e-300)/max(abs(z[h!=0])))
lq=20*np.log10(np.maximum(abs(q),1e-300)/max(abs(q[h!=0])))
z[lz < -100]=0; q[lq < -100]=0
coefficients=np.load(sys.argv[2])['G']
amp=coefficients[1:11]
G0=coefficients[0].real
Kmax=len(amp)
g=np.r_[np.conj(amp[::-1]),G0,amp]; N=2**int(np.ceil(np.log2(2*H+2)))
ax=np.arange(-N//2,N//2); n,m=np.meshgrid(ax,ax,indexing='ij')
Z=np.zeros(N,complex); Q=Z.copy(); Z[h+N//2]=z; Q[h+N//2]=q
k=-(n+m)/2; ok=((n+m)%2==0)&(abs(k)<=Kmax); G=np.zeros((N,N),complex); G[ok]=g[k[ok].astype(int)+Kmax]
base=Z[:,None]*Q[None,:]*G; reverse=Q[:,None]*Z[None,:]*G
inv=lambda c: np.fft.ifft2(np.fft.ifftshift(c)).real*N*N
rows=[]; cases={}
for M in range(2,9):
 target=2*np.pi*np.arange(M)/M; offsets=2*np.pi*np.arange(8)/8
 dist=abs(np.angle(np.exp(1j*(offsets[None,:]-target[:,None])))); p=np.zeros(M,int)
 for _ in range(M):
  # MATLAB min uses column-major ordering.
  ix=np.argmin(dist.ravel(order='F')); rr,j=np.unravel_index(ix,dist.shape,order='F')
  p[rr]=j; dist[rr,:]=np.inf; dist[:,j]=np.inf
 s=p[np.r_[0,np.arange(M-1,0,-1)]]
 for ideal in (0,1):
  A=((n-m)%M==0) if ideal else sum(np.exp(1j*(n*offsets[pr]+m*offsets[sr])) for pr,sr in zip(p,s))/M
  C=base*A; B=reverse*A; d1=inv(1j*m*C); d2=inv(1j*n*B); curl=d1-d2
  rms=np.sqrt(np.mean(curl**2)); norm=rms/max(np.sqrt(np.mean(d1**2+d2**2)),np.finfo(float).eps)
  pd=np.sum(abs(C[(n==m)&(n!=0)])**2); pc=np.sum(abs(C[n!=m])**2)
  # Independent periodic time averaging at several phase points, same retained signals.
  tau=2*np.pi*np.arange(2048)/2048; gt=G0+2*np.real(sum(a*np.exp(2j*(i+1)*tau) for i,a in enumerate(amp)))
  for a,b in [(0.13,.71),(1.1,2.3)]:
   delta=target if ideal else offsets[p]; epsilon=-target if ideal else offsets[s]
   direct=0.
   for dr,ds in zip(delta,epsilon):
    zz=np.real(np.exp(1j*(tau[:,None]+a+dr)*h)@z)
    qq=np.real(np.exp(1j*(tau[:,None]+b+ds)*h)@q)
    direct+=np.mean(gt*zz*qq)/M
   analytic=np.real(np.sum(C*np.exp(1j*(n*a+m*b))))
   assert abs(direct-analytic)<1e-9*max(1,np.sum(abs(C)))
  rows.append([M,ideal,pc/(pd+pc),rms,np.max(abs(curl)),norm])
  cases[f'M{M}_ideal{ideal}_curl']=curl
with (out/'sweep_summary.csv').open('w',newline='') as f:
 w=csv.writer(f);w.writerow(['M','IdealUniform','CrossPowerFraction','CurlRMS','CurlMax','NormalizedCurl']);w.writerows(rows)
positive=np.r_[G0,amp]
np.savetxt(out/'gain_harmonics.csv',np.c_[np.arange(Kmax+1),abs(positive),np.angle(positive)],delimiter=',',header='k,abs_Gk_S,phase_rad',comments='')
np.savetxt(out/'harmonics.csv',np.c_[h,abs(q),abs(z),lq,lz],delimiter=',',header='order,abs_Q_retained,abs_Z_retained,Q_dBc_raw,Z_dBc_raw',comments='')
np.savez(out/'independent_results.npz',h=h,q=q,z=z,g=g,summary=np.array(rows),**cases)
print(np.array(rows));print('retained maximum orders',max(abs(h[q!=0])),max(abs(h[z!=0])));print('Independent temporal averaging checks passed.')

