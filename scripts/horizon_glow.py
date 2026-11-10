"""Map every 48-pixel block of a plate-solved frame to altitude and azimuth, for reading horizon glow.

Usage: horizon_glow.py FRAME.NEF SOLVE_OUTPUT.txt OUTPUT_DIR LAT LON
SOLVE_OUTPUT.txt holds the solve_sky.py output line for the frame.
"""
import sys,re,math,json
import os
sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
sys.dont_write_bytecode=True
import numpy as np, rawpy
import solve_sky as ss
from scipy.spatial import cKDTree
path=sys.argv[1]; SOLVE=sys.argv[2]; S=sys.argv[3]; LAT,LON=float(sys.argv[4]),float(sys.argv[5]); ss.LAT,ss.LON=LAT,LON
n=re.search(r'DSC_(\d+)',path)[1]
os.makedirs(S,exist_ok=True)
line=[l for l in open(SOLVE) if 'DSC_%s'%n in l][0]
m=re.search(r'alt\s+([\d.-]+) az\s+([\d.-]+) roll\s+([\d.-]+)',line); alt0,az0,roll0=[float(v) for v in m.groups()]
ex=ss.exif(path); jd=ss.julian_day(*ex['utc'])
cells,_,_=ss.green_cells(path); h,w=cells.shape; c=(w/2.,h/2.); fc=16/(2*35.9/6016)
ra,dec,mag=ss.load_catalog(); alt,az=ss.star_altaz(ra,dec,jd); br=(mag<5.6)&(alt>2); v=ss.unit_vectors(alt[br],az[br])
src,_,_=ss.detect_sources(cells); tree=cKDTree(src)
best=None
for hd in (1.,-1.):
    x,y,d=ss.project(v,alt0,az0,roll0,fc,hd); ins=d>.3
    dist,_=tree.query(np.stack([x[ins]+c[0],y[ins]+c[1]],1),distance_upper_bound=14)
    s=int(np.isfinite(dist).sum())
    if best is None or s>best[0]: best=(s,hd)
hd=best[1]
right,up,view=ss.camera_basis(alt0,az0,roll0,hd)
r=rawpy.imread(path); im=r.raw_image_visible.astype(np.float32)-float(r.black_level_per_channel[0]); pat=r.raw_pattern
idx={}
for a in (0,1):
    for b in (0,1): idx.setdefault(chr(r.color_desc[pat[a][b]]),[]).append((a,b))
G=.5*(im[idx['G'][0][0]::2,idx['G'][0][1]::2]+im[idx['G'][1][0]::2,idx['G'][1][1]::2]); R=im[idx['R'][0][0]::2,idx['R'][0][1]::2]; B=im[idx['B'][0][0]::2,idx['B'][0][1]::2]
hh=min(G.shape[0],R.shape[0],B.shape[0]); ww=min(G.shape[1],R.shape[1],B.shape[1]); G,R,B=G[:hh,:ww],R[:hh,:ww],B[:hh,:ww]
bl=48; rows,cols=hh//bl,ww//bl
bm=lambda a: np.median(a[:rows*bl,:cols*bl].reshape(rows,bl,cols,bl),axis=(1,3))
g,rr,bb=bm(G),bm(R),bm(B)
yy,xx=np.mgrid[0:rows,0:cols]; dx=(xx+.5)*bl-c[0]; dy=(yy+.5)*bl-c[1]
th=np.arctan(np.hypot(dx,dy)/fc)
ray=view[None,None,:]*fc+right[None,None,:]*dx[...,None]-up[None,None,:]*dy[...,None]; ray/=np.linalg.norm(ray,axis=-1,keepdims=True)
balt=np.degrees(np.arcsin(ray[...,2])); baz=np.degrees(np.arctan2(ray[...,0],ray[...,1]))%360
gc=g/np.cos(th)**4
np.savez(S+'/%s.npz'%n,alt=balt,az=baz,g=gc,r=rr/np.cos(th)**4,b=bb/np.cos(th)**4,th=np.degrees(th),raw=g)
# preview
s=8; G2=G[::s,::s]; R2=R[::s,::s]; B2=B[::s,::s]
def st(a): return np.clip(a/900.,0,1)**0.45
rgb=np.stack([st(R2*1.9),st(G2),st(B2*1.6)],-1)
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
fig,ax=plt.subplots(figsize=(13,8.7)); ax.imshow(rgb,extent=(0,ww,hh,0),interpolation='bilinear')
X=(xx+.5)*bl; Y=(yy+.5)*bl
cs=ax.contour(X,Y,balt,levels=[0,5,10,20,30,45,60,75],colors='cyan',linewidths=.6); ax.clabel(cs,fmt='%d',fontsize=8)
a2=np.where(baz>180,baz-360,baz)
cz=ax.contour(X,Y,baz,levels=list(range(0,360,30)),colors='yellow',linewidths=.5); ax.clabel(cz,fmt='%d',fontsize=8)
ax.set_title('DSC_%s  alt %.1f az %.1f roll %.1f hd %+d stars %d'%(n,alt0,az0,roll0,hd,best[0])); ax.axis('off')
fig.tight_layout(); fig.savefig(S+'/%s.jpg'%n,dpi=95); print(n,alt0,az0,roll0,hd,best[0],'alt range',round(float(balt.min())),round(float(balt.max())),flush=True)
