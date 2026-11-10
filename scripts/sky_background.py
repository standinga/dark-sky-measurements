"""Sky background, colour ratios and star sharpness for every 30 s / ISO 800 NEF under a folder.

Usage: sky_background.py NEF_ROOT OUTPUT.json
"""
import sys,glob,os,subprocess,json
import numpy as np, rawpy
from scipy import ndimage as ndi
ROOT=sys.argv[1]
paths=sorted(glob.glob(ROOT+'/*/*.NEF')+glob.glob(ROOT+'/*/soft/*.NEF'))
out=subprocess.run(['exiftool','-T','-n','-FileName','-Directory','-DateTimeOriginal','-OffsetTimeOriginal','-ISO','-ExposureTime','-FNumber','-PitchAngle','-RollAngle']+paths,capture_output=True,text=True).stdout
rows=[]
for ln in out.strip().split('\n'):
    p=ln.split('\t')
    name,d,dt,off,iso,exp,fn,pitch,roll=p
    if float(iso)!=800 or abs(float(exp)-30)>1: continue
    if float(fn)>2.3: continue
    f=os.path.join(d,name)
    r=rawpy.imread(f); im=r.raw_image_visible.astype(np.float32)-float(r.black_level_per_channel[0]); pat=r.raw_pattern
    idx={}
    for a in (0,1):
        for b in (0,1):
            idx.setdefault(chr(r.color_desc[pat[a][b]]),[]).append((a,b))
    g=0.5*(im[idx['G'][0][0]::2,idx['G'][0][1]::2]+im[idx['G'][1][0]::2,idx['G'][1][1]::2])
    R=im[idx['R'][0][0]::2,idx['R'][0][1]::2]; B=im[idx['B'][0][0]::2,idx['B'][0][1]::2]
    h=min(g.shape[0],R.shape[0],B.shape[0]); w=min(g.shape[1],R.shape[1],B.shape[1]); g,R,B=g[:h,:w],R[:h,:w],B[:h,:w]
    block=64; rows_,cols_=h//block,w//block
    def bm(a): return np.median(a[:rows_*block,:cols_*block].reshape(rows_,block,cols_,block),axis=(1,3))
    bg,br,bb=bm(g),bm(R),bm(B)
    yy,xx=np.mgrid[0:rows_,0:cols_]; fc=float(16)/(2*35.9/6016)
    th=np.arctan(np.hypot((xx+.5)*block-w/2,(yy+.5)*block-h/2)/fc)
    res={}
    for c in (10,15):
        m=th<np.radians(c); v=bg[m]/np.cos(th[m])**4
        res[f's{c}']=float(np.percentile(v,25)); res[f'm{c}']=float(np.percentile(v,50))
    m=th<np.radians(15)
    res['rg']=float(np.median(br[m]/bg[m])); res['bg']=float(np.median(bb[m]/bg[m]))
    # sharpness: central crop, detect stars, second-moment width
    c=g[h//2-900:h//2+900,w//2-900:w//2+900]; bk=ndi.median_filter(c[::4,::4],size=15); bk=ndi.zoom(bk,4,order=1)[:c.shape[0],:c.shape[1]]
    d_=c[:bk.shape[0],:bk.shape[1]]-bk; sig=1.4826*np.median(np.abs(d_))
    lab,n=ndi.label(d_>6*sig)
    objs=ndi.find_objects(lab)[:5000]
    ws=[];pk=[]
    for i,sl in enumerate(objs):
        a=(lab[sl]==i+1).sum()
        if a<3: continue
        st=d_[max(sl[0].start-4,0):sl[0].stop+4,max(sl[1].start-4,0):sl[1].stop+4]; st=np.clip(st,0,None)
        pk.append(st.max())
        if st.max()>50*sig or st.max()<15*sig: continue
        y,x=np.indices(st.shape); t=st.sum(); xc=(x*st).sum()/t; yc=(y*st).sum()/t
        ws.append(np.sqrt(((((x-xc)**2+(y-yc)**2)*st).sum()/t)))
    res.update(n_det=int(n),n_w=len(ws),width=float(np.median(ws)) if ws else None)
    rows.append(dict(name=name,sess=os.path.relpath(d,ROOT),dt=dt,off=off,fn=float(fn),pitch=float(pitch),roll=float(roll),**res))
    print(name,rows[-1]['sess'],dt,f"s15 {res['s15']:.1f} s10 {res['s10']:.1f} rg {res['rg']:.3f} bg {res['bg']:.3f} det {n} wid {res['width'] if res['width'] is None else round(res['width'],2)}",flush=True)
json.dump(rows,open(sys.argv[2],'w'),indent=0)
