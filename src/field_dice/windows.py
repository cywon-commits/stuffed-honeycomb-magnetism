import numpy as np
from lswt import collinear_stable
def edge(sig,J,Jp,lo,hi,want,n=16,it=30):
    # bisection between lo (stable==want) and hi (stable!=want)
    for _ in range(it):
        m=(lo+hi)/2
        if collinear_stable(sig,m,J,Jp,n=n)==want: lo=m
        else: hi=m
    return (lo+hi)/2
def window(sig,J,Jp,hmax=8,n=16):
    hs=np.arange(0,hmax,0.05); st=[collinear_stable(sig,h,J,Jp,n=n) for h in hs]
    if not any(st): return None
    i=st.index(True); j=len(st)-1-st[::-1].index(True)
    a=hs[i] if i==0 else edge(sig,J,Jp,hs[i],hs[i-1],True,n)
    b=edge(sig,J,Jp,hs[j],hs[j+1],True,n)
    return a,b
print(" r    A-down UUD window      H-down UUD window   (S=1/2, J=1)")
for r in [0.25,0.5,0.75,0.9,1.0,1.1,1.25,1.5,1.75,2.0,2.5]:
    wa=window((-1,1,1),1,r); wh=window((1,1,-1),1,r)
    f=lambda w:'   none          ' if w is None else f'[{w[0]:.4f},{w[1]:.4f}]'
    print(f"{r:4.2f}  {f(wa)}   {f(wh)}")
