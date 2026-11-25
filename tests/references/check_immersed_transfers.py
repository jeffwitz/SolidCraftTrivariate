"""Check values and physical derivatives of the immersed reference transfers."""
import json
import numpy as np
from scipy.interpolate import BSpline
from immersed_multigrid import axis_transfer

if __name__=="__main__":
    records=[]
    for n,m,length in [(24,12,60.),(40,20,50.),(60,30,60.),(25,5,50.)]:
        p=axis_transfer(n,m,length)
        x=np.linspace(0,length,8*n+1)
        c=np.random.default_rng(20261125+n).normal(size=m+2)
        def spline(k,coeff):
            knots=np.r_[[0.,0.],np.linspace(0,length,k+1),[length,length]]
            return BSpline(knots,coeff,2)
        coarse=spline(m,c)
        fine=spline(n,p@c)
        error=[float(np.max(abs(fine(x,nu=j)-coarse(x,nu=j)))) for j in (0,1)]
        assert max(error)<1e-12
        records.append({"coarse":m,"fine":n,"length_mm":length,"max_value_error":error[0],"max_derivative_error":error[1]})
    print(json.dumps(records,indent=2))
