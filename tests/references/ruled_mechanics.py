"""Independent decimal volume moments and elastic work of the ideal swept plate."""
import argparse
from decimal import Decimal as D, localcontext
import json
from pathlib import Path

PI = D("3.14159265358979323846264338327950288419716939937510582097494")
CASES = [(8, 0, 8, 0, 0), (0, 8, 12, 4, -8), (-8, 1, 7, -16, 24)]


def measure(data, power, coefficient):
    a, b, height, _, _ = data
    h = D(height)
    real_change = D(a) / 8 - 1
    imag_change = D(b) / 8
    quadratic = real_change**2 + imag_change**2
    area = 400 - 16*PI
    modulus = D(210000)*(1-D(".3")) / ((1+D(".3"))*(1-2*D(".3")))

    def moment(k):
        # Integral of z^k over the complete solid from its exact section area.
        return area*h**(k+1)*(1/D(k+1) + 2*real_change/(k+2) + quadratic/(k+3))

    p, c = power, D(coefficient)
    energy = modulus*(p*c)**2*moment(2*p-2)/2
    body_z = -modulus*c*p*(p-1)*moment(p-2) if p>=2 else D(0)
    body_work = -modulus*c*c*p*(p-1)*moment(2*p-2)
    bottom = -modulus*c*area if p==1 else D(0)
    top = modulus*c*p*h**(p-1)*area*((1+real_change)**2+imag_change**2)
    values = dict(volume_mm3=moment(0),l2_squared_mm5=c*c*moment(2*p),
                  h1_seminorm_squared_mm3=(p*c)**2*moment(2*p-2),
                  strain_energy_n_mm=energy,external_work_n_mm=2*energy,
                  body_force_z_n=body_z,bottom_force_z_n=bottom,top_force_z_n=top,
                  lateral_force_z_n=-body_z-bottom-top,natural_force_z_n=-body_z-bottom,
                  body_work_n_mm=body_work,natural_work_n_mm=2*energy-body_work)
    return {"recipe":data, "power":p, "coefficient":coefficient,
            **{key:str(value) for key,value in values.items()}}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output",type=Path)
    args=parser.parse_args()
    with localcontext() as context:
        context.prec=52
        cases=[measure(case,p,c) for case in CASES for p,c in [(1,"0.001"),(3,"0.0001")]]
    args.output.write_text(json.dumps({"precision_decimal_digits":52,"cases":cases},indent=2)+"\n")


if __name__=="__main__":
    main()
