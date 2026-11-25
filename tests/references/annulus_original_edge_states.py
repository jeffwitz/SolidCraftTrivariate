"""Verify all three gallery solids against the original annular charts."""
from fractions import Fraction as F
import json
from pathlib import Path

import numpy as np
from annulus_prepared_volumes import volume

ROOT = Path(__file__).resolve().parents[2] / 'docs/site/data/models'


def eval_section(row, radial, axial):
    """Evaluate the section on one biquadratic span in exact arithmetic."""
    def basis(u):
        return ((1-u)**2, 2*u*(1-u), u*u)
    return [sum(b*c*row[i][j][k] for i,b in enumerate(basis(radial))
                for j,c in enumerate(basis(axial))) for k in range(4)]


def oracle(s, t, span, state):
    """Nominal (radius, height, denominator) before angular revolution."""
    if span == 1:
        return 4-2*s, F(1,4) + F(15,2)*t, F(2)
    if span == 2:
        t = 1-t
    den = 1+t*t
    if state == 'fillet':
        z = (t*t/2 + s*(1-s)*t*(1-t)/2)/den
        d = (1-t)**2/(4*den)
    else:
        z = (t*(1-t)/4 + t*t/2)/den
        d = F(1,4)-z if state == 'chamfer' else 0
    if span == 2:
        z = 8-z
    return 4-2*s+(2*s-1)*d, z, den


def check(path, state, source):
    data = json.loads(path.read_text())
    degrees = data['degrees']
    knots = data['knots']
    assert degrees == [2,2,2]
    assert knots == [[0,0,0,1,1,1]]*2 + [[0,0,0,.25,.25,7.75,7.75,8,8,8]]
    blocks = data['blocks']
    assert len(blocks) == len(source) == 4
    for idx, (block, parent) in enumerate(zip(blocks, source, strict=True)):
        poles = block['xyzw']
        assert parent['id'] == block['id']
        assert len(poles) == 63
        assert poles[-21:] == blocks[(idx+1)%4]['xyzw'][:21]
        for angle in range(3):
            x,y,_,w = parent['controls'][4*angle]
            for span in range(3):
                row = [[poles[21*angle+7*i+2*span+j] for j in range(3)] for i in range(3)]
                # Three distinct stations per direction determine each homogeneous
                # biquadratic polynomial identically on the entire span.
                for s in (F(0), F(1,2), F(1)):
                    for t in (F(0), F(1,2), F(1)):
                        r,z,h = oracle(s,t,span,state)
                        assert eval_section(row,s,t) == [8*x*r*h,8*y*r*h,32*w*z*h,32*w*h]
    charts = [dict(**b,shape=[3,3,7],degrees=degrees,knots=knots) for b in blocks]
    return [sum(volume(c,q) for c in charts) for q in (12,16)]


if __name__ == '__main__':
    source = json.loads((ROOT/'full-annulus.model.json').read_text())['blocks']
    for state,target in [('sharp',96*np.pi),('chamfer',299.2367002544223),('fillet',300.5816074146435)]:
        suffix = 'fillet-singular' if state == 'fillet' else state
        v12,v16 = check(ROOT/f'annulus-four-{suffix}.json',state,source)
        assert max(abs(v12-target),abs(v16-target)) < 1e-10
        print(f'{state}: exact tensor identities, 4/4 interfaces, V12={v12:.12f} V16={v16:.12f}')
