# Continuous certificate for the stored annular pair

The coefficient oracle checks both degree triples, all 24 homogeneous control points, the integer heights and the incidence. For this family the middle weight is the binary64 value $q$ of $\sqrt{1/2}$. All bounds below apply to the real-valued mathematical map defined by those stored numbers on both closed parameter cubes. The binary64 evaluator is checked separately.

Put $s=1-u$, $t=u$, $X=s^2+2qst$, $Y=t^2+2qst$ and $D=s^2+2qst+t^2$. The planar position is $(2-v)(X,Y)/D$. Since $D$ is a Bernstein combination of $(1,q,1)$, $q\le D\le 1$. Direct expansion gives

$$
X^2+Y^2-D^2=(4q^2-2)s^2t^2,\qquad
N=XY'-YX'=2q(s^2+t^2)+2st\ge q.
$$

The radius error relative to $2-v$ is at most $2|4q^2-2|/(16q^2)$ mm by rationalizing the difference of radii and using $st\le 1/4$. Both $X$ and $Y$ lie in $[0,1]$, so the polar angle derivative $N/(X^2+Y^2)$ is at least $q/2$. The endpoints give angles zero and $\pi/2$. At fixed angle the radius strictly decreases with $v$. Its slope magnitude is at least $1/\sqrt{2}>7/10$ because $X+Y=D+2qst\ge D$.

At each of the three nominal axial levels $n$, the represented elevation deviates from $n$ by $2st(\operatorname{round}_{64}(qn)-qn)/D$. Let $E$ be the maximum absolute product error at these levels divided by $2q$. Interpolation in $w$ bounds absolute elevation error by $E$ throughout each block and bounds the axial derivative from below by $h-2E$. The certificate checks the bounds against the fixed fidelity budget of $10^{-12}$ mm. It uses exact fractions for $q$, the stored products and all bound calculations.

With columns in $(u,v,w)$ order, the planar factor in the determinant is $(2-v)N/D^2\ge q$. The full determinant is at least $q(\min h-2E)>0$ mm³. The angle determines $u$, the radius determines $v$, and the elevation determines $w$. This proves injectivity on each closed cube. For every planar location, the two axial intervals have a common endpoint and disjoint interiors. Exact trace agreement and this interval ordering certify partition of the represented shape. The radial and axial error bounds specify its fidelity to the nominal quarter annulus. This result covers the fixed two-block template; the mandatory J1 corpus contains additional shapes.
