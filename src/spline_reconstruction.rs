//! Bounded four-prime rational reconstruction for exact jet elimination.
//! Unsigned intermediates are at most 244 bits. Only the reconstructed,
//! reduced i128 ratio leaves this module, after modular verification.
use crate::constraint_system::Ratio;
use std::cmp::Ordering;

type Result<T> = std::result::Result<T, &'static str>;
const OVERFLOW: &str = "indeterminate: field reconstruction arithmetic budget";
pub(crate) const PRIMES: [u64; 4] = [
    2_305_843_009_213_693_951,
    2_305_843_009_213_693_921,
    2_305_843_009_213_693_907,
    2_305_843_009_213_693_723,
];

#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
struct Wide([u64; 4]);

impl Wide {
    const ONE: Self = Self([1, 0, 0, 0]);
    fn cmp(self, b: Self) -> Ordering {
        self.0.iter().rev().cmp(b.0.iter().rev())
    }
    fn bits(self) -> usize {
        self.0
            .iter()
            .enumerate()
            .rev()
            .find(|(_, v)| **v != 0)
            .map_or(0, |(i, v)| 64 * i + 64 - v.leading_zeros() as usize)
    }
    fn add(self, b: Self) -> Result<Self> {
        let mut out = self;
        let mut carry = 0u128;
        for (x, y) in out.0.iter_mut().zip(b.0) {
            let s = u128::from(*x) + u128::from(y) + carry;
            *x = s as u64;
            carry = s >> 64;
        }
        if carry != 0 {
            return Err(OVERFLOW);
        }
        Ok(out)
    }
    fn sub(self, b: Self) -> Self {
        debug_assert!(!self.cmp(b).is_lt());
        let mut out = self;
        let mut borrow = false;
        for (x, y) in out.0.iter_mut().zip(b.0) {
            let (a, b1) = x.overflowing_sub(y);
            let (r, b2) = a.overflowing_sub(u64::from(borrow));
            *x = r;
            borrow = b1 || b2;
        }
        out
    }
    fn mul(self, b: Self) -> Result<Self> {
        let mut result = Self::default();
        for (i, &x) in self.0.iter().enumerate() {
            if x == 0 {
                continue;
            }
            let mut carry = 0u128;
            for (j, &y) in b.0.iter().enumerate() {
                if i + j >= 4 {
                    if y != 0 || carry != 0 {
                        return Err(OVERFLOW);
                    }
                    continue;
                }
                let s = u128::from(x) * u128::from(y) + u128::from(result.0[i + j]) + carry;
                result.0[i + j] = s as u64;
                carry = s >> 64;
            }
            if carry != 0 {
                return Err(OVERFLOW);
            }
        }
        Ok(result)
    }
    fn modulo(self, p: u64) -> u64 {
        self.0.iter().rev().fold(0u128, |r, &limb| {
            ((r << 64) + u128::from(limb)) % u128::from(p)
        }) as u64
    }
    // Binary long division, used only by the small Euclidean recovery after CRT.
    fn div_rem(self, b: Self) -> (Self, Self) {
        debug_assert!(b.bits() > 0);
        let mut q = Self::default();
        let mut rem = Self::default();
        for i in (0..self.bits()).rev() {
            let mut carry = (self.0[i / 64] >> (i % 64)) & 1;
            for x in &mut rem.0 {
                let next = *x >> 63;
                *x = (*x << 1) | carry;
                carry = next;
            }
            if !rem.cmp(b).is_lt() {
                rem = rem.sub(b);
                q.0[i / 64] |= 1u64 << (i % 64);
            }
        }
        (q, rem)
    }
    fn as_small(self) -> Result<i128> {
        if self.bits() > 127 {
            return Err(OVERFLOW);
        }
        Ok((i128::from(self.0[1]) << 64) + i128::from(self.0[0]))
    }
}
fn product(a: u64, b: u64, p: u64) -> u64 {
    (u128::from(a) * u128::from(b) % u128::from(p)) as u64
}
fn inverse(mut a: u64, p: u64) -> u64 {
    let mut k = p - 2;
    let mut r = 1;
    while k != 0 {
        if k & 1 != 0 {
            r = product(r, a, p);
        }
        a = product(a, a, p);
        k >>= 1;
    }
    r
}
pub(crate) fn recover(residues: [u64; 4]) -> Result<Ratio> {
    let mut modulus = Wide::ONE;
    let mut x = Wide::default();
    for (&r, p) in residues.iter().zip(PRIMES) {
        let delta = (r + p - x.modulo(p)) % p;
        let digit = product(delta, inverse(modulus.modulo(p), p), p);
        x = x.add(modulus.mul(Wide([digit, 0, 0, 0]))?)?;
        modulus = modulus.mul(Wide([p, 0, 0, 0]))?;
    }
    // With both magnitudes < 2^121 the candidate is unique since
    // 2 * 2^242 is smaller than the product of the four primes.
    let (mut m, mut r, mut t, mut u) = (modulus, x, Wide::default(), Wide::ONE);
    let mut negative = false;
    while r.bits() > 121 {
        let (q, next) = m.div_rem(r);
        (m, r) = (r, next);
        (t, u) = (u, t.add(q.mul(u)?)?);
        negative = !negative;
    }
    if u.bits() > 121 {
        return Err(OVERFLOW);
    }
    let v = Ratio::new(
        r.as_small()?,
        if negative {
            -u.as_small()?
        } else {
            u.as_small()?
        },
    )?;
    for (&original, p) in residues.iter().zip(PRIMES) {
        let n = v.numerator().rem_euclid(i128::from(p)) as u64;
        let d = v.denominator().rem_euclid(i128::from(p)) as u64;
        if d == 0 || product(original, d, p) != n {
            return Err("indeterminate: invalid four-prime reconstruction");
        }
    }
    Ok(v)
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn recovers_large_signed_fraction() {
        for (n, d) in [
            (0i128, 1i128),
            (1, 1),
            (-17, 3),
            (-57784330547159423007526818729, 181802432877686770735109425),
        ] {
            let values = PRIMES.map(|p| {
                product(
                    n.rem_euclid(i128::from(p)) as u64,
                    inverse((d % i128::from(p)) as u64, p),
                    p,
                )
            });
            assert_eq!(recover(values).unwrap(), Ratio::new(n, d).unwrap());
        }
    }
}
