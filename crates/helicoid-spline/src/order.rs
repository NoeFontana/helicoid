//! The order `K` and its exact tables (`0115`, `NUMERICS.md` §21).
//!
//! Every table is integer numerators over one integer denominator, computed by `const fn` at
//! compile time. Each number is asserted below `2^24`, so `S::lit` takes it exactly at binary32 and
//! binary64 and each entry costs one correctly rounded division (`0115` decision 2).

/// The order `K` (degree `K - 1`) of a spline, as a type.
#[derive(Clone, Copy, Debug)]
pub struct Order<const K: usize>;

mod sealed {
    pub trait Sealed {}
}

/// An order whose tables exist: `K = 2..=6` (`0113`, `0115` decision 4).
pub trait ValidOrder: sealed::Sealed {}

macro_rules! valid {
    ($($k:literal)*) => {$(
        impl sealed::Sealed for Order<$k> {}
        impl ValidOrder for Order<$k> {}
    )*};
}
valid!(2 3 4 5 6);

/// The exactness bound of `0115`: below it, an integer is exact in binary32.
const EXACT: i64 = 1 << 24;

const fn factorial(n: usize) -> i64 {
    let mut f = 1;
    let mut i = 2;
    while i <= n {
        f *= i as i64;
        i += 1;
    }
    f
}

const fn binomial(n: usize, k: usize) -> i64 {
    if k > n {
        return 0;
    }
    factorial(n) / (factorial(k) * factorial(n - k))
}

/// `b^e`, with `0^0 = 1` as the basis formula needs.
const fn power(b: i64, e: usize) -> i64 {
    let mut p = 1;
    let mut i = 0;
    while i < e {
        p *= b;
        i += 1;
    }
    p
}

const fn gcd(a: i64, b: i64) -> i64 {
    let (mut a, mut b) = (a.abs(), b.abs());
    while b != 0 {
        let t = a % b;
        a = b;
        b = t;
    }
    a
}

/// `n_{s,k}`: the non-cumulative basis numerators over `(K-1)!` (`NUMERICS.md` §21).
const fn basis<const K: usize>() -> [[i64; K]; K] {
    let mut n = [[0; K]; K];
    let mut s = 0;
    while s < K {
        let mut k = 0;
        while k < K {
            let mut sum = 0;
            let mut l = s;
            while l < K {
                let sign = if (l - s) % 2 == 0 { 1 } else { -1 };
                sum += sign * binomial(K, l - s) * power((K - 1 - l) as i64, K - 1 - k);
                l += 1;
            }
            n[s][k] = binomial(K - 1, k) * sum;
            k += 1;
        }
        s += 1;
    }
    n
}

/// `ñ_{j,k} = Σ_{s ≥ j} n_{s,k}`: the cumulative numerators.
const fn cumulative<const K: usize>() -> [[i64; K]; K] {
    let n = basis::<K>();
    let mut c = [[0; K]; K];
    let mut k = 0;
    while k < K {
        let mut acc = 0;
        let mut j = K;
        while j > 0 {
            j -= 1;
            acc += n[j][k];
            c[j][k] = acc;
        }
        k += 1;
    }
    c
}

/// `p! / (p - r)!`, the coefficient `r` derivatives put on `u^p`.
const fn falling(p: usize, r: usize) -> i64 {
    factorial(p) / factorial(p - r)
}

/// The Gram matrix `G^(r)` at `h = 1` as numerators over one common denominator (`NUMERICS.md`
/// §21). Each entry is formed over the denominator `L ((K-1)!)^2`, `L` the least common multiple of
/// every `p + q - 2r + 1`, and the whole table is then reduced by the gcd of its entries and the
/// denominator, which leaves the least common denominator.
const fn gram<const K: usize>(r: usize) -> ([[i64; K]; K], i64) {
    let n = basis::<K>();
    let f = factorial(K - 1);
    let mut lcm = 1;
    let mut m = 1;
    while m <= 2 * (K - 1 - r) + 1 {
        lcm = lcm / gcd(lcm, m as i64) * m as i64;
        m += 1;
    }
    let den = lcm * f * f;
    let mut g = [[0; K]; K];
    let mut common = den;
    let mut s = 0;
    while s < K {
        let mut t = 0;
        while t < K {
            let mut sum = 0;
            let mut p = r;
            while p < K {
                let mut q = r;
                while q < K {
                    sum += n[s][p]
                        * n[t][q]
                        * falling(p, r)
                        * falling(q, r)
                        * (lcm / (p + q - 2 * r + 1) as i64);
                    q += 1;
                }
                p += 1;
            }
            g[s][t] = sum;
            common = gcd(common, sum);
            t += 1;
        }
        s += 1;
    }
    let mut s = 0;
    while s < K {
        let mut t = 0;
        while t < K {
            g[s][t] /= common;
            t += 1;
        }
        s += 1;
    }
    (g, den / common)
}

const fn below_exact<const K: usize>(t: &[[i64; K]; K], scale: i64) -> bool {
    let mut j = 0;
    while j < K {
        let mut k = 0;
        while k < K {
            // `scale` bounds the derivative rows' `k (k - 1)` factor.
            if (t[j][k] * scale).abs() >= EXACT {
                return false;
            }
            k += 1;
        }
        j += 1;
    }
    true
}

/// The tables of one order.
pub(crate) struct Tables<const K: usize>;

impl<const K: usize> Tables<K> {
    /// `ñ_{j,k}`, row `j`, power `k`.
    pub(crate) const CUMULATIVE: [[i64; K]; K] = cumulative::<K>();
    /// `(K - 1)!`, the basis denominator.
    pub(crate) const FACTORIAL: i64 = factorial(K - 1);
    /// `G^(1)` and `G^(2)` at `h = 1`: numerators and their common denominator.
    pub(crate) const GRAM: [([[i64; K]; K], i64); 2] = [gram::<K>(1), Self::gram2()];

    /// `G^(2)` exists from `K = 3`; at `K = 2` it is the zero matrix over 1.
    const fn gram2() -> ([[i64; K]; K], i64) {
        if K >= 3 {
            gram::<K>(2)
        } else {
            ([[0; K]; K], 1)
        }
    }

    /// Evaluated wherever a table is read, so an order whose numbers stop being exact at binary32
    /// fails to build rather than rounding twice (`0115` decision 4).
    pub(crate) const EXACT: () = {
        let k1 = (K - 1) as i64;
        // `k (k - 1) <= (K - 1)^2` bounds the second-derivative row's factor.
        assert!(below_exact(
            &Self::CUMULATIVE,
            if k1 > 1 { k1 * k1 } else { 1 }
        ));
        assert!(Self::FACTORIAL < EXACT);
        let mut r = 0;
        while r < 2 {
            assert!(below_exact(&Self::GRAM[r].0, 1) && Self::GRAM[r].1 < EXACT);
            r += 1;
        }
    };
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_k4_rows_are_the_published_ones() {
        // `NUMERICS.md` §21's check: `(K-1)! M̃` at K = 4.
        assert_eq!(
            Tables::<4>::CUMULATIVE,
            [[6, 0, 0, 0], [5, 3, -3, 1], [1, 3, 3, -2], [0, 0, 0, 1]]
        );
        assert_eq!(Tables::<4>::FACTORIAL, 6);
    }

    #[test]
    fn the_gram_check_values_hold() {
        // G^(1)_00 = 1/20 at K = 4, G^(2)_00 = 1/252 at K = 6 (`NUMERICS.md` §21).
        let (g, d) = Tables::<4>::GRAM[0];
        assert_eq!(g[0][0] * 20, d);
        let (g, d) = Tables::<6>::GRAM[1];
        assert_eq!(g[0][0] * 252, d);
    }

    #[test]
    fn the_denominators_are_0115s() {
        // `0115`'s table: the reduced common denominators.
        assert_eq!(Tables::<6>::GRAM[0].1, 362_880);
        assert_eq!(Tables::<6>::GRAM[1].1, 5_040);
        assert_eq!(Tables::<5>::GRAM[0].1, 5_040);
        assert_eq!(Tables::<3>::GRAM[1].1, 1);
    }
}
