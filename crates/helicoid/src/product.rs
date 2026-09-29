//! Products of groups (`docs/PHASE3.md` §7): `Product<A, B>` with the pair of the factors'
//! tangents and the block-diagonal Jacobian [`ProductJac`].
//!
//! Every operation is the factors' operation, componentwise; every Jacobian row of `NUMERICS.md`
//! §2.3 is the block-diagonal composition of the factors' rows, because the tangent of a product
//! is the direct sum of the factors' tangents and `Exp`, `Log`, `Ad`, `ad` and the Jacobians of
//! `Exp` respect it. The dense order is `A` then `B` (`NUMERICS.md` §1).

use crate::side::Side;
use crate::traits::{tie_dof, Jac, LieGroup, Tangent};
use core::array;
use core::ops::Mul;
use helicoid_linalg::{Blend, Matrix, Real, StridedMut, Vector};

/// The product group `A x B`: composition is componentwise, spelled `a * b` (`docs/API.md` R1).
///
/// Its tangent is the pair `(A::Tangent, B::Tangent)`, its Jacobian is
/// [`ProductJac<A::Jac, B::Jac>`](ProductJac) and its `DOF` is `A::DOF + B::DOF`. Nesting is
/// nesting: `Product<Product<A, B>, C>` has tangent `((TA, TB), TC)` and the dense order `A`, `B`,
/// `C`. Composition is componentwise, so `Product<SO3<S>, Rn<S, 3>>` is not `SE(3)`: its
/// geodesic is the `tf2` slerp + lerp (`NUMERICS.md` §10).
///
/// There is no `Add` or `Sub`, on the group or on its tangent; this does not compile:
///
/// ```compile_fail,E0369
/// use helicoid::{Product, Rn};
/// use helicoid_linalg::Vector;
/// let a = Product(Rn(Vector([1.0_f64; 2])), Rn(Vector([1.0_f64; 3])));
/// let _ = a + a;
/// ```
///
/// Positive control: composition is `Mul`.
///
/// ```
/// use helicoid::{Product, Rn};
/// use helicoid_linalg::Vector;
/// let a = Product(Rn(Vector([1.0_f64; 2])), Rn(Vector([1.0_f64; 3])));
/// let b = a * a;
/// assert_eq!((b.0 .0 .0, b.1 .0 .0), ([2.0; 2], [2.0; 3]));
/// ```
#[derive(Clone, Copy, Debug)]
pub struct Product<A, B>(pub A, pub B);

/// A Jacobian of a [`Product`]: the block-diagonal matrix `diag(J_A, J_B)`.
///
/// It is closed under [`mul`](Jac::mul), [`inverse`](Jac::inverse) and [`neg`](Jac::neg) because
/// each factor's is, and a dense matrix leaves it only through [`write_dense`](Jac::write_dense)
/// and [`sandwich`](Jac::sandwich) (`docs/API.md` R5). The parameters are the factors' Jacobian
/// types, `Product<A, B>::Jac = ProductJac<A::Jac, B::Jac>`. There is no `PartialEq`, `Add` or
/// `Sub`.
#[derive(Clone, Copy, Debug)]
pub struct ProductJac<JA, JB>(pub JA, pub JB);

impl<A, B> Mul for Product<A, B>
where
    A: Mul<Output = A>,
    B: Mul<Output = B>,
{
    type Output = Self;
    #[inline]
    fn mul(self, o: Self) -> Self {
        Self(self.0 * o.0, self.1 * o.1)
    }
}

impl<S: Real, A: Blend<S>, B: Blend<S>> Blend<S> for Product<A, B> {
    #[inline]
    fn blend(m: S::Mask, t: Self, f: Self) -> Self {
        Self(A::blend(m, t.0, f.0), B::blend(m, t.1, f.1))
    }
}

impl<S: Real, JA: Blend<S>, JB: Blend<S>> Blend<S> for ProductJac<JA, JB> {
    #[inline]
    fn blend(m: S::Mask, t: Self, f: Self) -> Self {
        Self(JA::blend(m, t.0, f.0), JB::blend(m, t.1, f.1))
    }
}

// The tangent of a product is the pair; `Blend` for pairs is `helicoid-linalg`'s.
impl<S: Real, TA: Tangent<S>, TB: Tangent<S>> Tangent<S> for (TA, TB) {
    const DOF: usize = TA::DOF + TB::DOF;
    #[inline]
    fn zero() -> Self {
        (TA::zero(), TB::zero())
    }
    #[inline]
    fn add(&self, o: &Self) -> Self {
        (self.0.add(&o.0), self.1.add(&o.1))
    }
    #[inline]
    fn sub(&self, o: &Self) -> Self {
        (self.0.sub(&o.0), self.1.sub(&o.1))
    }
    #[inline]
    fn neg(&self) -> Self {
        (self.0.neg(), self.1.neg())
    }
    #[inline]
    fn scale(&self, k: S) -> Self {
        (self.0.scale(k), self.1.scale(k))
    }
    /// The sum of the two factors' inner products, each summed in index order: two partial
    /// sums, not one flat sum, since a buffer of `DOF` scalars cannot be sized here (`0005`).
    #[inline]
    fn dot(&self, o: &Self) -> S {
        self.0.dot(&o.0) + self.1.dot(&o.1)
    }
    #[inline]
    fn write_dense(&self, out: &mut [S]) {
        debug_assert!(
            out.len() == <Self as Tangent<S>>::DOF,
            "Tangent::write_dense: wrong length"
        );
        let (a, b) = out.split_at_mut(TA::DOF.min(out.len()));
        self.0.write_dense(a);
        self.1.write_dense(b);
    }
    #[inline]
    fn read_dense(src: &[S]) -> Self {
        debug_assert!(
            src.len() == <Self as Tangent<S>>::DOF,
            "Tangent::read_dense: wrong length"
        );
        let (a, b) = src.split_at(TA::DOF.min(src.len()));
        (TA::read_dense(a), TB::read_dense(b))
    }
}

/// `diag(J_A, J_B)` acting on a `D`-vector split at `TA::DOF`.
fn apply_flat<S, TA, TB, JA, JB, const D: usize>(
    j: &ProductJac<JA, JB>,
    v: &Vector<S, D>,
) -> Vector<S, D>
where
    S: Real,
    TA: Tangent<S>,
    TB: Tangent<S>,
    JA: Jac<S, TA>,
    JB: Jac<S, TB>,
{
    let mut out = [S::zero(); D];
    let n = TA::DOF.min(D);
    let (va, vb) = v.0.split_at(n);
    let (oa, ob) = out.split_at_mut(n);
    j.0.apply(&TA::read_dense(va)).write_dense(oa);
    j.1.apply(&TB::read_dense(vb)).write_dense(ob);
    Vector(out)
}

impl<S, TA, TB, JA, JB> Jac<S, (TA, TB)> for ProductJac<JA, JB>
where
    S: Real,
    TA: Tangent<S>,
    TB: Tangent<S>,
    JA: Jac<S, TA>,
    JB: Jac<S, TB>,
{
    #[inline]
    fn identity() -> Self {
        Self(JA::identity(), JB::identity())
    }
    #[inline]
    fn mul(&self, o: &Self) -> Self {
        Self(self.0.mul(&o.0), self.1.mul(&o.1))
    }
    /// `diag(J_A⁻¹, J_B⁻¹)`.
    ///
    /// # Domain
    ///
    /// Each factor is invertible, as its own `inverse` states.
    #[inline]
    fn inverse(&self) -> Self {
        Self(self.0.inverse(), self.1.inverse())
    }
    #[inline]
    fn neg(&self) -> Self {
        Self(self.0.neg(), self.1.neg())
    }
    #[inline]
    fn apply(&self, t: &(TA, TB)) -> (TA, TB) {
        (self.0.apply(&t.0), self.1.apply(&t.1))
    }
    #[inline]
    fn apply_transpose(&self, t: &(TA, TB)) -> (TA, TB) {
        (self.0.apply_transpose(&t.0), self.1.apply_transpose(&t.1))
    }
    /// Writes `diag(J_A, J_B)`, the two off-diagonal blocks as `+0`.
    ///
    /// # Domain
    ///
    /// `out` is `DOF x DOF`, checked by `debug_assert!`. A smaller view panics in the strided
    /// access, the one documented panic class (D11).
    #[inline]
    fn write_dense(&self, out: &mut StridedMut<'_, S>) {
        let (da, dof) = (TA::DOF, <(TA, TB) as Tangent<S>>::DOF);
        debug_assert!(
            out.rows() == dof && out.cols() == dof,
            "Jac::write_dense: the view is not DOF x DOF"
        );
        for c in 0..dof {
            for r in 0..dof {
                if (r < da) != (c < da) {
                    out.set(r, c, S::zero());
                }
            }
        }
        self.0.write_dense(&mut out.block(0, 0, da, da));
        self.1
            .write_dense(&mut out.block(da, da, dof - da, dof - da));
    }
    /// `J Σ Jᵀ` from the factors' `apply` alone: `M = J Σ` column by column, then `(J Mᵀ)ᵀ`, so
    /// `Σ` need not be symmetric and the result is not symmetrized. The factors' own `sandwich`
    /// takes a square block of a covariance whose size is `TA::DOF`, and stable Rust cannot
    /// build one from a generic `D`.
    ///
    /// `D` must be `DOF`, asserted at monomorphization (`cargo build` and `cargo test`, not
    /// `cargo check`):
    ///
    /// ```compile_fail,E0080
    /// use helicoid::{Jac, ProductJac, RnJac};
    /// use helicoid_linalg::Matrix;
    /// let j: ProductJac<RnJac<f64, 2>, RnJac<f64, 3>> = Jac::identity();
    /// let _ = j.sandwich::<6>(&Matrix::<f64, 6, 6>::identity());
    /// ```
    ///
    /// Positive control, so the failure above comes from the assertion and not from an import:
    ///
    /// ```
    /// use helicoid::{Jac, ProductJac, RnJac};
    /// use helicoid_linalg::Matrix;
    /// let j: ProductJac<RnJac<f64, 2>, RnJac<f64, 3>> = Jac::identity();
    /// let _ = j.sandwich::<5>(&Matrix::<f64, 5, 5>::identity());
    /// ```
    fn sandwich<const D: usize>(&self, cov: &Matrix<S, D, D>) -> Matrix<S, D, D> {
        const { assert!(D == <(TA, TB) as Tangent<S>>::DOF) };
        let m = Matrix::from_cols(array::from_fn(|c| {
            apply_flat::<S, TA, TB, JA, JB, D>(self, &cov.col(c))
        }));
        // The columns of `Mᵀ` are the rows of `M`, and `(J Mᵀ)ᵀ = M Jᵀ`.
        Matrix::from_cols(array::from_fn(|r| {
            apply_flat::<S, TA, TB, JA, JB, D>(self, &m.row(r))
        }))
        .transpose()
    }
}

// Every method, the provided ones included, is the factors' method, so a factor that overrides a
// provided method is what the product uses (no factor does yet) and each row of `NUMERICS.md`
// §2.3 is the factors' rows side by side.
impl<S: Real, A: LieGroup<S>, B: LieGroup<S>> LieGroup<S> for Product<A, B> {
    type Tangent = (A::Tangent, B::Tangent);
    type Jac = ProductJac<A::Jac, B::Jac>;
    const DOF: usize = <A as LieGroup<S>>::DOF + <B as LieGroup<S>>::DOF;

    #[inline]
    fn identity() -> Self {
        const { tie_dof::<S, Self>() };
        Self(A::identity(), B::identity())
    }
    #[inline]
    fn inverse(&self) -> Self {
        Self(self.0.inverse(), self.1.inverse())
    }
    #[inline]
    fn exp(tau: &Self::Tangent) -> Self {
        Self(A::exp(&tau.0), B::exp(&tau.1))
    }
    #[inline]
    fn log(&self) -> Self::Tangent {
        (self.0.log(), self.1.log())
    }
    #[inline]
    fn adjoint(&self) -> Self::Jac {
        ProductJac(self.0.adjoint(), self.1.adjoint())
    }
    #[inline]
    fn ad(tau: &Self::Tangent) -> Self::Jac {
        ProductJac(A::ad(&tau.0), B::ad(&tau.1))
    }
    #[inline]
    fn jr(tau: &Self::Tangent) -> Self::Jac {
        ProductJac(A::jr(&tau.0), B::jr(&tau.1))
    }
    /// # Domain
    ///
    /// Each factor's `jr_inv` domain.
    #[inline]
    fn jr_inv(tau: &Self::Tangent) -> Self::Jac {
        ProductJac(A::jr_inv(&tau.0), B::jr_inv(&tau.1))
    }
    #[inline]
    fn jl(tau: &Self::Tangent) -> Self::Jac {
        ProductJac(A::jl(&tau.0), B::jl(&tau.1))
    }
    #[inline]
    fn jl_inv(tau: &Self::Tangent) -> Self::Jac {
        ProductJac(A::jl_inv(&tau.0), B::jl_inv(&tau.1))
    }
    #[inline]
    fn rplus(&self, tau: &Self::Tangent) -> Self {
        Self(self.0.rplus(&tau.0), self.1.rplus(&tau.1))
    }
    #[inline]
    fn lplus(&self, tau: &Self::Tangent) -> Self {
        Self(self.0.lplus(&tau.0), self.1.lplus(&tau.1))
    }
    #[inline]
    fn rminus(&self, base: &Self) -> Self::Tangent {
        (self.0.rminus(&base.0), self.1.rminus(&base.1))
    }
    #[inline]
    fn lminus(&self, base: &Self) -> Self::Tangent {
        (self.0.lminus(&base.0), self.1.lminus(&base.1))
    }
    #[inline]
    fn rplus_jacobians(&self, tau: &Self::Tangent) -> (Self::Jac, Self::Jac) {
        let (a, b) = (
            self.0.rplus_jacobians(&tau.0),
            self.1.rplus_jacobians(&tau.1),
        );
        (ProductJac(a.0, b.0), ProductJac(a.1, b.1))
    }
    #[inline]
    fn lplus_jacobians(&self, tau: &Self::Tangent) -> (Self::Jac, Self::Jac) {
        let (a, b) = (
            self.0.lplus_jacobians(&tau.0),
            self.1.lplus_jacobians(&tau.1),
        );
        (ProductJac(a.0, b.0), ProductJac(a.1, b.1))
    }
    #[inline]
    fn rminus_jacobians(&self, base: &Self) -> (Self::Jac, Self::Jac) {
        let (a, b) = (
            self.0.rminus_jacobians(&base.0),
            self.1.rminus_jacobians(&base.1),
        );
        (ProductJac(a.0, b.0), ProductJac(a.1, b.1))
    }
    #[inline]
    fn lminus_jacobians(&self, base: &Self) -> (Self::Jac, Self::Jac) {
        let (a, b) = (
            self.0.lminus_jacobians(&base.0),
            self.1.lminus_jacobians(&base.1),
        );
        (ProductJac(a.0, b.0), ProductJac(a.1, b.1))
    }
    #[inline]
    fn compose_jacobians<Sd: Side>(&self, rhs: &Self) -> (Self::Jac, Self::Jac) {
        let (a, b) = (
            self.0.compose_jacobians::<Sd>(&rhs.0),
            self.1.compose_jacobians::<Sd>(&rhs.1),
        );
        (ProductJac(a.0, b.0), ProductJac(a.1, b.1))
    }
    #[inline]
    fn inverse_jacobian<Sd: Side>(&self) -> Self::Jac {
        ProductJac(
            self.0.inverse_jacobian::<Sd>(),
            self.1.inverse_jacobian::<Sd>(),
        )
    }
}
