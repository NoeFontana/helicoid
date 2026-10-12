//! `helicoid:unit-r`: the shipped `Ad` of SE_N(3) and Sim(3) with one variable changed, the
//! rotation block `R(q)` in Eigen's unit form `(1 − 2‖u‖²) I + 2uuᵀ + 2w[u]×` in place of
//! `Quat::to_matrix`'s homogeneous `(w² − ‖u‖²) I + …`. It is the twin `0076` decision 1 measured
//! by editing the library, kept as a subject so `cargo xtask compare` can re-measure it (`0078`).
//!
//! Assembled here from the public API, since the library's `SEn3Jac` and `Sim3Jac` blocks are
//! crate-private (`0028`): the element parsed as [`Helicoid`] parses it, the blocks of
//! `SEn3::adjoint` and `Sim3::adjoint` formed in the same operation order, written dense as
//! `Jac::write_dense` writes them. That the assembly changes only `R` is pinned:
//! with `Quat::to_matrix` it is [`Helicoid`]'s answer to the bit on every record
//! (`the_assembly_with_the_shipped_r_is_the_shipped_ad`).
//!
//! Planted, as `helicoid:host-std` is: a twin is not a candidate, and no bar reads its rows.

use helicoid::Quat;
use helicoid_linalg::{Mat3, Matrix, Precision, Real, Vec3, Vector};

use super::{sen3_of, sim3_of, Helicoid, Sen3, Sim3Id};
use crate::conformance::corpus::Record;
use crate::conformance::subject::{Output, Registered, Subject};

/// `R(q)` as Eigen's `Quaternion::toRotationMatrix` rounds it: the diagonal `1 − (2y² + 2z²)`,
/// the off-diagonal entries `Quat::to_matrix`'s, from the same nine products.
pub(crate) fn unit_form<S: Real>(q: &Quat<S>) -> Mat3<S> {
    let (w, x, y, z) = (q.w, q.x, q.y, q.z);
    let two = S::lit(2.0);
    let (tw, tx, ty, tz) = (two * w, two * x, two * y, two * z);
    let (xx, yy, zz) = (tx * x, ty * y, tz * z);
    let (xy, xz, yz) = (tx * y, tx * z, ty * z);
    let (wx, wy, wz) = (tw * x, tw * y, tw * z);
    let one = S::one();
    Mat3::from_cols([
        Vector([one - (yy + zz), xy + wz, xz - wy]),
        Vector([xy - wz, one - (xx + zz), yz + wx]),
        Vector([xz + wy, yz - wx, one - (xx + yy)]),
    ])
}

/// `[v]_× B` in the library's `so3::hat_mul` order, which is crate-private.
fn hat_mul<S: Real>(v: Vec3<S>, b: &Mat3<S>) -> Mat3<S> {
    let [x, y, z] = v.0;
    Matrix::from_cols(core::array::from_fn(|c| {
        let (b0, b1, b2) = (b.get(0, c), b.get(1, c), b.get(2, c));
        Vector([-z * b1 + y * b2, z * b0 + -x * b2, -y * b0 + x * b1])
    }))
}

/// A dense `d × d` column-major buffer, zero but for the blocks [`Dense::put`] writes.
struct Dense<S> {
    d: usize,
    buf: Vec<S>,
}

impl<S: Real + Into<f64>> Dense<S> {
    fn new(d: usize) -> Self {
        Self {
            d,
            buf: vec![S::zero(); d * d],
        }
    }

    fn put(&mut self, r0: usize, c0: usize, m: &Mat3<S>) {
        for c in 0..3 {
            for r in 0..3 {
                self.buf[(c0 + c) * self.d + r0 + r] = m.get(r, c);
            }
        }
    }

    fn out(self) -> Output {
        let v = self.buf.into_iter().map(Into::into).collect();
        Output::from([("Ad".to_string(), v)])
    }
}

/// `Ad_X = R + ε [x_i]_× R` (`NUMERICS.md` §5.2), `R` from `r_of`.
fn sen3_ad<S: Real + Into<f64> + From<f64>, const N: usize>(
    record: &Record,
    r_of: fn(&Quat<S>) -> Mat3<S>,
) -> Output {
    let Some(x) = sen3_of::<S, N>(record, "q", "x") else {
        return Output::new();
    };
    let (rot, cols) = x.parts();
    let r = r_of(&rot.quat());
    let mut dense = Dense::new(3 + 3 * N);
    dense.put(0, 0, &r);
    for (i, v) in cols.into_iter().enumerate() {
        dense.put(3 + 3 * i, 3 + 3 * i, &r);
        dense.put(3 + 3 * i, 0, &hat_mul(v, &r));
    }
    dense.out()
}

/// `Ad_X = [R 0 0; [t]_× R  sR  −t; 0 0 1]` (SM.9), `R` from `r_of`.
fn sim3_ad<S: Real + Into<f64>>(record: &Record, r_of: fn(&Quat<S>) -> Mat3<S>) -> Output {
    let Some(x) = sim3_of::<S>(record) else {
        return Output::new();
    };
    let (rot, t, _) = x.parts();
    let r = r_of(&rot.quat());
    let mut dense = Dense::new(7);
    dense.put(0, 0, &r);
    dense.put(3, 0, &hat_mul(t, &r));
    dense.put(3, 3, &r.scale(x.scale()));
    for (i, f) in (-t).0.into_iter().enumerate() {
        dense.buf[6 * 7 + 3 + i] = f;
    }
    dense.buf[6 * 7 + 6] = S::one();
    dense.out()
}

/// The four `Ad` ids at `precision`, `R` from `r64` or `r32`; nothing for any other id, and the
/// SE_N(3) ids at binary64 alone, as [`Helicoid`] answers them.
fn ad(
    fn_id: &str,
    record: &Record,
    precision: Precision,
    r64: fn(&Quat<f64>) -> Mat3<f64>,
    r32: fn(&Quat<f32>) -> Mat3<f32>,
) -> Output {
    if let Some(Sim3Id::Ad) = Sim3Id::of_fn(fn_id) {
        return match precision {
            Precision::F64 => sim3_ad::<f64>(record, r64),
            Precision::F32 => sim3_ad::<f32>(record, r32),
        };
    }
    match (Sen3::of_fn(fn_id), precision) {
        (Some((Sen3::Ad, 1)), Precision::F64) => sen3_ad::<f64, 1>(record, r64),
        (Some((Sen3::Ad, 2)), Precision::F64) => sen3_ad::<f64, 2>(record, r64),
        (Some((Sen3::Ad, 3)), Precision::F64) => sen3_ad::<f64, 3>(record, r64),
        _ => Output::new(),
    }
}

/// The ids the twin answers: the four whose answer reads `R(q)` as a matrix.
const IDS: [&str; 4] = ["sen3_ad_n1", "sen3_ad_n2", "sen3_ad_n3", "sim3_ad"];

pub(crate) struct UnitR;

impl Subject for UnitR {
    fn name(&self) -> &str {
        "helicoid:unit-r"
    }

    fn supports(&self, fn_id: &str) -> bool {
        IDS.contains(&fn_id) && Helicoid.supports(fn_id)
    }

    fn eval(&self, fn_id: &str, record: &Record, precision: Precision) -> Output {
        ad(fn_id, record, precision, unit_form, unit_form)
    }
}

pub(crate) fn registered() -> Registered {
    let version = format!("{}@unit-r", env!("CARGO_PKG_VERSION"));
    Registered {
        version: version.clone(),
        version_f32: version,
        no_f32: None,
        subject: Box::new(UnitR),
        planted: true,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::conformance::{corpus, corpus_dir};

    fn bits(o: &Output) -> Vec<(String, Vec<u64>)> {
        o.iter()
            .map(|(k, v)| (k.clone(), v.iter().map(|x| x.to_bits()).collect()))
            .collect()
    }

    /// The control: every variable but `R` is the shipped program's.
    #[test]
    fn the_assembly_with_the_shipped_r_is_the_shipped_ad() -> Result<(), String> {
        let dir = corpus_dir()?;
        let (mut checked, mut moved) = (0, 0);
        for entry in corpus::manifest(&dir)?
            .into_iter()
            .filter(|e| IDS.contains(&e.fn_id.as_str()))
        {
            for rec in corpus::read(&dir, &entry)? {
                let p = match rec.is_f32_stratum() {
                    true => Precision::F32,
                    false => Precision::F64,
                };
                let shipped = Helicoid.eval(&entry.fn_id, &rec, p);
                let twin = ad(&entry.fn_id, &rec, p, Quat::to_matrix, Quat::to_matrix);
                assert!(!shipped.is_empty(), "{} {}", entry.fn_id, rec.id);
                assert_eq!(bits(&twin), bits(&shipped), "{} {}", entry.fn_id, rec.id);
                checked += 1;
                let unit = UnitR.eval(&entry.fn_id, &rec, p);
                moved += usize::from(bits(&unit) != bits(&shipped));
            }
        }
        // Every record of the four ids, and the variable is a variable on most of them.
        assert_eq!(checked, 3 * 324 + 1008);
        assert!(moved > checked / 2, "{moved} of {checked}");
        Ok(())
    }

    #[test]
    fn the_two_forms_are_one_matrix_on_an_exact_unit_quaternion() {
        // A quarter turn about z and the identity: every product is exact, so both forms are.
        let h = 0.5f64.sqrt();
        for q in [
            Quat {
                w: 1.0,
                x: 0.0,
                y: 0.0,
                z: 0.0,
            },
            Quat {
                w: 0.0,
                x: 0.0,
                y: 0.0,
                z: 1.0,
            },
            Quat {
                w: 0.5,
                x: 0.5,
                y: 0.5,
                z: 0.5,
            },
            Quat {
                w: h,
                x: 0.0,
                y: 0.0,
                z: h,
            },
        ] {
            let (u, m) = (unit_form(&q), q.to_matrix());
            for (r, c) in (0..3).flat_map(|r| (0..3).map(move |c| (r, c))) {
                assert!(
                    (u.get(r, c) - m.get(r, c)).abs() <= 4.0 * f64::EPSILON,
                    "{q:?}"
                );
            }
        }
    }
}
