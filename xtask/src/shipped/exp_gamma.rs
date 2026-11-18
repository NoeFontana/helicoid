//! `so3::{exp_gamma_apply, exp_gamma_apply_jacobians}` replayed over the committed records of
//! `so3_gamma2` and `so3_exp` (`0119`): each output is the separate calls' to the bit, at `f64` on
//! the binary64 strata, at `f32` on the `@f32` ones, and on `Dual<_, 1>` at both, every lane
//! compared. Neither function is a corpus id of its own: what it promises is that no record of
//! the three it fuses moves, and this is that promise read on their records.

#[cfg(test)]
mod tests {
    use helicoid::{so3, LieGroup, SO3Tangent, SO3};
    use helicoid_linalg::{Dual, Mat3, Precision, Real, Vec3, Vector};

    use crate::conformance::corpus::{self, Record};
    use crate::conformance::corpus_dir;
    use crate::shipped::scalars;

    /// Every lane's bits: one for a float, the value's and then the derivative's for `Dual`.
    trait Lanes: Real {
        fn lanes(self, out: &mut Vec<u64>);
        /// `x` exactly, and for `Dual` the derivative lane `dx`.
        fn seed(x: Self::Inner, dx: f64) -> Self;
        type Inner: Real;
    }

    impl Lanes for f64 {
        type Inner = f64;
        fn lanes(self, out: &mut Vec<u64>) {
            out.push(self.to_bits());
        }
        fn seed(x: f64, _: f64) -> Self {
            x
        }
    }

    impl Lanes for f32 {
        type Inner = f32;
        fn lanes(self, out: &mut Vec<u64>) {
            out.push(u64::from(self.to_bits()));
        }
        fn seed(x: f32, _: f64) -> Self {
            x
        }
    }

    impl<S: Lanes<Inner = S>> Lanes for Dual<S, 1> {
        type Inner = S;
        fn lanes(self, out: &mut Vec<u64>) {
            self.v.lanes(out);
            self.d[0].lanes(out);
        }
        fn seed(x: S, dx: f64) -> Self {
            let dx = match S::PRECISION {
                Precision::F64 => S::lit(dx),
                Precision::F32 => S::lit(f64::from(dx as f32)),
            };
            Dual { v: x, d: [dx] }
        }
    }

    fn of_vec<S: Lanes>(v: Vec3<S>, out: &mut Vec<u64>) {
        v.0.iter().for_each(|x| x.lanes(out));
    }

    fn of_mat<S: Lanes>(m: &Mat3<S>, out: &mut Vec<u64>) {
        (0..9).for_each(|k| m.get(k / 3, k % 3).lanes(out));
    }

    fn of_rot<S: Lanes>(r: &SO3<S>, out: &mut Vec<u64>) {
        let q = r.quat();
        [q.w, q.x, q.y, q.z].iter().for_each(|x| x.lanes(out));
    }

    fn of_triple<S: Lanes>(t: &so3::GammaJacobians<S>, out: &mut Vec<u64>) {
        of_vec(t.0, out);
        of_mat(&t.1, out);
        of_mat(&t.2, out);
    }

    /// Whether both fused functions are the separate calls at `(φ, v)`, to the bit.
    fn fused_is_separate<S: Lanes>(phi: &SO3Tangent<S>, v: Vec3<S>) -> bool {
        let (mut got, mut want) = (Vec::new(), Vec::new());
        let (r, g1v, g2v) = so3::exp_gamma_apply(phi, v);
        of_rot(&r, &mut got);
        of_vec(g1v, &mut got);
        of_vec(g2v, &mut got);
        let (r, t1, t2) = so3::exp_gamma_apply_jacobians(phi, v);
        of_rot(&r, &mut got);
        of_triple(&t1, &mut got);
        of_triple(&t2, &mut got);
        let r = SO3::exp(phi);
        of_rot(&r, &mut want);
        of_vec(so3::gamma1(phi) * v, &mut want);
        of_vec(so3::gamma2(phi) * v, &mut want);
        of_rot(&r, &mut want);
        of_triple(&so3::gamma_apply_jacobians::<1, S>(phi, v), &mut want);
        of_triple(&so3::gamma_apply_jacobians::<2, S>(phi, v), &mut want);
        got == want
    }

    /// `v` along `φ`, across it, and neither, each exact at binary32.
    fn vs(phi: [f64; 3]) -> [[f64; 3]; 3] {
        let [x, y, z] = phi;
        [[x, y, z], [-y, x, 0.0], [1.0, -2.0, 0.5]]
    }

    /// The records of `fn_id` at `precision` that `S` holds exactly, each at three `v`, plain and on
    /// `Dual` seeded along a fixed direction; the count checked, `Err` naming the first that parts.
    fn replay<S: Lanes<Inner = S>>(fn_id: &str, precision: Precision) -> Result<usize, String> {
        const DIR: [f64; 3] = [0.6, -0.48, 0.64];
        let dir = corpus_dir()?;
        let entry = corpus::manifest(&dir)?
            .into_iter()
            .find(|e| e.fn_id == fn_id)
            .ok_or(format!("no {fn_id} in the manifest"))?;
        let records: Vec<Record> = corpus::read(&dir, &entry)?
            .into_iter()
            .filter(|r| r.is_f32_stratum() == (precision == Precision::F32))
            .collect();
        let mut n = 0;
        for rec in &records {
            let Some(phi) = scalars::<S>(rec, "phi").filter(|p| p.len() == 3) else {
                return Err(format!("{fn_id} {}: no φ", rec.id));
            };
            let wide = phi.iter().map(|x| x.value_f64()).collect::<Vec<_>>();
            for v in vs([wide[0], wide[1], wide[2]]) {
                let v = v.map(|x| match S::PRECISION {
                    Precision::F64 => S::lit(x),
                    Precision::F32 => S::lit(f64::from(x as f32)),
                });
                let plain = SO3Tangent {
                    phi: Vector([phi[0], phi[1], phi[2]]),
                };
                let dual = SO3Tangent {
                    phi: Vector(core::array::from_fn(|i| {
                        <Dual<S, 1> as Lanes>::seed(phi[i], DIR[i])
                    })),
                };
                let dv = Vector(v.map(|x| <Dual<S, 1> as Lanes>::seed(x, 0.0)));
                if !fused_is_separate(&plain, Vector(v)) || !fused_is_separate(&dual, dv) {
                    return Err(format!("{fn_id} {} ({}) at {wide:?}", rec.id, rec.stratum));
                }
                n += 1;
            }
        }
        Ok(n)
    }

    #[test]
    fn exp_gamma_apply_is_the_separate_calls_on_the_corpus() -> Result<(), String> {
        for fn_id in ["so3_gamma2", "so3_exp"] {
            let wide = replay::<f64>(fn_id, Precision::F64)?;
            assert!(wide > 0, "{fn_id}: no binary64 record");
            if fn_id == "so3_gamma2" {
                let narrow = replay::<f32>(fn_id, Precision::F32)?;
                assert!(narrow > 0, "{fn_id}: no binary32 record");
            }
        }
        Ok(())
    }
}
