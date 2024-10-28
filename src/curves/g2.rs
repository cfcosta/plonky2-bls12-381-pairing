use ark_bls12_381::{Fq, Fq2, G2Affine};
use ark_ec::{short_weierstrass::SWCurveConfig, AffineRepr};
use ark_ff::BitIteratorBE;
use num::One;
use plonky2::{
    field::extension::Extendable, hash::hash_types::RichField,
    plonk::circuit_builder::CircuitBuilder,
};

use crate::{
    fields::{fq2_target::Fq2Target, fq_target::FqTarget},
    native::miller_loop::G2Projective,
    utils::constants::BLS_X,
};

#[derive(Clone, Debug)]
pub struct G2AffineTarget<F: RichField + Extendable<D>, const D: usize> {
    x: Fq2Target<F, D>,
    y: Fq2Target<F, D>,
    infinity: bool,
}

impl<F: RichField + Extendable<D>, const D: usize> G2AffineTarget<F, D> {
    pub fn constant(builder: &mut CircuitBuilder<F, D>, g2: G2Affine) -> Self {
        Self {
            x: Fq2Target::constant(builder, g2.x().unwrap().clone()),
            y: Fq2Target::constant(builder, g2.y().unwrap().clone()),
            infinity: false,
        }
    }

    fn xy(&self) -> Option<(&self::Fq2Target<F, D>, &self::Fq2Target<F, D>)> {
        (!self.infinity).then(|| (&self.x, &self.y))
    }
}

#[derive(Clone, Debug)]
struct G2ProjectiveTarget<F: RichField + Extendable<D>, const D: usize> {
    x: Fq2Target<F, D>,
    y: Fq2Target<F, D>,
    z: Fq2Target<F, D>,
}

impl<F: RichField + Extendable<D>, const D: usize> G2ProjectiveTarget<F, D> {
    pub fn constant(builder: &mut CircuitBuilder<F, D>, new_g2_projective: &G2Projective) -> Self {
        Self {
            x: Fq2Target::constant(builder, new_g2_projective.x),
            y: Fq2Target::constant(builder, new_g2_projective.y),
            z: Fq2Target::constant(builder, new_g2_projective.z),
        }
    }

    pub fn connect(builder: &mut CircuitBuilder<F, D>, lhs: &Self, rhs: &Self) {
        Fq2Target::connect(builder, &lhs.x, &rhs.x);
        Fq2Target::connect(builder, &lhs.y, &rhs.y);
        Fq2Target::connect(builder, &lhs.z, &rhs.z);
    }

    fn double_in_place(
        &mut self,
        builder: &mut CircuitBuilder<F, D>,
        two_inv: &FqTarget<F, D>,
    ) -> EllCoeffTarget<F, D> {
        let a = self.x.mul(builder, &self.y);
        let a = a.mul_assign_by_fp(builder, two_inv.clone());
        let b = self.y.simple_square(builder);
        let c = self.z.simple_square(builder);
        let coeff_b = ark_bls12_381::g2::Config::COEFF_B;
        let coeff_b = Fq2Target::constant(builder, coeff_b);
        let c_double = c.double(builder);
        let c_double = c_double.add(builder, &c);
        let e = coeff_b.mul(builder, &c_double);
        let e_double = e.double(builder);
        let f = e_double.add(builder, &e);
        let g = b.add(builder, &f);
        let g = g.mul_assign_by_fp(builder, two_inv.clone());
        let y_add_z = self.y.add(builder, &self.z);
        let y_z_squared = y_add_z.simple_square(builder);
        let b_add_c = b.add(builder, &c);
        let h = y_z_squared.sub(builder, &b_add_c);
        let i = e.sub(builder, &b);
        let j = self.x.simple_square(builder);
        let e_square = e.simple_square(builder);

        let b_sub_f = b.sub(builder, &f);

        self.x = a.mul(builder, &b_sub_f);
        let g_square = g.simple_square(builder);
        let e_sq_double = e_square.double(builder);
        let e_sq_double = e_sq_double.add(builder, &e_square);
        self.y = g_square.sub(builder, &e_sq_double);
        self.z = b.mul(builder, &h);

        let j_double = j.double(builder);
        let j_triple = j_double.add(builder, &j);

        (i, j_triple, h.neg(builder))
    }

    fn add_in_place(
        &mut self,
        builder: &mut CircuitBuilder<F, D>,
        q: &G2AffineTarget<F, D>,
    ) -> EllCoeffTarget<F, D> {
        let (qx, qy) = q.xy().unwrap();
        let qy_z = qy.mul(builder, &self.z);
        let theta = self.y.sub(builder, &qy_z);
        let qx_z = qx.mul(builder, &self.z);
        let lambda = self.x.sub(builder, &qx_z);
        let c = theta.simple_square(builder);
        let d = lambda.simple_square(builder);
        let e = lambda.mul(builder, &d);
        let f = self.z.mul(builder, &c);
        let g = self.x.mul(builder, &d);
        let g_double = g.double(builder);
        let e_add_f = e.add(builder, &f);
        let h = e_add_f.sub(builder, &g_double);
        self.x = lambda.mul(builder, &h);
        let g_sub_h = g.sub(builder, &h);
        let e_mul_y = e.mul(builder, &self.y);
        let theta_mul_g_sub_h = theta.mul(builder, &g_sub_h);
        self.y = theta_mul_g_sub_h.sub(builder, &e_mul_y);
        self.z = self.z.mul(builder, &e);
        let lambda_qy = lambda.mul(builder, &qy);
        let theta_qx = theta.mul(builder, &qx);
        let j = theta_qx.sub(builder, &lambda_qy);

        (j, theta.neg(builder), lambda)
    }
}

#[derive(Clone, Debug)]
pub struct G2PreparedTarget<F: RichField + Extendable<D>, const D: usize> {
    /// Stores the coefficients of the line evaluations as calculated in
    /// <https://eprint.iacr.org/2013/722.pdf>
    pub ell_coeffs: Vec<EllCoeffTarget<F, D>>,
    pub infinity: bool,
}

pub(crate) type EllCoeffTarget<F, const D: usize> =
    (Fq2Target<F, D>, Fq2Target<F, D>, Fq2Target<F, D>);

impl<F: RichField + Extendable<D>, const D: usize> G2PreparedTarget<F, D> {
    pub fn is_zero(&self) -> bool {
        self.infinity
    }

    pub fn from(builder: &mut CircuitBuilder<F, D>, q: G2AffineTarget<F, D>) -> Self {
        let one = FqTarget::constant(builder, Fq::one()); // Fq::two
        let two = one.add(builder, &one);
        let two_inv = two.inv(builder);
        let zero = G2PreparedTarget {
            ell_coeffs: vec![],
            infinity: true,
        };

        q.xy().map_or(zero, |(q_x, q_y)| {
            let mut ell_coeffs = vec![];
            let mut r = G2ProjectiveTarget {
                x: q_x.clone(),
                y: q_y.clone(),
                z: Fq2Target::constant(builder, Fq2::one()),
            };

            for i in BitIteratorBE::new([BLS_X]).skip(1) {
                ell_coeffs.push(r.double_in_place(builder, &two_inv));

                if i {
                    ell_coeffs.push(r.add_in_place(builder, &q));
                }
            }

            Self {
                ell_coeffs,
                infinity: false,
            }
        })
    }
}

#[cfg(test)]
mod tests {
    use ark_bls12_381::{G1Affine, G2Affine};
    use ark_ec::AffineRepr;
    use ark_ff::UniformRand;
    use plonky2::{
        field::goldilocks_field::GoldilocksField,
        iop::witness::PartialWitness,
        plonk::{
            circuit_builder::CircuitBuilder, circuit_data::CircuitConfig,
            config::PoseidonGoldilocksConfig,
        },
    };

    use crate::{
        curves::g2::G2AffineTarget, fields::fq_target::FqTarget, native::miller_loop::G2Projective,
    };

    use super::G2ProjectiveTarget;

    type F = GoldilocksField;
    type C = PoseidonGoldilocksConfig;
    const D: usize = 2;

    #[test]
    fn test_g2_projective_double_in_place() {
        let config = CircuitConfig::standard_recursion_config();
        let mut builder = CircuitBuilder::<F, D>::new(config);
        let rng = &mut rand::thread_rng();
        let rand_g1_affine = G1Affine::rand(rng);
        let rand_x = rand_g1_affine.x().unwrap();
        let rand_x_t = FqTarget::constant(&mut builder, *rand_x);

        let rand_g2 = G2Projective::random();
        let mut rand_g2_t = G2ProjectiveTarget::constant(&mut builder, &rand_g2);
        let mut m_rand_g2 = rand_g2;
        m_rand_g2.double_in_place(rand_x);
        rand_g2_t.double_in_place(&mut builder, &rand_x_t);

        let expected_g2 = G2ProjectiveTarget::constant(&mut builder, &m_rand_g2);

        G2ProjectiveTarget::connect(&mut builder, &expected_g2, &rand_g2_t);

        let pw = PartialWitness::new();
        let data = builder.build::<C>();
        let _proof = data.prove(pw);
    }

    #[test]
    fn test_g2_projective_add_in_place() {
        let config = CircuitConfig::standard_recursion_config();
        let mut builder = CircuitBuilder::<F, D>::new(config);
        let rng = &mut rand::thread_rng();
        let rand_g2_affine = G2Affine::rand(rng);
        let rand_g2_affine_t = G2AffineTarget::constant(&mut builder, rand_g2_affine);

        let rand_g2 = G2Projective::random();
        let mut rand_g2_t = G2ProjectiveTarget::constant(&mut builder, &rand_g2);
        let mut m_rand_g2 = rand_g2;
        m_rand_g2.add_in_place(&rand_g2_affine);
        rand_g2_t.add_in_place(&mut builder, &rand_g2_affine_t);

        let expected_g2 = G2ProjectiveTarget::constant(&mut builder, &m_rand_g2);

        G2ProjectiveTarget::connect(&mut builder, &expected_g2, &rand_g2_t);
        let pw = PartialWitness::new();
        let data = builder.build::<C>();
        let _proof = data.prove(pw);
    }
}
