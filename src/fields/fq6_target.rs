use ark_bls12_381::{Fq, Fq6};
use ark_ff::Field;
use itertools::Itertools;
use num::BigUint;
use plonky2::{
    field::extension::Extendable,
    hash::hash_types::RichField,
    iop::{
        generator::{GeneratedValues, SimpleGenerator},
        target::{BoolTarget, Target},
        witness::PartitionWitness,
    },
    plonk::circuit_builder::CircuitBuilder,
    util::serialization::Buffer,
};
use plonky2_ecdsa::gadgets::{
    biguint::{GeneratedValuesBigUint, WitnessBigUint},
    nonnative::CircuitBuilderNonNative,
};

use super::{fq2_target::Fq2Target, fq_target::FqTarget};
use crate::utils::{helpers::from_biguint_to_fq, my_fq6::MyFq6};

#[derive(Debug, Clone)]
pub struct Fq6Target<F: RichField + Extendable<D>, const D: usize> {
    pub coeffs: [FqTarget<F, D>; 6],
}

impl<F: RichField + Extendable<D>, const D: usize> Fq6Target<F, D> {
    pub fn empty(builder: &mut CircuitBuilder<F, D>) -> Self {
        let coeffs = [(); 6]
            .iter()
            .map(|_| FqTarget::empty(builder))
            .collect_vec()
            .try_into()
            .unwrap();
        Fq6Target { coeffs }
    }

    pub fn new(coeffs: Vec<FqTarget<F, D>>) -> Self {
        Fq6Target {
            coeffs: coeffs.try_into().unwrap(),
        }
    }

    pub fn connect(builder: &mut CircuitBuilder<F, D>, lhs: &Self, rhs: &Self) {
        for i in 0..6 {
            builder.connect_nonnative(&lhs.coeffs[i].target, &rhs.coeffs[i].target);
        }
    }

    pub fn select(
        builder: &mut CircuitBuilder<F, D>,
        a: &Self,
        b: &Self,
        flag: &BoolTarget,
    ) -> Self {
        let selected = a
            .coeffs
            .iter()
            .zip(b.coeffs.iter())
            .map(|(a, b)| FqTarget::select(builder, a, b, flag))
            .collect_vec();

        Self {
            coeffs: selected.try_into().unwrap(),
        }
    }

    pub fn constant(builder: &mut CircuitBuilder<F, D>, c: Fq6) -> Self {
        let c: MyFq6 = c.into();
        let coeffs = c
            .coeffs
            .iter()
            .map(|x| FqTarget::constant(builder, *x))
            .collect_vec()
            .try_into()
            .unwrap();
        Self { coeffs }
    }

    // this method fails if self is zero
    pub fn inv(&self, builder: &mut CircuitBuilder<F, D>) -> Self {
        let inv = Self::empty(builder);
        builder.add_simple_generator(Fq6InverseGenerator::<F, D> {
            x: self.clone(),
            inv: inv.clone(),
        });
        let one = Self::constant(builder, Fq6::ONE);
        let x_mul_inv = self.mul(builder, &inv);
        Self::connect(builder, &x_mul_inv, &one);
        inv
    }

    pub fn add(&self, builder: &mut CircuitBuilder<F, D>, rhs: &Self) -> Self {
        let coeffs = self
            .coeffs
            .iter()
            .enumerate()
            .map(|(i, x)| x.add(builder, &rhs.coeffs[i]))
            .collect_vec()
            .try_into()
            .unwrap();
        Fq6Target { coeffs }
    }

    pub fn neg(&self, builder: &mut CircuitBuilder<F, D>) -> Self {
        let coeffs = self
            .coeffs
            .iter()
            .map(|x| x.neg(builder))
            .collect_vec()
            .try_into()
            .unwrap();
        Fq6Target { coeffs }
    }

    pub fn sub(&self, builder: &mut CircuitBuilder<F, D>, rhs: &Self) -> Self {
        let coeffs = self
            .coeffs
            .iter()
            .enumerate()
            .map(|(i, x)| x.sub(builder, &rhs.coeffs[i]))
            .collect_vec()
            .try_into()
            .unwrap();
        Fq6Target { coeffs }
    }

    pub fn mul(&self, builder: &mut CircuitBuilder<F, D>, rhs: &Self) -> Self {
        let a = self;
        let b = rhs;
        let mut a0b0_coeffs: Vec<FqTarget<F, D>> = Vec::with_capacity(5);
        let mut a0b1_coeffs: Vec<FqTarget<F, D>> = Vec::with_capacity(5);
        let mut a1b0_coeffs: Vec<FqTarget<F, D>> = Vec::with_capacity(5);
        let mut a1b1_coeffs: Vec<FqTarget<F, D>> = Vec::with_capacity(5);
        for i in 0..3 {
            for j in 0..3 {
                let coeff00 = a.coeffs[i].mul(builder, &b.coeffs[j]);
                let coeff01 = a.coeffs[i].mul(builder, &b.coeffs[j + 3]);
                let coeff10 = a.coeffs[i + 3].mul(builder, &b.coeffs[j]);
                let coeff11 = a.coeffs[i + 3].mul(builder, &b.coeffs[j + 3]);
                if i + j < a0b0_coeffs.len() {
                    a0b0_coeffs[i + j] = a0b0_coeffs[i + j].add(builder, &coeff00);
                    a0b1_coeffs[i + j] = a0b1_coeffs[i + j].add(builder, &coeff01);
                    a1b0_coeffs[i + j] = a1b0_coeffs[i + j].add(builder, &coeff10);
                    a1b1_coeffs[i + j] = a1b1_coeffs[i + j].add(builder, &coeff11);
                } else {
                    a0b0_coeffs.push(coeff00);
                    a0b1_coeffs.push(coeff01);
                    a1b0_coeffs.push(coeff10);
                    a1b1_coeffs.push(coeff11);
                }
            }
        }

        let mut a0b0_minus_a1b1: Vec<FqTarget<F, D>> = Vec::with_capacity(5);
        let mut a0b1_plus_a1b0: Vec<FqTarget<F, D>> = Vec::with_capacity(5);
        for i in 0..5 {
            let a0b0_minus_a1b1_entry = a0b0_coeffs[i].sub(builder, &a1b1_coeffs[i]);
            let a0b1_plus_a1b0_entry = a0b1_coeffs[i].add(builder, &a1b0_coeffs[i]);
            a0b0_minus_a1b1.push(a0b0_minus_a1b1_entry);
            a0b1_plus_a1b0.push(a0b1_plus_a1b0_entry);
        }

        let const_one = FqTarget::constant(builder, Fq::from(1));
        let mut out_coeffs: Vec<FqTarget<F, D>> = Vec::with_capacity(6);
        for i in 0..3 {
            if i < 2 {
                let term0 = a0b0_minus_a1b1[i].clone();
                let term1 = a0b0_minus_a1b1[i + 3].mul(builder, &const_one);
                let term2 = a0b1_plus_a1b0[i + 3].neg(builder);
                let term0_plus_term1 = term0.add(builder, &term1);
                let coeff = term0_plus_term1.add(builder, &term2);
                out_coeffs.push(coeff);
            } else {
                out_coeffs.push(a0b0_minus_a1b1[i].clone());
            }
        }
        for i in 0..3 {
            if i < 2 {
                let term0 = a0b1_plus_a1b0[i].clone();
                let term1 = a0b0_minus_a1b1[i + 3].clone();
                let term2 = a0b1_plus_a1b0[i + 3].mul(builder, &const_one);
                let term0_plus_term1 = term0.add(builder, &term1);
                let coeff = term0_plus_term1.add(builder, &term2);
                out_coeffs.push(coeff);
            } else {
                out_coeffs.push(a0b1_plus_a1b0[i].clone());
            }
        }
        Self {
            coeffs: out_coeffs.try_into().unwrap(),
        }
    }

    pub fn mul_by_01(
        self,
        builder: &mut CircuitBuilder<F, D>,
        c0: &Fq2Target<F, D>,
        c1: &Fq2Target<F, D>,
    ) -> Self {
        let fq6_c00 = &self.coeffs[0];
        let fq6_c10 = &self.coeffs[1];
        let fq6_c20 = &self.coeffs[2];
        let fq6_c01 = &self.coeffs[3];
        let fq6_c11 = &self.coeffs[4];
        let fq6_c21 = &self.coeffs[5];

        let fq6_c0 = Fq2Target::new(vec![fq6_c00.clone(), fq6_c01.clone()]);
        let fq6_c1 = Fq2Target::new(vec![fq6_c10.clone(), fq6_c11.clone()]);
        let fq6_c2 = Fq2Target::new(vec![fq6_c20.clone(), fq6_c21.clone()]);

        let a_a = fq6_c0.mul(builder, c0);
        let b_b = fq6_c1.mul(builder, c1);

        let t1 = fq6_c2.mul(builder, c1);
        let t1 = t1.mul_by_nonresidue(builder);
        let t1 = t1.add(builder, &a_a);

        let c0_add_c1 = c0.add(builder, c1);
        let fq6_c0_add_fq6_c1 = fq6_c0.add(builder, &fq6_c1);
        let t2 = c0_add_c1.mul(builder, &fq6_c0_add_fq6_c1);
        let t2 = t2.sub(builder, &a_a);
        let t2 = t2.sub(builder, &b_b);

        let t3 = fq6_c2.mul(builder, c0);
        let t3 = t3.add(builder, &b_b);

        Self::new(vec![
            t1.coeffs[0].clone(),
            t2.coeffs[0].clone(),
            t3.coeffs[0].clone(),
            t1.coeffs[1].clone(),
            t2.coeffs[1].clone(),
            t3.coeffs[1].clone(),
        ])
    }

    pub fn mul_by_1(self, builder: &mut CircuitBuilder<F, D>, c1: &Fq2Target<F, D>) -> Self {
        let fq6_c00 = &self.coeffs[0];
        let fq6_c10 = &self.coeffs[1];
        let fq6_c20 = &self.coeffs[2];
        let fq6_c01 = &self.coeffs[3];
        let fq6_c11 = &self.coeffs[4];
        let fq6_c21 = &self.coeffs[5];

        let fq6_c0 = Fq2Target::new(vec![fq6_c00.clone(), fq6_c01.clone()]);
        let fq6_c1 = Fq2Target::new(vec![fq6_c10.clone(), fq6_c11.clone()]);
        let fq6_c2 = Fq2Target::new(vec![fq6_c20.clone(), fq6_c21.clone()]);

        let b_b = fq6_c1.clone();
        let b_b = b_b.mul(builder, c1);

        let t1 = c1;
        let tmp = fq6_c1.clone();
        let tmp = tmp.add(builder, &fq6_c2);

        let t1 = t1.mul(builder, &tmp);
        let t1 = t1.sub(builder, &b_b);
        let t1 = t1.mul_by_nonresidue(builder);

        let t2 = c1;
        let tmp = fq6_c0;
        let tmp = tmp.add(builder, &fq6_c1);

        let t2 = t2.mul(builder, &tmp);
        let t2 = t2.sub(builder, &b_b);

        Self::new(vec![
            t1.coeffs[0].clone(),
            t2.coeffs[0].clone(),
            b_b.coeffs[0].clone(),
            t1.coeffs[1].clone(),
            t2.coeffs[1].clone(),
            b_b.coeffs[1].clone(),
        ])
    }

    /// Multiply by quadratic nonresidue v.
    pub fn mul_by_nonresidue(&self, builder: &mut CircuitBuilder<F, D>) -> Self {
        // Given a + bv + cv^2, this produces
        //     av + bv^2 + cv^3
        // but because v^3 = u + 1, we have
        //     c(u + 1) + av + v^2

        let fq6_c00 = &self.coeffs[0];
        let fq6_c10 = &self.coeffs[1];
        let fq6_c20 = &self.coeffs[2];
        let fq6_c01 = &self.coeffs[3];
        let fq6_c11 = &self.coeffs[4];
        let fq6_c21 = &self.coeffs[5];

        let fq6_c0 = Fq2Target::new(vec![fq6_c00.clone(), fq6_c01.clone()]);
        let fq6_c1 = Fq2Target::new(vec![fq6_c10.clone(), fq6_c11.clone()]);
        let fq6_c2 = Fq2Target::new(vec![fq6_c20.clone(), fq6_c21.clone()]);

        let c0 = fq6_c2.mul_by_nonresidue(builder);
        let c1 = fq6_c0;
        let c2 = fq6_c1;

        Self::new(vec![
            c0.coeffs[0].clone(),
            c1.coeffs[0].clone(),
            c2.coeffs[0].clone(),
            c0.coeffs[1].clone(),
            c1.coeffs[1].clone(),
            c2.coeffs[1].clone(),
        ])
    }

    pub fn conditional_mul(
        &self,
        builder: &mut CircuitBuilder<F, D>,
        x: &Self,
        flag: &BoolTarget,
    ) -> Self {
        let muled = self.mul(builder, x);
        Self::select(builder, &muled, self, flag)
    }

    // pub fn div(&self, builder: &mut CircuitBuilder<F, D>, other: &Self) -> Self {
    // }

    // pub fn inv(&self, builder: &mut CircuitBuilder<F, D>) -> Self {
    // }

    // pub fn confugate(&self, builder: &mut CircuitBuilder<F, D>) -> Self {
    // }
}

#[derive(Debug)]
struct Fq6InverseGenerator<F: RichField + Extendable<D>, const D: usize> {
    x: Fq6Target<F, D>,
    inv: Fq6Target<F, D>,
}

impl<F: RichField + Extendable<D>, const D: usize> SimpleGenerator<F, D>
    for Fq6InverseGenerator<F, D>
{
    fn dependencies(&self) -> Vec<Target> {
        self.x
            .coeffs
            .iter()
            .flat_map(|coeff| coeff.target.value.limbs.iter().map(|&l| l.0))
            .collect_vec()
    }

    fn run_once(&self, witness: &PartitionWitness<F>, out_buffer: &mut GeneratedValues<F>) {
        let coeffs: Vec<Fq> = self
            .x
            .coeffs
            .iter()
            .map(|x| from_biguint_to_fq(witness.get_biguint_target(x.target.value.clone())))
            .collect_vec();
        let x = MyFq6 {
            coeffs: coeffs.try_into().unwrap(),
        };
        let x: Fq6 = x.into();
        let inv_x: Fq6 = x.inverse().unwrap();
        let inv_x: MyFq6 = inv_x.into();
        let inv_x_biguint: Vec<BigUint> = inv_x
            .coeffs
            .iter()
            .cloned()
            .map(|coeff| coeff.into())
            .collect_vec();

        for i in 0..6 {
            out_buffer.set_biguint_target(&self.inv.coeffs[i].target.value, &inv_x_biguint[i]);
        }
    }

    fn id(&self) -> std::string::String {
        "Fq6InverseGenerator".to_string()
    }

    fn serialize(
        &self,
        _dst: &mut Vec<u8>,
        _common_data: &plonky2::plonk::circuit_data::CommonCircuitData<F, D>,
    ) -> plonky2::util::serialization::IoResult<()> {
        todo!()
    }

    fn deserialize(
        _src: &mut Buffer,
        _common_data: &plonky2::plonk::circuit_data::CommonCircuitData<F, D>,
    ) -> plonky2::util::serialization::IoResult<Self>
    where
        Self: Sized,
    {
        todo!()
    }
}

#[cfg(test)]
mod tests {
    use ark_bls12_381::{Fq12Config, Fq2, Fq6};
    use ark_ff::{Field, Fp12Config, UniformRand};
    use plonky2::{
        field::goldilocks_field::GoldilocksField,
        iop::witness::PartialWitness,
        plonk::{
            circuit_builder::CircuitBuilder, circuit_data::CircuitConfig,
            config::PoseidonGoldilocksConfig,
        },
    };

    use super::Fq6Target;
    use crate::fields::fq2_target::Fq2Target;

    type F = GoldilocksField;
    type C = PoseidonGoldilocksConfig;
    const D: usize = 2;

    #[test]
    fn test_fq6_inv_circuit() {
        let rng = &mut rand::thread_rng();
        let x: Fq6 = Fq6::rand(rng);
        let inv_x_expected = x.inverse().unwrap();

        let config = CircuitConfig::wide_ecc_config();
        let mut builder = CircuitBuilder::<F, D>::new(config);
        let x_t = Fq6Target::constant(&mut builder, x);
        let inv_x_t = x_t.inv(&mut builder);
        let inv_x_expected_t = Fq6Target::constant(&mut builder, inv_x_expected);

        Fq6Target::connect(&mut builder, &inv_x_t, &inv_x_expected_t);

        let pw = PartialWitness::new();
        let data = builder.build::<C>();
        let _proof = data.prove(pw);
    }

    #[test]
    fn test_mul_by_01() {
        let rng = &mut rand::thread_rng();
        let x: Fq6 = Fq6::rand(rng);
        let c0: Fq2 = Fq2::rand(rng);
        let c1: Fq2 = Fq2::rand(rng);

        let config = CircuitConfig::wide_ecc_config();
        let mut builder = CircuitBuilder::<F, D>::new(config);
        let x_t = Fq6Target::constant(&mut builder, x);
        let x_c0 = Fq2Target::constant(&mut builder, c0);
        let x_c1 = Fq2Target::constant(&mut builder, c1);
        let mut x = x;
        x.mul_by_01(&c0, &c1);

        let x_t = x_t.mul_by_01(&mut builder, &x_c0, &x_c1);

        let x_expected_t = Fq6Target::constant(&mut builder, x);

        Fq6Target::connect(&mut builder, &x_t, &x_expected_t);

        let pw = PartialWitness::new();
        let data = builder.build::<C>();
        let _proof = data.prove(pw);
    }

    #[test]
    fn test_mul_by_1() {
        let rng = &mut rand::thread_rng();
        let x: Fq6 = Fq6::rand(rng);
        let c1: Fq2 = Fq2::rand(rng);

        let config = CircuitConfig::wide_ecc_config();
        let mut builder = CircuitBuilder::<F, D>::new(config);
        let x_t = Fq6Target::constant(&mut builder, x);
        let x_c1 = Fq2Target::constant(&mut builder, c1);
        let mut x = x;
        x.mul_by_1(&c1);

        let x_t = x_t.mul_by_1(&mut builder, &x_c1);

        let x_expected_t = Fq6Target::constant(&mut builder, x);

        Fq6Target::connect(&mut builder, &x_t, &x_expected_t);

        let pw = PartialWitness::new();
        let data = builder.build::<C>();
        let _proof = data.prove(pw);
    }

    #[test]
    fn test_mul_fq6_by_nonresidue_in_place() {
        let rng = &mut rand::thread_rng();
        let x: Fq6 = Fq6::rand(rng);

        let config = CircuitConfig::wide_ecc_config();
        let mut builder = CircuitBuilder::<F, D>::new(config);
        let x_t = Fq6Target::constant(&mut builder, x);
        let mut x = x;
        Fq12Config::mul_fp6_by_nonresidue_in_place(&mut x);

        let x_t = x_t.mul_by_nonresidue(&mut builder);

        let x_expected_t = Fq6Target::constant(&mut builder, x);

        Fq6Target::connect(&mut builder, &x_t, &x_expected_t);

        let pw = PartialWitness::new();
        let data = builder.build::<C>();
        let _proof = data.prove(pw);
    }
}
