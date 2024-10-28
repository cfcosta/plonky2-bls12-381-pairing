use ark_bls12_381::{Fq, Fq12};
use ark_ff::Field;
use itertools::Itertools;
use num_bigint::BigUint;
use plonky2::{
    field::extension::Extendable,
    hash::hash_types::RichField,
    iop::{
        generator::{GeneratedValues, SimpleGenerator},
        target::{BoolTarget, Target},
        witness::{PartitionWitness, WitnessWrite},
    },
    plonk::circuit_builder::CircuitBuilder,
    util::serialization::Buffer,
};
use plonky2_ecdsa::gadgets::{
    biguint::{GeneratedValuesBigUint, WitnessBigUint},
    nonnative::CircuitBuilderNonNative,
};

use super::{fq2_target::Fq2Target, fq6_target::Fq6Target, fq_target::FqTarget};
use crate::utils::helpers::{from_biguint_to_fq, MyFq12};

#[derive(Debug, Clone)]
pub struct Fq12Target<F: RichField + Extendable<D>, const D: usize> {
    pub coeffs: [FqTarget<F, D>; 12],
}

impl<F: RichField + Extendable<D>, const D: usize> Fq12Target<F, D> {
    pub fn empty(builder: &mut CircuitBuilder<F, D>) -> Self {
        let coeffs = [(); 12]
            .iter()
            .map(|_| FqTarget::empty(builder))
            .collect_vec()
            .try_into()
            .unwrap();
        Fq12Target { coeffs }
    }

    pub fn new(coeffs: Vec<FqTarget<F, D>>) -> Self {
        Fq12Target {
            coeffs: coeffs.try_into().unwrap(),
        }
    }

    pub fn connect(builder: &mut CircuitBuilder<F, D>, lhs: &Self, rhs: &Self) {
        for i in 0..12 {
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

    pub fn constant(builder: &mut CircuitBuilder<F, D>, c: Fq12) -> Self {
        let c: MyFq12 = c.into();
        let coeffs = c
            .coeffs
            .iter()
            .map(|x| FqTarget::constant(builder, *x))
            .collect_vec()
            .try_into()
            .unwrap();
        Self { coeffs }
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
        Fq12Target { coeffs }
    }

    pub fn neg(&self, builder: &mut CircuitBuilder<F, D>) -> Self {
        let coeffs = self
            .coeffs
            .iter()
            .map(|x| x.neg(builder))
            .collect_vec()
            .try_into()
            .unwrap();
        Fq12Target { coeffs }
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
        Fq12Target { coeffs }
    }

    pub fn mul(&self, builder: &mut CircuitBuilder<F, D>, rhs: &Self) -> Self {
        let a = self;
        let b = rhs;
        let mut a0b0_coeffs: Vec<FqTarget<F, D>> = Vec::with_capacity(11);
        let mut a0b1_coeffs: Vec<FqTarget<F, D>> = Vec::with_capacity(11);
        let mut a1b0_coeffs: Vec<FqTarget<F, D>> = Vec::with_capacity(11);
        let mut a1b1_coeffs: Vec<FqTarget<F, D>> = Vec::with_capacity(11);
        for i in 0..6 {
            for j in 0..6 {
                let coeff00 = a.coeffs[i].mul(builder, &b.coeffs[j]);
                let coeff01 = a.coeffs[i].mul(builder, &b.coeffs[j + 6]);
                let coeff10 = a.coeffs[i + 6].mul(builder, &b.coeffs[j]);
                let coeff11 = a.coeffs[i + 6].mul(builder, &b.coeffs[j + 6]);
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

        let mut a0b0_minus_a1b1: Vec<FqTarget<F, D>> = Vec::with_capacity(11);
        let mut a0b1_plus_a1b0: Vec<FqTarget<F, D>> = Vec::with_capacity(11);
        for i in 0..11 {
            let a0b0_minus_a1b1_entry = a0b0_coeffs[i].sub(builder, &a1b1_coeffs[i]);
            let a0b1_plus_a1b0_entry = a0b1_coeffs[i].add(builder, &a1b0_coeffs[i]);
            a0b0_minus_a1b1.push(a0b0_minus_a1b1_entry);
            a0b1_plus_a1b0.push(a0b1_plus_a1b0_entry);
        }

        let const_one = FqTarget::constant(builder, Fq::from(1));
        let mut out_coeffs: Vec<FqTarget<F, D>> = Vec::with_capacity(12);
        for i in 0..6 {
            if i < 5 {
                // let coeff: Fq = a0b0_minus_a1b1[i] + Fq::from(1) * a0b0_minus_a1b1[i + 6]
                //     - a0b1_plus_a1b0[i + 6];
                let term0 = a0b0_minus_a1b1[i].clone();
                let term1 = a0b0_minus_a1b1[i + 6].mul(builder, &const_one);
                let term2 = a0b1_plus_a1b0[i + 6].neg(builder);
                let term0_plus_term1 = term0.add(builder, &term1);
                let coeff = term0_plus_term1.add(builder, &term2);
                out_coeffs.push(coeff);
            } else {
                out_coeffs.push(a0b0_minus_a1b1[i].clone());
            }
        }
        for i in 0..6 {
            if i < 5 {
                // let coeff: Fq = a0b1_plus_a1b0[i]
                //     + a0b0_minus_a1b1[i + 6]
                //     + Fq::from(1) * a0b1_plus_a1b0[i + 6];
                let term0 = a0b1_plus_a1b0[i].clone();
                let term1 = a0b0_minus_a1b1[i + 6].clone();
                let term2 = a0b1_plus_a1b0[i + 6].mul(builder, &const_one);
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

    pub fn convert_to_fq6(&self) -> (Fq6Target<F, D>, Fq6Target<F, D>) {
        //                   Fq12
        //                  /    \
        //                 /      \
        //                /        \
        //               /          \
        //              /            \
        //             /              \
        //            /                \
        //           c0                c1
        //        /  |  \            /  |  \
        //       /   |   \          /   |   \
        //      c0   c1   c2       c0   c1   c2
        //   c0 c1 c0 c1 c0 c1  c0 c1  c0 c1 c0 c1

        // c000 = 0 c100 = 1 c010 = 2 c110 = 3 c020 = 4 c120 = 5
        // c001 = 6 c101 = 7 c011 = 8 c111 = 9 c021 = 10 c121 = 11

        let c000 = &self.coeffs[0];
        let c001 = &self.coeffs[6];
        let c010 = &self.coeffs[2];
        let c011 = &self.coeffs[8];
        let c020 = &self.coeffs[4];
        let c021 = &self.coeffs[10];
        let c100 = &self.coeffs[1];
        let c101 = &self.coeffs[7];
        let c110 = &self.coeffs[3];
        let c111 = &self.coeffs[9];
        let c120 = &self.coeffs[5];
        let c121 = &self.coeffs[11];

        (
            Fq6Target::new(vec![
                c000.clone(),
                c010.clone(),
                c020.clone(),
                c001.clone(),
                c011.clone(),
                c021.clone(),
            ]),
            Fq6Target::new(vec![
                c100.clone(),
                c110.clone(),
                c120.clone(),
                c101.clone(),
                c111.clone(),
                c121.clone(),
            ]),
        )
    }

    pub fn from_fq6(c0: Fq6Target<F, D>, c1: Fq6Target<F, D>) -> Self {
        Self::new(vec![
            c0.coeffs[0].clone(),
            c1.coeffs[0].clone(),
            c0.coeffs[1].clone(),
            c1.coeffs[1].clone(),
            c0.coeffs[2].clone(),
            c1.coeffs[2].clone(),
            c0.coeffs[3].clone(),
            c1.coeffs[3].clone(),
            c0.coeffs[4].clone(),
            c1.coeffs[4].clone(),
            c0.coeffs[5].clone(),
            c1.coeffs[5].clone(),
        ])
    }

    pub fn mul_by_014(
        &self,
        builder: &mut CircuitBuilder<F, D>,
        c0: &Fq2Target<F, D>,
        c1: &Fq2Target<F, D>,
        c4: &Fq2Target<F, D>,
    ) -> Self {
        let (self_c0, self_c1) = self.convert_to_fq6();

        let aa = self_c0.clone();
        let aa = aa.mul_by_01(builder, c0, c1);
        let bb = self_c1.clone();
        let bb = bb.mul_by_1(builder, c4);
        let o = c1;
        let o = o.add(builder, c4);
        let self_c1 = self_c1.add(builder, &self_c0);
        let self_c1 = self_c1.mul_by_01(builder, c0, &o);
        let self_c1 = self_c1.sub(builder, &aa);
        let self_c1 = self_c1.sub(builder, &bb);
        let self_c0 = bb;
        let self_c0 = self_c0.mul_by_nonresidue(builder);
        let self_c0 = self_c0.add(builder, &aa);

        Self::new(vec![
            self_c0.coeffs[0].clone(),
            self_c1.coeffs[0].clone(),
            self_c0.coeffs[1].clone(),
            self_c1.coeffs[1].clone(),
            self_c0.coeffs[2].clone(),
            self_c1.coeffs[2].clone(),
            self_c0.coeffs[3].clone(),
            self_c1.coeffs[3].clone(),
            self_c0.coeffs[4].clone(),
            self_c1.coeffs[4].clone(),
            self_c0.coeffs[5].clone(),
            self_c1.coeffs[5].clone(),
        ])
    }

    pub fn div(&self, builder: &mut CircuitBuilder<F, D>, other: &Self) -> Self {
        let inv = other.inv(builder);
        self.mul(builder, &inv)
    }

    pub fn inv(&self, builder: &mut CircuitBuilder<F, D>) -> Self {
        let inv = Self::empty(builder);
        builder.add_simple_generator(Fq12InverseGenerator::<F, D> {
            x: self.clone(),
            inv: inv.clone(),
        });
        let one = Self::constant(builder, Fq12::ONE);
        let x_mul_inv = self.mul(builder, &inv);
        Self::connect(builder, &x_mul_inv, &one);
        inv
    }

    pub fn conjugate(&self, builder: &mut CircuitBuilder<F, D>) -> Self {
        let mut coeffs = self.coeffs.clone();
        coeffs[1] = coeffs[1].neg(builder);
        coeffs[3] = coeffs[3].neg(builder);
        coeffs[5] = coeffs[5].neg(builder);
        coeffs[7] = coeffs[7].neg(builder);
        coeffs[9] = coeffs[9].neg(builder);
        coeffs[11] = coeffs[11].neg(builder);
        Self { coeffs }
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

    /// Iterates over the entire iterator, multiplying all the elements
    pub fn multiply_elements(
        builder: &mut CircuitBuilder<F, D>,
        iter: impl Iterator<Item = Self>,
    ) -> Option<Self> {
        let mut result: Option<Self> = None;

        for item in iter {
            match result {
                Some(val) => {
                    result = Some(val.mul(builder, &item));
                }
                None => {
                    result = Some(item);
                }
            }
        }

        result
    }
}

#[derive(Debug)]
struct Fq12InverseGenerator<F: RichField + Extendable<D>, const D: usize> {
    x: Fq12Target<F, D>,
    inv: Fq12Target<F, D>,
}

impl<F: RichField + Extendable<D>, const D: usize> SimpleGenerator<F, D>
    for Fq12InverseGenerator<F, D>
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
        let x = MyFq12 {
            coeffs: coeffs.try_into().unwrap(),
        };
        let x: Fq12 = x.into();
        let inv_x: Fq12 = x.inverse().unwrap();
        let inv_x: MyFq12 = inv_x.into();
        let inv_x_biguint: Vec<BigUint> = inv_x
            .coeffs
            .iter()
            .cloned()
            .map(|coeff| coeff.into())
            .collect_vec();

        for i in 0..12 {
            out_buffer.set_biguint_target(&self.inv.coeffs[i].target.value, &inv_x_biguint[i]);
        }
    }

    fn id(&self) -> std::string::String {
        "Fq12InverseGenerator".to_string()
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

impl<F: RichField + Extendable<D>, const D: usize> Fq12Target<F, D> {
    pub fn to_vec(&self) -> Vec<Target> {
        self.coeffs.iter().flat_map(|c| c.to_vec()).collect()
    }

    pub fn from_vec(builder: &mut CircuitBuilder<F, D>, input: &[Target]) -> Self {
        let num_limbs = 12;
        assert_eq!(input.len(), 12 * num_limbs);
        let coeffs = input
            .iter()
            .cloned()
            .chunks(num_limbs)
            .into_iter()
            .map(|chunk| FqTarget::from_vec(builder, &chunk.collect_vec()))
            .collect_vec();
        Fq12Target {
            coeffs: coeffs.try_into().unwrap(),
        }
    }

    pub fn set_witness<W: WitnessWrite<F>>(&self, pw: &mut W, value: &Fq12) {
        let my_value: MyFq12 = (*value).into();
        self.coeffs
            .iter()
            .cloned()
            .zip(my_value.coeffs)
            .map(|(c_t, c)| c_t.set_witness(pw, &c))
            .for_each(drop);
    }
}

#[cfg(test)]
mod tests {
    use ark_bls12_381::{Fq, Fq12, Fq2};
    use ark_ff::Field;
    use ark_std::UniformRand;
    use num_bigint::BigUint;
    use plonky2::{
        field::goldilocks_field::GoldilocksField,
        iop::witness::PartialWitness,
        plonk::{
            circuit_builder::CircuitBuilder, circuit_data::CircuitConfig,
            config::PoseidonGoldilocksConfig,
        },
    };

    use super::{from_biguint_to_fq, Fq12Target};
    use crate::fields::{fq2_target::Fq2Target, fq6_target::Fq6Target};

    type F = GoldilocksField;
    type C = PoseidonGoldilocksConfig;
    const D: usize = 2;

    #[test]
    fn test_from_to_vec() {
        let rng = &mut rand::thread_rng();
        let a = Fq12::rand(rng);
        let config = CircuitConfig::wide_ecc_config();
        let mut builder = CircuitBuilder::<F, D>::new(config);
        let a_t = Fq12Target::constant(&mut builder, a);

        let a_vec = a_t.to_vec();
        let restored_a_t = Fq12Target::from_vec(&mut builder, &a_vec);

        Fq12Target::connect(&mut builder, &a_t, &restored_a_t);

        let pw = PartialWitness::new();
        let data = builder.build::<C>();
        let _proof = data.prove(pw);
    }

    #[test]
    fn test_from_biguint_to_fq() {
        let rng = &mut rand::thread_rng();
        let x = Fq::rand(rng);
        let x_biguint: BigUint = x.into();
        let converted_x = from_biguint_to_fq(x_biguint);
        assert_eq!(x, converted_x);
    }

    #[test]
    fn test_fq12_mul_circuit() {
        let rng = &mut rand::thread_rng();
        let a = Fq12::rand(rng);
        let b = Fq12::rand(rng);
        let c_expected = a * b;

        let config = CircuitConfig::wide_ecc_config();
        let mut builder = CircuitBuilder::<F, D>::new(config);
        let a_t = Fq12Target::constant(&mut builder, a);
        let b_t = Fq12Target::constant(&mut builder, b);
        let c_t = a_t.mul(&mut builder, &b_t);
        let c_expected_t = Fq12Target::constant(&mut builder, c_expected);

        Fq12Target::connect(&mut builder, &c_expected_t, &c_t);

        let pw = PartialWitness::new();
        let data = builder.build::<C>();
        let _proof = data.prove(pw);
    }

    #[test]
    fn test_fq12_inv_circuit() {
        let rng = &mut rand::thread_rng();
        let x: Fq12 = Fq12::rand(rng);
        let inv_x_expected = x.inverse().unwrap();

        let config = CircuitConfig::wide_ecc_config();
        let mut builder = CircuitBuilder::<F, D>::new(config);
        let x_t = Fq12Target::constant(&mut builder, x);
        let inv_x_t = x_t.inv(&mut builder);
        let inv_x_expected_t = Fq12Target::constant(&mut builder, inv_x_expected);

        Fq12Target::connect(&mut builder, &inv_x_t, &inv_x_expected_t);

        let pw = PartialWitness::new();
        let data = builder.build::<C>();
        let _proof = data.prove(pw);
    }

    #[test]
    fn test_mul_by_014() {
        let rng = &mut rand::thread_rng();
        let x: Fq12 = Fq12::rand(rng);
        let c0: Fq2 = Fq2::rand(rng);
        let c1: Fq2 = Fq2::rand(rng);
        let c4: Fq2 = Fq2::rand(rng);
        let config = CircuitConfig::wide_ecc_config();
        let mut builder = CircuitBuilder::<F, D>::new(config);
        let c0_t = Fq2Target::constant(&mut builder, c0);
        let c1_t = Fq2Target::constant(&mut builder, c1);
        let c4_t = Fq2Target::constant(&mut builder, c4);
        let r_t = Fq12Target::constant(&mut builder, x);
        let mut r_expected = x;
        r_expected.mul_by_014(&c0, &c1, &c4);
        let r_t = r_t.mul_by_014(&mut builder, &c0_t, &c1_t, &c4_t);

        let r_expected_t = Fq12Target::constant(&mut builder, r_expected);

        Fq12Target::connect(&mut builder, &r_t, &r_expected_t);

        let pw = PartialWitness::new();
        let data = builder.build::<C>();
        let _proof = data.prove(pw);
    }

    #[test]
    fn test_from_fq6_to_fq12_casting() {
        let rng = &mut rand::thread_rng();
        let x: Fq12 = Fq12::rand(rng);

        let config = CircuitConfig::wide_ecc_config();
        let mut builder = CircuitBuilder::<F, D>::new(config);
        let c0 = Fq6Target::constant(&mut builder, x.c0);
        let c1 = Fq6Target::constant(&mut builder, x.c1);
        let r_t = Fq12Target::from_fq6(c0, c1);
        let r_expected_t = Fq12Target::constant(&mut builder, x);

        Fq12Target::connect(&mut builder, &r_t, &r_expected_t);

        let pw = PartialWitness::new();
        let data = builder.build::<C>();
        let _proof = data.prove(pw);
    }

    #[test]
    fn test_convert_to_fq6() {
        let rng = &mut rand::thread_rng();
        let x: Fq12 = Fq12::rand(rng);
        let config = CircuitConfig::wide_ecc_config();
        let mut builder = CircuitBuilder::<F, D>::new(config);

        let r_t = Fq12Target::constant(&mut builder, x);

        let (c0_t, c1_t) = r_t.convert_to_fq6();
        let c0_expected_t = Fq6Target::constant(&mut builder, x.c0);
        let c1_expected_t = Fq6Target::constant(&mut builder, x.c1);

        Fq6Target::connect(&mut builder, &c0_t, &c0_expected_t);
        Fq6Target::connect(&mut builder, &c1_t, &c1_expected_t);

        let pw = PartialWitness::new();
        let data = builder.build::<C>();
        let _proof = data.prove(pw);
    }

    #[test]
    fn test_conjugate() {
        let rng = &mut rand::thread_rng();
        let x: Fq12 = Fq12::rand(rng);
        let config = CircuitConfig::wide_ecc_config();
        let mut builder = CircuitBuilder::<F, D>::new(config);
        let x_t = Fq12Target::constant(&mut builder, x);

        let mut x = x;
        x.conjugate_in_place();
        let x_t_conjugate = x_t.conjugate(&mut builder);
        let x_expected_t = Fq12Target::constant(&mut builder, x);
        Fq12Target::connect(&mut builder, &x_t_conjugate, &x_expected_t);

        let pw = PartialWitness::new();
        let data = builder.build::<C>();
        let _proof = data.prove(pw);
    }

    #[test]
    fn test_multiply_elements() {
        let rng = &mut rand::thread_rng();
        let x: Fq12 = Fq12::rand(rng);
        let y: Fq12 = Fq12::rand(rng);
        let z: Fq12 = Fq12::rand(rng);
        let config = CircuitConfig::wide_ecc_config();
        let mut builder = CircuitBuilder::<F, D>::new(config);
        let x_t = Fq12Target::constant(&mut builder, x);
        let y_t = Fq12Target::constant(&mut builder, y);
        let z_t = Fq12Target::constant(&mut builder, z);

        let result = x * y * z;
        let result_t = Fq12Target::constant(&mut builder, result);

        let expected_result_array = [x_t, y_t, z_t].into_iter();
        let expected_result =
            Fq12Target::multiply_elements(&mut builder, expected_result_array.clone()).unwrap();

        Fq12Target::connect(&mut builder, &result_t, &expected_result);

        let pw = PartialWitness::new();
        let data = builder.build::<C>();
        let _proof = data.prove(pw);
    }
}
