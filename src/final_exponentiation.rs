use plonky2::{
    field::extension::Extendable,
    hash::hash_types::RichField,
    plonk::circuit_builder::CircuitBuilder,
};

use crate::{
    fields::fq12_target::Fq12Target,
    final_exponentiation_helpers::{frobenius_map, pow_target},
    utils::constants::BLS_X,
};

pub fn easy_part<F: RichField + Extendable<D>, const D: usize>(
    builder: &mut CircuitBuilder<F, D>,
    a: &Fq12Target<F, D>,
) -> Fq12Target<F, D> {
    let f1 = a.conjugate(builder);
    let f2 = f1.div(builder, a);
    let f3 = frobenius_map(builder, &f2, 2);

    f3.mul(builder, &f2)
}

pub fn hard_part_target<F: RichField + Extendable<D>, const D: usize>(
    builder: &mut CircuitBuilder<F, D>,
    r: Fq12Target<F, D>,
) -> Fq12Target<F, D> {
    let mut y0 = r.mul(builder, &r); // optimize square
    let mut y1 = pow_target(builder, r.clone(), vec![BLS_X]);
    y1 = y1.conjugate(builder);
    let mut y2 = r.clone();
    y2 = y2.conjugate(builder);

    y1 = y1.mul(builder, &y2);
    y2 = pow_target(builder, y1.clone(), vec![BLS_X]);
    y2 = y2.conjugate(builder);
    y1 = y1.conjugate(builder);
    y1 = y1.mul(builder, &y2);
    y2 = pow_target(builder, y1.clone(), vec![BLS_X]);
    y2 = y2.conjugate(builder);
    y1 = frobenius_map(builder, &y1, 1);
    y1 = y1.mul(builder, &y2);
    let r = r.mul(builder, &y0);
    y0 = pow_target(builder, y1.clone(), vec![BLS_X]);
    y0 = y0.conjugate(builder);
    y2 = pow_target(builder, y0, vec![BLS_X]);
    y2 = y2.conjugate(builder);
    y0 = y1.clone();
    y0 = frobenius_map(builder, &y0, 2);
    y1 = y1.conjugate(builder);
    y1 = y1.mul(builder, &y2);
    y1 = y1.mul(builder, &y0);

    r.mul(builder, &y1)
}

pub fn final_exponentiation<F: RichField + Extendable<D>, const D: usize>(
    builder: &mut CircuitBuilder<F, D>,
    a: Fq12Target<F, D>,
) -> Fq12Target<F, D> {
    let f0 = easy_part(builder, &a);

    hard_part_target(builder, f0)
}

#[cfg(test)]
mod tests {
    use std::time::Instant;

    use ark_bls12_381::Fq12;
    use ark_ec::pairing::{MillerLoopOutput, Pairing};
    use ark_ff::UniformRand;
    use plonky2::{
        field::goldilocks_field::GoldilocksField,
        iop::witness::PartialWitness,
        plonk::{
            circuit_builder::CircuitBuilder,
            circuit_data::CircuitConfig,
            config::PoseidonGoldilocksConfig,
        },
    };

    use crate::{fields::fq12_target::Fq12Target, final_exponentiation::final_exponentiation};

    type F = GoldilocksField;
    type C = PoseidonGoldilocksConfig;
    const D: usize = 2;

    #[test]
    fn test_final_exponentiation_circuit() {
        let config = CircuitConfig::standard_recursion_config();
        let mut builder = CircuitBuilder::<F, D>::new(config);
        let rng = &mut rand::thread_rng();
        let x = Fq12::rand(rng);
        let input_t = Fq12Target::constant(&mut builder, x);

        let output = final_exponentiation::<F, D>(&mut builder, input_t);
        let output_expected = ark_bls12_381::Bls12_381::final_exponentiation(MillerLoopOutput(x))
            .unwrap()
            .0;
        let output_expected_t = Fq12Target::constant(&mut builder, output_expected.into());

        Fq12Target::connect(&mut builder, &output, &output_expected_t);

        let now = Instant::now();
        let pw = PartialWitness::new();
        let data = builder.build::<C>();
        let _proof = data.prove(pw);
        println!("time: {:?}", now.elapsed());
    }
}
