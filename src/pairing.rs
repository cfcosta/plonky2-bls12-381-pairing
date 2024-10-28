use plonky2::{
    field::extension::Extendable, hash::hash_types::RichField,
    plonk::circuit_builder::CircuitBuilder,
};

use crate::{
    curves::{g1::G1PreparedTarget, g2::G2PreparedTarget},
    fields::fq12_target::Fq12Target,
    final_exponentiation::final_exponentiation,
    miller_loop::multi_miller_loop,
};

pub fn pairing<F: RichField + Extendable<D>, const D: usize>(
    builder: &mut CircuitBuilder<F, D>,
    a: impl IntoIterator<Item = impl Into<G1PreparedTarget<F, D>>>,
    b: impl IntoIterator<Item = impl Into<G2PreparedTarget<F, D>>>,
) -> Fq12Target<F, D> {
    let f = multi_miller_loop(builder, a, b);
    final_exponentiation::<F, D>(builder, f)
}

#[cfg(test)]
mod tests {
    use ark_bls12_381::{G1Affine, G2Affine};
    use ark_ec::pairing::Pairing;
    use ark_ff::UniformRand;
    use plonky2::{
        field::goldilocks_field::GoldilocksField,
        iop::witness::PartialWitness,
        plonk::{
            circuit_builder::CircuitBuilder, circuit_data::CircuitConfig,
            config::PoseidonGoldilocksConfig,
        },
    };

    use super::pairing;
    use crate::{
        curves::{
            g1::{G1AffineTarget, G1PreparedTarget},
            g2::{G2AffineTarget, G2PreparedTarget},
        },
        fields::fq12_target::Fq12Target,
    };

    type F = GoldilocksField;
    type C = PoseidonGoldilocksConfig;
    const D: usize = 2;

    #[test]
    fn test_pairing_circuit() {
        let mut config = CircuitConfig::standard_recursion_config();
        config.num_wires = 300;
        let mut builder = CircuitBuilder::<F, D>::new(config);
        let rng = &mut rand::thread_rng();
        let p = G1Affine::rand(rng);
        let q = G2Affine::rand(rng);

        let result_expected = ark_bls12_381::Bls12_381::pairing(p, q).0;

        let p_prepared_t = [G1PreparedTarget(G1AffineTarget::constant(&mut builder, p))];
        let q_t = G2AffineTarget::constant(&mut builder, q);
        let q_prepared_t = [G2PreparedTarget::from(&mut builder, q_t)];
        let result_circuit = pairing(&mut builder, p_prepared_t, q_prepared_t);

        let result_expected_t = Fq12Target::constant(&mut builder, result_expected);

        Fq12Target::connect(&mut builder, &result_expected_t, &result_circuit);

        let pw = PartialWitness::<F>::new();
        let data = builder.build::<C>();
        let _proof = data.prove(pw);
    }
}
