use std::ops::Div;

use ark_bls12_381::{Fq, Fq2};
use ark_ff::Field;
use ark_std::Zero;
use num::BigUint;
use num_traits::One;
use plonky2::{
    field::extension::Extendable, hash::hash_types::RichField,
    plonk::circuit_builder::CircuitBuilder,
};

use crate::fields::{fq12_target::Fq12Target, fq2_target::Fq2Target};

pub fn frobenius_map<F: RichField + Extendable<D>, const D: usize>(
    builder: &mut CircuitBuilder<F, D>,
    a: &Fq12Target<F, D>,
    power: usize,
) -> Fq12Target<F, D> {
    let neg_one: BigUint = Fq::from(-1).into();
    let modulus = neg_one + BigUint::from(1u64);
    assert_eq!(modulus.clone() % 4u64, BigUint::from(3u64));
    assert_eq!(modulus % 6u64, BigUint::from(1u64));
    let pow = power % 12;

    let mut out_fp2 = Vec::with_capacity(6);
    for i in 0..6 {
        let frob_coeff = frob_coeffs(pow).pow([i as u64]);
        let mut a_fp2 = Fq2Target {
            coeffs: [a.coeffs[i].clone(), a.coeffs[i + 6].clone()],
        };
        if pow % 2 != 0 {
            a_fp2 = a_fp2.conjugate(builder);
        }
        if frob_coeff == Fq2::one() {
            out_fp2.push(a_fp2);
        } else if frob_coeff.c1 == Fq::zero() {
            let frob_fixed = Fq2::new(frob_coeff.c0, Fq::zero());
            let frob_fixed_t = Fq2Target::constant(builder, frob_fixed);
            let out_nocarry = a_fp2.mul(builder, &frob_fixed_t);
            out_fp2.push(out_nocarry);
        } else {
            let frob_fixed = Fq2::new(frob_coeff.c0, frob_coeff.c1);
            let frob_fixed_t = Fq2Target::constant(builder, frob_fixed);
            let out_nocarry = a_fp2.mul(builder, &frob_fixed_t);
            out_fp2.push(out_nocarry);
        }
    }
    let out_coeffs = out_fp2
        .iter()
        .map(|x| x.coeffs[0].clone())
        .chain(out_fp2.iter().map(|x| x.coeffs[1].clone()))
        .collect::<Vec<_>>();

    Fq12Target {
        coeffs: out_coeffs.try_into().unwrap(),
    }
}

pub fn pow_target<F: RichField + Extendable<D>, const D: usize>(
    builder: &mut CircuitBuilder<F, D>,
    a: Fq12Target<F, D>,
    exp: Vec<u64>,
) -> Fq12Target<F, D> {
    let mut res = a.clone();
    let mut is_started = false;
    let naf = get_bls_naf(exp);

    for &z in naf.iter().rev() {
        if is_started {
            res = res.mul(builder, &res);
        }

        if z != 0 {
            assert!(z == 1 || z == -1);
            if is_started {
                res = res.mul(builder, &a);
            } else {
                assert_eq!(z, 1);
                is_started = true;
            }
        }
    }

    res
}

pub fn get_bls_naf(mut exp: Vec<u64>) -> Vec<i8> {
    // https://en.wikipedia.org/wiki/Non-adjacent_form
    // NAF for exp:
    let mut naf: Vec<i8> = Vec::with_capacity(64 * exp.len());
    let len = exp.len();

    // generate the NAF for exp
    for idx in 0..len {
        let mut e: u64 = exp[idx];
        for _ in 0..64 {
            if e & 1 == 1 {
                let z = 2i8 - (e % 4) as i8;
                // Excluded since our constant in NAF form doesn't contain negative ones
                // if z == -1 {
                //     e += 1;
                // }
                naf.push(z);
            } else {
                naf.push(0);
            }
            // Moving this outside the if and else statements since we are not checking if z == -1
            e /= 2;
        }
        if e != 0 {
            println!("enters e != 0");
            assert_eq!(e, 1);
            let mut j = idx + 1;
            while j < exp.len() && exp[j] == u64::MAX {
                exp[j] = 0;
                j += 1;
            }
            if j < exp.len() {
                exp[j] += 1;
            } else {
                exp.push(1);
            }
        }
    }
    if exp.len() != len {
        println!("enters exp.len() != len");
        assert_eq!(len, exp.len() + 1);
        assert!(exp[len] == 1);
        naf.push(1);
    }

    naf
}

pub fn frob_coeffs(index: usize) -> Fq2 {
    let neg_one: BigUint = Fq::from(-1).into();
    let modulus = neg_one + 1u64;

    let num: BigUint = modulus.pow(index as u32) - 1u64;
    let k: BigUint = num.div(6u64);

    let c = Fq2::new(Fq::from(1), Fq::one());
    c.pow(k.to_u64_digits())
}
