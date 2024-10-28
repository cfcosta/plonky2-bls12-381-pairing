use ark_bls12_381::{Fq, Fq2};
use ark_ff::{Field, Fp2};
use itertools::Itertools;
use num::BigUint;
use num::{One, Zero};
use std::ops::Div;

use crate::utils::helpers::MyFq12;

pub fn frobenius_map_native(a: MyFq12, power: usize) -> MyFq12 {
    let neg_one: BigUint = Fq::from(-1).into();
    let modulus = neg_one + BigUint::from(1u64);
    assert_eq!(modulus.clone() % 4u64, BigUint::from(3u64));
    assert_eq!(modulus % 6u64, BigUint::from(1u64));
    let pow = power % 12;

    let mut out_fp2 = Vec::with_capacity(6);

    for i in 0..6 {
        let frob_coeff = frob_coeffs(pow).pow([i as u64]);
        let mut a_fp2 = Fq2::new(a.coeffs[i], a.coeffs[i + 6]);
        if pow % 2 != 0 {
            a_fp2 = conjugate_fp2(a_fp2);
        }

        if frob_coeff == Fq2::one() {
            out_fp2.push(a_fp2);
        } else if frob_coeff.c1 == Fq::zero() {
            let frob_fixed = Fq2::new(frob_coeff.c0, Fq::zero());
            let out_nocarry = a_fp2 * frob_fixed;
            out_fp2.push(out_nocarry);
        } else {
            let frob_fixed = Fq2::new(frob_coeff.c0, frob_coeff.c1);
            out_fp2.push(a_fp2 * frob_fixed);
        }
    }

    let out_coeffs = out_fp2
        .iter()
        .map(|x| x.c0)
        .chain(out_fp2.iter().map(|x| x.c1))
        .collect_vec();

    MyFq12 {
        coeffs: out_coeffs.try_into().unwrap(),
    }
}

pub fn frob_coeffs(index: usize) -> Fq2 {
    let neg_one: BigUint = Fq::from(-1).into();
    let modulus = neg_one + 1u64;

    let num: BigUint = modulus.pow(index as u32) - 1u64;
    let k: BigUint = num.div(6u64);

    let c = Fq2::new(Fq::from(1), Fq::one());
    c.pow(k.to_u64_digits())
}

pub fn pow_native(a: MyFq12, exp: Vec<u64>) -> MyFq12 {
    let mut res = a;
    let mut is_started = false;
    let naf = get_naf(exp);

    for &z in naf.iter().rev() {
        if is_started {
            res = res * res;
        }

        if z != 0 {
            assert!(z == 1 || z == -1);
            if is_started {
                res = res * a;
            } else {
                assert_eq!(z, 1);
                is_started = true;
            }
        }
    }

    res
}

pub fn get_naf(mut exp: Vec<u64>) -> Vec<i8> {
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

pub fn conjugate_fp12(a: MyFq12) -> MyFq12 {
    let coeffs: Vec<Fq> = a
        .coeffs
        .iter()
        .enumerate()
        .map(|(i, c)| if i % 2 == 0 { *c } else { -*c })
        .collect();
    MyFq12 {
        coeffs: coeffs.try_into().unwrap(),
    }
}

pub fn conjugate_fp2(x: Fq2) -> Fq2 {
    Fp2 {
        c0: x.c0,
        c1: -x.c1,
    }
}
