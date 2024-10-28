use crate::utils::{constants::BLS_X, helpers::MyFq12};
use ark_bls12_381::Fq12;
use ark_ff::{CyclotomicMultSubgroup, Field};

use super::fin_exp_helpers::{conjugate_fp12, frobenius_map_native, pow_native};

// out = in^{ (q^6 - 1)*(q^2 + 1) }
pub fn easy_part_native<'v>(a: MyFq12) -> MyFq12 {
    let f1 = conjugate_fp12(a);
    let f2 = {
        let f1_fp12: Fq12 = f1.into();
        let a_fp12: Fq12 = a.into();
        let divided = f1_fp12 / a_fp12;
        divided.into()
    };
    let f3 = frobenius_map_native(f2, 2);
    
    f3 * f2
}

pub fn ark_easy_part(f: Fq12) -> Option<Fq12> {
    // f1 = r.cyclotomic_inverse_in_place() = f^(p^6)
    let f = f;
    let mut f1 = f;
    f1.cyclotomic_inverse_in_place();

    f.inverse().map(|mut f2| {
        // f2 = f^(-1);
        // r = f^(p^6 - 1)
        let mut r = f1 * &f2;

        // f2 = f^(p^6 - 1)
        f2 = r;
        // r = f^((p^6 - 1)(p^2))
        r.frobenius_map_in_place(2);

        // r = f^((p^6 - 1)(p^2) + (p^6 - 1))
        // r = f^((p^6 - 1)(p^2 + 1))
        r *= &f2;
        r
    })
}

pub fn hard_part_native(r: Fq12) -> Fq12 {
    let mut y0 = r.square();
    let mut y1: Fq12 = pow_native(r.into(), vec![BLS_X]).into();
    y1 = *y1.conjugate_in_place();
    let mut y2 = r;
    y2 = *y2.conjugate_in_place();

    y1 *= y2;
    y2 = pow_native(y1.into(), vec![BLS_X]).into();
    y2 = *y2.conjugate_in_place();
    y1 = *y1.conjugate_in_place();
    y1 *= y2;
    y2 = pow_native(y1.into(), vec![BLS_X]).into();
    y2 = *y2.conjugate_in_place();
    y1 = frobenius_map_native(y1.into(), 1).into();
    y1 *= y2;
    let r = r * y0;
    y0 = pow_native(y1.into(), vec![BLS_X]).into();
    y0 = *y0.conjugate_in_place();
    y2 = pow_native(y0.into(), vec![BLS_X]).into();
    y2 = *y2.conjugate_in_place();
    y0 = y1;
    y0 = frobenius_map_native(y0.into(), 2).into();
    y1 = *y1.conjugate_in_place();
    y1 *= y2;
    y1 *= y0;
    
    r * y1
}

pub fn final_exponentiation(a: MyFq12) -> Fq12 {
    let f0 = easy_part_native(a);
    
    hard_part_native(f0.into())
}

#[cfg(test)]
mod tests {
    use ark_bls12_381::{Fq12, Fr};
    use ark_ec::pairing::{MillerLoopOutput, Pairing};
    use ark_ff::{Field, UniformRand};

    use crate::native::final_exponentiation::{
        ark_easy_part, final_exponentiation, hard_part_native,
    };

    #[test]
    fn test_native_final_exponentiation() {
        let rng = &mut rand::thread_rng();
        let x = Fq12::rand(rng);
        let ark_final_exp_result =
            ark_bls12_381::Bls12_381::final_exponentiation(MillerLoopOutput(x))
                .unwrap()
                .0;
        let final_exp_native_result = final_exponentiation(x.into());

        assert_eq!(ark_final_exp_result, final_exp_native_result);
    }

    #[test]
    fn test_ark_final_exponentiation() {
        use ark_bls12_381::{Fq, Fq12, Fr};
        use ark_ff::{Field, UniformRand};
        use num::BigUint;
        let rng = &mut rand::thread_rng();
        let rand_x = Fq12::rand(rng);
        let ark_final_exponentiation =
            ark_bls12_381::Bls12_381::final_exponentiation(MillerLoopOutput(rand_x))
                .unwrap()
                .0;

        use ark_ff::PrimeField;
        let p: BigUint = Fq::MODULUS.into();
        let r: BigUint = Fr::MODULUS.into();
        let exponent = (p.pow(12) - 1u32) / r;
        let fixed_large_exponent = rand_x.pow(&exponent.to_u64_digits());

        assert_eq!(ark_final_exponentiation, fixed_large_exponent);
    }

    #[test]
    fn test_ark_easy_part() {
        use ark_bls12_381::{Fq, Fq12};
        use ark_ff::{Field, UniformRand};
        use num::BigUint;
        let rng = &mut rand::thread_rng();
        let rand_x = Fq12::rand(rng);
        let ark_easy_part = ark_easy_part(rand_x).unwrap();

        // f_easy = x^((p^6 - 1)*(p^2 + 1))
        use ark_ff::PrimeField;
        let p: BigUint = Fq::MODULUS.into();
        let exponent = (p.pow(6) - 1u32) * (p.pow(2) + 1u32);
        let fixed_large_exponent = rand_x.pow(&exponent.to_u64_digits());

        // Passes
        assert_eq!(ark_easy_part, fixed_large_exponent);
    }

    #[test]
    fn test_ark_hard_part() {
        use ark_bls12_381::{Fq, Fq12};
        use ark_ff::UniformRand;
        use num::BigUint;
        let rng = &mut rand::thread_rng();
        let rand_x = Fq12::rand(rng);
        let ark_hard_part = hard_part_native(rand_x);

        // (p^12 - 1) / r = f_easy * ((p^4 - p^2 + 1) / r)
        use ark_ff::PrimeField;
        let p: BigUint = Fq::MODULUS.into();
        let r: BigUint = Fr::MODULUS.into();

        let fixed_large_exponent = (p.pow(12) - 1u32) / r.clone();
        let easy_part = (p.pow(6) - 1u32) * (p.pow(2) + 1u32);
        let hard_part = (p.pow(4) - p.pow(2) + 1u64) / r;
        let expected_hard_part = fixed_large_exponent.clone() / easy_part.clone();
        let expected_final_exponentiation = hard_part.clone() * easy_part;

        // Passes
        assert_eq!(expected_final_exponentiation, fixed_large_exponent);

        // Passes
        assert_eq!(expected_hard_part, hard_part);

        // Fails
        let expected_result = rand_x.pow(&hard_part.to_u64_digits());
        assert_eq!(ark_hard_part, expected_result);
    }
}
