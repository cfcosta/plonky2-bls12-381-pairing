use ark_bls12_381::G1Affine;
use ark_ec::AffineRepr;
use plonky2::{
    field::extension::Extendable,
    hash::hash_types::RichField,
    plonk::circuit_builder::CircuitBuilder,
};

use crate::fields::fq_target::FqTarget;

#[derive(Clone, Debug)]
pub struct G1AffineTarget<F: RichField + Extendable<D>, const D: usize> {
    pub x: FqTarget<F, D>,
    pub y: FqTarget<F, D>,
    pub infinity: bool,
}

impl<F: RichField + Extendable<D>, const D: usize> G1AffineTarget<F, D> {
    pub fn is_zero(&self) -> bool {
        self.infinity
    }

    pub fn xy(&self) -> Option<(&self::FqTarget<F, D>, &self::FqTarget<F, D>)> {
        (!self.infinity).then_some((&self.x, &self.y))
    }

    pub fn constant(builder: &mut CircuitBuilder<F, D>, g1: G1Affine) -> Self {
        Self {
            x: FqTarget::constant(builder, g1.x().unwrap()),
            y: FqTarget::constant(builder, g1.y().unwrap()),
            infinity: false,
        }
    }

    pub fn connect(builder: &mut CircuitBuilder<F, D>, lhs: &Self, rhs: &Self) {
        // The only cases in which two points are equal are
        // 1. infinity is set on both
        // 2. infinity is not set on both, and their coordinates are equal
        assert!((lhs.infinity & rhs.infinity) | (!lhs.infinity) & (!rhs.infinity));
        FqTarget::connect(builder, &lhs.x, &rhs.x);
        FqTarget::connect(builder, &lhs.y, &rhs.y);
    }
}

#[derive(Clone, Debug)]
pub struct G1PreparedTarget<F: RichField + Extendable<D>, const D: usize>(pub G1AffineTarget<F, D>);
