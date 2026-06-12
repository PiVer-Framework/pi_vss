use std::ops::Mul;

use blake3::Hasher;
use curve25519_dalek::{RistrettoPoint, Scalar, ristretto::CompressedRistretto};
use num_bigint::BigUint;
use rayon::prelude::*;
use zeroize::Zeroize;

use crate::{
    Q,
    error::{
        Error,
        ErrorKind::{CountMismatch, PointDecompressionError},
    },
};

pub fn mod_pow(i: usize, t: usize) -> Scalar {
    let i_m = BigUint::from(i);
    let t_m = BigUint::from(t);
    let modulus = BigUint::from_bytes_le(&Q);

    let bytes = i_m.modpow(&t_m, &modulus).to_bytes_le();
    let mut tmp_bytes: [u8; 32] = [0u8; 32];

    for i in 0..bytes.len() {
        tmp_bytes[i] = bytes[i];
    }

    Scalar::from_bytes_mod_order(tmp_bytes)
}

pub fn pointwise_op_in_place(
    op: fn(Scalar, Scalar) -> Scalar,
    a: &mut Vec<Scalar>,
    b: &Vec<Scalar>,
) {
    a.par_iter_mut().zip(b.par_iter()).for_each(|(a_i, b_i)| {
        *a_i = op(*a_i, *b_i);
    });
}

pub fn compute_d_from_hash_commitments(
    hasher: &mut Hasher,
    buf: &mut [u8; 64],
    commitments: &[[u8; 64]],
) -> Scalar {
    commitments.iter().for_each(|c| {
        hasher.update(c);
    });

    hasher.finalize_xof().fill(buf);
    hasher.reset();

    let d = Scalar::from_bytes_mod_order_wide(buf);
    buf.zeroize();
    d
}

pub fn compute_d_from_point_commitments(
    hasher: &mut Hasher,
    buf: &mut [u8; 64],
    commitments: &[CompressedRistretto],
) -> Scalar {
    commitments.iter().for_each(|c| {
        hasher.update(c.as_bytes());
    });

    hasher.finalize_xof().fill(buf);
    hasher.reset();

    let d = Scalar::from_bytes_mod_order_wide(buf);
    buf.zeroize();
    d
}

pub fn compute_d_powers_from_hash_commitments(
    hasher: &mut Hasher,
    buf: &mut [u8; 64],
    commitments: &[[u8; 64]],
    k: usize,
) -> Vec<Scalar> {
    let d = compute_d_from_hash_commitments(hasher, buf, commitments);

    compute_d_powers(k, &d)
}

pub fn compute_d_powers_from_point_commitments(
    hasher: &mut Hasher,
    buf: &mut [u8; 64],
    commitments: &[CompressedRistretto],
    k: usize,
) -> Vec<Scalar> {
    let d = compute_d_from_point_commitments(hasher, buf, commitments);

    compute_d_powers(k, &d)
}

pub fn compute_d_powers(k: usize, d: &Scalar) -> Vec<Scalar> {
    let mut d_vals: Vec<Scalar> = Vec::with_capacity(k);
    // [d^1,
    d_vals.push(*d);

    // d^2, d^3, ... d^k]
    for i in 1..k {
        d_vals.push(d_vals[i - 1] * d);
    }
    //
    d_vals
}

pub fn precompute_lambda(n: usize, t: usize) -> Vec<Scalar> {
    (1..=n)
        .into_par_iter()
        .map(|i| {
            let zq_i = Scalar::from(i as u64);
            let mut lambda_i = Scalar::ONE;
            for j in 1..=(t + 1) {
                if j != i {
                    let zq_j = Scalar::from(j as u64);

                    lambda_i *= zq_j * ((zq_j - zq_i).invert());
                }
            }
            lambda_i
        })
        .collect()
}

pub fn compute_lagrange_bases(qualified_set: &Vec<usize>) -> Vec<Scalar> {
    qualified_set
        .par_iter()
        .map(|i| compute_lagrange_basis(*i, qualified_set))
        .collect()
}

pub fn compute_lagrange_basis(i: usize, qualified_set: &Vec<usize>) -> Scalar {
    let zq_i = Scalar::from(i as u64);

    qualified_set
        .par_iter()
        .map(|j| {
            if i == *j {
                Scalar::ONE
            } else {
                let zq_j = Scalar::from(*j as u64);
                zq_j * ((zq_j - zq_i).invert())
            }
        })
        .reduce(|| Scalar::ONE, Scalar::mul)
}

pub fn decompress_ristretto_point(
    compressed_point: CompressedRistretto,
) -> Result<RistrettoPoint, Error> {
    match compressed_point.decompress() {
        Some(decompressed_point) => Ok(decompressed_point),
        None => Err(Error::from_kind(PointDecompressionError(format!(
            "{compressed_point:?}",
        )))),
    }
}

pub fn batch_decompress_ristretto_points(
    compressed_points: &[CompressedRistretto],
) -> Result<Vec<RistrettoPoint>, Error> {
    compressed_points
        .par_iter()
        .map(|compressed_point| decompress_ristretto_point(*compressed_point))
        .collect()
}
pub fn batch_decompress_batched_ristretto_points(
    // vec[vec[_;k]; n]
    batch_compressed_points: &Vec<Vec<CompressedRistretto>>,
) -> Result<Vec<Vec<RistrettoPoint>>, Error> {
    batch_compressed_points
        .par_iter()
        .map(|compressed_points| batch_decompress_ristretto_points(compressed_points))
        .collect()
}

pub fn ingest_public_keys(
    n: usize,
    own_public_key: &RistrettoPoint,
    own_index: usize,
    public_keys: &[CompressedRistretto],
) -> Result<Vec<RistrettoPoint>, Error> {
    if public_keys.len() == n - 1 {
        match batch_decompress_ristretto_points(public_keys) {
            Ok(mut pks) => {
                pks.insert(own_index - 1, *own_public_key);
                Ok(pks)
            }
            Err(x) => Err(x),
        }
    } else {
        Err(CountMismatch(n, "parties", public_keys.len(), "public_keys").into())
    }
}
