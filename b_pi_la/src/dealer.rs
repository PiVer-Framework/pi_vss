use common::{
    error::{Error, ErrorKind::CountMismatch},
    polynomial::Polynomial,
    secret_sharing::generate_shares_batched,
    utils::{batch_decompress_ristretto_points, compute_d_powers_from_hash_commitments},
};
use rand::{CryptoRng, Rng};

use blake3::Hasher;
use curve25519_dalek::{RistrettoPoint, Scalar, ristretto::CompressedRistretto};
use rayon::prelude::*;

pub struct Dealer {
    pub t: usize,
    pub public_keys: Vec<RistrettoPoint>,
    pub(crate) secrets: Option<Vec<Scalar>>,
}

impl Dealer {
    pub fn new(n: usize, t: usize, public_keys: &[CompressedRistretto]) -> Result<Self, Error> {
        if public_keys.len() != n {
            return Err(CountMismatch(n, "parties", public_keys.len(), "public keys").into());
        }
        match batch_decompress_ristretto_points(public_keys) {
            Ok(pks) => Ok(Self {
                t,
                public_keys: pks,
                secrets: None,
            }),
            Err(x) => Err(x),
        }
    }

    pub fn t(&self) -> usize {
        self.t
    }

    pub fn get_pk0(&self) -> &RistrettoPoint {
        &self.public_keys[0]
    }

    pub fn publish_f0(&self) -> Vec<Scalar> {
        self.secrets.clone().unwrap()
    }

    pub fn deal_secrets<R>(
        &mut self,
        rng: &mut R,
        hasher: &mut Hasher,
        buf: &mut [u8; 64],
        x_pows: &Vec<Vec<Scalar>>,
        secrets: &Vec<Scalar>,
    ) -> (Vec<Vec<Scalar>>, (Vec<[u8; 64]>, Polynomial))
    where
        R: CryptoRng + Rng,
    {
        // number of secrets to share
        let k = secrets.len();

        let (mut f_polynomials, f_evals) =
            generate_shares_batched(self.public_keys.len(), self.t, x_pows, secrets);

        let mut c_buf = vec![[0u8; 64]; self.public_keys.len()];

        let z = self.generate_proof(
            rng,
            hasher,
            buf,
            &mut c_buf,
            x_pows,
            k,
            &mut f_polynomials,
            &f_evals,
        );

        (f_evals, (c_buf, z))
    }

    pub fn generate_proof<R>(
        &self,
        rng: &mut R,
        hasher: &mut Hasher,
        buf: &mut [u8; 64],
        c_buf: &mut Vec<[u8; 64]>,
        x_pows: &Vec<Vec<Scalar>>,
        k: usize,
        f_polynomials: &Vec<Polynomial>,
        f_evals: &Vec<Vec<Scalar>>,
    ) -> Polynomial
    where
        R: CryptoRng,
    {
        let mut r = Polynomial::sample(self.t, rng);

        let r_evals = r.evaluate_range_precomp(x_pows, 1, self.public_keys.len());

        c_buf
            .par_iter_mut()
            .zip(f_evals.par_iter().zip(r_evals.par_iter()))
            .for_each_init(
                || Hasher::new(),
                |l_hasher, (l_buf, (fi, ri))| {
                    fi.iter().for_each(|fi_k| {
                        l_hasher.update(fi_k.as_bytes());
                    });

                    l_hasher.update(ri.as_bytes());

                    l_hasher.finalize_xof().fill(l_buf);
                    l_hasher.reset();
                },
            );

        // [d, d^2, ..., d^k]
        let d_vals = compute_d_powers_from_hash_commitments(hasher, buf, &c_buf, k);

        // z == r += ( ∑ d_j * f_j )
        r.compute_z(f_polynomials, &d_vals);

        r
    }
}
