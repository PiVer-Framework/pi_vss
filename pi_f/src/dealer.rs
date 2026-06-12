use common::{
    error::{Error, ErrorKind::CountMismatch},
    polynomial::Polynomial,
    secret_sharing::generate_shares,
    utils::{batch_decompress_ristretto_points, compute_d_from_point_commitments},
};
use rand::{CryptoRng, Rng};
use rayon::prelude::*;

use blake3::Hasher;
use curve25519_dalek::{RistrettoPoint, Scalar, ristretto::CompressedRistretto};

pub struct Dealer {
    pub t: usize,
    pub g1: RistrettoPoint,
    pub g2: RistrettoPoint,
    pub public_keys: Vec<RistrettoPoint>,
    pub(crate) secret: Option<Scalar>,
}

impl Dealer {
    pub fn new(
        g1: RistrettoPoint,
        g2: RistrettoPoint,
        n: usize,
        t: usize,
        public_keys: &[CompressedRistretto],
    ) -> Result<Self, Error> {
        if public_keys.len() != n {
            return Err(CountMismatch(n, "parties", public_keys.len(), "public keys").into());
        }
        match batch_decompress_ristretto_points(public_keys) {
            Ok(pks) => Ok(Self {
                t,
                public_keys: pks,
                secret: None,
                g1: g1.clone(),
                g2: g2.clone(),
            }),
            Err(x) => Err(x),
        }
    }

    pub fn t(&self) -> usize {
        self.t
    }

    pub fn deal_secret<R>(
        &mut self,
        rng: &mut R,
        hasher: &mut Hasher,
        buf: &mut [u8; 64],
        x_pows: &Vec<Vec<Scalar>>,
        secret: &Scalar,
    ) -> (Vec<Scalar>, (Vec<CompressedRistretto>, Polynomial))
    where
        R: CryptoRng + Rng,
    {
        let (f_polynomial, f_evals) =
            generate_shares(rng, self.public_keys.len(), self.t, x_pows, secret);

        let mut c_buf: Vec<CompressedRistretto> = Vec::with_capacity(self.public_keys.len());

        let z = self.generate_proof(rng, hasher, buf, &mut c_buf, x_pows, f_polynomial, &f_evals);

        (f_evals, (c_buf, z))
    }

    pub fn generate_proof<R>(
        &self,
        rng: &mut R,
        hasher: &mut Hasher,
        buf: &mut [u8; 64],
        c_buf: &mut Vec<CompressedRistretto>,
        x_pows: &Vec<Vec<Scalar>>,
        f_polynomial: Polynomial,
        f_evals: &Vec<Scalar>,
    ) -> Polynomial
    where
        R: CryptoRng,
    {
        let mut r = Polynomial::sample(self.t, rng);
        let r_evals = r.evaluate_range_precomp(x_pows, 1, self.public_keys.len());

        f_evals
            .par_iter()
            .zip(r_evals.par_iter())
            .map(|(fi, ri)| (self.g1 * fi + self.g2 * ri).compress())
            .collect_into_vec(c_buf);

        let d = compute_d_from_point_commitments(hasher, buf, &c_buf);

        // z == r +=  d * f
        r.compute_z(&[f_polynomial], &[d]);

        r
    }

    pub fn get_pk0(&self) -> &RistrettoPoint {
        &self.public_keys[0]
    }

    pub fn publish_f0(&self) -> Scalar {
        self.secret.unwrap()
    }
}
