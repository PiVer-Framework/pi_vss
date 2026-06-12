use common::{
    error::{Error, ErrorKind::CountMismatch},
    polynomial::Polynomial,
    secret_sharing::generate_shares_batched,
    utils::batch_decompress_ristretto_points,
};
use rand::{CryptoRng, Rng};

use curve25519_dalek::{RistrettoPoint, Scalar, ristretto::CompressedRistretto, traits::Identity};

use rayon::prelude::*;

pub struct Dealer {
    pub t: usize,
    // [g1...gk]
    pub g: Vec<RistrettoPoint>,
    pub g0: RistrettoPoint,
    pub public_keys: Vec<RistrettoPoint>,
    pub(crate) secret: Option<Scalar>,
}

impl Dealer {
    pub fn new(
        g: Vec<RistrettoPoint>,
        g0: RistrettoPoint,
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
                g: g.clone(),
                g0: g0.clone(),
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
        x_pows: &Vec<Vec<Scalar>>,
        secrets: &Vec<Scalar>,
    ) -> (Vec<Vec<Scalar>>, (Vec<Scalar>, Vec<CompressedRistretto>))
    where
        R: CryptoRng + Rng,
    {
        let (f_polynomials, f_evals) =
            generate_shares_batched(self.public_keys.len(), self.t, x_pows, secrets);

        let mut c_buf: Vec<CompressedRistretto> = vec![CompressedRistretto::identity(); self.t + 1];

        let r_evals = self.generate_proof(rng, &mut c_buf, x_pows, &f_polynomials);
        (f_evals, (r_evals, c_buf))
    }

    pub fn generate_proof<R>(
        &self,
        rng: &mut R,
        c_buf: &mut Vec<CompressedRistretto>,
        x_pows: &Vec<Vec<Scalar>>,
        f_polynomials: &Vec<Polynomial>,
    ) -> Vec<Scalar>
    where
        R: CryptoRng,
    {
        let r = Polynomial::sample(self.t, rng);
        let r_evals = r.evaluate_range_precomp(x_pows, 1, self.public_keys.len());

        c_buf
            .par_iter_mut()
            .zip(r.coef_ref().par_iter())
            .enumerate()
            .for_each(|(t, (c, r_coef))| {
                *c = f_polynomials
                    .par_iter()
                    .zip(self.g.par_iter())
                    .map(|(fk, gk)| gk * fk.coef_at_unchecked(t))
                    .reduce(|| self.g0 * r_coef, |acc, prod| acc + prod)
                    .compress()
            });

        r_evals
    }

    pub fn get_pk0(&self) -> &RistrettoPoint {
        &self.public_keys[0]
    }

    pub fn publish_f0(&self) -> Scalar {
        self.secret.unwrap()
    }
}
