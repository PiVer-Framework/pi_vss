use common::{
    error::{Error, ErrorKind::CountMismatch},
    polynomial::Polynomial,
    secret_sharing::generate_shares,
};

use blake3::Hasher;
use common::utils::batch_decompress_ristretto_points;
use curve25519_dalek::{RistrettoPoint, Scalar, ristretto::CompressedRistretto};

use rand::{CryptoRng, Rng};
use rayon::prelude::*;
use zeroize::Zeroize;

pub struct Dealer {
    t: usize,
    public_keys: Vec<RistrettoPoint>,
    pub(crate) secret: Option<Scalar>,
}

impl Dealer {
    pub fn new(n: usize, t: usize, public_keys: &[CompressedRistretto]) -> Result<Self, Error> {
        if public_keys.len() != n {
            return Err(CountMismatch(n, "parties", public_keys.len(), "public keys").into());
        }
        match batch_decompress_ristretto_points(public_keys) {
            Ok(pks) => Ok(Self {
                t,
                public_keys: pks.par_iter().map(|pk| *pk).collect(),
                secret: None,
            }),
            Err(x) => Err(x),
        }
    }

    pub fn deal_secret<R>(
        &mut self,
        rng: &mut R,
        hasher: &mut Hasher,
        buf: &mut [u8; 64],
        x_pows: &Vec<Vec<Scalar>>,
        secret: &Scalar,
    ) -> (Vec<CompressedRistretto>, (Scalar, Polynomial))
    where
        R: CryptoRng + Rng,
    {
        self.secret = Some(*secret);
        let (f_polynomial, f_evals) = self.generate_encrypted_shares(rng, x_pows, secret);

        let (d, z) = self.generate_proof(rng, hasher, buf, x_pows, f_polynomial, &f_evals);

        (f_evals, (d, z))
    }

    pub fn generate_encrypted_shares<R>(
        &self,
        rng: &mut R,
        x_pows: &Vec<Vec<Scalar>>,
        secret: &Scalar,
    ) -> (Polynomial, Vec<CompressedRistretto>)
    where
        R: CryptoRng,
    {
        let (f_polynomial, f_evals) =
            generate_shares(rng, self.public_keys.len(), self.t, x_pows, secret);

        let encrypted_shares = f_evals
            .par_iter()
            .zip(self.public_keys.par_iter())
            .map(|(f_eval, pub_key)| (f_eval * pub_key).compress())
            .collect();

        (f_polynomial, encrypted_shares)
    }

    pub fn generate_proof<R>(
        &self,
        rng: &mut R,
        hasher: &mut Hasher,
        buf: &mut [u8; 64],
        x_pows: &Vec<Vec<Scalar>>,
        f_polynomial: Polynomial,
        f_evals: &Vec<CompressedRistretto>,
    ) -> (Scalar, Polynomial)
    where
        R: CryptoRng,
    {
        let mut r = Polynomial::sample(self.t, rng);
        let r_evals = r.evaluate_range_precomp(x_pows, 1, self.public_keys.len());

        let encrypted_r_evals: Vec<CompressedRistretto> = r_evals
            .par_iter()
            .zip(self.public_keys.par_iter())
            .map(|(f_eval, pub_key)| (f_eval * pub_key).compress())
            .collect();

        f_evals
            .iter()
            .chain(encrypted_r_evals.iter())
            .for_each(|x| {
                hasher.update(x.as_bytes());
            });

        hasher.finalize_xof().fill(buf);

        let d = Scalar::from_bytes_mod_order_wide(buf);

        hasher.reset();
        buf.zeroize();

        r.compute_z(&[f_polynomial], &[d]);

        (d, r)
    }
}
