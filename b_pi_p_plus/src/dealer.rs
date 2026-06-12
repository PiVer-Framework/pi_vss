use common::{
    error::{Error, ErrorKind::CountMismatch},
    polynomial::Polynomial,
    random::random_scalars,
    secret_sharing::generate_shares_batched,
    utils::{batch_decompress_ristretto_points, compute_d_powers_from_point_commitments},
};
use rand::{CryptoRng, Rng};
use rayon::prelude::*;

use blake3::Hasher;
use curve25519_dalek::{RistrettoPoint, Scalar, ristretto::CompressedRistretto};
use zeroize::Zeroize;

pub struct Dealer {
    pub t: usize,
    // [g1...gk]
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
        secrets: &Vec<Scalar>,
    ) -> (
        Vec<Vec<Scalar>>,
        (Vec<Scalar>, Vec<CompressedRistretto>, Polynomial),
    )
    where
        R: CryptoRng + Rng,
    {
        let k = secrets.len();

        let (f_polynomials, f_evals) =
            generate_shares_batched(self.public_keys.len(), self.t, x_pows, secrets);

        let mut c_buf: Vec<CompressedRistretto> = Vec::with_capacity(self.public_keys.len());

        let (g, z) = self.generate_proof(
            rng,
            hasher,
            buf,
            &mut c_buf,
            x_pows,
            k,
            &f_polynomials,
            &f_evals,
        );
        (f_evals, (g, c_buf, z))
    }

    pub fn generate_proof<R>(
        &self,
        rng: &mut R,
        hasher: &mut Hasher,
        buf: &mut [u8; 64],
        c_buf: &mut Vec<CompressedRistretto>,
        x_pows: &Vec<Vec<Scalar>>,
        k: usize,
        f_polynomials: &Vec<Polynomial>,
        f_evals: &Vec<Vec<Scalar>>,
    ) -> (Vec<Scalar>, Polynomial)
    where
        R: CryptoRng,
    {
        let mut r = Polynomial::sample(self.t, rng);
        let r_evals = r.evaluate_range_precomp(x_pows, 1, self.public_keys.len());

        let g: Vec<Scalar> = random_scalars(rng, self.public_keys.len());

        f_evals
            .par_iter()
            .zip(r_evals.par_iter())
            .zip(g.par_iter())
            .map_init(
                || (Hasher::new(), [0u8; 64]),
                |(l_hasher, l_buf), ((fi, ri), gi)| {
                    fi.iter().for_each(|fi_k| {
                        l_hasher.update(fi_k.as_bytes());
                    });
                    l_hasher.update(ri.as_bytes());

                    l_hasher.finalize_xof().fill(l_buf);
                    l_hasher.reset();

                    let h = Scalar::from_bytes_mod_order_wide(l_buf);
                    l_buf.zeroize();

                    ((self.g1 * h) + (self.g2 * gi)).compress()
                },
            )
            .collect_into_vec(c_buf);

        // [d, d^2, ..., d^k]
        let d_vals = compute_d_powers_from_point_commitments(hasher, buf, &c_buf, k);

        // z == r +=  d * f
        // if self.g1 == self.g2 * d {
        //     panic!("g1 == g2^d");
        // } else {
        r.compute_z(f_polynomials, &d_vals);
        (g, r)
        // }
    }

    pub fn get_pk0(&self) -> &RistrettoPoint {
        &self.public_keys[0]
    }

    pub fn publish_f0(&self) -> Scalar {
        self.secret.unwrap()
    }
}
