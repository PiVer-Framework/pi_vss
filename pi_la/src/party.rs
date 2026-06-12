use blake3::Hasher;
use curve25519_dalek::{RistrettoPoint, Scalar, ristretto::CompressedRistretto};
use rand::{CryptoRng, Rng};
use zeroize::Zeroize;

use common::{
    error::{
        Error,
        ErrorKind::{CountMismatch, InvalidPararmeterSet, InvalidProof, UninitializedValue},
    },
    polynomial::Polynomial,
    random::random_scalar,
    utils::compute_d_from_hash_commitments,
};
use rayon::prelude::*;

#[derive(Clone)]
pub struct Party {
    pub private_key: Scalar,
    pub public_key: (CompressedRistretto, RistrettoPoint),
    pub index: usize,
    pub n: usize,
    pub t: usize,
    pub public_keys: Option<Vec<RistrettoPoint>>,
    pub dealer_proof: Option<(Vec<[u8; 64]>, Polynomial)>,
    pub validated_shares: Vec<usize>,
    pub share: Option<Scalar>,
    pub d: Option<Scalar>,
    pub shares: Option<Vec<Scalar>>,
    pub qualified_set: Option<Vec<(usize, Scalar)>>,
}

impl Party {
    pub fn new<R>(
        g: &RistrettoPoint,
        rng: &mut R,
        n: usize,
        t: usize,
        index: usize,
    ) -> Result<Self, Error>
    where
        R: CryptoRng + Rng,
    {
        let private_key = random_scalar(rng);
        let public_key = g * &private_key;

        if index <= n && t < n && t as f32 == ((n - 1) as f32 / 2.0).floor() {
            Ok(Self {
                private_key,
                public_key: (public_key.compress(), public_key),
                index,
                n,
                t,
                dealer_proof: None,
                share: None,
                public_keys: None,
                validated_shares: vec![],
                d: None,
                shares: None,
                qualified_set: None,
            })
        } else {
            Err(InvalidPararmeterSet(n, t as isize, index).into())
        }
    }

    pub fn ingest_share(&mut self, share: &Scalar) {
        self.share = Some(share.clone());
    }

    pub fn ingest_dealer_proof(
        &mut self,
        proof: (&Vec<[u8; 64]>, &Polynomial),
    ) -> Result<(), Error> {
        if proof.1.len() != self.t + 1 {
            Err(InvalidProof(format!("z len: {}, t: {}", proof.1.len(), self.t + 1)).into())
        } else if proof.0.len() != self.n {
            Err(InvalidProof(format!("c_vals len: {}, n: {}", proof.0.len(), self.n)).into())
        } else {
            self.dealer_proof = Some((proof.0.clone(), proof.1.clone()));
            Ok(())
        }
    }

    pub fn verify_share(
        &self,
        hasher: &mut Hasher,
        buf: &mut [u8; 64],
        x_pows: &Vec<Vec<Scalar>>,
    ) -> Result<bool, Error> {
        match &self.dealer_proof {
            Some((c_vals, z)) => match &self.share {
                Some(fi) => {
                    let d = compute_d_from_hash_commitments(hasher, buf, c_vals);

                    hasher.update(fi.as_bytes());
                    hasher.update(
                        Polynomial::compute_r_eval(
                            &z.evaluate_precomp(x_pows, self.index),
                            &[*fi],
                            &[d],
                        )
                        .as_bytes(),
                    );

                    hasher.finalize_xof().fill(buf);
                    hasher.reset();

                    let check_bit = &c_vals[self.index - 1] == buf;
                    buf.zeroize();

                    Ok(check_bit)
                }
                None => Err(UninitializedValue("party.share").into()),
            },
            None => Err(UninitializedValue("party.dealer_proof").into()),
        }
    }

    pub fn verify_shares(
        &mut self,
        hasher: &mut Hasher,
        buf: &mut [u8; 64],
        x_pows: &Vec<Vec<Scalar>>,
    ) -> Result<bool, Error> {
        match &self.dealer_proof {
            Some((cvals, z)) => match &self.shares {
                Some(shares) => {
                    let d = compute_d_from_hash_commitments(hasher, buf, cvals);
                    let z_evals = z.evaluate_range_precomp(x_pows, 1, self.n);

                    self.validated_shares = shares
                        .par_iter()
                        .zip(z_evals.par_iter())
                        .enumerate()
                        .map_init(
                            || (Hasher::new(), [0u8; 64]),
                            |(l_hasher, l_buf), (i, (fi, zi))| {
                                l_hasher.update(fi.as_bytes());
                                l_hasher.update(
                                    (Polynomial::compute_r_eval(zi, &[*fi], &[d])).as_bytes(),
                                );

                                l_hasher.finalize_xof().fill(l_buf);
                                l_hasher.reset();

                                let check_bit = cvals[i] == *l_buf;
                                l_buf.zeroize();
                                if check_bit { Some(i) } else { None }
                            },
                        )
                        .filter(Option::is_some)
                        .map(|res| res.unwrap())
                        .collect();
                    Ok(self.validated_shares.len() > self.t)
                }
                None => Err(UninitializedValue("party.share").into()),
            },
            None => Err(UninitializedValue("party.dealer_proof").into()),
        }
    }

    pub fn ingest_shares(&mut self, shares: &Vec<Scalar>) -> Result<(), Error> {
        if shares.len() == self.n {
            self.shares = Some(shares.clone());
            Ok(())
        } else {
            Err(CountMismatch(self.n, "parties", shares.len(), "ingestable shares").into())
        }
    }
}

pub fn generate_parties<R>(g: &RistrettoPoint, rng: &mut R, n: usize, t: usize) -> Vec<Party>
where
    R: CryptoRng + Rng,
{
    (1..=n)
        .map(|i| Party::new(g, rng, n, t, i).unwrap())
        .collect()
}
