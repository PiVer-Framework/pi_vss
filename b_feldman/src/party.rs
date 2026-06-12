use curve25519_dalek::{RistrettoPoint, Scalar, ristretto::CompressedRistretto, traits::Identity};
use rand::{CryptoRng, Rng};

use common::{
    error::{
        Error,
        ErrorKind::{
            CountMismatch, InvalidPararmeterSet, InvalidProof, PointDecompressionError,
            UninitializedValue,
        },
    },
    random::random_scalar,
    utils::mod_pow,
};
use rayon::prelude::*;

#[derive(Clone)]
pub struct Party {
    pub g: Vec<RistrettoPoint>,
    pub private_key: Scalar,
    pub public_key: (CompressedRistretto, RistrettoPoint),
    pub index: usize,
    pub n: usize,
    pub t: usize,
    pub public_keys: Option<Vec<RistrettoPoint>>,
    pub dealer_proof: Option<(Vec<CompressedRistretto>, Vec<RistrettoPoint>)>,
    pub validated_shares: Vec<usize>,
    pub share: Option<Vec<Scalar>>,
    pub d: Option<Scalar>,
    pub shares: Option<Vec<Vec<Scalar>>>,
    pub qualified_set: Option<Vec<(usize, Vec<Scalar>)>>,
}

impl Party {
    pub fn new<R>(
        generator: &RistrettoPoint,
        g: Vec<RistrettoPoint>,
        rng: &mut R,
        n: usize,
        t: usize,
        index: usize,
    ) -> Result<Self, Error>
    where
        R: CryptoRng + Rng,
    {
        let private_key = random_scalar(rng);
        let public_key = generator * &private_key;

        if index <= n && t < n && t as f32 == ((n - 1) as f32 / 2.0).floor() {
            Ok(Self {
                g: g.clone(),
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

    pub fn ingest_share(&mut self, share: &Vec<Scalar>) {
        assert!(share.len() == self.g.len());
        self.share = Some(share.clone());
    }

    pub fn ingest_dealer_proof(&mut self, proof: &Vec<CompressedRistretto>) -> Result<(), Error> {
        if proof.len() != self.t + 1 {
            Err(InvalidProof(format!("c_vals len: {}, t: {}", proof.len(), self.t + 1)).into())
        } else {
            let mut decompressed_c_vals = Vec::with_capacity(self.n);

            for c_i in proof {
                match c_i.decompress() {
                    Some(c) => decompressed_c_vals.push(c),
                    None => {
                        return Err(Error::from_kind(PointDecompressionError(format!(
                            "{c_i:?}",
                        ))));
                    }
                }
            }
            self.dealer_proof = Some((proof.clone(), decompressed_c_vals));
            Ok(())
        }
    }

    pub fn verify_share(&self) -> Result<bool, Error> {
        match &self.dealer_proof {
            Some((_, cvals)) => match &self.share {
                Some(fi) => {
                    let a: RistrettoPoint = fi
                        .par_iter()
                        .zip(self.g.par_iter())
                        .map(|(fik, gk)| fik * gk)
                        .sum();

                    let b: RistrettoPoint = cvals
                        .par_iter()
                        .enumerate()
                        // .map(|(t, c)| c * Scalar::from(self.index.pow(t as u32) as u64))
                        .map(|(t, c)| c * mod_pow(self.index, t))
                        .sum();

                    Ok(a == b)
                }
                None => Err(UninitializedValue("party.share").into()),
            },
            None => Err(UninitializedValue("party.dealer_proof").into()),
        }
    }

    pub fn verify_shares(&mut self) -> Result<bool, Error> {
        match &self.dealer_proof {
            Some((_, cvals)) => match &self.shares {
                Some(shares) => {
                    self.validated_shares = shares
                        .iter()
                        .enumerate()
                        .map(|(i, fi)| {
                            let a = fi
                                .par_iter()
                                .zip(self.g.par_iter())
                                .map(|(fik, gk)| fik * gk)
                                .reduce(|| RistrettoPoint::identity(), |acc, prod| acc + prod);

                            let b = cvals
                                .par_iter()
                                .enumerate()
                                // .map(|(t, c)| c * Scalar::from((i + 1).pow(t as u32) as u64))
                                .map(|(t, c)| c * mod_pow(i + 1, t))
                                .reduce(|| RistrettoPoint::identity(), |acc, prod| acc + prod);

                            if a == b { Some(i) } else { None }
                        })
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

    pub fn ingest_shares(&mut self, shares: &Vec<Vec<Scalar>>) -> Result<(), Error> {
        if shares.len() == self.n && shares[0].len() == self.g.len() {
            self.shares = Some(shares.clone());
            Ok(())
        } else {
            Err(CountMismatch(self.n, "parties", shares.len(), "ingestable shares").into())
        }
    }
}

pub fn generate_parties<R>(
    generator: &RistrettoPoint,
    g: &Vec<RistrettoPoint>,
    rng: &mut R,
    n: usize,
    t: usize,
) -> Vec<Party>
where
    R: CryptoRng + Rng,
{
    (1..=n)
        .map(|i| Party::new(generator, g.clone(), rng, n, t, i).unwrap())
        .collect()
}
