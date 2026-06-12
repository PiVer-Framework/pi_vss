use blake3::Hasher;
use curve25519_dalek::{RistrettoPoint, Scalar, ristretto::CompressedRistretto};
use rand::{CryptoRng, Rng};
use zeroize::Zeroize;

use common::{
    error::{
        Error,
        ErrorKind::{
            CountMismatch, InvalidPararmeterSet, InvalidProof, PointDecompressionError,
            UninitializedValue,
        },
    },
    polynomial::Polynomial,
    random::random_scalar,
    utils::compute_d_powers_from_point_commitments,
};
use rayon::prelude::*;

#[derive(Clone)]
pub struct Party {
    pub g1: RistrettoPoint,
    pub g2: RistrettoPoint,
    pub private_key: Scalar,
    pub public_key: (CompressedRistretto, RistrettoPoint),
    pub index: usize,
    pub n: usize,
    pub t: usize,
    pub public_keys: Option<Vec<RistrettoPoint>>,
    pub dealer_proof: Option<(Vec<CompressedRistretto>, Vec<RistrettoPoint>, Polynomial)>,
    pub validated_shares: Vec<usize>,
    // (fi, gamma_i)
    pub share: Option<(Vec<Scalar>, Scalar)>,
    pub d: Option<Scalar>,
    pub shares: Option<Vec<(Vec<Scalar>, Scalar)>>,
    pub qualified_set: Option<Vec<(usize, Vec<Scalar>)>>,
}

impl Party {
    pub fn new<R>(
        g: &RistrettoPoint,
        g1: RistrettoPoint,
        g2: RistrettoPoint,
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
                g1: g1.clone(),
                g2: g2.clone(),
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

    pub fn ingest_share(&mut self, share: (&Vec<Scalar>, &Scalar)) {
        self.share = Some((share.0.clone(), share.1.clone()));
    }

    pub fn ingest_dealer_proof(
        &mut self,
        proof: (&Vec<CompressedRistretto>, &Polynomial),
    ) -> Result<(), Error> {
        if proof.1.len() != self.t + 1 {
            Err(InvalidProof(format!("z len: {}, t: {}", proof.1.len(), self.t + 1)).into())
        } else if proof.0.len() != self.n {
            Err(InvalidProof(format!("c_vals len: {}, n: {}", proof.0.len(), self.n)).into())
        } else {
            let mut decompressed_c_vals = Vec::with_capacity(self.n);

            for c_i in proof.0 {
                match c_i.decompress() {
                    Some(c) => decompressed_c_vals.push(c),
                    None => {
                        return Err(Error::from_kind(PointDecompressionError(format!(
                            "{c_i:?}",
                        ))));
                    }
                }
            }
            self.dealer_proof = Some((proof.0.clone(), decompressed_c_vals, proof.1.clone()));
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
            Some((compressed_cvals, cvals, z)) => match &self.share {
                Some((fi, gi)) => {
                    let k = fi.len();

                    let d_vals =
                        compute_d_powers_from_point_commitments(hasher, buf, &compressed_cvals, k);

                    let zi = z.evaluate_precomp(x_pows, self.index);

                    let expected_c = cvals[self.index - 1];

                    fi.iter().for_each(|fi_k| {
                        hasher.update(fi_k.as_bytes());
                    });

                    hasher.update(Polynomial::compute_r_eval(&zi, &fi, &d_vals).as_bytes());

                    hasher.finalize_xof().fill(buf);
                    hasher.reset();

                    let h = Scalar::from_bytes_mod_order_wide(buf);
                    buf.zeroize();

                    let c = (self.g1 * h) + (self.g2 * gi);

                    Ok(expected_c == c)
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
            Some((compressed_cvals, cvals, z)) => match &self.shares {
                Some(shares) => {
                    let k = shares[0].0.len();

                    let d_vals =
                        compute_d_powers_from_point_commitments(hasher, buf, &compressed_cvals, k);
                    let z_evals = z.evaluate_range_precomp(x_pows, 1, self.n);

                    self.validated_shares = shares
                        .par_iter()
                        .zip(z_evals.par_iter())
                        .enumerate()
                        .map_init(
                            || (Hasher::new(), [0u8; 64]),
                            |(l_hasher, l_buf), (i, ((fi, gi), zi))| {
                                fi.iter().for_each(|fi_k| {
                                    l_hasher.update(fi_k.as_bytes());
                                });

                                l_hasher.update(
                                    Polynomial::compute_r_eval(&zi, &fi, &d_vals).as_bytes(),
                                );

                                l_hasher.finalize_xof().fill(l_buf);
                                l_hasher.reset();

                                let h = Scalar::from_bytes_mod_order_wide(l_buf);
                                l_buf.zeroize();

                                if cvals[i] == ((self.g1 * h) + (self.g2 * gi)) {
                                    Some(i)
                                } else {
                                    None
                                }
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

    pub fn ingest_shares(
        &mut self,
        shares: (&Vec<Vec<Scalar>>, &Vec<Scalar>),
    ) -> Result<(), Error> {
        if shares.0.len() == self.n && shares.1.len() == self.n {
            self.shares = Some(
                shares
                    .0
                    .clone()
                    .into_iter()
                    .zip(shares.1.clone().into_iter())
                    .collect(),
            );
            Ok(())
        } else {
            Err(CountMismatch(self.n, "parties", shares.0.len(), "ingestable shares").into())
        }
    }
}

pub fn generate_parties<R>(
    g: &RistrettoPoint,
    g1: &RistrettoPoint,
    g2: &RistrettoPoint,
    rng: &mut R,
    n: usize,
    t: usize,
) -> Vec<Party>
where
    R: CryptoRng + Rng,
{
    (1..=n)
        .map(|i| Party::new(g, g1.clone(), g2.clone(), rng, n, t, i).unwrap())
        .collect()
}
