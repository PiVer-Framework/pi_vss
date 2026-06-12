use blake3::Hasher;
use curve25519_dalek::{RistrettoPoint, Scalar, ristretto::CompressedRistretto, traits::Identity};
use rand::{CryptoRng, Rng};
use zeroize::Zeroize;

use common::{
    error::{
        Error,
        ErrorKind::{CountMismatch, InvalidPararmeterSet, InvalidProof, UninitializedValue},
    },
    polynomial::Polynomial,
    random::random_scalar,
    utils::{batch_decompress_batched_ristretto_points, compute_d_powers},
};
use rayon::prelude::*;

#[derive(Clone)]
pub struct Party {
    pub g: RistrettoPoint,
    pub private_key: Scalar,
    pub public_key: (CompressedRistretto, RistrettoPoint),
    pub index: usize,
    pub n: usize,
    pub t: usize,
    pub public_keys: Option<Vec<RistrettoPoint>>,
    pub dealer_proof: Option<(Scalar, Polynomial)>,
    pub validated_shares: Vec<usize>,
    pub encrypted_share: Option<Vec<RistrettoPoint>>,
    pub decrypted_share: Option<Vec<RistrettoPoint>>,
    pub encrypted_shares: Option<(Vec<Vec<CompressedRistretto>>, Vec<Vec<RistrettoPoint>>)>,
    pub decrypted_shares: Option<Vec<Vec<RistrettoPoint>>>,
    pub share_proof: Option<(Scalar, Scalar)>,
    pub share_proofs: Option<Vec<(Scalar, Scalar)>>,
    pub qualified_set: Option<Vec<(usize, Vec<RistrettoPoint>)>>,
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
                g: g.clone(),
                private_key,
                public_key: (public_key.compress(), public_key),
                index,
                n,
                t,
                dealer_proof: None,
                public_keys: None,
                validated_shares: vec![],
                encrypted_shares: None,
                decrypted_shares: None,
                encrypted_share: None,
                decrypted_share: None,
                share_proof: None,
                share_proofs: None,
                qualified_set: None,
            })
        } else {
            Err(InvalidPararmeterSet(n, t as isize, index).into())
        }
    }

    pub fn ingest_encrypted_share(&mut self, share: &Vec<CompressedRistretto>) {
        self.encrypted_share = Some(
            share
                .par_iter()
                .map(|share_k| share_k.decompress().unwrap())
                .collect(),
        );
    }

    pub fn ingest_dealer_proof(&mut self, proof: (&Scalar, &Polynomial)) -> Result<(), Error> {
        if proof.1.len() != self.t + 1 {
            Err(InvalidProof(format!("z len: {}, t: {}", proof.1.len(), self.t)).into())
        } else {
            self.dealer_proof = Some((proof.0.clone(), proof.1.clone()));
            Ok(())
        }
    }

    pub fn verify_encrypted_shares(
        &mut self,
        hasher: &mut Hasher,
        buf: &mut [u8; 64],
        x_pows: &Vec<Vec<Scalar>>,
    ) -> Result<bool, Error> {
        match &self.dealer_proof {
            Some((d, z)) => match (&self.encrypted_shares, &self.public_keys) {
                (Some(encrypted_shares), Some(public_keys)) => {
                    hasher.reset();
                    buf.zeroize();
                    let k = encrypted_shares.0[0].len();

                    let d_vals = compute_d_powers(k, d);

                    let z_evals = z.evaluate_range_precomp(x_pows, 1, public_keys.len());

                    let suite: Vec<CompressedRistretto> = z_evals
                        .par_iter()
                        .zip(public_keys.par_iter().zip(encrypted_shares.1.par_iter()))
                        .map(|(z_eval, (public_key, encrypted_shares_i))| {
                            ((z_eval * public_key)
                                - d_vals
                                    .par_iter()
                                    .zip(encrypted_shares_i)
                                    .map(|(d_val, encrypted_shares_i_k)| {
                                        encrypted_shares_i_k * d_val
                                    })
                                    .reduce(|| RistrettoPoint::identity(), |acc, x| acc + x))
                            // .fold(RistrettoPoint::identity(), |acc, x| acc + x))
                            .compress()
                        })
                        .collect();

                    encrypted_shares
                        .0
                        .iter()
                        .flatten()
                        .chain(suite.iter())
                        .for_each(|x| {
                            hasher.update(x.as_bytes());
                        });

                    hasher.finalize_xof().fill(buf);
                    hasher.reset();

                    let d_comp = Scalar::from_bytes_mod_order_wide(buf);
                    buf.zeroize();

                    Ok(*d == d_comp)
                }
                (Some(_), None) => Err(UninitializedValue("party.public_keys").into()),
                (None, Some(_)) => Err(UninitializedValue("party.encrypted_shares").into()),
                (None, None) => {
                    Err(UninitializedValue("party.{encrypted_shares, public_keys}").into())
                }
            },
            None => Err(UninitializedValue("party.dealer_proof").into()),
        }
    }

    pub fn decrypt_shares(&mut self) -> Result<(), Error> {
        let inv_private_key = self.private_key.invert();
        match &self.encrypted_share {
            Some(encrypted_share) => {
                self.decrypted_share = Some(
                    encrypted_share
                        .par_iter()
                        .map(|enc_share| enc_share * inv_private_key)
                        .collect(),
                );
                Ok(())
            }
            None => Err(UninitializedValue("party.encrypted_share").into()),
        }
    }

    pub fn dleq_share<R>(
        &mut self,
        g: &RistrettoPoint,
        rng: &mut R,
        hasher: &mut Hasher,
        buf: &mut [u8; 64],
    ) -> Result<(), Error>
    where
        R: CryptoRng + Rng,
    {
        match (&self.decrypted_share, &self.encrypted_share) {
            (Some(decrypted_shares), Some(encrypted_shares)) => {
                let r = random_scalar(rng);
                // g^si
                hasher.update(self.public_key.0.as_bytes());
                // g^r
                hasher.update((g * &r).compress().as_bytes());

                decrypted_shares.iter().zip(encrypted_shares).for_each(
                    |(decrypted_share, encrypted_share)| {
                        // g^(si * fi)
                        hasher.update(encrypted_share.compress().as_bytes());
                        // g^(fi * r)
                        hasher.update((decrypted_share * r).compress().as_bytes());
                    },
                );

                hasher.finalize_xof().fill(buf);

                let d = Scalar::from_bytes_mod_order_wide(buf);
                let z = r + d * self.private_key;
                hasher.reset();
                buf.zeroize();
                self.share_proof = Some((d, z));

                Ok(())
            }
            (None, Some(_)) => Err(UninitializedValue("party.decrypted_share").into()),
            (Some(_), None) => Err(UninitializedValue("party.encrypted_shares").into()),
            (None, None) => {
                Err(UninitializedValue("party.{decrypted_share, encrypted_shares}").into())
            }
        }
    }

    pub fn ingest_encrypted_shares(
        &mut self,
        encrypted_shares: &Vec<Vec<CompressedRistretto>>,
    ) -> Result<(), Error> {
        if encrypted_shares.len() == self.n {
            match batch_decompress_batched_ristretto_points(encrypted_shares) {
                Ok(enc_shares) => {
                    self.encrypted_share = Some(enc_shares[self.index - 1].clone());
                    self.encrypted_shares = Some((encrypted_shares.to_vec(), enc_shares));
                    Ok(())
                }
                Err(x) => Err(x),
            }
        } else {
            Err(CountMismatch(
                self.n,
                "parties",
                encrypted_shares.len(),
                "encrypted shares",
            )
            .into())
        }
    }

    pub fn ingest_decrypted_shares_and_proofs(
        &mut self,
        decrypted_shares: &Vec<Vec<CompressedRistretto>>,
        proofs: Vec<(Scalar, Scalar)>,
    ) -> Result<(), Error> {
        if decrypted_shares.len() == self.n - 1 {
            if proofs.len() == decrypted_shares.len() {
                match batch_decompress_batched_ristretto_points(decrypted_shares) {
                    Ok(mut dec_shares) => match (&self.decrypted_share, &self.share_proof) {
                        (Some(own_dec_share), Some(own_proof)) => {
                            dec_shares.insert(self.index - 1, own_dec_share.clone());
                            self.decrypted_shares = Some(dec_shares);
                            let mut proofs = proofs;
                            proofs.insert(self.index - 1, own_proof.clone());
                            self.share_proofs = Some(proofs);
                            Ok(())
                        }
                        (None, Some(_)) => Err(UninitializedValue("party.decrypted_share").into()),
                        (Some(_), None) => Err(UninitializedValue("party.share_proof").into()),
                        (None, None) => {
                            Err(UninitializedValue("party.{decrypted_share, share_proof}").into())
                        }
                    },
                    Err(x) => Err(x),
                }
            } else {
                Err(CountMismatch(self.n, "parties", proofs.len(), "proofs").into())
            }
        } else {
            Err(CountMismatch(
                self.n,
                "parties",
                decrypted_shares.len(),
                "decrypted shares",
            )
            .into())
        }
    }

    pub fn verify_decrypted_shares(&mut self) -> Result<bool, Error> {
        match (&self.public_keys, &self.encrypted_shares) {
            (Some(public_keys), Some(enc_shares)) => {
                match (&self.decrypted_shares, &self.share_proofs) {
                    (Some(dec_shares), Some(proofs)) => {
                        self.validated_shares = dec_shares
                            .par_iter()
                            .zip(
                                proofs
                                    .par_iter()
                                    .zip(public_keys.par_iter().zip(enc_shares.1.par_iter())),
                            )
                            .enumerate()
                            .map_init(|| (blake3::Hasher::new(), [0u8; 64]),
                        |(hasher, buf),(i, (dec_share, ((d, z), (public_key, enc_share))))| {
                            hasher.update(public_key.compress().as_bytes());

                            hasher.update(((self.g * z) -(public_key*d)).compress().as_bytes());
dec_share.iter().zip(enc_share).for_each(|(dec_share_k, enc_share_k)| {
                                    hasher.update(enc_share_k.compress().as_bytes());
                                    hasher.update(((dec_share_k * z)-(enc_share_k *d)).compress().as_bytes());
                                });

                                hasher.finalize_xof().fill(buf);

                                hasher.reset();
                                    let reconstructed_d = Scalar::from_bytes_mod_order_wide(buf);

                                    buf.zeroize();
                               if reconstructed_d == *d {
                                    Some(i)
                                } else {
                                    None
                                }
                            })
                            .filter(Option::is_some)
                            .map(|res| res.unwrap())
                            .collect();
                        Ok(self.validated_shares.len() > self.t)
                    }
                    (None, Some(_)) => Err(UninitializedValue("party.decrypted_shares").into()),
                    (Some(_), None) => Err(UninitializedValue("party.share_proofs").into()),
                    (None, None) => {
                        Err(UninitializedValue("party.{decrypted_shares, share_proofs}").into())
                    }
                }
            }
            (None, Some(_)) => Err(UninitializedValue("party.encrypted_shares").into()),
            (Some(_), None) => Err(UninitializedValue("party.public_keys").into()),
            (None, None) => Err(UninitializedValue("party.{public_keys, encrypted_shares}").into()),
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
