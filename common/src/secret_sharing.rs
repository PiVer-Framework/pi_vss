use curve25519_dalek::{RistrettoPoint, Scalar, ristretto::CompressedRistretto};
use rand::{CryptoRng, Rng, seq::SliceRandom};

use crate::{
    error::{
        Error,
        ErrorKind::{InsufficientShares, UninitializedValue},
    },
    polynomial::Polynomial,
};
use rayon::prelude::*;

pub fn generate_encrypted_shares_batched(
    t: usize,
    x_pows: &Vec<Vec<Scalar>>,
    public_keys: &[RistrettoPoint],
    secrets: &Vec<Scalar>,
) -> (Vec<Polynomial>, Vec<Vec<CompressedRistretto>>) {
    let (f_polynomials, f_evals) = generate_shares_batched(public_keys.len(), t, x_pows, secrets);

    let encrypted_shares = f_evals
        .par_iter()
        .zip(public_keys.par_iter())
        .map(|(fk, pub_key)| {
            fk.par_iter()
                .map(|f_eval| (f_eval * pub_key).compress())
                .collect()
        })
        .collect();

    (f_polynomials, encrypted_shares)
}

pub fn decrypt_share(private_key: &Scalar, encrypted_share: &RistrettoPoint) -> RistrettoPoint {
    private_key.invert() * encrypted_share
}

pub fn generate_encrypted_shares<R>(
    rng: &mut R,
    t: usize,
    x_pows: &Vec<Vec<Scalar>>,
    public_keys: &[RistrettoPoint],
    secret: &Scalar,
) -> (Polynomial, Vec<CompressedRistretto>)
where
    R: CryptoRng,
{
    let (f_polynomial, f_evals) = generate_shares(rng, public_keys.len(), t, x_pows, secret);

    let encrypted_shares = f_evals
        .par_iter()
        .zip(public_keys.par_iter())
        .map(|(fi, pub_key)| (fi * pub_key).compress())
        .collect();

    (f_polynomial, encrypted_shares)
}

pub fn generate_shares_batched(
    n: usize,
    t: usize,
    x_pows: &Vec<Vec<Scalar>>,
    secrets: &Vec<Scalar>,
) -> (Vec<Polynomial>, Vec<Vec<Scalar>>) {
    // This contains k * f_polynomial
    let f_polynomials = Polynomial::sample_n_set_f0(secrets.len(), t, secrets).unwrap();
    // evals is vec[vec[k]; n]
    let f_evals = Polynomial::evaluate_many_range_precomp(x_pows, &f_polynomials, 1, n);
    (f_polynomials, f_evals)
}

pub fn generate_shares<R>(
    rng: &mut R,
    n: usize,
    t: usize,
    x_pows: &Vec<Vec<Scalar>>,
    secret: &Scalar,
) -> (Polynomial, Vec<Scalar>)
where
    R: CryptoRng,
{
    let polynomial = Polynomial::sample_set_f0(t, rng, secret);
    let evals = polynomial.evaluate_range_precomp(x_pows, 1, n);
    (polynomial, evals)
}

pub fn select_qualified_set<R, T>(
    rng: &mut R,
    t: usize,
    shares: &Option<Vec<T>>,
    validated_shares: &Vec<usize>,
) -> Result<Vec<(usize, T)>, Error>
where
    T: Clone,
    R: CryptoRng + Rng,
{
    match shares {
        Some(shares) => {
            if shares.len() > t {
                let mut tmp = validated_shares.clone();
                tmp.shuffle(rng);

                Ok(tmp
                    .into_iter()
                    .take(t + 1)
                    .map(|x| (x + 1, shares[x].clone()))
                    .collect())
            } else {
                Err(InsufficientShares(validated_shares.len(), t).into())
            }
        }
        None => Err(UninitializedValue("party.{shares || decrypted_shares}").into()),
    }
}

pub fn reconstruct_secrets_exponent(
    qualified_set: &Option<Vec<(usize, Vec<RistrettoPoint>)>>,
    lambdas: &Vec<Scalar>,
) -> Result<Vec<RistrettoPoint>, Error> {
    match qualified_set {
        Some(qualified_set) => {
            let k = qualified_set[0].1.len();

            Ok((0..k)
                .into_par_iter()
                .map(|k| {
                    qualified_set
                        .par_iter()
                        .zip(lambdas.par_iter())
                        .map(|((_, poly_share), lambda)| lambda * poly_share[k])
                        .sum()
                })
                .collect())
        }
        None => Err(UninitializedValue("party.qualified_set").into()),
    }
}

pub fn reconstruct_secrets(
    qualified_set: &Option<Vec<(usize, Vec<Scalar>)>>,
    lambdas: &Vec<Scalar>,
) -> Result<Vec<Scalar>, Error> {
    match qualified_set {
        Some(qualified_set) => {
            let k = qualified_set[0].1.len();
            Ok((0..k)
                .into_par_iter()
                .map(|k| {
                    qualified_set
                        .par_iter()
                        .zip(lambdas.par_iter())
                        .map(|((_, poly_share), lambda)| lambda * poly_share[k])
                        .sum()
                })
                .collect())
        }
        None => Err(UninitializedValue("party.qualified_set").into()),
    }
}

pub fn reconstruct_secret(
    qualified_set: &Option<Vec<(usize, Scalar)>>,
    lambdas: &Vec<Scalar>,
) -> Result<Scalar, Error> {
    match qualified_set {
        Some(qualified_set) => Ok(qualified_set
            .par_iter()
            .zip(lambdas.par_iter())
            .map(|((_, decrypted_share), lambda)| lambda * decrypted_share)
            .sum()),
        None => Err(UninitializedValue("party.qualified_set").into()),
    }
}

pub fn reconstruct_secret_exponent(
    qualified_set: &Option<Vec<(usize, RistrettoPoint)>>,
    lambdas: &Vec<Scalar>,
) -> Result<RistrettoPoint, Error> {
    match qualified_set {
        Some(qualified_set) => Ok(qualified_set
            .par_iter()
            .zip(lambdas.par_iter())
            .map(|((_, decrypted_share), lambda)| lambda * decrypted_share)
            .sum()),
        None => Err(UninitializedValue("party.qualified_set").into()),
    }
}

#[cfg(test)]
mod test {
    use curve25519_dalek::RistrettoPoint;

    use crate::{
        precompute::gen_powers,
        random::{random_scalar, random_scalars},
        secret_sharing::{
            decrypt_share, generate_encrypted_shares, generate_encrypted_shares_batched,
            generate_shares, generate_shares_batched, reconstruct_secret,
            reconstruct_secret_exponent, reconstruct_secrets, reconstruct_secrets_exponent,
            select_qualified_set,
        },
        utils::compute_lagrange_bases,
    };

    use rayon::prelude::*;

    #[test]
    fn gen_shares() {
        let n = 128;
        let t = 63;
        let mut rng = rand::rng();

        let x_pows = gen_powers(n, t);
        let secret = random_scalar(&mut rng);

        let shares = generate_shares(&mut rng, n, t, &x_pows, &secret);

        let qualified_set = select_qualified_set(
            &mut rng,
            t,
            &Some(shares.1),
            &(0..n).collect::<Vec<usize>>(),
        )
        .unwrap();

        let indices: Vec<usize> = qualified_set.iter().map(|(index, _)| *index).collect();

        let lagrange_bases = compute_lagrange_bases(&indices);

        let q = Some(qualified_set);

        assert_eq!(secret, reconstruct_secret(&q, &lagrange_bases).unwrap());
    }
    #[test]
    fn gen_shares_batch() {
        let n = 128;
        let t = 63;
        let k = 10;

        let mut rng = rand::rng();

        let x_pows = gen_powers(n, t);

        let secrets = random_scalars(&mut rng, k);
        let shares = generate_shares_batched(n, t, &x_pows, &secrets);

        let qualified_set = select_qualified_set(
            &mut rng,
            t,
            &Some(shares.1),
            &(0..n).collect::<Vec<usize>>(),
        )
        .unwrap();

        let indices: Vec<usize> = qualified_set.iter().map(|(index, _)| *index).collect();

        let lagrange_bases = compute_lagrange_bases(&indices);

        let q = Some(qualified_set);

        assert_eq!(secrets, reconstruct_secrets(&q, &lagrange_bases).unwrap())
    }
    #[test]
    fn gen_encrypted_shares() {
        let n = 128;
        let t = 63;
        let mut rng = rand::rng();
        let x_pows = gen_powers(n, t);

        let private_keys = random_scalars(&mut rng, n);

        let public_keys: Vec<RistrettoPoint> = private_keys
            .par_iter()
            .map(|private_key| RistrettoPoint::mul_base(private_key))
            .collect();

        let secret = random_scalar(&mut rng);

        let (f, encrypted_shares) =
            generate_encrypted_shares(&mut rng, t, &x_pows, &public_keys, &secret);

        assert_eq!(f.coef_at_unchecked(0), &secret);

        let decrypted_shares: Vec<RistrettoPoint> = encrypted_shares
            .par_iter()
            .zip(private_keys.par_iter())
            .map(|(encrypted_share, private_key)| {
                decrypt_share(private_key, &encrypted_share.decompress().unwrap())
            })
            .collect();

        let qualified_set = select_qualified_set(
            &mut rng,
            t,
            &Some(decrypted_shares),
            &(0..n).collect::<Vec<usize>>(),
        )
        .unwrap();

        let indices: Vec<usize> = qualified_set.iter().map(|(index, _)| *index).collect();

        let lagrange_bases = compute_lagrange_bases(&indices);

        let q = Some(qualified_set);

        let secret_exp = RistrettoPoint::mul_base(&secret);

        assert_eq!(
            secret_exp,
            reconstruct_secret_exponent(&q, &lagrange_bases).unwrap()
        );
    }
    #[test]
    fn gen_encrypted_shares_batch() {
        let n = 128;
        let t = 63;
        let k = 10;

        let mut rng = rand::rng();
        let x_pows = gen_powers(n, t);

        let private_keys = random_scalars(&mut rng, n);

        let public_keys: Vec<RistrettoPoint> = private_keys
            .par_iter()
            .map(|private_key| RistrettoPoint::mul_base(private_key))
            .collect();

        let secrets = random_scalars(&mut rng, k);

        let (_fk, encrypted_shares) =
            generate_encrypted_shares_batched(t, &x_pows, &public_keys, &secrets);

        let decrypted_shares: Vec<Vec<RistrettoPoint>> = encrypted_shares
            .par_iter()
            .zip(private_keys.par_iter())
            .map(|(encrypted_shares_i, private_key)| {
                encrypted_shares_i
                    .par_iter()
                    .map(|encrypted_share| {
                        decrypt_share(private_key, &encrypted_share.decompress().unwrap())
                    })
                    .collect()
            })
            .collect();

        let qualified_set = select_qualified_set(
            &mut rng,
            t,
            &Some(decrypted_shares),
            &(0..n).collect::<Vec<usize>>(),
        )
        .unwrap();

        let indices: Vec<usize> = qualified_set.iter().map(|(index, _)| *index).collect();

        let lagrange_bases = compute_lagrange_bases(&indices);

        let q = Some(qualified_set);

        let secret_exps: Vec<RistrettoPoint> = secrets
            .par_iter()
            .map(|secret| RistrettoPoint::mul_base(secret))
            .collect();

        assert_eq!(
            secret_exps,
            reconstruct_secrets_exponent(&q, &lagrange_bases).unwrap()
        );
    }
}
