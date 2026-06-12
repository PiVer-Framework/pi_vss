use curve25519_dalek::{RistrettoPoint, Scalar};
use rand::*;

use zeroize::Zeroize;

pub fn random_scalar<R>(rng: &mut R) -> Scalar
where
    R: CryptoRng + Rng,
{
    let mut bytes = [0u8; 64];
    rng.fill_bytes(&mut bytes);
    let scalar = Scalar::from_bytes_mod_order_wide(&bytes);
    bytes.zeroize();
    scalar
}

pub fn random_scalars<R>(rng: &mut R, n: usize) -> Vec<Scalar>
where
    R: CryptoRng + Rng,
{
    let mut bytes = [0u8; 64];
    (0..n)
        .map(|_| {
            rng.fill_bytes(&mut bytes);
            let scalar = Scalar::from_bytes_mod_order_wide(&bytes);
            bytes.zeroize();
            scalar
        })
        .collect()
}

pub fn random_point<R>(rng: &mut R) -> RistrettoPoint
where
    R: CryptoRng + Rng,
{
    let mut bytes = [0u8; 64];
    rng.fill_bytes(&mut bytes);
    let point = RistrettoPoint::from_uniform_bytes(&bytes);
    bytes.zeroize();
    point
}

pub fn random_points<R>(rng: &mut R, n: usize) -> Vec<RistrettoPoint>
where
    R: CryptoRng + Rng,
{
    let mut bytes = [0u8; 64];
    (0..n)
        .map(|_| {
            rng.fill_bytes(&mut bytes);
            let point = RistrettoPoint::from_uniform_bytes(&bytes);
            bytes.zeroize();
            point
        })
        .collect()
}

#[cfg(test)]
mod test {
    use curve25519_dalek::Scalar;
    use rand::*;
    use rayon::prelude::*;

    #[test]
    fn test_rand() {
        let mut rng = rand::rng();
        let n = 10;

        let mut bytes = vec![0u8; 64 * n];
        rng.fill_bytes(&mut bytes);

        let v: Vec<Scalar> = bytes
            .as_chunks::<64>()
            .0
            .into_par_iter()
            .map(|chunk| Scalar::from_bytes_mod_order_wide(chunk))
            .collect();

        println!("{:?}", v);
    }
}
