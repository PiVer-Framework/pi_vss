use curve25519_dalek::{RistrettoPoint, Scalar, ristretto::CompressedRistretto};

use rand::{CryptoRng, Rng};
use rayon::prelude::*;

use crate::{
    error::{Error, ErrorKind::CountMismatch},
    random::{random_scalar, random_scalars},
    utils::pointwise_op_in_place,
};

#[derive(Clone)]
pub struct Polynomial {
    pub coefficients: Vec<Scalar>,
}

impl Polynomial {
    pub fn len(&self) -> usize {
        self.coefficients.len()
    }

    pub fn coef_ref<'a>(&'a self) -> &'a Vec<Scalar> {
        &self.coefficients
    }

    pub fn coef_mut<'a>(&'a mut self) -> &'a mut Vec<Scalar> {
        &mut self.coefficients
    }

    pub fn coef_at(&self, index: usize) -> Option<Scalar> {
        if index < self.coefficients.len() - 1 {
            Some(self.coefficients[index])
        } else {
            None
        }
    }
    pub fn coef_at_unchecked<'a>(&'a self, index: usize) -> &'a Scalar {
        &self.coefficients[index]
    }
    pub fn from_coefficients(coefs: Vec<Scalar>) -> Self {
        Self {
            coefficients: coefs,
        }
    }
    pub fn sample<R>(degree: usize, rng: &mut R) -> Self
    where
        R: CryptoRng + Rng,
    {
        Polynomial::from_coefficients(random_scalars(rng, degree + 1))
    }

    pub fn sample_set_f0<R>(degree: usize, rng: &mut R, f0: &Scalar) -> Self
    where
        R: CryptoRng + Rng,
    {
        let mut coefs: Vec<Scalar> = (0..=degree).map(|_| random_scalar(rng)).collect();
        coefs[0] = *f0;

        Polynomial {
            coefficients: coefs,
        }
    }

    pub fn sample_two<R>(degree: usize, rng: &mut R) -> (Self, Self)
    where
        R: CryptoRng + Rng,
    {
        let (coefs_1, coefs_2) = (0..=degree)
            .map(|_| (random_scalar(rng), random_scalar(rng)))
            .collect();

        (
            Polynomial {
                coefficients: coefs_1,
            },
            Polynomial {
                coefficients: coefs_2,
            },
        )
    }

    pub fn sample_n(n: usize, degree: usize) -> Vec<Self> {
        (0..n)
            .into_par_iter()
            .map_init(|| rand::rng(), |mut rng, _| Self::sample(degree, &mut rng))
            .collect()
    }

    pub fn sample_n_set_f0(
        n: usize,
        degree: usize,
        f0_vals: &Vec<Scalar>,
    ) -> Result<Vec<Self>, Error> {
        match f0_vals.len() == n {
            true => Ok((0..n)
                .into_par_iter()
                .zip(f0_vals)
                .map_init(
                    || rand::rng(),
                    |mut rng, (_, f0)| Self::sample_set_f0(degree, &mut rng, &f0),
                )
                .collect()),
            false => Err(CountMismatch(n, "degree", f0_vals.len(), "f0 values").into()),
        }
    }

    pub fn sample_two_set_f0<R>(degree: usize, f0: &Scalar, rng: &mut R) -> (Self, Self)
    where
        R: CryptoRng + Rng,
    {
        let (mut coefs_1, coefs_2): (Vec<Scalar>, Vec<Scalar>) = (0..=degree)
            .map(|_| (random_scalar(rng), random_scalar(rng)))
            .collect();

        coefs_1[0] = *f0;

        (
            Polynomial {
                coefficients: coefs_1,
            },
            Polynomial {
                coefficients: coefs_2,
            },
        )
    }

    // assuming all polynomials are of same degree, panics otherwise
    pub fn evaluate_many_range(polynomials: &[Self], from: usize, to: usize) -> Vec<Vec<Scalar>> {
        (from..=to)
            .into_par_iter()
            .map(|i| {
                let mut x_powers: Vec<Scalar> = vec![Scalar::ONE, Scalar::from(i as u64)];

                for i in 2..polynomials[0].coefficients.len() {
                    x_powers.push(x_powers[1] * x_powers[i - 1]);
                }

                polynomials
                    .par_iter()
                    .map(|polynomial| {
                        polynomial
                            .coefficients
                            .iter()
                            .zip(&x_powers)
                            .map(|(coef, x_pow)| coef * x_pow)
                            .sum()
                    })
                    .collect()
            })
            .collect()
    }

    // assuming all polynomials are of same degree, panics otherwise
    pub fn evaluate_many_range_precomp(
        x_powers: &Vec<Vec<Scalar>>,
        polynomials: &[Self],
        from: usize,
        to: usize,
    ) -> Vec<Vec<Scalar>> {
        (from..=to)
            .into_par_iter()
            .map(|i| {
                polynomials
                    .par_iter()
                    .map(|polynomial| {
                        polynomial
                            .coefficients
                            .iter()
                            .zip(&x_powers[i])
                            .map(|(coef, x_pow)| coef * x_pow)
                            .sum()
                    })
                    .collect()
            })
            .collect()
    }

    // assuming all polynomials are of same degree, panics otherwise
    pub fn evaluate_range_precomp(
        &self,
        x_powers: &Vec<Vec<Scalar>>,
        from: usize,
        to: usize,
    ) -> Vec<Scalar> {
        (from..=to)
            .into_par_iter()
            .map(|i| {
                self.coefficients
                    .par_iter()
                    .zip(&x_powers[i])
                    .map(|(coef, x_pow)| coef * x_pow)
                    .sum()
            })
            .collect()
    }

    // assuming all polynomials are of same degree, panics otherwise
    pub fn evaluate_range(&self, from: usize, to: usize) -> Vec<Scalar> {
        (from..=to)
            .into_par_iter()
            .map(|i| {
                let mut x_powers: Vec<Scalar> = vec![Scalar::ONE, Scalar::from(i as u64)];

                for i in 2..self.coefficients.len() {
                    x_powers.push(x_powers[1] * x_powers[i - 1]);
                }

                self.coefficients
                    .par_iter()
                    .zip(x_powers)
                    .map(|(coef, x_pow)| coef * x_pow)
                    .sum()
            })
            .collect()
    }

    pub fn evaluate_precomp(&self, x_powers: &Vec<Vec<Scalar>>, x: usize) -> Scalar {
        self.coefficients
            .par_iter()
            .zip(&x_powers[x])
            .map(|(coef, x_pow)| coef * x_pow)
            .sum()
    }

    pub fn evaluate(&self, x: usize) -> Scalar {
        let mut x_powers: Vec<Scalar> = vec![Scalar::ONE, Scalar::from(x as u64)];

        for i in 2..self.coefficients.len() {
            x_powers.push(x_powers[1] * x_powers[i - 1]);
        }

        self.coefficients
            .par_iter()
            .zip(x_powers)
            .map(|(coef, x_pow)| coef * x_pow)
            .sum()
    }

    pub fn evaluate_two_range(
        &self,
        other: &Self,
        from: usize,
        to: usize,
    ) -> (Vec<Scalar>, Vec<Scalar>) {
        (from..=to)
            .into_par_iter()
            .map(|i| {
                let mut x_powers: Vec<Scalar> = vec![Scalar::ONE, Scalar::from(i as u64)];

                for i in 2..self.coefficients.len() {
                    x_powers.push(x_powers[1] * x_powers[i - 1]);
                }

                self.coefficients
                    .iter()
                    .zip(other.coefficients.iter())
                    .zip(x_powers)
                    .fold(
                        (Scalar::ZERO, Scalar::ZERO),
                        |(acc_f, acc_r), ((coef_f, coef_r), x_pow)| {
                            ((acc_f + (coef_f * x_pow)), (acc_r + (coef_r * x_pow)))
                        },
                    )
            })
            .unzip()
    }

    pub fn evaluate_multiply(
        &self,
        points: &Vec<RistrettoPoint>,
        x_0: usize,
    ) -> (Vec<Scalar>, Vec<RistrettoPoint>) {
        points
            .par_iter()
            .enumerate()
            .map(|(i, point)| {
                // (i = 0 => 1st party (x = x_0))
                // for pvss x_0 == 1; for ppvss x_0 == 0
                // i is the index of the party
                // i+1 here means start evaluating at x=1
                // x_powers[0] <= constant term multiplier
                // x_powers[1] = x(^1)
                let mut x_powers: Vec<Scalar> = vec![Scalar::ONE, Scalar::from((i + x_0) as u64)];

                // x_powers[2] = x_powers[1] * x_powers[2-1] == x_powers[1] * x_powers[1] = x * x = x^2
                // x_powers[3] = x_powers[1] * x_powers[3-1] == x_powers[1] * x_powers[2] = x * x^2 == x^3
                // ...
                // i+1 here means start evaluating at till x = 32
                for j in 2..self.coefficients.len() {
                    x_powers.push(x_powers[1] * x_powers[j - 1]);
                }

                let f_val = self
                    .coefficients
                    .iter()
                    .zip(x_powers)
                    .fold(Scalar::ZERO, |acc_f, (coef_f, x_pow)| {
                        acc_f + coef_f * x_pow
                    });

                (f_val, (f_val * point))
            })
            .unzip()
    }

    pub fn evaluate_multiply_two_ppvss(
        &self,
        other: &Self,
        points: &Vec<RistrettoPoint>,
    ) -> (Vec<CompressedRistretto>, Vec<CompressedRistretto>) {
        points
            .par_iter()
            .enumerate()
            .map(|(i, point)| {
                // i is the index of the party
                // i here means start evaluating at x=0
                // x_powers[0] <= constant term multiplier
                // x_powers[1] = x(^1)

                let mut x_powers: Vec<Scalar> = Vec::with_capacity(self.coefficients.len());
                x_powers.push(Scalar::ONE);
                x_powers.push(Scalar::from((i) as u64));

                // x_powers[2] = x_powers[1] * x_powers[2-1] == x_powers[1] * x_powers[1] = x * x = x^2
                // x_powers[3] = x_powers[1] * x_powers[3-1] == x_powers[1] * x_powers[2] = x * x^2 == x^3
                // ...
                // i+1 here means start evaluating at till x = 32
                for j in 2..self.coefficients.len() {
                    x_powers.push(x_powers[1] * x_powers[j - 1]);
                }

                let (f_val, r_val) = self
                    .coefficients
                    .par_iter()
                    .zip(other.coefficients.par_iter())
                    .zip(x_powers)
                    .map(|((coef_f, coef_r), x_pow)| (coef_f * x_pow, coef_r * x_pow))
                    .reduce(
                        || (Scalar::ZERO, Scalar::ZERO),
                        |(mut acc_f, mut acc_r), (fval, rval)| {
                            acc_f += fval;
                            acc_r += rval;
                            (acc_f, acc_r)
                        },
                    );

                ((f_val * point).compress(), (r_val * point).compress())
            })
            .collect()
    }

    pub fn evaluate_multiply_two(
        &self,
        other: &Self,
        points: &Vec<RistrettoPoint>,
    ) -> (Vec<CompressedRistretto>, Vec<CompressedRistretto>) {
        points
            .par_iter()
            .enumerate()
            .map(|(i, point)| {
                // i is the index of the party
                // i+1 here means start evaluating at x=1
                // x_powers[0] <= constant term multiplier
                // x_powers[1] = x(^1)
                let mut x_powers: Vec<Scalar> = vec![Scalar::ONE, Scalar::from((i + 1) as u64)];

                // x_powers[2] = x_powers[1] * x_powers[2-1] == x_powers[1] * x_powers[1] = x * x = x^2
                // x_powers[3] = x_powers[1] * x_powers[3-1] == x_powers[1] * x_powers[2] = x * x^2 == x^3
                // ...
                // i+1 here means start evaluating at till x = 32
                for j in 2..self.coefficients.len() {
                    x_powers.push(x_powers[1] * x_powers[j - 1]);
                }

                let (f_val, r_val) = self
                    .coefficients
                    .iter()
                    .zip(other.coefficients.iter())
                    .zip(x_powers)
                    .fold(
                        (Scalar::ZERO, Scalar::ZERO),
                        |(acc_f, acc_r), ((coef_f, coef_r), x_pow)| {
                            (acc_f + coef_f * x_pow, acc_r + coef_r * x_pow)
                        },
                    );

                ((f_val * point).compress(), (r_val * point).compress())
            })
            .collect()
    }

    pub fn sum(&self, p: &Polynomial) -> Self {
        Self {
            coefficients: self
                .coefficients
                .par_iter()
                .zip(p.coefficients.par_iter())
                .map(|(a, b)| a + b)
                .collect(),
        }
    }

    pub fn op_in_place(&mut self, op: fn(Scalar, Scalar) -> Scalar, p2: &Self) {
        pointwise_op_in_place(op, self.coef_mut(), p2.coef_ref());
    }

    pub fn op_many_in_place(&mut self, op: fn(Scalar, Scalar) -> Scalar, ps: &[Polynomial]) {
        ps.iter().for_each(|p| {
            pointwise_op_in_place(op, self.coef_mut(), p.coef_ref());
        });
    }

    pub fn fold_op(op: fn(Scalar, Scalar) -> Scalar, polynomials: &Vec<Self>) -> Self {
        let mut accumulator = polynomials[0].clone();
        accumulator.op_many_in_place(op, &polynomials[1..polynomials.len()]);

        accumulator
    }

    pub fn fold_op_into(&mut self, op: fn(Scalar, Scalar) -> Scalar, polynomials: &Vec<Self>) {
        self.op_many_in_place(op, &polynomials);
    }

    pub fn coef_op(&self, f: fn(Scalar, Scalar) -> Scalar, x: &Scalar) -> Self {
        Self {
            coefficients: self
                .coefficients
                .par_iter()
                .map(|coef| f(*coef, *x))
                .collect(),
        }
    }

    pub fn mul_with_point_compress(&self, point: &RistrettoPoint) -> Vec<CompressedRistretto> {
        self.coefficients
            .par_iter()
            .map(|coef| (coef * point).compress())
            .collect()
    }

    pub fn coef_op_in_place(&mut self, f: fn(Scalar, Scalar) -> Scalar, x: &Scalar) {
        self.coefficients
            .par_iter_mut()
            .for_each(|coef| *coef = f(*coef, *x));
    }

    pub fn mul_sum(&mut self, mul_val: &Scalar, p2: &Self) {
        self.coefficients
            .par_iter_mut()
            .zip(p2.coefficients.par_iter())
            .for_each(|(p1_coef, p2_coef)| {
                *p1_coef *= mul_val;
                *p1_coef += p2_coef
            });
    }
    pub fn mul_many_sum(&mut self, mul_val: &Scalar, p2: &Self) {
        self.coefficients
            .par_iter_mut()
            .zip(p2.coefficients.par_iter())
            .for_each(|(p1_coef, p2_coef)| {
                *p1_coef *= mul_val;
                *p1_coef += p2_coef
            });
    }

    // The input here is &mut r(x), &[f1...fk] , &[d1...dk]
    // z = r + ( ∑ d_j * f_j )
    pub fn compute_z(&mut self, f_polynomials: &[Self], d_vals: &[Scalar]) {
        self.coef_mut()
            .into_par_iter()
            .enumerate()
            .for_each(|(i, r_coef)| {
                *r_coef += f_polynomials
                    .par_iter()
                    .zip(d_vals.par_iter())
                    .map(|(f_k, d_k)| f_k.coef_ref()[i] * d_k)
                    .sum::<Scalar>();
            });
    }

    // The input here is &z(x), &[f1(x)...fk(x)] , &[d1...dk]
    pub fn compute_r_eval(z_eval: &Scalar, f_evals: &[Scalar], d_vals: &[Scalar]) -> Scalar {
        z_eval
            - f_evals
                .par_iter()
                .zip(d_vals.par_iter())
                .map(|(f_k, d_k)| f_k * d_k)
                .sum::<Scalar>()
    }
}

impl std::fmt::Display for Polynomial {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{}",
            self.coefficients
                .par_iter()
                .map(|coef| format!("{:x?}", coef.as_bytes()))
                .collect::<Vec<String>>()
                .join(",")
        )
    }
}

#[cfg(test)]
mod test {
    use curve25519_dalek::Scalar;

    use crate::{polynomial::Polynomial, random::random_scalar};

    #[test]
    fn test_thing() {
        // f1(x) = 13x^3 + 2x^2 + 7x + 128
        let f1 = Polynomial {
            coefficients: vec![128u8.into(), 7u8.into(), 2u8.into(), 13u8.into()],
        };

        // f2(x) = 81x^3 + 7x^2 + 153x + 32
        let f2 = Polynomial {
            coefficients: vec![32u8.into(), 153u8.into(), 7u8.into(), 81u8.into()],
        };

        // r(x) = 7x^3 + 15x^2 + 81x + 2
        let r = Polynomial {
            coefficients: vec![2u8.into(), 81u8.into(), 15u8.into(), 7u8.into()],
        };

        let f1_at_5 = f1.evaluate(5);
        let f2_at_5 = f2.evaluate(5);
        let r_at_5 = r.evaluate(5);

        assert_eq!(f1_at_5, Scalar::from(1838u16));
        assert_eq!(f2_at_5, Scalar::from(11097u16));
        assert_eq!(r_at_5, Scalar::from(1657u16));

        let d1 = Scalar::from(127u8);
        let d2 = Scalar::from(17u8);

        // Proof step
        let mut r_test_z1 = r.clone();
        r_test_z1.compute_z(&[f1.clone()], &[d1]);

        let z1_at_5 = r_test_z1.evaluate(5);
        assert_eq!(z1_at_5, Scalar::from(235083u32));

        let mut r_test_z2 = r.clone();
        r_test_z2.compute_z(&[f1.clone(), f2.clone()], &[d1, d2]);

        let z2_at_5 = r_test_z2.evaluate(5);
        assert_eq!(z2_at_5, Scalar::from(423732u32));

        // Verification step below

        let potential_r_5_1 = Polynomial::compute_r_eval(&z1_at_5, &[f1_at_5], &[d1]);

        let potential_r_5_2 = Polynomial::compute_r_eval(&z2_at_5, &[f1_at_5, f2_at_5], &[d1, d2]);

        assert_eq!(potential_r_5_1, r_at_5);
        assert_eq!(potential_r_5_2, r_at_5);
    }

    #[test]
    fn test_thing_big() {
        let mut rng = rand::rng();
        let polynomials = Polynomial::sample_n(3, 10);

        let (f1, f2, r) = (
            polynomials[0].clone(),
            polynomials[1].clone(),
            polynomials[2].clone(),
        );

        let f1_at_5 = f1.evaluate(5);
        let f2_at_5 = f2.evaluate(5);
        let r_at_5 = r.evaluate(5);

        let d1 = random_scalar(&mut rng);
        let d2 = random_scalar(&mut rng);

        // Proof step
        let mut r_test_z1 = r.clone();
        r_test_z1.compute_z(&[f1.clone()], &[d1]);

        let z1_at_5 = r_test_z1.evaluate(5);

        let mut r_test_z2 = r.clone();
        r_test_z2.compute_z(&[f1.clone(), f2.clone()], &[d1, d2]);

        let z2_at_5 = r_test_z2.evaluate(5);

        // Verification step below

        let potential_r_5_1 = Polynomial::compute_r_eval(&z1_at_5, &[f1_at_5], &[d1]);

        let potential_r_5_2 = Polynomial::compute_r_eval(&z2_at_5, &[f1_at_5, f2_at_5], &[d1, d2]);

        assert_eq!(potential_r_5_1, r_at_5);
        assert_eq!(potential_r_5_2, r_at_5);
    }
}
