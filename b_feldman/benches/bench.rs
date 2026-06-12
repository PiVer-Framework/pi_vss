use b_feldman::{dealer::Dealer, party::generate_parties};
use criterion::{BatchSize, Criterion, criterion_group, criterion_main};
use curve25519_dalek::{RistrettoPoint, ristretto::CompressedRistretto, traits::Identity};

use common::{
    BENCH_K, BENCH_N_T,
    precompute::gen_powers,
    random::{random_point, random_points, random_scalars},
    secret_sharing::generate_shares_batched,
    utils::ingest_public_keys,
};

fn vss(c: &mut Criterion) {
    for (n, t) in BENCH_N_T {
        let mut rng = rand::rng();

        let generator: RistrettoPoint = random_point(&mut rng);

        let xpows = gen_powers(n, t);
        for k in 1..=16 {
            let g: Vec<RistrettoPoint> = random_points(&mut rng, k);
            let mut parties = generate_parties(&generator, &g, &mut rng, n, t);

            let public_keys: Vec<CompressedRistretto> =
                parties.iter().map(|party| party.public_key.0).collect();

            let mut dealer = Dealer::new(g, n, t, &public_keys).unwrap();

            for party in &mut parties {
                let public_keys: Vec<CompressedRistretto> = public_keys
                    .iter()
                    .filter(|pk| &party.public_key.0 != *pk)
                    .copied()
                    .collect();

                party.public_keys = Some(
                    ingest_public_keys(n, &party.public_key.1, party.index, &public_keys).unwrap(),
                );
            }

            let secrets = random_scalars(&mut rng, k);

            let (f_polynomials, _) = generate_shares_batched(n, t, &xpows, &secrets);

            c.bench_function(
                &format!(
                    "(k: {}, n: {}, t: {}) | B_Feldman VSS | Dealer: Generate Proof",
                    k, n, t
                ),
                |b| {
                    b.iter_batched(
                        || vec![CompressedRistretto::identity(); t + 1],
                        |mut c_buf| dealer.generate_proof(&mut c_buf, &f_polynomials),
                        BatchSize::PerIteration,
                    )
                },
            );

            let (shares, c_vals) = dealer.deal_secret(&xpows, &secrets);

            let p = &mut parties[0];
            p.ingest_dealer_proof(&c_vals).unwrap();

            p.ingest_shares(&shares).unwrap();

            c.bench_function(
                &format!(
                    "(k: {}, n: {}, t: {}) | B_Feldman VSS | Party: Verify Shares",
                    k, n, t
                ),
                |b| {
                    b.iter(|| {
                        assert!(p.verify_shares().unwrap());
                    })
                },
            );
        }
    }
}

criterion_group!(benches, vss);
criterion_main!(benches);
