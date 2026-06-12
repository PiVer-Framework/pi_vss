use b_pedersen::{dealer::Dealer, party::generate_parties};
use criterion::{BatchSize, Criterion, criterion_group, criterion_main};
use curve25519_dalek::{RistrettoPoint, ristretto::CompressedRistretto, traits::Identity};

use common::{
    BENCH_K, BENCH_N_T, Q,
    precompute::gen_powers,
    random::{random_point, random_points, random_scalars},
    secret_sharing::generate_shares_batched,
    utils::ingest_public_keys,
};

fn vss(c: &mut Criterion) {
    for (n, t) in [(1024, 511)] {
        // for (n, t) in BENCH_N_T {
        let mut rng = rand::rng();

        let generator: RistrettoPoint = random_point(&mut rng);
        let g2: RistrettoPoint = random_point(&mut rng);

        let xpows = gen_powers(n, t);
        for k in 1..=16 {
            let g: Vec<RistrettoPoint> = random_points(&mut rng, k);
            let mut parties = generate_parties(&generator, &g, &g2, &mut rng, n, t);

            let public_keys: Vec<CompressedRistretto> =
                parties.iter().map(|party| party.public_key.0).collect();

            let mut dealer = Dealer::new(g, g2, n, t, &public_keys).unwrap();

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
                    "(k: {}, n: {}, t: {}) | B_Pedersen VSS | Dealer: Generate Proof",
                    k, n, t
                ),
                |b| {
                    b.iter_batched(
                        || vec![CompressedRistretto::identity(); dealer.public_keys.len()],
                        |mut c_buf| {
                            dealer.generate_proof(&mut rng, &mut c_buf, &xpows, &f_polynomials)
                        },
                        BatchSize::PerIteration,
                    )
                },
            );

            let (shares, (r_evals, c_vals)) = dealer.deal_secret(&mut rng, &xpows, &secrets);

            let p = &mut parties[0];
            p.ingest_dealer_proof(&c_vals).unwrap();

            p.ingest_shares((&shares, &r_evals)).unwrap();

            c.bench_function(
                &format!(
                    "(k: {}, n: {}, t: {}) | B_Pedersen VSS | Party: Verify Shares",
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
