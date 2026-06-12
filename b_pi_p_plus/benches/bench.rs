use b_pi_p_plus::{dealer::Dealer, party::generate_parties};
use criterion::{BatchSize, Criterion, criterion_group, criterion_main};
use curve25519_dalek::{RistrettoPoint, ristretto::CompressedRistretto, traits::Identity};

use common::{
    BENCH_K, BENCH_N_T,
    precompute::gen_powers,
    random::{random_point, random_scalars},
    secret_sharing::generate_shares_batched,
    utils::ingest_public_keys,
};

fn pvss(c: &mut Criterion) {
    for (n, t) in BENCH_N_T {
        // for (n, t) in [(256, 127)] {
        let mut rng = rand::rng();
        let mut hasher = blake3::Hasher::new();
        let mut buf: [u8; 64] = [0u8; 64];

        let g: RistrettoPoint = random_point(&mut rng);
        let g1: RistrettoPoint = random_point(&mut rng);
        let g2: RistrettoPoint = random_point(&mut rng);

        let xpows = gen_powers(n, t);

        let mut parties = generate_parties(&g, &g1, &g2, &mut rng, n, t);

        let public_keys: Vec<CompressedRistretto> =
            parties.iter().map(|party| party.public_key.0).collect();

        let mut dealer = Dealer::new(g1, g2, n, t, &public_keys).unwrap();

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

        // for k in BENCH_K {
        for k in [10000] {
            let secrets = random_scalars(&mut rng, k);

            let (f_polynomials, f_evals) = generate_shares_batched(n, t, &xpows, &secrets);

            c.bench_function(
                &format!(
                    "(k: {}, n: {}, t: {}) | B_Pi_P+ VSS | Dealer: Generate Proof",
                    k, n, t
                ),
                |b| {
                    b.iter_batched(
                        || {
                            (
                                blake3::Hasher::new(),
                                [0u8; 64],
                                vec![CompressedRistretto::identity(); dealer.public_keys.len()],
                            )
                        },
                        |(mut hasher, mut buf, mut c_buf)| {
                            dealer.generate_proof(
                                &mut rng,
                                &mut hasher,
                                &mut buf,
                                &mut c_buf,
                                &xpows,
                                k,
                                &f_polynomials,
                                &f_evals,
                            )
                        },
                        BatchSize::PerIteration,
                    )
                },
            );

            let (shares, (g, c_vals, z)) =
                dealer.deal_secret(&mut rng, &mut hasher, &mut buf, &xpows, &secrets);

            let p = &mut parties[0];

            p.ingest_shares((&shares, &g)).unwrap();
            p.ingest_dealer_proof((&c_vals, &z)).unwrap();

            assert!(
                p.verify_shares(&mut hasher, &mut buf, &xpows).unwrap(),
                "share verification failure"
            );

            c.bench_function(
                &format!(
                    "(k: {}, n: {}, t: {}) | B_Pi_P+ VSS | Party: Verify Shares",
                    k, n, t
                ),
                |b| {
                    b.iter_batched(
                        || (blake3::Hasher::new(), [0u8; 64]),
                        |(mut hasher, mut buf)| {
                            assert!(p.verify_shares(&mut hasher, &mut buf, &xpows).unwrap())
                        },
                        BatchSize::PerIteration,
                    )
                },
            );
        }
    }
}

criterion_group!(benches, pvss);
criterion_main!(benches);
