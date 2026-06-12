use b_pi_s::{dealer::Dealer, party::generate_parties};
use criterion::{BatchSize, Criterion, criterion_group, criterion_main};
use curve25519_dalek::{RistrettoPoint, Scalar, ristretto::CompressedRistretto};

use blake3::Hasher;
use common::{
    BENCH_K, BENCH_N_T,
    precompute::gen_powers,
    random::{random_point, random_scalars},
    secret_sharing::generate_encrypted_shares_batched,
    utils::ingest_public_keys,
};
use rayon::prelude::*;

fn pvss(c: &mut Criterion) {
    // for (n, t) in BENCH_N_T {
    for (n, t) in [(16, 7)] {
        let mut rng = rand::rng();
        let mut hasher = Hasher::new();
        let mut buf: [u8; 64] = [0u8; 64];

        let g: RistrettoPoint = random_point(&mut rng);
        let xpows = gen_powers(n, t);

        let mut parties = generate_parties(&g, &mut rng, n, t);

        let public_keys: Vec<CompressedRistretto> =
            parties.iter().map(|party| party.public_key.0).collect();

        let mut dealer = Dealer::new(n, t, &public_keys).unwrap();

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

        for k in [1] {
            let secrets = random_scalars(&mut rng, k);

            let (f_polynomials, f_evals) = generate_encrypted_shares_batched(
                t,
                &xpows,
                &parties[0].public_keys.as_ref().unwrap(),
                &secrets,
            );

            c.bench_function(
                &format!(
                    "(k: {}, n: {}, t: {}) | B_Pi_S PVSS | Dealer: Generate Proof",
                    k, n, t
                ),
                |b| {
                    b.iter_batched(
                        || (blake3::Hasher::new(), [0u8; 64]),
                        |(mut hasher, mut buf)| {
                            dealer.generate_proof(
                                &mut rng,
                                &mut hasher,
                                &mut buf,
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

            let (shares, (c_vals, z)) =
                dealer.deal_secrets(&mut rng, &mut hasher, &mut buf, &xpows, &secrets);

            for p in &mut parties {
                p.ingest_dealer_proof((&c_vals, &z)).unwrap();

                p.ingest_encrypted_shares(&shares).unwrap();
            }

            c.bench_function(
                &format!(
                    "(k: {}, n: {}, t: {}) | B_Pi_S PVSS | Party: Verify Shares",
                    k, n, t
                ),
                |b| {
                    b.iter_batched(
                        || (blake3::Hasher::new(), [0u8; 64]),
                        |(mut hasher, mut buf)| {
                            assert!(
                                parties[0]
                                    .verify_encrypted_shares(&mut hasher, &mut buf, &xpows)
                                    .unwrap()
                            )
                        },
                        BatchSize::PerIteration,
                    )
                },
            );
            c.bench_function(
                &format!(
                    "(k: {}, n: {}, t: {}) | B_Pi_S PVSS | Party: Decrypt Shares",
                    k, n, t
                ),
                |b| b.iter(|| parties[0].decrypt_shares().unwrap()),
            );

            let (mut decrypted_shares, mut share_proofs): (
                Vec<Vec<CompressedRistretto>>,
                Vec<(Scalar, Scalar)>,
            ) = parties
                .iter_mut()
                .map(|p| {
                    p.decrypt_shares().unwrap();
                    p.dleq_share(&g, &mut rng, &mut hasher, &mut buf).unwrap();

                    (
                        p.decrypted_share
                            .clone()
                            .unwrap()
                            .par_iter()
                            .map(|ds| ds.compress())
                            .collect(),
                        p.share_proof.clone().unwrap(),
                    )
                })
                .collect();

            decrypted_shares.remove(parties[0].index - 1);
            share_proofs.remove(parties[0].index - 1);
            parties[0]
                .ingest_decrypted_shares_and_proofs(&decrypted_shares, share_proofs)
                .unwrap();

            c.bench_function(
                &format!(
                    "(k: {}, n: {}, t: {}) | B_Pi_S PVSS | Party: Verify Decrypted Shares",
                    k, n, t
                ),
                |b| b.iter(|| parties[0].verify_decrypted_shares().unwrap()),
            );
        }
    }
}

criterion_group!(benches, pvss);
criterion_main!(benches);
