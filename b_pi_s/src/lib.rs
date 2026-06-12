pub mod dealer;
pub mod party;

#[cfg(test)]

mod tests {
    use curve25519_dalek::{RistrettoPoint, Scalar, ristretto::CompressedRistretto};

    use crate::{dealer::Dealer, party::generate_parties};

    use common::{
        precompute::gen_powers,
        random::{random_point, random_scalars},
        secret_sharing::{reconstruct_secrets_exponent, select_qualified_set},
        utils::{compute_lagrange_bases, ingest_public_keys},
    };
    use rayon::prelude::*;

    #[test]
    fn end_to_end() {
        const N: usize = 16;
        const T: usize = 7;
        const K: usize = 3;

        let mut rng = rand::rng();
        let mut hasher = blake3::Hasher::new();
        let mut buf = [0u8; 64];

        let g: RistrettoPoint = random_point(&mut rng);

        // let xpows = XPowTable::from_params("../table.json", N, T);
        let xpows = gen_powers(N, T);

        let mut parties = generate_parties(&g, &mut rng, N, T);

        let public_keys: Vec<CompressedRistretto> =
            parties.iter().map(|party| party.public_key.0).collect();

        let mut dealer = Dealer::new(N, T, &public_keys).unwrap();

        for party in &mut parties {
            let public_keys: Vec<CompressedRistretto> = public_keys
                .iter()
                .filter(|pk| &party.public_key.0 != *pk)
                .copied()
                .collect();

            party.public_keys = Some(
                ingest_public_keys(N, &party.public_key.1, party.index, &public_keys).unwrap(),
            );
        }

        let secrets = random_scalars(&mut rng, K);

        let (shares, (c_vals, z)) =
            dealer.deal_secrets(&mut rng, &mut hasher, &mut buf, &xpows, &secrets);

        for p in &mut parties {
            p.ingest_dealer_proof((&c_vals, &z)).unwrap();

            p.ingest_encrypted_shares(&shares).unwrap();

            let verif_result = p
                .verify_encrypted_shares(&mut hasher, &mut buf, &xpows)
                .unwrap();

            assert!(verif_result, "share verification failure");
        }

        let (decrypted_shares, share_proofs): (
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

        for p in &mut parties {
            let (mut decrypted_shares, mut share_proofs) =
                (decrypted_shares.clone(), share_proofs.clone());

            decrypted_shares.remove(p.index - 1);
            share_proofs.remove(p.index - 1);
            p.ingest_decrypted_shares_and_proofs(&decrypted_shares, share_proofs)
                .unwrap();

            assert!(p.verify_decrypted_shares().unwrap());

            p.qualified_set = Some(
                select_qualified_set(&mut rng, p.t, &p.decrypted_shares, &p.validated_shares)
                    .unwrap(),
            );

            let indices: Vec<usize> = p
                .qualified_set
                .as_ref()
                .unwrap()
                .iter()
                .map(|(index, _)| *index)
                .collect();

            let lagrange_bases = compute_lagrange_bases(&indices);

            let sec = reconstruct_secrets_exponent(&p.qualified_set, &lagrange_bases).unwrap();

            sec.iter()
                .zip(secrets.iter())
                .for_each(|(secret, dealer_secret)| assert_eq!(g * dealer_secret, *secret));
        }
    }
}
