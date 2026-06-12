pub mod dealer;
pub mod party;

#[cfg(test)]

mod tests {
    use curve25519_dalek::{RistrettoPoint, ristretto::CompressedRistretto};

    use crate::{dealer::Dealer, party::generate_parties};

    use common::{
        precompute::gen_powers,
        random::{random_point, random_scalars},
        secret_sharing::{reconstruct_secrets, select_qualified_set},
        utils::{compute_lagrange_bases, ingest_public_keys},
    };

    #[test]
    fn end_to_end() {
        const N: usize = 128;
        const T: usize = 63;
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

            p.ingest_share(&shares[p.index - 1]);

            assert!(
                p.verify_share(&mut hasher, &mut buf, &xpows).unwrap(),
                "share verification failure"
            );

            p.ingest_shares(&shares).unwrap();

            let verif_result = p.verify_shares(&mut hasher, &mut buf, &xpows).unwrap();

            assert!(verif_result, "share verification failure");

            p.qualified_set =
                Some(select_qualified_set(&mut rng, p.t, &p.shares, &p.validated_shares).unwrap());

            let indices: Vec<usize> = p
                .qualified_set
                .as_ref()
                .unwrap()
                .iter()
                .map(|(index, _)| *index)
                .collect();

            let lagrange_bases = compute_lagrange_bases(&indices);

            let sec = reconstruct_secrets(&p.qualified_set, &lagrange_bases).unwrap();

            assert!(secrets == sec, "Invalid Reconstructed Secret");
        }
    }
    #[test]
    fn beaver_mult() {
        const N: usize = 8;
        const T: usize = 3;
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

        // a, b
        let mut secrets = random_scalars(&mut rng, 2);
        // c = a * b
        secrets.push(secrets[0] * secrets[1]);

        let (shares, (c_vals, z)) =
            dealer.deal_secrets(&mut rng, &mut hasher, &mut buf, &xpows, &secrets);

        for p in &mut parties {
            p.ingest_dealer_proof((&c_vals, &z)).unwrap();

            p.ingest_share(&shares[p.index - 1]);

            assert!(
                p.verify_share(&mut hasher, &mut buf, &xpows).unwrap(),
                "share verification failure"
            );

            p.ingest_shares(&shares).unwrap();

            let verif_result = p.verify_shares(&mut hasher, &mut buf, &xpows).unwrap();

            assert!(verif_result, "share verification failure");

            p.qualified_set =
                Some(select_qualified_set(&mut rng, p.t, &p.shares, &p.validated_shares).unwrap());

            let indices: Vec<usize> = p
                .qualified_set
                .as_ref()
                .unwrap()
                .iter()
                .map(|(index, _)| *index)
                .collect();

            let lagrange_bases = compute_lagrange_bases(&indices);

            let sec = reconstruct_secrets(&p.qualified_set, &lagrange_bases).unwrap();

            assert!(secrets == sec, "Invalid Reconstructed Secret");
        }
    }
}
