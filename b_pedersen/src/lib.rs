pub mod dealer;
pub mod party;

#[cfg(test)]

mod tests {
    use curve25519_dalek::{RistrettoPoint, ristretto::CompressedRistretto};

    use crate::{dealer::Dealer, party::generate_parties};

    use common::{
        precompute::gen_powers,
        random::{random_point, random_points, random_scalars},
        secret_sharing::{reconstruct_secrets, select_qualified_set},
        utils::{compute_lagrange_bases, ingest_public_keys},
    };

    #[test]
    fn end_to_end() {
        const N: usize = 1024;
        const T: usize = 511;
        const K: usize = 16;

        let mut rng = rand::rng();

        let generator: RistrettoPoint = random_point(&mut rng);
        let g: Vec<RistrettoPoint> = random_points(&mut rng, K);

        let g0: RistrettoPoint = random_point(&mut rng);

        let xpows = gen_powers(N, T);

        let mut parties = generate_parties(&generator, &g, &g0, &mut rng, N, T);

        let public_keys: Vec<CompressedRistretto> =
            parties.iter().map(|party| party.public_key.0).collect();

        let mut dealer = Dealer::new(g, g0, N, T, &public_keys).unwrap();

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

        let (shares, (r_evals, c_vals)) = dealer.deal_secret(&mut rng, &xpows, &secrets);

        for p in &mut parties {
            p.ingest_dealer_proof(&c_vals).unwrap();

            p.ingest_share((&shares[p.index - 1], &r_evals[p.index - 1]));
            assert!(
                p.verify_share().unwrap(),
                "individual share verification failure"
            );
            println!("pass own share: {}", p.index);

            p.ingest_shares((&shares, &r_evals)).unwrap();

            assert!(p.verify_shares().unwrap(), "share verification failure");

            let party_shares = p
                .shares
                .as_ref()
                .unwrap()
                .iter()
                .map(|(p_share, _)| p_share.clone())
                .collect();

            p.qualified_set = Some(
                select_qualified_set(&mut rng, p.t, &Some(party_shares), &p.validated_shares)
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

            let sec = reconstruct_secrets(&p.qualified_set, &lagrange_bases).unwrap();

            assert!(secrets == sec, "Invalid Reconstructed Secret");
        }
    }
}
