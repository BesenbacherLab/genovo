//! Reproducible random numbers.
//!
//! Every (purpose, region, mutation type) gets its own random number stream,
//! derived from the `--seed` value. Results therefore depend only on the seed,
//! not on the order in which regions are processed (which varies between runs
//! because they are kept in hash maps), and running a single region with
//! `--id` gives the same numbers for that region as a run over all regions.

use rand::SeedableRng;
use rand_chacha::ChaCha8Rng;

use crate::MutationType;

/// 64-bit FNV-1a hash: fixed, unlike the standard library's hasher.
fn fnv1a(mut hash: u64, bytes: &[u8]) -> u64 {
    for &b in bytes {
        hash ^= u64::from(b);
        hash = hash.wrapping_mul(0x0000_0100_0000_01b3);
    }
    hash
}

/// The random number stream for `purpose` (e.g. "sample") in one region and
/// mutation type.
pub fn stream_rng(seed: u64, purpose: &str, region: &str, mutation_type: MutationType) -> ChaCha8Rng {
    let mut h = 0xcbf2_9ce4_8422_2325;
    h = fnv1a(h, purpose.as_bytes());
    h = fnv1a(h, &[0]);
    h = fnv1a(h, region.as_bytes());
    h = fnv1a(h, &[0, mutation_type as u8]);
    let mut key = [0u8; 32];
    key[..8].copy_from_slice(&seed.to_le_bytes());
    key[8..16].copy_from_slice(&h.to_le_bytes());
    ChaCha8Rng::from_seed(key)
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::Rng;

    const PINNED: [u32; 4] = [1903280384, 2051161785, 2927687284, 2031500480];

    #[test]
    fn streams_are_reproducible_and_distinct() {
        let draw = |seed, purpose, region, t| -> Vec<u32> {
            let mut r = stream_rng(seed, purpose, region, t);
            (0..4).map(|_| r.gen()).collect()
        };
        let a = draw(1, "sample", "ENST1", MutationType::Missense);
        assert_eq!(a, draw(1, "sample", "ENST1", MutationType::Missense));
        assert_ne!(a, draw(2, "sample", "ENST1", MutationType::Missense));
        assert_ne!(a, draw(1, "lof", "ENST1", MutationType::Missense));
        assert_ne!(a, draw(1, "sample", "ENST2", MutationType::Missense));
        assert_ne!(a, draw(1, "sample", "ENST1", MutationType::Nonsense));
        // The streams must not change between versions (old results must stay
        // reproducible).
        assert_eq!(a, PINNED);
    }
}
