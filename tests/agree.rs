use linkcell::{knearest, knearest_brute, OrthoBox};

fn almost(a: f64, b: f64) -> bool {
    (a - b).abs() <= 1e-12 * (1.0 + a.abs().max(b.abs()))
}

#[test]
fn rejects_zero_k() {
    let b = OrthoBox::new(10.0, 10.0, 10.0).unwrap();
    let xyz = [[0.0, 0.0, 0.0]];
    assert!(knearest(&xyz, &b, 0, None, None).is_err());
}

#[test]
fn rejects_bad_box() {
    assert!(OrthoBox::new(0.0, 1.0, 1.0).is_err());
}

#[test]
fn two_points_k1() {
    let b = OrthoBox::new(10.0, 10.0, 10.0).unwrap();
    let xyz = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]];
    let rows = knearest(&xyz, &b, 1, None, Some(2.0)).unwrap();
    assert_eq!(rows[0].indices, vec![1]);
    assert_eq!(rows[1].indices, vec![0]);
    assert!(almost(rows[0].dist2[0], 1.0));
}

#[test]
fn periodic_image_is_nearer() {
    let b = OrthoBox::new(10.0, 10.0, 10.0).unwrap();
    let xyz = [[0.2, 0.0, 0.0], [9.4, 0.0, 0.0]];
    let rows = knearest(&xyz, &b, 1, None, Some(2.0)).unwrap();
    assert_eq!(rows[0].indices, vec![1]);
    // 0.2 + (10-9.4) = 0.8, not 9.2
    assert!(almost(rows[0].dist2[0], 0.8 * 0.8));
}

#[test]
fn mask_drops_sources_and_candidates() {
    let b = OrthoBox::new(10.0, 10.0, 10.0).unwrap();
    let xyz = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
    let mask = [true, false, true];
    let rows = knearest(&xyz, &b, 1, Some(&mask), None).unwrap();
    assert!(rows[1].indices.is_empty());
    assert_eq!(rows[0].indices, vec![2]);
    assert_eq!(rows[2].indices, vec![0]);
}

#[test]
fn agrees_with_brute_force_on_a_random_cell() {
    let b = OrthoBox::new(12.0, 11.0, 13.0).unwrap();
    let mut xyz = Vec::new();
    let mut s: u64 = 1;
    for _ in 0..64 {
        s = s.wrapping_mul(6364136223846793005).wrapping_add(1);
        let x = ((s >> 11) as f64 / (1u64 << 53) as f64) * 12.0;
        s = s.wrapping_mul(6364136223846793005).wrapping_add(1);
        let y = ((s >> 11) as f64 / (1u64 << 53) as f64) * 11.0;
        s = s.wrapping_mul(6364136223846793005).wrapping_add(1);
        let z = ((s >> 11) as f64 / (1u64 << 53) as f64) * 13.0;
        xyz.push([x, y, z]);
    }
    let cell = knearest(&xyz, &b, 4, None, Some(2.5)).unwrap();
    let brute = knearest_brute(&xyz, &b, 4, None).unwrap();
    for i in 0..xyz.len() {
        assert_eq!(cell[i].indices, brute[i].indices, "source {i}");
        for (a, c) in cell[i].dist2.iter().zip(brute[i].dist2.iter()) {
            assert!(almost(*a, *c), "{a} vs {c}");
        }
    }
}
