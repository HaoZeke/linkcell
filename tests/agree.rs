use linkcell::{knearest, knearest_brute, Cell};

fn almost(a: f64, b: f64) -> bool {
    (a - b).abs() <= 1e-12 * (1.0 + a.abs().max(b.abs()))
}

#[test]
fn rejects_zero_k() {
    let b = Cell::ortho(10.0, 10.0, 10.0).unwrap();
    let xyz = [[0.0, 0.0, 0.0]];
    assert!(knearest(&xyz, &b, 0, None, None).is_err());
}

#[test]
fn rejects_bad_box() {
    assert!(Cell::ortho(0.0, 1.0, 1.0).is_err());
    assert!(Cell::from_vectors(
        [1.0, 0.0, 0.0],
        [2.0, 0.0, 0.0],
        [0.0, 0.0, 1.0],
        [0.0, 0.0, 0.0]
    )
    .is_err());
}

#[test]
fn two_points_k1() {
    let b = Cell::ortho(10.0, 10.0, 10.0).unwrap();
    let xyz = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]];
    let rows = knearest(&xyz, &b, 1, None, Some(2.0)).unwrap();
    assert_eq!(rows[0].indices, vec![1]);
    assert_eq!(rows[1].indices, vec![0]);
    assert!(almost(rows[0].dist2[0], 1.0));
}

#[test]
fn periodic_image_is_nearer() {
    let b = Cell::ortho(10.0, 10.0, 10.0).unwrap();
    let xyz = [[0.2, 0.0, 0.0], [9.4, 0.0, 0.0]];
    let rows = knearest(&xyz, &b, 1, None, Some(2.0)).unwrap();
    assert_eq!(rows[0].indices, vec![1]);
    // 0.2 + (10-9.4) = 0.8, not 9.2
    assert!(almost(rows[0].dist2[0], 0.8 * 0.8));
}

#[test]
fn mask_drops_sources_and_candidates() {
    let b = Cell::ortho(10.0, 10.0, 10.0).unwrap();
    let xyz = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
    let mask = [true, false, true];
    let rows = knearest(&xyz, &b, 1, Some(&mask), None).unwrap();
    assert!(rows[1].indices.is_empty());
    assert_eq!(rows[0].indices, vec![2]);
    assert_eq!(rows[2].indices, vec![0]);
}

#[test]
fn agrees_with_brute_force_on_a_random_cell() {
    let b = Cell::ortho(12.0, 11.0, 13.0).unwrap();
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

fn hex_prism() -> Cell {
    // a = (10,0,0), b = (5, 5*sqrt(3), 0), c = (0,0,10)
    Cell::from_vectors(
        [10.0, 0.0, 0.0],
        [5.0, 8.660254037844386, 0.0],
        [0.0, 0.0, 10.0],
        [0.0, 0.0, 0.0],
    )
    .unwrap()
}

#[test]
fn triclinic_image_beats_the_cartesian_far_point() {
    let b = hex_prism();
    // Near a and near a+b, which are 10 apart along a, but the
    // minimum image across -b + a is shorter than the raw vector.
    let xyz = [[0.2, 0.1, 1.0], [9.7, 0.1, 1.0]];
    let rows = knearest(&xyz, &b, 1, None, Some(2.0)).unwrap();
    assert_eq!(rows[0].indices, vec![1]);
    let brute = knearest_brute(&xyz, &b, 1, None).unwrap();
    assert_eq!(rows[0].indices, brute[0].indices);
    assert!(almost(rows[0].dist2[0], brute[0].dist2[0]));
    assert!(rows[0].dist2[0] < 1.0);
}

#[test]
fn agrees_with_brute_force_on_a_sheared_cell() {
    let b = hex_prism();
    let mut xyz = Vec::new();
    let mut s: u64 = 7;
    for _ in 0..48 {
        s = s.wrapping_mul(6364136223846793005).wrapping_add(1);
        let u = (s >> 11) as f64 / (1u64 << 53) as f64;
        s = s.wrapping_mul(6364136223846793005).wrapping_add(1);
        let v = (s >> 11) as f64 / (1u64 << 53) as f64;
        s = s.wrapping_mul(6364136223846793005).wrapping_add(1);
        let w = (s >> 11) as f64 / (1u64 << 53) as f64;
        let r = [
            10.0 * u + 5.0 * v,
            8.660254037844386 * v,
            10.0 * w,
        ];
        xyz.push(r);
    }
    let cell = knearest(&xyz, &b, 4, None, Some(2.0)).unwrap();
    let brute = knearest_brute(&xyz, &b, 4, None).unwrap();
    for i in 0..xyz.len() {
        assert_eq!(cell[i].indices, brute[i].indices, "source {i}");
        for (a, c) in cell[i].dist2.iter().zip(brute[i].dist2.iter()) {
            assert!(almost(*a, *c), "{a} vs {c} at {i}");
        }
    }
}
