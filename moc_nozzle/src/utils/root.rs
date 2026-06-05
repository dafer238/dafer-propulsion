pub fn bisection<F: Fn(f64) -> f64>(f: F, mut a: f64, mut b: f64) -> f64 {
    for _ in 0..50 {
        let mid = 0.5 * (a + b);
        if f(mid) * f(a) < 0.0 {
            b = mid;
        } else {
            a = mid;
        }
    }
    0.5 * (a + b)
}
