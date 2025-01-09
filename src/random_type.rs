#[derive(Debug, Clone, Copy)]
pub enum RandomType {
    Uniform,
    Fixed(f64)
}

impl RandomType {
    pub fn random(&self) -> f64 {
        match self {
            RandomType::Uniform => rand::random::<f64>(),
            RandomType::Fixed(value) => *value
        }
    }
}