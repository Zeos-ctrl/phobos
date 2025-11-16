#[derive(Debug, Clone, Copy)]
pub struct Complex {
    pub real: f64,
    pub imag: f64,
}

impl Complex {
    pub fn new(real: f64, imag: f64) -> Self {
        Complex { real, imag }
    }

    pub fn add(&self, other: &Complex) -> Complex {
        Complex {
            real: self.real + other.real,
            imag: self.imag + other.imag
        }
    }

    pub fn subtract(&self, other: &Complex) -> Complex {
        Complex {
            real: self.real - other.real,
            imag: self.imag - other.imag
        }
    }

    pub fn multiply(&self, other: &Complex) -> Complex {
        Complex {
            real: (self.real * other.real) - (self.imag * other.imag),
            imag: (self.real * other.imag) + (self.imag * other.real)
        }
    }

    pub fn scale(&self, scalar: f64) -> Complex {
        Complex {
            real: self.real * scalar,
            imag: self.imag * scalar
        }
    }

    pub fn magnitude_squared(&self) -> f64 {
        self.real * self.real + self.imag * self.imag
    }

    pub fn conjugate(&self) -> Complex {
        Complex {
            real: self.real,
            imag: -self.imag
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_new() {
        let c = Complex::new(3.0, 4.0);
        assert_eq!(c.real, 3.0);
        assert_eq!(c.imag, 4.0);
    }

    #[test]
    fn test_addition() {
        let a = Complex::new(3.0, 4.0);
        let b = Complex::new(1.0, 2.0);
        let c = a.add(&b);
        assert_eq!(c.real, 4.0);
        assert_eq!(c.imag, 6.0);
    }

    #[test]
    fn test_multiplication() {
        let a = Complex::new(2.0, 3.0);
        let b = Complex::new(4.0, 5.0);
        let c = a.multiply(&b);
        assert_eq!(c.real, -7.0);
        assert_eq!(c.imag, 22.0);
    }

    #[test]
    fn test_scale() {
        let a = Complex::new(2.0, 3.0);
        let b = a.scale(2.0);
        assert_eq!(b.real, 4.0);
        assert_eq!(b.imag, 6.0);
    }

    #[test]
    fn test_magnitude_squared() {
        let a = Complex::new(3.0, 4.0);
        let b = a.magnitude_squared();
        assert_eq!(b, 25.0);
    }

    #[test]
    fn test_conjugate() {
        let a = Complex::new(3.0, 4.0);
        let b = a.conjugate();
        assert_eq!(b.real, 3.0);
        assert_eq!(b.imag, -4.0);
    }

}
