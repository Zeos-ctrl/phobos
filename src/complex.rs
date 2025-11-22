/// Complex number representation, these are the representation
/// of Qubits within the circuit.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Complex {
    /// Real part to the number
    pub real: f64,
    /// Imaginary part of the number
    pub imag: f64,
}


impl Complex {
    /// Creates a new complex number with the specified real and
    /// imaginary counterparts.
    ///
    /// # Arguments
    /// * `real` - The real part of the number
    /// * `imag` - The imaginary part of the number
    ///
    /// # Examples
    /// ```
    /// use phobos::Complex;
    ///
    /// let a = Complex::new(1.0, 0.0);
    /// ```
    pub fn new(real: f64, imag: f64) -> Self {
        Complex { real, imag }
    }

    /// Addition operation between two complex numbers.
    /// num1 = a + bi
    /// num2 = c + di 
    /// result = (a + c) + (b + d)i
    ///
    /// # Arguments
    /// * `other` - The complex number to be added
    ///
    /// # Examples
    /// ```
    /// use phobos::Complex;
    ///
    /// let a = Complex::new(1.0, 0.0);
    /// let b = Complex::new(0.0, 1.0);
    /// let c = a.add(&b);
    /// ```
    pub fn add(&self, other: &Complex) -> Complex {
        Complex {
            real: self.real + other.real,
            imag: self.imag + other.imag
        }
    }

    /// Subtraction operation between two complex numbers.
    /// num1 = a + bi
    /// num2 = c + di 
    /// result = (a - c) + (b - d)i
    ///
    /// # Arguments
    /// * `other` - The complex number to be added
    ///
    /// # Examples
    /// ```
    /// use phobos::Complex;
    ///
    /// let a = Complex::new(1.0, 0.0);
    /// let b = Complex::new(0.0, 1.0);
    /// let c = a.subtract(&b);
    /// ```
    pub fn subtract(&self, other: &Complex) -> Complex {
        Complex {
            real: self.real - other.real,
            imag: self.imag - other.imag
        }
    }

    /// Multiplication operation between two complex numbers.
    /// num1 = a + bi
    /// num2 = c + di 
    /// result = (ac - bd) + (ad + bc)i
    ///
    /// # Arguments
    /// * `other` - The complex number to be added
    ///
    /// # Examples
    /// ```
    /// use phobos::Complex;
    ///
    /// let a = Complex::new(1.0, 0.0);
    /// let b = Complex::new(0.0, 1.0);
    /// let c = a.multiply(&b);
    /// ```
    pub fn multiply(&self, other: &Complex) -> Complex {
        Complex {
            real: (self.real * other.real) - (self.imag * other.imag),
            imag: (self.real * other.imag) + (self.imag * other.real)
        }
    }

    /// Scale the complex number by a scalar value.
    /// num1 = a + bi
    /// scalar = x
    /// result = (x * a) + (x * b)i
    ///
    /// # Arguments
    /// * `scalar` - The scalar value
    ///
    /// # Examples
    /// ```
    /// use phobos::Complex;
    ///
    /// let a = Complex::new(2.0, 2.0);
    /// let b = a.scale(2.0);
    /// ```
    pub fn scale(&self, scalar: f64) -> Complex {
        Complex {
            real: self.real * scalar,
            imag: self.imag * scalar
        }
    }

    /// Calculate the magnitude squared of the complex number.
    /// Used to get the probability of the quantum state vector.
    /// num1 = a + bi
    /// result = |a|^2 + (|b|^2)i
    ///
    /// # Examples
    /// ```
    /// use phobos::Complex;
    ///
    /// let a = Complex::new(1.0, 0.0);
    /// let probability = a.magnitude_squared();
    /// ```
    pub fn magnitude_squared(&self) -> f64 {
        self.real * self.real + self.imag * self.imag
    }

    /// Calculate the complex conjugate (flip the sign on the imaginary part)
    /// num1 = a + bi
    /// result = a - bi
    ///
    /// # Examples
    /// ```
    /// use phobos::Complex;
    ///
    /// let a = Complex::new(1.0, 0.0);
    /// let b = a.conjugate();
    /// ```
    pub fn conjugate(&self) -> Complex {
        Complex {
            real: self.real,
            imag: -self.imag
        }
    }

    /// Compute the complex exponential e^(a+bi)
    /// Using Euler's formula: e^(a+bi) = e^a * (cos(b) + i*sin(b))
    pub fn exp(&self) -> Complex {
        let exp_real = self.real.exp();  // e^a
        Complex {
            real: exp_real * self.imag.cos(),   // e^a * cos(b)
            imag: exp_real * self.imag.sin()    // e^a * sin(b)
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
