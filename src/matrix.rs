use crate::Complex;
use std::ops::Index;
use std::f64;

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Matrix2 {
    data: [[Complex; 2]; 2]  // 2D array
}

impl Matrix2 {
    pub fn new(m00: Complex, m01: Complex, m10: Complex, m11: Complex) -> Self {
        Matrix2 { 
            data: [
                [m00, m01],  // Row 0
                [m10, m11]   // Row 1
            ]
        }
    }

    pub fn zeros() -> Self {
        Matrix2 {
            data : [
                [Complex::new(0.0,0.0), Complex::new(0.0, 0.0)],
                [Complex::new(0.0,0.0), Complex::new(0.0, 0.0)]
            ]
        }
    } 

    pub fn multiply(&self, other: &Matrix2) -> Self {
        Matrix2 {
            data: [
                [
                    self[(0,0)].multiply(&other[(0,0)]).add(&self[(0,1)].multiply(&other[(1,0)])),
                    self[(0,0)].multiply(&other[(0,1)]).add(&self[(0,1)].multiply(&other[(1,1)])),
                ],
                [
                    self[(1,0)].multiply(&other[(0,0)]).add(&self[(1,1)].multiply(&other[(1,0)])),
                    self[(1,0)].multiply(&other[(0,1)]).add(&self[(1,1)].multiply(&other[(1,1)])),
                ]
            ]
        }
    }

    pub fn add(&self, other: &Matrix2) -> Self {
        Matrix2 {
            data: [
                [
                    self[(0,0)].add(&other[(0,0)]),
                    self[(0,1)].add(&other[(0,1)]),
                ],
                [
                    self[(1,0)].add(&other[(1,0)]),
                    self[(1,1)].add(&other[(1,1)]),
                ]
            ]
        }
    }

    pub fn scale(&self, scalar: &Complex) -> Self {
        Matrix2 {
            data: [
                [
                    self[(0,0)].multiply(scalar),
                    self[(0,1)].multiply(scalar),
                ],
                [
                    self[(1,0)].multiply(scalar),
                    self[(1,1)].multiply(scalar),
                ]
            ]
        }
    }

    pub fn conjugate_transpose(&self) -> Self {
        Matrix2 {
            data: [
                [
                    self[(0,0)].conjugate(),
                    self[(1,0)].conjugate(),
                ],
                [
                    self[(0,1)].conjugate(),
                    self[(1,1)].conjugate(),
                ]
            ]
        }
    }

    pub fn eigenvalues_hermitian(&self) -> (f64, f64) {
        let a = self[(0,0)].real;  // Diagonal elements are real
        let d = self[(1,1)].real;
        let b = self[(0,1)];       // This is a Complex number

        let modifier = (a - d).powi(2)/4.0 + b.magnitude_squared();
        
        let lamda1 = (a + d)/2.0 + modifier.sqrt();
        let lamda2 = (a + d)/2.0 - modifier.sqrt();

        (lamda1, lamda2)
    }

    pub fn eigenvectors_hermitian(&self, lambda1: f64, lambda2: f64) -> (Vector2, Vector2) {
        let a = self[(0,0)].real;
        let b = self[(0,1)];
        
        // Check if b is essentially zero (diagonal matrix)
        if b.magnitude_squared() < 1e-10 {
            // Special case: diagonal matrix
            // Eigenvectors are just the basis vectors
            let v1 = Vector2::new(Complex::new(0.0, 0.0), Complex::new(1.0, 0.0));
            let v2 = Vector2::new(Complex::new(1.0, 0.0), Complex::new(0.0, 0.0));
            (v1, v2)
        } else {
            // General case: use the formula
            let v1 = Vector2::new(b, Complex::new(lambda1 - a, 0.0)).normalize();
            let v2 = Vector2::new(b, Complex::new(lambda2 - a, 0.0)).normalize();
            (v1, v2)
        }
    }
}

impl Index<(usize, usize)> for Matrix2 {
    type Output = Complex;
    
    fn index(&self, idx: (usize, usize)) -> &Self::Output {
        let row = idx.0;
        let col = idx.1;

        &self.data[row][col]
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Vector2 {
    pub data: [Complex; 2]
}

impl Vector2 {
    pub fn new(v0: Complex, v1: Complex) -> Self {
        Vector2 { data: [v0, v1] }
    }

    pub fn normalize(&self) -> Self {
        let mag_squared = self.data[0].magnitude_squared() + self.data[1].magnitude_squared();
        let magnitude = mag_squared.sqrt();  // Get actual magnitude
        
        let v0_normalized = self.data[0].scale(1.0/magnitude);
        let v1_normalized = self.data[1].scale(1.0/magnitude);
        
        Vector2 { data: [v0_normalized, v1_normalized] }
    }

    pub fn outer_product(&self) -> Matrix2 {
        Matrix2::new(
            self.data[0].multiply(&self.data[0].conjugate()),
            self.data[0].multiply(&self.data[1].conjugate()),
            self.data[1].multiply(&self.data[0].conjugate()),
            self.data[1].multiply(&self.data[1].conjugate()),
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_new() {
        let m00 = Complex::new(1.0, 0.0);
        let m01 = Complex::new(2.0, 0.0);
        let m10 = Complex::new(3.0, 0.0);
        let m11 = Complex::new(4.0, 0.0);
        let m = Matrix2::new(m00, m01, m10, m11);

        assert_eq!(m.data[0][0].real, 1.0);
        assert_eq!(m.data[0][0].imag, 0.0);
    }

    #[test]
    fn test_multiply() {
        // Matrix A = [[1, 2], [3, 4]]
        let a = Matrix2::new(
            Complex::new(1.0, 0.0), Complex::new(2.0, 0.0),
            Complex::new(3.0, 0.0), Complex::new(4.0, 0.0)
        );
        
        // Matrix B = [[5, 6], [7, 8]]
        let b = Matrix2::new(
            Complex::new(5.0, 0.0), Complex::new(6.0, 0.0),
            Complex::new(7.0, 0.0), Complex::new(8.0, 0.0)
        );
        
        let c = a.multiply(&b);

        println!("A[(0,0)] = {}, A[(0,1)] = {}", a[(0,0)].real, a[(0,1)].real);
        println!("B[(0,0)] = {}, B[(1,0)] = {}", b[(0,0)].real, b[(1,0)].real);
        
        assert_eq!(c[(0,0)].real, 19.0);
        assert_eq!(c[(0,1)].real, 22.0);
        assert_eq!(c[(1,0)].real, 43.0);
        assert_eq!(c[(1,1)].real, 50.0);
    }

    
    // Matrix2 Tests
    
    #[test]
    fn test_zeros() {
        let m = Matrix2::zeros();

        assert_eq!(m[(0,0)], Complex::new(0.0, 0.0));
        assert_eq!(m[(0,1)], Complex::new(0.0, 0.0));
        assert_eq!(m[(1,0)], Complex::new(0.0, 0.0));
        assert_eq!(m[(1,1)], Complex::new(0.0, 0.0));
    }
    
    #[test]
    fn test_add() {
        // Matrix A = [[1+2i, 3+4i], [5+6i, 7+8i]]
        let a = Matrix2::new(
            Complex::new(1.0, 2.0), Complex::new(3.0, 4.0),
            Complex::new(5.0, 6.0), Complex::new(7.0, 8.0)
        );
        // Matrix B = [[1, 1], [1, 1]] (all ones)
        let b = Matrix2::new(
            Complex::new(1.0, 0.0), Complex::new(1.0, 0.0),
            Complex::new(1.0, 0.0), Complex::new(1.0, 0.0)
        );
        
        let c = a.add(&b);

        assert_eq!(c[(0,0)], Complex::new(2.0, 2.0));
        assert_eq!(c[(0,1)], Complex::new(4.0, 4.0));
        assert_eq!(c[(1,0)], Complex::new(6.0, 6.0));
        assert_eq!(c[(1,1)], Complex::new(8.0, 8.0));
    }
    
    #[test]
    fn test_scale() {
        // Matrix M = [[1, 2], [3, 4]]
        let m = Matrix2::new(
            Complex::new(1.0, 0.0), Complex::new(2.0, 0.0),
            Complex::new(3.0, 0.0), Complex::new(4.0, 0.0)
        );
        // Scalar s = 2+0i
        let scalar = Complex::new(2.0, 0.0);
        
        let result = m.scale(&scalar);

        assert_eq!(result[(0,0)], Complex::new(2.0, 0.0));
        assert_eq!(result[(0,1)], Complex::new(4.0, 0.0));
        assert_eq!(result[(1,0)], Complex::new(6.0, 0.0));
        assert_eq!(result[(1,1)], Complex::new(8.0, 0.0));
    }
    
    #[test]
    fn test_conjugate_transpose() {
        // Matrix M = [[1+2i, 3+4i], [5+6i, 7+8i]]
        let m = Matrix2::new(
            Complex::new(1.0, 2.0), Complex::new(3.0, 4.0),
            Complex::new(5.0, 6.0), Complex::new(7.0, 8.0)
        );
        
        let m_dagger = m.conjugate_transpose();
        
        assert_eq!(m_dagger[(0,0)], Complex::new(1.0, -2.0));
        assert_eq!(m_dagger[(0,1)], Complex::new(5.0, -6.0));
        assert_eq!(m_dagger[(1,0)], Complex::new(3.0, -4.0));
        assert_eq!(m_dagger[(1,1)], Complex::new(7.0, -8.0));
    }

    // Vector2 Tests
    
    #[test]
    fn test_normalize() {
        // Vector v = [3, 4] (real components)
        let v = Vector2::new(
            Complex::new(3.0, 0.0),
            Complex::new(4.0, 0.0)
        );
        
        let v_normalized = v.normalize();

        assert!((v_normalized.data[0].real - 0.6).abs() < 1e-10);
        assert!((v_normalized.data[1].real - 0.8).abs() < 1e-10);
        
        let mag_sq = v_normalized.data[0].magnitude_squared() + 
                     v_normalized.data[1].magnitude_squared();

        assert!((mag_sq - 1.0).abs() < 1e-10);
    }
    
    #[test]
    fn test_outer_product() {
        // Vector v = [1, 0]
        let v = Vector2::new(
            Complex::new(1.0, 0.0),
            Complex::new(0.0, 0.0)
        );
        
        let proj = v.outer_product();
        
        assert_eq!(proj[(0,0)].real, 1.0);
        assert_eq!(proj[(0,1)].real, 0.0);
        assert_eq!(proj[(1,0)].real, 0.0);
        assert_eq!(proj[(1,1)].real, 0.0);
    }

    // Eigendecomposition Tests
    
    #[test]
    fn test_eigenvalues_diagonal() {
        // Diagonal matrix H = [[2, 0], [0, 5]]
        let h = Matrix2::new(
            Complex::new(2.0, 0.0), Complex::new(0.0, 0.0),
            Complex::new(0.0, 0.0), Complex::new(5.0, 0.0)
        );
        
        let (lambda1, lambda2) = h.eigenvalues_hermitian();
        
        assert!((lambda1 - 5.0).abs() < 1e-10);
        assert!((lambda2 - 2.0).abs() < 1e-10);
    }
    
    #[test]
    fn test_pauli_z_eigenvalues() {
        // Pauli-Z: H = [[1, 0], [0, -1]]
        let h = Matrix2::new(
            Complex::new(1.0, 0.0), Complex::new(0.0, 0.0),
            Complex::new(0.0, 0.0), Complex::new(-1.0, 0.0)
        );
        
        let (lambda1, lambda2) = h.eigenvalues_hermitian();
        
        assert!((lambda1 - 1.0).abs() < 1e-10);
        assert!((lambda2 + 1.0).abs() < 1e-10);
    }
    
    #[test]
    fn test_eigenvectors_diagonal() {
        let h = Matrix2::new(
            Complex::new(2.0, 0.0), Complex::new(0.0, 0.0),
            Complex::new(0.0, 0.0), Complex::new(5.0, 0.0)
        );
        
        let (lambda1, lambda2) = h.eigenvalues_hermitian();
        let (v1, v2) = h.eigenvectors_hermitian(lambda1, lambda2);

        println!("lambda1 = {}, lambda2 = {}", lambda1, lambda2);
        println!("v1 = [{}, {}]", v1.data[0].real, v1.data[1].real);
        println!("v2 = [{}, {}]", v2.data[0].real, v2.data[1].real);
        
        assert!((v1.data[0].real - 0.0).abs() < 1e-10);
        assert!((v1.data[1].real - 1.0).abs() < 1e-10);
        assert!((v2.data[0].real - 1.0).abs() < 1e-10);
        assert!((v2.data[1].real - 0.0).abs() < 1e-10);
    }

    #[test]
    fn test_eigenvalues() {
        // H = [3  1]
        //     [1  3]
        let h = Matrix2::new(
            Complex::new(3.0, 0.0), Complex::new(1.0, 0.0),
            Complex::new(1.0, 0.0), Complex::new(3.0, 0.0)
        );
        
        let (lambda1, lambda2) = h.eigenvalues_hermitian();
        
        // From our hand calculation: λ₁ = 4, λ₂ = 2
        assert!((lambda1 - 4.0).abs() < 1e-10);
        assert!((lambda2 - 2.0).abs() < 1e-10);
    }

    #[test]
    fn test_eigenvectors() {
        // H = [3  1]
        //     [1  3]
        let h = Matrix2::new(
            Complex::new(3.0, 0.0), Complex::new(1.0, 0.0),
            Complex::new(1.0, 0.0), Complex::new(3.0, 0.0)
        );
        
        let (lambda1, lambda2) = h.eigenvalues_hermitian();
        let (v1, _) = h.eigenvectors_hermitian(lambda1, lambda2);
        
        // v1 should be approximately [1/√2, 1/√2] for λ₁ = 4
        // v2 should be approximately [1/√2, -1/√2] for λ₂ = 2
        
        let sqrt2 = 2.0_f64.sqrt();
        assert!((v1.data[0].real - 1.0/sqrt2).abs() < 1e-10);
        assert!((v1.data[1].real - 1.0/sqrt2).abs() < 1e-10);
    }

}
