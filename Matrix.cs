//Keaton Harlow
using System;
namespace MathLibrary {
    public class Matrix {
        // Creating the matrix
        private double[,] m, p, l, u;
        private static Random rnd = new Random();
        private int pswap;
        private bool LU = false;
        
        // Creates a matrix(1, 1)
        public Matrix() {
            m = new double[1, 1];
        }
        
        // Creates a matrix(r, c)
        public Matrix(int r, int c) {
            m = new double[r, c];
        }
        
        // Fill the matrix with random numbers
        public void setRand(double low, double high) {
            for(int i = 0; i < m.GetLength(0); i++) {
                for(int j = 0; j < m.GetLength(1); j++) {
                    m[i, j] = rnd.NextDouble() * (high - low) + low;
                }
            }
        }
        
        // Fill the matrix with an identity matrix
        public void identMat() {
            int size = (m.GetLength(0) > m.GetLength(1)) ? m.GetLength(0) : m.GetLength(1);
            for(int i = 0; i < size; i++) m[i, i] = 1.0;
        }
        
        // Returns the number of rows in the matrix
        public int getRow() {
            return m.GetLength(0);
        }
        
        // Returns the number of columns in the matrix
        public int getCol() {
            return m.GetLength(1);
        }
        
        // Sets a certain element to the value sent.
        // (value, row value, column value)
        public void setElem(double n, int r, int c) {
            m[r, c] = n;
        }
        
        // Set the matrix equal to the received matrix
        public int setMat(double[,] mat) {
            if(mat.GetLength(0) == m.GetLength(0) && mat.GetLength(1) == m.GetLength(1)) {
                m = mat;
                return 1;
            }
            
            else return 0;
        }
        
        // Returns the element from the matrix at location (r, c)
        public double getElem(int r, int c) {
            return m[r, c];
        }
        
        // Returns the matrix
        public double[,] getMat() {
            return m;
        }
        
        // Print the matrix
        public void print() {
            for(int i = 0; i < m.GetLength(0); i++) {
                for(int j = 0; j < m.GetLength(1); j++) {
                    Console.Write("{0:F2} ", m[i, j]);
                }
                Console.WriteLine();
            }
        }
        
        // For internal use only, troubleshooting
        private void print(double[, ] mat) {
            for(int i = 0; i < mat.GetLength(0); i++) {
                for(int j = 0; j < mat.GetLength(1); j++) {
                    Console.Write("{0:F2} ", mat[i, j]);
                }
                Console.WriteLine();
            }
        }
        
        // LU decomposition
        public void LUDecomp() {
            int r = m.GetLength(0), c = m.GetLength(1);
            if(r == c) {
                LU = true;
                pswap = 0;
                int imax = 0;
                double max = 0.0, tmpd = 0.0;
                p = new double[r, c];
                l = new double[r, c];
                u = new double[r, c];
                
                // Create a temporary matrix to manipulate
                double[,] tmpm = new Double[r, c];
                for(int i = 0; i < r; i++)
                    for(int j = 0; j < c; j++)
                        tmpm[i, j] = m[i, j];
                
                // Set p and l to be an identity matrix
                for(int i = 0; i < r; i++) {
                    p[i, i] = 1;
                    l[i, i] = 1;
                }
                
                // Row swap the matrix
                for(int j = 0; j < c - 1; j++) {
                    // Find the row with the highest value
                    max = Math.Abs(tmpm[j, j]);
                    for(int i = j; i < r; i++) {
                        if(tmpm[i, j] == 1) {
                            imax = i;
                            break;
                        }
                        if(Math.Abs(tmpm[i, j]) > max) {
                            max = Math.Abs(tmpm[i, j]);
                            imax = i;
                        }
                    }
                    if(imax != j) {
                        // Row swap tmpm and p
                        for(int k = 0; k < c; k++) {
                            // Row swap tmpm
                            tmpd = tmpm[j, k];
                            tmpm[j, k] = tmpm[imax, k];
                            tmpm[imax, k] = tmpd;
                            
                            // Row swap p
                            tmpd = p[j, k];
                            p[j, k] = p[imax, k];
                            p[imax, k] = tmpd;
                        }
                        pswap++;
                    }
                }
                
                // Compute the L and U matrices
                // Goes through the columns        
                for(int j = 0; j < c; j++) {

                    // For the u matrix
                    for(int i = 0; i < j + 1; i++) {
                        
                        // For each element in the u matrix
                        u[i, j] = tmpm[i, j];
                        for(int k = 0; k < i; k++) {
                            u[i, j] -= l[i, k] * u[k, j];
                        }
                    }
                    
                    // For the l matrix
                    for(int i = j + 1; i < r; i++) {
                        
                        // For each element in the l matrix
                        l[i, j] = tmpm[i, j];
                        for(int k = 0; k < j; k++) {
                            l[i, j] -= l[i, k] * u[k, j];
                        }
                        l[i, j] /= u[j, j];
                    }
                }
            }
        }
        
        // Finding the determinant using the LUDecomp
        public double LUDeterminant() {
            double d = 0.0;
            double sign = 0.0;
            if(!LU) {
                this.LUDecomp();
            }
            d = u[0, 0];
            for(int i = 1; i < u.GetLength(0); i++) {
                d *= u[i, i];
            }
            sign = (pswap % 2 == 0) ? 1.0 : -1.0;
            return sign * d;    
        }
        
        // Calls the recursive determinant function 
        public double determinant() {
            return determinant(m);
        }
        
        // Calculates and returns the determinant of a recieved matrix
        // Very inefficient. O(n!)
        private double determinant(double[,] mat) {
            double d = 0;
            int mod = 1;
            // Check if square
            if(mat.GetLength(0) != mat.GetLength(1)) {
                return 0;
            }
            // If only one element, return element
            if(mat.GetLength(0) == 1) {
                return mat[0, 0];
            }
            // If only 2x2, return determinant
            if(mat.GetLength(0) == 2) {
                return (mat[0, 0] * mat[1, 1] - mat[1, 0] * mat[0, 1]);
            }
            double[, ] tmp = new double[mat.GetLength(0) - 1, mat.GetLength(1) - 1];
            for(int i = 0; i < mat.GetLength(0); i++) {
                for(int j = 0; j < mat.GetLength(0) - 1; j++) {
                    for(int k = 0; k < mat.GetLength(1) - 1; k++) {
                        if(k < i) tmp[j, k] = mat[j+1, k];
                        else tmp[j, k] = mat[j+1, k+1];
                    }
                }
                if(i % 2 == 0) mod = 1;
                else mod = -1;
                d += mod * mat[0, i] * determinant(tmp);
            }
            return d;
        }
        
        // Checks that the dimensions of the matrices are equal
        // then adds the two matrices together. 
        // If the dimensions are not equal, returns a blank 1x1
        public static Matrix operator+ (Matrix A, Matrix B) {
            if(A.getRow() == B.getRow() && A.getCol() == B.getCol()) {
                Matrix Sum = new Matrix(A.getRow(), A.getCol());
                for(int i = 0; i < Sum.getRow(); i++) {
                    for(int j = 0; j < Sum.getCol(); j++) {
                        Sum.setElem(A.getElem(i, j) + B.getElem(i, j), i, j);
                    }
                }
                return Sum;
            }
            
            else{
                Console.WriteLine("These Matrices can not be added.");
                return new Matrix();
            }
        }
        
        // Checks that the dimensions of the matrices are equal
        // then subtracts the two matrices together. 
        // If the dimensions are not equal, returns a blank 1x1
        public static Matrix operator- (Matrix A, Matrix B) {
            if(A.getRow() == B.getRow() && A.getCol() == B.getCol()) {
                Matrix Sum = new Matrix(A.getRow(), A.getCol());
                for(int i = 0; i < Sum.getRow(); i++) {
                    for(int j = 0; j < Sum.getCol(); j++) {
                        Sum.setElem(A.getElem(i, j) - B.getElem(i, j), i, j);
                    }
                }
                return Sum;
            }
            
            else{
                Console.WriteLine("These Matrices can not be added.");
                return new Matrix();
            }
        }
        
        // Checks that the matrices can be multiplied. Multiplies
        // matrices.
        // If the matrices can't be multiplied, returns a blank 1x1
        public static Matrix operator* (Matrix A, Matrix B) {
            if(A.getCol() == B.getRow()) {
                int x = A.getRow(), y = A.getCol(), z = B.getCol();
                Matrix Product = new Matrix(x, z);
                for(int i = 0; i < x; i++) {
                    for(int j = 0; j < z; j++) {
                        for(int k = 0; k < y; k++) {
                            Product.setElem(Product.getElem(i, j) + (A.getElem(i, k) * B.getElem(k, j)), i, j);
                        }
                    }
                }
                return Product;
            }
            
            else {
                Console.WriteLine("These Matrices can not be multiplied.");
                return new Matrix();
            }
        }
    }
}