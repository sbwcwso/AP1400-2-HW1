#include "hw1.h"

namespace algebra
{
    /** Create a n x m matrix with all elements equal to zero. */
    Matrix zeros(size_t n, size_t m)
    {
        return Matrix(n, vector<double>(m, 0.0));
    }

    /** create a n x m matrix with all elements equal to one. */
    Matrix ones(size_t n, size_t m)
    {
        return Matrix(n, vector<double>(m, 1.0));
    }

    /** create a n x m matrix with all elements a random number between min and max. */
    Matrix random(size_t n, size_t m, double min, double max)
    {
        if (min > max)
        {
            throw logic_error("Minimun number cannot be greater than maximun value.");
        }
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<> dis(min, max);
        Matrix matrix(n, vector<double>(m));
        for (vector<double> &row : matrix)
        {
            for (double &col : row)
            {
                col = dis(gen);
            }
        }
        return matrix;
    }

    /** show the matrix, each element of the matrix should have
     *  exactly 3 decimal places. */
    void show(const Matrix &matrix)
    {
        for (vector<double> row : matrix)
        {
            for (double col : row)
            {
                cout << setprecision(3) << col << "\t";
            }
            cout << endl;
        }
    }

    /**  multiplies the matrix into the constant scalar c*/
    Matrix multiply(const Matrix &matrix, double c)
    {
        if (matrix.empty())
        {
            return Matrix{};
        }
        size_t m = matrix.size();
        size_t n = matrix[0].size();
        Matrix result = zeros(m, n);
        for (int i = 0; i < m; i++)
        {
            for (int j = 0; j < n; j++)
            {
                result[i][j] = c * matrix[i][j];
            }
        }
        return result;
    }

    /**  multiplies the matrix1 into matrix2. (this is not an element-wise multiplication) */
    Matrix multiply(const Matrix &matrix1, const Matrix &matrix2)
    {
        if (matrix1.empty() && matrix2.empty())
        {
            return Matrix{};
        }

        if (matrix1.empty() || matrix2.empty())
        {
            throw logic_error("matrices with wrong dimensions cannot be multiplied.");
        }

        size_t m1 = matrix1.size();
        size_t m2 = matrix2.size();
        size_t n1 = matrix1[0].size();
        size_t n2 = matrix2[0].size();
        if (n1 != m2)
        {
            throw logic_error("matrices with wrong dimensions cannot be multiplied.");
        }
        Matrix result = zeros(m1, n2);
        for (int i = 0; i < m1; i++)
        {
            for (int j = 0; j < n2; j++)
            {
                for (int k = 0; k < n1; k++)
                {
                    result[i][j] += matrix1[i][k] * matrix2[k][j];
                }
            }
        }
        return result;
    }

    /** adds the constant number c to every element of matrix. */
    Matrix sum(const Matrix &matrix, double c)
    {
        if (matrix.empty())
        {
            return Matrix{};
        }
        size_t m = matrix.size();
        size_t n = matrix[0].size();
        Matrix result = zeros(m, n);
        for (int i = 0; i < m; i++)
        {
            for (int j = 0; j < n; j++)
            {
                result[i][j] = matrix[i][j] + c;
            }
        }
        return result;
    }

    /** adds 2 matrices to each other. */
    Matrix sum(const Matrix &matrix1, const Matrix &matrix2)
    {
        if (matrix1.empty() && matrix2.empty())
        {
            return Matrix{};
        }

        size_t m1 = matrix1.size();
        size_t m2 = matrix2.size();
        if (m1 != m2)
        {
            throw logic_error("Matrices with wrong dimensions cannot be summed.");
        }

        size_t n1 = matrix1[0].size();
        size_t n2 = matrix2[0].size();
        if (n1 != n2)
        {
            throw logic_error("Matrices with wrong dimensions cannot be summed.");
        }

        Matrix result = zeros(m1, n1);
        for (int i = 0; i < m1; i++)
        {
            for (int j = 0; j < n1; j++)
            {
                result[i][j] = matrix1[i][j] + matrix2[i][j];
            }
        }
        return result;
    }

    /** generate the transpose matrix of the input matrix */
    Matrix transpose(const Matrix &matrix)
    {
        if (matrix.empty())
        {
            return Matrix{};
        }
        size_t m = matrix.size();
        size_t n = matrix[0].size();

        Matrix result = zeros(n, m);
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < m; j++)
            {
                result[i][j] = matrix[j][i];
            }
        }
        return result;
    }

    /** create the minor of the input matrix with respect
     * to nth row and mth column. */
    Matrix minor(const Matrix &matrix, size_t n, size_t m)
    {
        size_t rows = matrix.size();
        // if (rows == 0) {
        //     throw logic_error("input matrix is empty!");
        // }
        size_t columns = matrix[0].size();
        // if (rows != columns) {
        //     throw logic_error("input matrix should be a square matrix!");
        // }

        // if (n >= rows || m >= columns) {
        //     throw logic_error("the index out of range!");
        // }

        Matrix minor_matix = zeros(rows - 1, columns - 1);
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < m; j++)
            {
                minor_matix[i][j] = matrix[i][j];
            }
            for (int j = m + 1; j < columns; j++)
            {
                minor_matix[i][j - 1] = matrix[i][j];
            }
        }

        for (int i = n + 1; i < rows; i++)
        {
            for (int j = 0; j < m; j++)
            {
                minor_matix[i - 1][j] = matrix[i][j];
            }
            for (int j = m + 1; j < columns; j++)
            {
                minor_matix[i - 1][j - 1] = matrix[i][j];
            }
        }

        return minor_matix;
    }

    /** calculates the determinant of the input matrixa */
    double determinant(const Matrix &matrix)
    {
        size_t rows = matrix.size();
        if (rows == 0)
        {
            return 1.0;
        }
        size_t columns = matrix[0].size();
        if (rows != columns)
        {
            throw logic_error("non-square matrices have no determinant");
        }

        double result = 0;
        double sign = 1;
        for (int j = 0; j < columns; j++)
        {
            result += sign * matrix[0][j] * determinant(minor(matrix, 0, j));
            sign = -sign;
        }
        return result;
    }

    /** generates the matrix's inverse. */
    Matrix inverse(const Matrix &matrix)
    {
        size_t rows = matrix.size();
        if (rows == 0)
        {
            return Matrix{};
        }
        size_t columns = matrix[0].size();
        if (rows != columns)
        {
            throw logic_error("non-square matrices have no inverse");
        }
        double determinant_matrix = determinant(matrix);
        if (determinant_matrix == 0)
        {
            throw logic_error("singular matrices have no inverse");
        }

        Matrix adj = zeros(rows, rows);
        double sign_row = 1;
        for (size_t i = 0; i < rows; i++)
        {
            double sign_column = sign_row;
            for (size_t j = 0; j < rows; j++)
            {
                adj[i][j] = sign_column * determinant(minor(matrix, j, i));
                sign_column = -sign_column;
            }
            sign_row = -sign_row;
        }

        return multiply(adj, 1.0 / determinant_matrix);
    }

    /** concatenate matrix1 and matrix2 along the specified axis. (axis=0: on top of each other | axis=1: alongside each other).*/
    Matrix concatenate(const Matrix &matrix1, const Matrix &matrix2, int axis)
    {
        if (matrix1.empty() || matrix2.empty())
        {
            throw logic_error("the input matrix is empty");
        }
        size_t row1 = matrix1.size();
        size_t row2 = matrix2.size();
        size_t col1 = matrix1[0].size();
        size_t col2 = matrix2[0].size();
        if (axis == 0)
        {
            if (col1 != col2)
                throw logic_error("matrices with wrong dimensions cannot be concatenated");
            Matrix result = zeros(row1 + row2, col1);
            for (size_t i = 0; i < row1; i++)
                for (size_t j = 0; j < col1; j++)
                    result[i][j] = matrix1[i][j];

            for (size_t i = 0; i < row2; i++)
                for (size_t j = 0; j < col2; j++)
                    result[i + row1][j] = matrix2[i][j];
            return result;
        }
        else if (axis == 1)
        {
            if (row1 != row2)
                throw logic_error("matrices with wrong dimensions cannot be concatenated");
            Matrix result = zeros(row1, col1 + col2);
            for (size_t i = 0; i < row1; i++)
                for (size_t j = 0; j < col1; j++)
                    result[i][j] = matrix1[i][j];

            for (size_t i = 0; i < row2; i++)
                for (size_t j = 0; j < col2; j++)
                    result[i][j + col1] = matrix2[i][j];
            return result;
        }
        else
            throw logic_error("axis must be 0 or 1");
    }

    /** elementary row operation (ERO) */
    /** swaps r1th row with r2th.*/
    Matrix ero_swap(const Matrix &matrix, size_t r1, size_t r2)
    {
        size_t row = matrix.size();
        if (r1 >= row || r2 >= row)
            throw logic_error("r1 or r2 inputs are out of range.");

        if (matrix.empty())
            return Matrix{};

        size_t col = matrix[0].size();
        Matrix result = zeros(row, col);
        for (size_t i = 0; i < row; i++)
        {
            int ii = i;
            if (i == r1)
                ii = r2;
            else if (i == r2)
                ii = r1;
            for (size_t j = 0; j < col; j++)
                result[ii][j] = matrix[i][j];
        }
        return result;
    }

    /** multiplies every element in rth row with constant number c.*/
    Matrix ero_multiply(const Matrix &matrix, size_t r, double c)
    {
        size_t row = matrix.size();
        if (r >= row)
            throw logic_error("r input is out of range.");

        if (matrix.empty())
            return Matrix{};

        size_t col = matrix[0].size();
        Matrix result = zeros(row, col);
        for (size_t i = 0; i < row; i++)
            for (size_t j = 0; j < col; j++)
                if (i == r)
                    result[i][j] = c * matrix[i][j];
                else
                    result[i][j] = matrix[i][j];
        return result;
    }

    /** adds r1th x c into r2th row.*/
    Matrix ero_sum(const Matrix &matrix, size_t r1, double c, size_t r2)
    {
        size_t row = matrix.size();
        if (r1 >= row || r2 >= row)
            throw logic_error("r1 or r2 inputs are out of range.");

        if (matrix.empty())
            return Matrix{};

        size_t col = matrix[0].size();
        Matrix result = zeros(row, col);

        for (size_t i = 0; i < row; i++)
            for (size_t j = 0; j < col; j++)
                if (i == r2)
                    result[i][j] = c * matrix[r1][j] + matrix[i][j];
                else
                    result[i][j] = matrix[i][j];
        return result;
    }

    /** calculate the upper triangular form of the matrix using the ERO operations.*/
    Matrix upper_triangular(const Matrix &matrix)
    {
        if (matrix.empty())
            return Matrix{};

        size_t row = matrix.size();
        size_t col = matrix[0].size();

        if (row != col)
        {
            throw logic_error("non-square matrices have no upper triangular form");
        }

        Matrix result = zeros(row, col);

        for (size_t i = 0; i < row; i++)
            for (size_t j = 0; j < row; j++)
                result[i][j] = matrix[i][j];

        for (size_t i = 0; i < row - 1; i++)
        {
            if (result[i][i] == 0)
            {
                bool found = false;
                for (size_t j = i + 1; j < row; j++)
                    if (result[j][i] != 0)
                    {
                        result = ero_swap(result, i, j);
                        found = true;
                        break;
                    }
                if (!found)
                    continue;
            }
            for (size_t j = i + 1; j < row; j++)
            {
                double c = -result[j][i] / result[i][i];
                result = ero_sum(result, i, c, j);
            }
        }

        return result;
    }
}