/* Experimental illustration of the relativity of information.
/* The program tests whether, whenever the result is computed for an input,
   it is computed for all possible inputs simultaneously,
   in an "encrypted" form that can be "decrypted".
   The program needs two input text files, which will be created if absent,
   and saves another text file containing the results on the computer. */

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>

typedef std::vector<int> Vector;
typedef std::vector<Vector> Matrix;

// Converts a vector to a string
std::string ToString(const Matrix& V)
{
    std::string str = "(";
    for (size_t i = 0; i < V.size(); i++)
    {
        str += std::to_string(V[i][0]);
        str += (i + 1 == V.size() ? ")" : ",");
    }
    return str;
}

// Prints a matrix
void PrintMatrix(const Matrix& M)
{
    for (size_t row = 0; row < M.size(); row++)
    {
        for (size_t col = 0; col < M[row].size(); col++)
            std::cout << M[row][col] << "\t";

        std::cout << std::endl;
    }
}

// Computes the product of two matrices
Matrix Multiply(const Matrix& A, const Matrix& B)
{
    Matrix AxB(A.size(), Vector(B[0].size(), 0));

    for (size_t row = 0; row < AxB.size(); row++)
    {
        for (size_t col = 0; col < AxB[row].size(); col++)
        {
            for (size_t i = 0; i < A[row].size(); i++)
                AxB[row][col] += A[row][i] * B[i][col];
        }
    }
    
    return AxB;
}

// Computes the transposed of a matrix
Matrix Transposed(const Matrix& M)
{
    Matrix M_transposed(M[0].size(), Vector(M.size(), 0));

    for (size_t row = 0; row < M_transposed.size(); row++)
    {
        for (size_t col = 0; col < M_transposed[row].size(); col++)
            M_transposed[row][col] = M[col][row];
    }

    return M_transposed;
}

// Computes -Vector
Vector operator-(const Vector& V)
{
    Vector minusV;
    for (size_t i = 0; i < V.size(); i++)
        minusV.push_back(-V[i]);

    return minusV;
}

// Computes -Matrix
Matrix operator-(const Matrix& M)
{
    Matrix minusM;
    for (size_t row = 0; row < M.size(); row++)
        minusM.push_back(-M[row]);

    return minusM;
}

// Computes the U matrix
Matrix ComputeMatrixU(const Matrix& v)
{
    // Initialize matrix U with 0
    Matrix U(v.size() + 1, Vector(v.size() + 1, 0));

    // Fill U with 1 on the main diagonal
    for (size_t i = 0; i <= v.size(); i++)
        U[i][i] = 1;

    // Fill the first row of U with the elements of v
    for (size_t i = 0; i < v.size(); i++)
        U[0][i + 1] = v[i][0];

    return U;
}

// Computes the inverse of the U matrix
Matrix ComputeMatrixU_inverse(const Matrix& v)
{
    return ComputeMatrixU(-v);
}

// Converts partition matrix to extended S matrix
Matrix PartitionToSignMatrix1(const Matrix& w)
{
    Matrix SignMatrix1(w.size() + 1, Vector(w.size() + 1, 0));
    SignMatrix1[0][0] = 1;
    for (size_t i = 0; i < w.size(); i++)
        SignMatrix1[i + 1][i + 1] = w[i][0];

    return SignMatrix1;
}

// Converts a string of bits to a partition matrix
Matrix BitStringToPartition(const std::string& strBits)
{
    Matrix w(strBits.size(), Vector(1, 0));
    for (size_t i = 0; i < strBits.size(); i++)
    {
        if (strBits[i] == '1')
            w[i][0] = +1;
        else
            w[i][0] = -1;
    }

    return w;
}

// Converts a number to a string of bits
std::string DecimalToBitString(unsigned int numDec, int size)
{
    std::string strBits;
    for (int i = size-1; i >= 0; i--)
    {
        int k = numDec >> i;
        if (k & 1)
            strBits.append({ '1' });
        else
            strBits.append({ '0' });
    }

    return strBits;
}

// Computes the matrix that transforms the result for w0 to that for w
Matrix ComputeTransformationMatrix(const Matrix& v, const Matrix& w0, const Matrix& w)
{
    const Matrix U = ComputeMatrixU(v);
    const Matrix U_inverse = ComputeMatrixU_inverse(v);
    const Matrix S10 = PartitionToSignMatrix1(w0);
    const Matrix S1 = PartitionToSignMatrix1(w);
    const Matrix S = Multiply(S10, S1);
    return Multiply(U, Multiply(S, U_inverse));
}

// Computes U applied to w
Matrix ComputePartitionEvaluator(const Matrix& v, const Matrix& w)
{
    Matrix U = ComputeMatrixU(v);
    Matrix w1 = w;
    w1.insert(w1.begin(), { 0 });
    return Multiply(U, w1);
}

// Verifies if the partition w of v is fair
bool IsPartitionFair(const Matrix& v, const Matrix& w)
{
    Matrix Uw1 = ComputePartitionEvaluator(v, w);
    return Uw1[0][0] == 0;
}

// Reads the list of numbers to be partitioned
Matrix ReadListOfNumbers()
{
    Vector v;
    std::fstream fileNumbers("list of numbers.txt", std::ios_base::in);
    if (fileNumbers.good())
    {
        int number;
        while (fileNumbers >> number)
            v.push_back(number);

        return Transposed({ v });
    }
    else
    {
        // If the file doesn't exist, create it and fill it
        std::ofstream newFile("list of numbers.txt");
        newFile << "5 5 5 10 10 25";
        newFile.close();
        return ReadListOfNumbers();
    }
}

// Reads the partition into the vector w, interpreting '0' as -1 and '1' as +1
Matrix ReadPartition()
{
    std::fstream filePartition("partition.txt", std::ios_base::in);
    if (filePartition.good())
    {
        std::string line;
        std::getline(filePartition, line);
        return BitStringToPartition(line);
    }
    else
    {
        // If the file doesn't exist, create it and fill it
        std::ofstream newFile("partition.txt");
        newFile << "000000";
        newFile.close();
        return ReadPartition();
    }
}

// Test a partition and returns the result in a string
bool TestPartition(const Matrix& evaluator_w0, const Matrix& v, const Matrix& w0, const std::string& strBits, std::string& strResult)
{
    const Matrix w = BitStringToPartition(strBits);
    Matrix evaluator_w_direct = ComputePartitionEvaluator(v, w);
    Matrix R = ComputeTransformationMatrix(v, w0, w);
    Matrix evaluator_w_transformation = Multiply(R, evaluator_w0);

    strResult += strBits + "\t";
    strResult += std::to_string(evaluator_w_transformation[0][0]) + "\t";
    strResult += std::to_string(evaluator_w_direct[0][0]) + "\t";

    bool bIsCorrect = true;
    if (evaluator_w_transformation[0][0] == evaluator_w_direct[0][0])
        strResult += "Yes\t";
    else
    {
        strResult += "No\t";
        bIsCorrect = false;
    }

    if (evaluator_w_direct[0][0] == 0)
        strResult += "Yes";
    else
        strResult += "No";

    return bIsCorrect;
}

// Tests whether all partitions are computed simultaneously with w0
bool TestAllPartitions(const Matrix& evaluator_w0, const Matrix& v, const Matrix& w0)
{
    // Create and open a file
    const std::string fileTestResultsName("test results.csv");
    std::ofstream fileTestResults(fileTestResultsName);
    fileTestResults << "Partition\tTransformed\tDirect\tCorrect\tFair" << std::endl;

    const unsigned int maxDec = (int)std::pow(2, v.size());
    bool bIsCorrect = true;
    for (unsigned int decNum = 0; decNum < maxDec; decNum++)
    {
        const std::string strBits = DecimalToBitString(decNum, (int)v.size());
        std::string strResult;
        if (!TestPartition(evaluator_w0, v, w0, strBits, strResult))
            bIsCorrect = false;

        fileTestResults << strResult << std::endl;
    }
    fileTestResults.close();

    if (bIsCorrect)
    {
        std::cout << "All alternatives were found! ";
        std::cout << "The results are saved in the file \"" << fileTestResultsName << "\". ";
        std::cout << "It can be opened as a spreadsheet or as a text file.\n";
    }
    else
        std::cout << "Verification failed.\n";

    return bIsCorrect;
}

// Tests if a user defined partition is computed simultaneously with w0
bool TestUserDefinedPartition(const Matrix& evaluator_w0, const Matrix& v, const Matrix& w0)
{
    std::cout << "Please input " << v.size() << " bits (with values \'0\' or \'1\')\n";

    std::string strBits;

    while (strBits.size() < v.size())
    {
        std::string strInput;
        std::cin >> strInput;
        for (size_t i = 0; i < strInput.size(); i++)
        {
            if (strInput[i] == '0' || strInput[i] == '1')
                strBits.append({ strInput[i] });
        }
    }

    if (strBits.size() > v.size())
        strBits.erase(v.size(), std::string::npos);

    std::string strResult;
    const bool bIsCorrect = TestPartition(evaluator_w0, v, w0, strBits, strResult);
    if (bIsCorrect)
    {
        std::cout << "Partition\tTransformed\tDirect\tCorrect\tFair\n";
        std::cout << strResult << "\n";
        std::cout << "Verification successfull!\n";
    }
    else
        std::cout << "Verification failed.\n";

    return bIsCorrect;
}

// main function of the program
int main()
{
    Matrix v = ReadListOfNumbers();
    std::cout << "v = " << ToString(v) << std::endl;
    Matrix w0 = ReadPartition();
    std::cout << "w0 = " << ToString(w0) << std::endl;

    Matrix evaluator_w0 = ComputePartitionEvaluator(v, w0);

    std::cout << "Help:\n";
    std::cout << "Press \"i\" to input a partition to test.\n";
    std::cout << "Press \"a\" to test that all possible partitions are verified simultaneously.\n";
    std::cout << "Press \"x\" to exit.\n";

    char key = 0;
    while(key != 'x')
    {
        std::cin.get(key);
        key = std::tolower(key);
        switch (key)
        {
        case 'a':
            TestAllPartitions(evaluator_w0, v, w0);
            break;
        case 'i':
            TestUserDefinedPartition(evaluator_w0, v, w0);
            break;
        }
    }

    return 0;
}
