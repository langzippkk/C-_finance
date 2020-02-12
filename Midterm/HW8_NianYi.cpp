//#include <vector>
//#include <iostream>
//#include <iomanip>
//#include <cmath>
//#include <fstream>
//#include <ctime>
//#include <cstdlib>
//#include <algorithm>
//#include <complex>
//#include "C:\Users\nian yi\Documents\Visual Studio 2013\newmat11\newmatap.h"           
//#include "C:\Users\nian yi\Documents\Visual Studio 2013\newmat11\newmat.h"
//#include "C:\Users\nian yi\Documents\Visual Studio 2013\newmat11\newmatio.h" 
//using namespace std;
//class Squaring{
//public:
//	int num_rows, exp;
//	Matrix B, A_store;
//	Matrix repeated_squaring(Matrix A, int exponent, int no_rows)
//	{
//		num_rows = no_rows;
//		exp = exponent;
//		A_store = A;
//		if (0 == exponent)
//			return IdentityMatrix(no_rows);
//		if ((exponent % 2) == 1)
//			return (A*repeated_squaring(A*A, (exponent - 1) / 2, no_rows));
//		else
//			return (repeated_squaring(A*A, exponent / 2, no_rows));
//	}
//	Matrix brute_multiplication(Matrix A, int exponent, int no_rows)
//	{
//		Matrix C(no_rows, no_rows);
//		C = A;
//		for (int i = 1; i < exponent; i++)
//			C = A*C;
//		return C;
//	}
//	void squaring_output(Matrix A, int exponent, int no_rows)
//	{
//		const clock_t begin_time = clock();
//		float diff;
//		B = repeated_squaring(A, exponent, no_rows);
//		clock_t time_after = clock();
//		diff = ((float)time_after - (float)begin_time);
//		cout << "it took " << (float)diff / CLOCKS_PER_SEC << " second to complete" << endl;
//		cout << "The exponent = " << exponent << endl;
//		cout << "The number of rows and columns of the matrix = " << no_rows << endl;
//		for (int i = 1; i <= no_rows; i++)
//		{
//			for (int j = 1; j <= no_rows; j++)
//				cout << B(i, j) << " ";
//			cout << endl;
//		}
//	}
//	void brute_output(Matrix A, int exponent, int no_rows)
//	{
//		const clock_t begin_time = clock();
//		float diff;
//		B = brute_multiplication(A, exponent, no_rows);
//		clock_t time_after = clock();
//		diff = ((float)time_after - (float)begin_time);
//		cout << time_after << endl;
//		cout << begin_time << endl;
//		cout << "it took " <<setprecision(6)<<(float) diff / CLOCKS_PER_SEC << " second to complete" << endl;
//		cout << "The exponent = " << exponent << endl;
//		cout << "The number of rows and columns of the matrix = " << no_rows << endl;
//		for (int i = 1; i <= no_rows; i++)
//		{
//			for (int j = 1; j <= no_rows; j++)
//				cout << B(i, j) << " ";
//			cout << endl;
//		}
//	}
//	void plot(Matrix A, int no_rows)
//	{
//		float diff;
//		clock_t begin_time, time_after;
//		ofstream output_file("file_matrix");
//		output_file << "brute:" << endl;
//		for (int i = 10000; i <= 11000; i++)
//		{
//			clock_t begin_time = clock();
//			B = brute_multiplication(A, i, no_rows);
//			clock_t time_after = clock();
//			diff = ((float)time_after - (float)begin_time)/CLOCKS_PER_SEC;
//			output_file << (float)diff << endl;
//		}
//		//output_file << "square:" << endl;
//		/*for (int i = 0; i < 10500; i++)
//		{
//			begin_time = clock();
//			B = repeated_squaring(A, i, no_rows);
//			time_after = clock();
//			diff = ((float)time_after - (float)begin_time) / CLOCKS_PER_SEC;
//			output_file << (float)diff << endl;
//		}*/
//	}
//};
//int main(){
//	Matrix A(5,5),B(5,5);
//	for(int i = 1; i <= 5; i++)
//	for (int j = 1; j <= 5; j++)
//		A(i, j) = (i + j) ;
//	/*int n = 3;
//	B = brute_multiplication(A, 11, 3);
//	for (int i = 1; i <= n; i++)
//	{
//		for (int j = 1; j <= n; j++)
//			cout << B(i, j) << " ";
//		cout << endl;
//	}*/
//	Squaring x;
//	//x.repeated_squaring(A, 11, 3);
//	//x.brute_output(A, 95, 5);
//	x.plot(A, 5);
//	system("pause");
//}
//// Public member function that reads the incomplete puzzle we are not doing any checks on the input puzzle -- that is, we are assuming they are indeed valid
///*void run_the_filter(int argc, char * const argv[])
//{
//sscanf(argv[1], "%d", &no_of_row);
//sscanf(argv[2], "%d", &exponent_k);
//
//std::cout << "No. of rows is" << argv[3] << std::endl;
//std::cout << "k = " << no_of_data_points << endl;
//std::cout << "result: " << no_of_terms << endl;
//std::cout << "time needed " << argv[4] << std::endl;
//
//ColumnVector data(no_of_data_points), filtered_data(no_of_data_points);
//
//// get ticker data
//get_data(argv[3], data);
//
//// filter the ticker data
//filter_the_data(data, filtered_data, no_of_terms);
//
//// write the filtered data
//write_data(argv[4], filtered_data);
//}*/
//
//
//
//
//
//
//
//
