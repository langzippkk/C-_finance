#ifndef sudoku
#define sudoku

#include <vector>
#include <fstream>
#include <iostream>
using std::vector;
using namespace std;
class Sudoku
{
	// Private
	int puzzle[9][9];
	int Counter = 0;
	// Private member function that checks if the named row is valid
	bool row_valid(int row)
	{
		for (int i = 0; i < 8; i++)
		{
			for (int j = i + 1; j < 9; j++)
			{

				{
					if (puzzle[row][i] == puzzle[row][j] && (puzzle[row][i] != 0 || puzzle[row][j] != 0))
					return false;
				}
			}
			// write code that checks if "row" is valid
		}
		return true;
	}

	// Private member function that checks if the named column is valid
	bool col_valid(int col)
	{
		for (int i = 0; i < 8; i++)
		{
			for (int j = i + 1; j < 9; j++)
			{


				if (puzzle[i][col] == puzzle[j][col] && (puzzle[i][col] != 0 || puzzle[j][col] != 0))
					return false;

			}
			// check validity of "col" 
		}
		return true;
	}

	// Private member function that checks if the named 3x3 block is valid
	bool block_valid(int row, int col)
	{
		// check 3 x 3 block validity 
		int a = (row / 3 * 3);
		int b = (col / 3 * 3);
		for (int i = a; i < (a + 3); i++){
			for (int j = b; j< (b + 3); j++)
			{
				if ((i == row && j == col))
				{
					continue;
				}
				if (puzzle[i][j] == puzzle[row][col])
				{
					return false;
				}
			}
		}
		return true;
	}

public:
	// Public member function that reads the incomplete puzzle
	// we are not doing any checks on the input puzzle -- that is,
	// we are assuming they are indeed valid
	void read_puzzle(int argc, char * const argv[])
	{

		int value_just_read_from_file;
		vector <int> P; //vector of numbers

		ifstream input_file("input4"); // The input file name is "input1"
		// It contains the numbers P_1 P_2 ... P_M

		if (input_file.is_open()) // If Data exists in the local directory
		{
			while (input_file >> value_just_read_from_file) // As long as your are not at the "end-of-file" 
			{
				P.push_back(value_just_read_from_file);
			}
		}
		int k = 0;
		for (int i = 0; i < 9; i++)
		{
			for (int j = 0; j < 9; j++)
			{
				puzzle[i][j] = P[k];
				k++;
			}

		}
	}
	/*	cout << "Vector P has " << P.size() << " entries, and they are:" << endl;
	for (int i = 0; i < P.size(); i++)
	cout << P[i] << "\t";
	cout << endl;*/
	// write code that reads the input puzzle using the 
	// guidelines of figure 23 of the bootcamp material


	// Public member function that prints the puzzle when called
	void print_puzzle()
	{
		std::cout << std::endl << "Board Position" << std::endl;
		cout << endl;
		for (int i = 0; i < 9; i++)
		{
			for (int j = 0; j < 9; j++)
			{
				// check if we have a legitimate integer between 1 and 9
				if ((puzzle[i][j] >= 1) && (puzzle[i][j] <= 9))
				{
					// printing initial value of the puzzle with some formatting
					std::cout << puzzle[i][j] << " ";
				}
				else {
					// printing initial value of the puzzle with some formatting
					std::cout << "X ";
				}
			}
			std::cout << std::endl;
		}
	}

	// Public member function that (recursively) implements the brute-force 
	// search for possible solutions to the incomplete Sudoku puzzle
	/*bool next(int x, int y)
	{
	if (x == 8 && y == 8) return true;
	else if (x == 8) return Solve(0, y + 1);
	else return Solve(x + 1, y);
	}*/
	// I discussed this checkposition method with XianZhe Xu.
	bool checkposition(int &row, int &col){
		for (row = 0; row < 9; row++)
		{
			for (col = 0; col < 9; col++)
			{
				if (puzzle[row][col] == 0)
					return true;
			}
		}
		return false;
	}

	bool Solve(int row, int col)
	{
		// this part of the code identifies the row and col number of the 
		// first incomplete (i.e. 0) entry in the puzzle.  If the puzzle has
		// no zeros, the variable row will be 9 => the puzzle is done, as 
		// each entry is row-, col- and block-valid...		

		if (!checkposition(row, col))
		{
			cout << endl;
			Counter++;
			cout << "the_" << Counter << " th Solution" << endl;
			print_puzzle();
			goto L1;
			return true;
		}
		for (int k = 1; k < 10; k++)
		{

			puzzle[row][col] = k;
			if ((row_valid(row)) && (col_valid(col)) && (block_valid(row, col)))
			{
				Solve(row, col);
			}
			puzzle[row][col] = 0;
		}
		return false;
		L1:
		system("pause");
		return true;
		// use the pseudo code of figure 3 of the description
	}


};

#endif