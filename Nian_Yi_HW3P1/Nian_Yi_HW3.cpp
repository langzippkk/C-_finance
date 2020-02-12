//
//// There is an ILP formulation for any 9 x 9 Sudoku Puzzle
//// See the handout for the formulation -- there are 729 variables
//// with a whole lot of constraints that is driven by the input
//// (which is the incomplete board-layout).  You can write this
//// (moderate-size) ILP automatically in C++, which is then solved
//// using the API for Lpsolve.
////#define WIN32
//#include <iostream>
//#include <iomanip>
//#include <cmath>
//#include <fstream>
//#include <cstdlib>
//
//#include "lp_lib.h"
//
//using namespace std;
//
//// Global variables for the problem 9*9*9 = 729
//// See my handout for further explanations
//lprec *lp;
//double solution[729];
//
//// This sets the generic ILP for the Sudoku Puzzle
//// It does not use any of the pre-determined board-positions,
//// that part is done after the input file is read.
//void set_sudoku_ilp()
//{
//	// setting the problem up: 9*9*9 = 729 variables
//	lp = make_lp(0, 729);
//
//	// This keeps the message reporting of lp_solve to a minimum
//	set_verbose(lp, 3);
//
//	{
//		for (int i = 1; i <= 9; i++)
//		{
//			for (int j = 1; j <= 9; j++)
//			{
//				// constraint that says each (i,j) entry in the 9x9 table
//				// should have one number in it
//				// creating a row with 1's at the right spot (the first
//				// entry has to be a zero; this is idiosynchratic of lpsolve)
//
//				double row[730];
//				for (int k = 0; k < 730; k++)
//					row[k] = 0;
//
//				for (int k = 1; k <= 9; k++)
//					row[(81 * (i - 1)) + (9 * (j - 1)) + k] = 1;
//
//				// adding the constraint
//				add_constraint(lp, row, EQ, 1);
//			}
//		}
//	}
//	//     check row;
//	{
//
//		for (int i = 1; i <= 9; i++)
//		{
//			for (int k = 1; k <= 9; k++)
//			{
//				double row[730];
//				for (int j = 0; j < 730; j++)
//					row[j] = 0;
//
//				for (int j = 1; j <= 9; j++)
//					row[(81 * (i - 1)) + (9 * (j - 1)) + k] = 1;
//				// adding the constraint
//				add_constraint(lp, row, EQ, 1);
//			}
//		}
//	}
//	//   check column;
//	{
//		for (int j = 1; j <= 9; j++)
//		{
//			for (int k = 1; k <= 9; k++)
//			{
//				double row[730];
//				for (int i = 0; i < 730; i++)
//					row[i] = 0;
//
//				for (int i = 1; i <= 9; i++)
//					row[(81 * (i - 1)) + (9 * (j - 1)) + k] = 1;
//				// adding the constraint
//				add_constraint(lp, row, EQ, 1);
//			}
//		}
//	//   making sure each number occurs once within each block
//	//	 Block 1
//		{
//			for (int k = 1; k <= 9; k++)
//			{
//				double row[730];
//				for (int i = 0; i < 730; i++)
//					row[i] = 0;
//				for (int i = 1; i <= 3; i++)
//				{
//					for (int j = 1; j <= 3; j++)
//					{
//						row[(81 * (i - 1)) + (9 * (j - 1)) + k] = 1;
//					}
//				}
//				add_constraint(lp, row, EQ, 1);
//			}
//			
//		}
//		// Block 2
//		{
//			for (int k = 1; k <= 9; k++)
//			{
//				double row[730];
//				for (int i = 0; i < 730; i++)
//					row[i] = 0;
//				for (int i = 1; i <= 3; i++)
//				{
//					for (int j = 4; j <= 6; j++)
//					{
//
//						row[(81 * (i - 1)) + (9 * (j - 1)) + k] = 1;
//					}
//				}
//				add_constraint(lp, row, EQ, 1);
//			}
//		}
//
//		// Block 3
//		{
//			for (int k = 1; k <= 9; k++)
//			{
//				double row[730];
//				for (int i = 0; i < 730; i++)
//					row[i] = 0;
//				for (int i = 1; i <= 3; i++)
//				{
//					for (int j = 7; j <= 9; j++)
//					{
//						row[(81 * (i - 1)) + (9 * (j - 1)) + k] = 1;
//
//					}
//				}
//				add_constraint(lp, row, EQ, 1);
//			}
//		}
//
//		// Block 4
//		{
//			for (int k = 1; k <= 9; k++)
//			{
//				double row[730];
//				for (int i = 0; i < 730; i++)
//					row[i] = 0;
//				for (int i = 4; i <= 6; i++)
//				{
//					for (int j = 1; j <= 3; j++)
//					{
//						row[(81 * (i - 1)) + (9 * (j - 1)) + k] = 1;
//
//					}
//				}
//				add_constraint(lp, row, EQ, 1);
//			}
//		}
//
//		// Block 5
//		{
//			for (int k = 1; k <= 9; k++)
//			{
//				double row[730];
//				for (int i = 0; i < 730; i++)
//					row[i] = 0;
//				for (int i = 4; i <= 6; i++)
//				{
//					for (int j = 4; j <= 6; j++)
//					{
//						row[(81 * (i - 1)) + (9 * (j - 1)) + k] = 1;
//					}
//				}
//				add_constraint(lp, row, EQ, 1);
//			}
//		}
//
//
//		// Block 6
//		{
//			for (int k = 1; k <= 9; k++)
//			{
//				double row[730];
//				for (int i = 0; i < 730; i++)
//					row[i] = 0;
//				for (int i = 4; i <= 6; i++)
//				{
//					for (int j = 7; j <= 9; j++)
//					{
//						row[(81 * (i - 1)) + (9 * (j - 1)) + k] = 1;
//						
//					}
//				}
//				add_constraint(lp, row, EQ, 1);
//			}
//		}
//
//		// Block 7
//		{
//			for (int k = 1; k <= 9; k++)
//			{
//				double row[730];
//				for (int i = 0; i < 730; i++)
//					row[i] = 0;
//				for (int i = 7; i <= 9; i++)
//				{
//					for (int j = 1; j <= 3; j++)
//					{
//						row[(81 * (i - 1)) + (9 * (j - 1)) + k] = 1;
//
//					}
//				}
//				add_constraint(lp, row, EQ, 1);
//			}
//		}
//
//
//		// Block 8
//		{
//			for (int k = 1; k <= 9; k++)
//			{
//				double row[730];
//				for (int i = 0; i < 730; i++)
//					row[i] = 0;
//				for (int i = 7; i <= 9; i++)
//				{
//					for (int j = 4; j <= 6; j++)
//					{
//
//						row[(81 * (i - 1)) + (9 * (j - 1)) + k] = 1;
//					}
//				}
//				add_constraint(lp, row, EQ, 1);
//			}
//		}
//
//
//		// Block 9
//		{
//			for (int k = 1; k <= 9; k++)
//			{
//				double row[730];
//				for (int i = 0; i < 730; i++)
//					row[i] = 0;
//				for (int i = 7; i <= 9; i++)
//				{
//					for (int j = 7; j <= 9; j++)
//					{
//
//						row[(81 * (i - 1)) + (9 * (j - 1)) + k] = 1;
//					}
//				}
//				add_constraint(lp, row, EQ, 1);
//			}
//		}
//
//		//// upper-bound each variable by 1
//		for (int i = 1; i <= 729; i++)
//		{
//			double row[730];
//			for (int j = 0; j < 730; j++)
//				row[j] = 0;
//			row[i] = 1;
//			add_constraint(lp, row, LE, 1);
//		}
//
//		// it does not matter what the objective function (why?)
//		// I am minimizing the sum of all variables.
//		{
//			double row[730];
//			for (int i = 1; i < 730; i++)
//				row[i] = 1;
//			row[0] = 0;
//
//			set_obj_fn(lp, row);
//		}
//
//		// set the variables to be integers
//		for (int i = 1; i <= 729; i++)
//			set_int(lp, i, TRUE);
//	}
//}
//
//	// This subroutine reads the incomplete board-information from the
//	// input file and sets the appropriate constraints before the ILP
//	// is solved.
//	void read_input_data(char* argv[])
//	{
//		// reading the input filename from commandline
//		//ifstream input_filename(argv[1]);
//		ifstream input_filename(argv[1]);
//
//		if (input_filename.is_open()) {
//			//	cout << "Input File Name: " << argv[1] << endl;
//			cout << endl << "Initial Board Position" << endl;
//			for (int i = 1; i <= 9; i++)
//			{
//				for (int j = 1; j <= 9; j++)
//				{
//					int value_just_read;
//					input_filename >> value_just_read;
//
//					if ((value_just_read >= 1) && (value_just_read <= 9))
//					{
//						// printing initial value of the puzzle with some formatting
//						cout << value_just_read << " ";
//						double row[730];
//						for (int k = 0; k <= 729;k++)
//						{
//							row[k] = 0;
//						}
//						row[(81 * (i - 1)) + (9 * (j - 1)) + value_just_read] = 1;
//						add_constraint(lp, row, EQ, 1);
//						// add appropriate constraints that bind the value of the
//						// appropriate variables based on the incomplete information
//						// that was read from the input file
//					}
//					else {
//						// printing initial value of the puzzle with some formatting
//						cout << "X ";
//					}
//				}
//				cout << endl;
//			}
//		}
//		else {
//			cout << "Input file missing" << endl;
//			exit(0);
//		}
//	}
//
//	// The ILP formulation is solved using the API for Lpsolve
//	// Pay attention to how the solution is interpretted...
//	void solve_the_puzzle()
//	{
//		int ret;
//		double row[730];
//		// solve the lp
//		ret = solve(lp);
//		// Check if you had a solution 
//		// (see online guide for the codes at http://lpsolve.sourceforge.net/5.0/ )
//		if (ret == 0)
//		{
//			// get the optimal assignment
//			get_variables(lp, solution);
//	
//			// print the solution
//
//
//			cout << endl << "Final Solution" << endl;
//			{
//				int puzzle[9][9];
//				for (int i = 0; i < 9; i++)
//				{
//					for (int j = 0; j < 9; j++)
//					{
//						puzzle[i][j] = 0;
//					}
//				}
//
//				row[0] = 0;
//				for (int k = 1; k <= 729; k++)
//					row[k] = 0;
//
//				for (int k = 1; k <= 729; k++)
//				row[k] = solution[k - 1];
//				
//
//				for (int i = 0; i <9; i++)
//				{
//					for (int j = 0; j <9; j++)
//					{
//						for (int k = 1; k <=9; k++)
//						{
//							if (row[81 * i + 9 * j + k] == 1)
//								puzzle[i][j] = k;
//
//						}
//					}
//				}
//
//				for (int i = 0; i < 9; i++)
//				{
//					for (int j = 0; j < 9; j++)
//					{
//						cout << puzzle[i][j];
//						cout << " ";
//					}
//
//					cout << endl;
//				}
//			}
//		}
//		else {
//			cout << "Check the input file... looks like there is no solution" << endl;
//		}
//	      delete_lp(lp);
//	}
//
//
//      int main(int argc, char* argv[])
//{
//	// formulate the non-input related part of the puzzle
//	set_sudoku_ilp();
//
//	// read the incomplete input information, and set appropriate constraints
//	read_input_data(argv);
//
//	// solve the puzzle and print the solution
//	solve_the_puzzle();
//
//	return(0);
//}
