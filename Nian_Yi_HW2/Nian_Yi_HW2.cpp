#include <iostream>
#include "Nian_Yi_HW2.h"

int main(int argc, char * const argv[])
{
	Sudoku x;
	x.read_puzzle(argc, argv);
	x.print_puzzle();
	x.Solve(0, 0);
	//x.alternate_Solve(0, 0);
	//x.print_puzzle();

	return 0;
}
