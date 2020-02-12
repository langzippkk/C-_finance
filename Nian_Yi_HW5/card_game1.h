//#include <vector>
//#include <iostream>
//#include <algorithm>
//using namespace std;
//
//#ifndef	CARDGAME_H
//#define CARDGAME_H
//using namespace std;
//class Construct
//{
//	//define a constructor
//private:
//	int total_cards;
//	//create a vector to memory the latest value
//	vector<float> V;
//public:
//	// create a constructor to initialize the vector
//	float value_of_game(int r, int b);
//	Construct(int r, int b);
//};
//inline float Construct::value_of_game(int r, int b)
//{
//	int N = r*total_cards + b - 1;
//
//
//	if (0 == r)
//		return ((float)b);
//	if (0 == b)
//		return (0);
//	if (V[N] != 0)
//		return V[N];
//
//
//	float term1 = ((float)r / (r + b)) * value_of_game(r - 1, b);
//	float term2 = ((float)b / (r + b)) * value_of_game(r, b - 1);
//	float result = max((term1 + term2), (float)(b - r));
//
//	V[N] = result;
//	return result;
//}
//Construct::Construct(int r, int b) :total_cards(r + b), V(r*total_cards + b, 0)
//{
//
//}
//#endif