////  main.cpp
////  Stable Marriage Problem
////
//#include <iostream>
//#include <vector>
//#include <fstream>
//#include <sstream>
//
//using namespace std;
//
//class stable_marriage_instance
//{
//	// Private
//	int no_of_couples;
//	vector <vector <int> > Preference_of_men;
//	vector <vector <int> > Preference_of_women;
//	vector <int> match_for_men;
//	vector <int> match_for_women;
//
//	// private member function: checks if anybody is free in boolean "my_array"
//	// returns the index of the first-person who is free in "my_array"
//	// if no one is free it returns a -1.
//	int anybody_free(vector <bool> my_array)
//	{
//		int temp = -1;
//		for (int i = 0; i < no_of_couples; i++)
//		{
//			if (my_array[i] == true)
//			{
//				temp = i;
//				break;
//			}
//		}
//		return temp;
//
//	}
//
//	// private member function: if index1 is ranked higher than index2
//	// in a preference array called "my_array" it returns "true"; otherwise
//	// it returns "false"
//	bool rank_check(vector <int> my_array, int index1, int index2)
//	{
//		int a, b;
//		for (int i = 0; i < my_array.size(); i++)
//		{
//			if (my_array[i] == index1)
//				a = i;
//			if (my_array[i] == index2)
//				b = i;
//		}
//			if (a < b)
//			{
//				return  true;
//			}
//			else
//				return false;
//		}
//
//
//	// fill the necessary code for this function
//
//	// private member function: implements the Gale-Shapley algorithm
//	void Gale_Shapley()
//	{
//		vector <bool> is_man_free;
//		vector <bool> is_woman_free;
//		vector <vector <bool> > has_this_man_proposed_to_this_woman;
//		int man_index, woman_index;
//
//		// initializing everything
//		for (int i = 0; i < no_of_couples; i++)
//		{
//			is_woman_free.push_back(true);
//			is_man_free.push_back(true);
//			match_for_men.push_back(0);
//			match_for_women.push_back(0);
//		}
//
//		// initialize the has_this_man_proposed_to_this_women matrix;
//		for (int i = 0; i < no_of_couples; i++)
//		{
//			vector<bool>solution;
//			has_this_man_proposed_to_this_woman.push_back(solution);
//		}
//		// do the necessary initialization here
//		for (int i = 0; i < no_of_couples; i++)
//		{
//			for (int j = 0; j < no_of_couples; j++)
//			{
//				has_this_man_proposed_to_this_woman[i].push_back(false);
//			}
//		}
//         
//		//push back J times of false's;
//
//
//		// Gale-Shapley Algorithm
//		while ((man_index = anybody_free(is_man_free)) >= 0)
//		{
//			for (int i = 0; i < no_of_couples; i++)
//			{
//				if (is_man_free[i] == true)
//				{
//					man_index = i;
//					break;
//				}
//			}
//		
//			//find the first women, who is unpaired, on the preference list;
//			for (int i = 0; i < no_of_couples; i++)
//			{
//				int w = Preference_of_men[man_index][i];
//				if (has_this_man_proposed_to_this_woman[man_index][w] == false)
//				{
//					woman_index = w;
//					break;
//				}
//			}
//
//			//  inspired by stackflow,Gale-Shapley algorithm;
//			//  https://stackoverflow.com/questions/12858734/c-implementation-of-gale-shapley-algorithm;
//			has_this_man_proposed_to_this_woman[man_index][woman_index] = true;
//
//			if (is_woman_free[woman_index])
//			{
//				is_man_free[man_index] = false;
//				is_woman_free[woman_index] = false;
//				match_for_men[man_index] = woman_index;
//				match_for_women[woman_index] = man_index;
//				//match them
//			}
//			else {
//				// find the competitor use loop;
//
//				/*	int index_competitor;
//					for (int i = 0; i < no_of_couples; i++)
//					{
//					if (match_for_men[i] == woman_index)
//					{
//					index_competitor = i;
//					break;
//					}*/
//				int index_competitor;
//				index_competitor = match_for_women[woman_index];
//				//compare these two man and match the better one;
//
//				if (rank_check(Preference_of_women[woman_index], man_index, index_competitor))
//				{
//					is_man_free[man_index] = false;
//					is_man_free[index_competitor] = true;
//					is_woman_free[woman_index] = false;
//					match_for_men[man_index] = woman_index;
//					match_for_women[woman_index] = man_index;
//				}
//				else
//				{
//					is_man_free[man_index] = true;
//					is_man_free[index_competitor] = false;
//					is_woman_free[woman_index] = false;
//				}
//			}
//		}
//	}
//	
//
//	// private member function: reads data
//	void read_data(int argc, const char * argv[])
//	{
//		// fill the necessary code here.  The first entry in the input
//		// file is the #couples, followed by the preferences of the men
//		// and preferences of the women.  Keep in mind all indices start
//		// from 0.
//		ifstream input_filename(argv[1]);
//		int value_read_from_file;
//		vector<int>P;
//
//		if
//			(input_filename.is_open())
//		{
//			cout << "input File Name:" << argv[1] << endl;
//			while (input_filename >> value_read_from_file)
//			{
//				if ((value_read_from_file >= 0) && (value_read_from_file <= 9))
//				{
//					P.push_back(value_read_from_file);
//				}
//			}
//		}
//		no_of_couples = P[0];
//		cout << endl << "Number of Couples is " << no_of_couples << endl;
//		for (int i = 0; i < no_of_couples; i++)
//		{
//			vector<int>Preference;
//			Preference_of_men.push_back(Preference);
//		}
//		for (int i = 0; i < no_of_couples; i++)
//		{
//			vector<int>Preference;
//			Preference_of_women.push_back(Preference);
//		}
//		static int count = 1;
//		for (int i = 0; i < no_of_couples;i++)
//		{
//			for (int j = 0; j < no_of_couples;j++)
//			{
//				Preference_of_men[i].push_back(P[count]);
//				count++;
//			}
//		}
//		cout << "----------------------" << endl;
//		for (int i = 0; i < no_of_couples;i++)
//		{
//			cout << "(" << i << ")";
//			for (int j = 0; j < no_of_couples;j++)
//			{
//				cout << Preference_of_men[i][j]<<" ";
//			}
//			cout << endl;
//		}
//		cout << endl;
//		cout << "Preference of women" << endl;
//		cout << "-------------------- - " << endl;
//		for (int i = 0; i < no_of_couples; i++)
//		{
//			for (int j = 0; j < no_of_couples; j++)
//			{
//				Preference_of_women[i].push_back(P[count]);
//				count++;
//			}
//		}
//		for (int i = 0; i < no_of_couples; i++)
//		{
//			cout << "(" << i << ")";
//			for (int j = 0; j < no_of_couples; j++)
//			{
//				cout << Preference_of_women[i][j]<<" ";
//			}
//			cout << endl;
//		}
//			}
//
//	// private member function: print solution
//	void print_soln()
//	{
//			// write the appropriate code here
//			cout << endl;
//			cout<<"Matching for men"<<endl;
//			for (int i = 0; i < no_of_couples; i++)
//				cout << "Man: " << i << " -> " << "Woman: " << match_for_men[i] << endl;
//			cout << endl;
//			cout << "Matching for women" << endl;
//			for (int i = 0; i < no_of_couples; i++)
//				cout << "Woman: " << i << " -> " << "man: " << match_for_women[i] << endl;
//		}
//		// write the appropriate code here
//
//
//public:
//	void solve_it(int argc, const char * argv[])
//	{
//		read_data(argc, argv);
//
//		Gale_Shapley();
//
//		print_soln();
//	}
//};
//
//int main(int argc, const char * argv[])
//{
//	stable_marriage_instance x;
//	x.solve_it(argc, argv);
//}
//
