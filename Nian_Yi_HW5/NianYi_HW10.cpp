//American & European option

#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <cstdlib>
#include <vector>
using namespace std;

double risk_free_rate, strike_price, up_factor, uptick_prob, down_prob;
double initial_stock_price, expiration_time, volatility, R;

int no_of_divisions;
vector<vector<float>> V;
vector<vector<float>> M;
vector<vector<float>> EURO_CALL;
vector<vector<float>> EURO_PUT;




double max(double a, double b) {
	return (b < a) ? a : b;
}

double american_call_option(int k, int i, double current_stock_price) {

	if (k == no_of_divisions)
	{
		V[k][no_of_divisions + i] = max(0.0, (current_stock_price - strike_price));
		return max(0.0, (current_stock_price - strike_price));
	}
	else if (V[k][no_of_divisions + i] >= 0.0)
	{
		return V[k][no_of_divisions + i];
	}
	else
	{
		V[k][no_of_divisions + i] = max((current_stock_price - strike_price),
			(uptick_prob*american_call_option(k + 1, i + 1, current_stock_price*up_factor) +
			down_prob*american_call_option(k + 1, i - 1, current_stock_price / up_factor) +
			(1 - uptick_prob - down_prob)*american_call_option(k + 1, i, current_stock_price)) / R);
		return V[k][no_of_divisions + i];
	}
}

double american_put_option(int k, int i, double current_stock_price) {
	if (k == no_of_divisions)
	{
		M[k][no_of_divisions + i] = max(0.0, (strike_price - current_stock_price));
		return max(0.0, (strike_price - current_stock_price));
	}
	else if (M[k][no_of_divisions + i] >= 0.0)
	{
		return M[k][no_of_divisions + i];
	}
	else
	{
		M[k][no_of_divisions + i] = max((strike_price - current_stock_price),
			(uptick_prob*american_put_option(k + 1, i + 1, current_stock_price*up_factor) +
			down_prob*american_put_option(k + 1, i - 1, current_stock_price / up_factor) +
			(1 - uptick_prob - down_prob)*american_put_option(k + 1, i, current_stock_price)) / R);
		return M[k][no_of_divisions + i];
	}

}

/*****************************************************************************/
//European
double european_call_option(int k, int i) {
	if (k == no_of_divisions)
	{
		EURO_CALL[k][no_of_divisions + i] = max(0.0, (initial_stock_price*pow(up_factor, ((float)i)) - strike_price));
		return max(0.0, (initial_stock_price*pow(up_factor, ((float)i)) - strike_price));
	}
	else if (EURO_CALL[k][no_of_divisions + i] >= 0.0)
	{
		return EURO_CALL[k][no_of_divisions + i];
	}
	else
	{
		EURO_CALL[k][no_of_divisions + i] = (uptick_prob*european_call_option(k + 1, i + 1) +
			down_prob*european_call_option(k + 1, i - 1) +
			(1 - uptick_prob - down_prob)*european_call_option(k + 1, i)) / R;


		return EURO_CALL[k][no_of_divisions + i];
	}
}
double european_put_option(int k, int i) {
	if (k == no_of_divisions)
	{
		EURO_PUT[k][no_of_divisions + i] = max(0.0, (strike_price - (initial_stock_price*pow(up_factor, ((float)i)))));
		return max(0.0, (strike_price - (initial_stock_price*pow(up_factor, ((float)i)))));
	}
	else if (EURO_PUT[k][no_of_divisions + i] >= 0.0)
	{
		return EURO_PUT[k][no_of_divisions + i];
	}
	else
	{
		EURO_PUT[k][no_of_divisions + i] = (uptick_prob*european_put_option(k + 1, i + 1) +
			down_prob*european_put_option(k + 1, i - 1) +
			(1 - uptick_prob - down_prob)*european_put_option(k + 1, i)) / R;
		return EURO_PUT[k][no_of_divisions + i];

	}
}


/*****************************************************************************/

int main(int argc, char* argv[])
{

	sscanf_s(argv[1], "%lf", &expiration_time);
	sscanf_s(argv[2], "%d", &no_of_divisions);
	sscanf_s(argv[3], "%lf", &risk_free_rate);
	sscanf_s(argv[4], "%lf", &volatility);
	sscanf_s(argv[5], "%lf", &initial_stock_price);
	sscanf_s(argv[6], "%lf", &strike_price);

	int dimension = 2 * no_of_divisions + 1;
	vector<float> C;
	vector<float> D;
	vector<float> E;
	vector<float> F;

	for (int i = 0; i <= no_of_divisions; i++)
	{
		for (int j = 0; j <= dimension; j++)
		{
			C.push_back(-1.0);
		}
		V.push_back(C);
		C.clear();
	}

	for (int i = 0; i <= no_of_divisions; i++)
	{
		for (int j = 0; j <= dimension; j++)
		{
			D.push_back(-1.0);
		}
		M.push_back(D);
		D.clear();
	}

	for (int i = 0; i <= no_of_divisions; i++)
	{
		for (int j = 0; j <= dimension; j++)
		{
			E.push_back(-1.0);
		}
		EURO_CALL.push_back(E);
		E.clear();
	}

	for (int i = 0; i <= no_of_divisions; i++)
	{
		for (int j = 0; j <= dimension; j++)
		{
			F.push_back(-1.0);
		}
		EURO_PUT.push_back(F);
		F.clear();
	}

	/****************************************************************************/
	up_factor = exp(volatility*sqrt(2 * expiration_time / ((float)no_of_divisions)));
	R = exp(risk_free_rate*expiration_time / ((float)no_of_divisions));
	uptick_prob = pow((sqrt(R) - (1 / sqrt(up_factor))) / (sqrt(up_factor) - (1 / sqrt(up_factor))), 2);
	down_prob = pow((sqrt(up_factor) - sqrt(R)) / (sqrt(up_factor) - (1 / sqrt(up_factor))), 2);
	cout << "Recursive Binomial American-Asian Option Pricing" << endl;
	cout << "Expiration Time (Years) = " << expiration_time << endl;
	cout << "Number of Divisions = " << no_of_divisions << endl;
	cout << "Risk Free Interest Rate = " << risk_free_rate << endl;
	cout << "Volatility (%age of stock value) = " << volatility * 100 << endl;
	cout << "Initial Stock Price = " << initial_stock_price << endl;
	cout << "Strike Price = " << strike_price << endl;
	cout << "--------------------------------------" << endl;
	cout << "R = " << R << endl;
	cout << "Up Factor = " << up_factor << endl;
	cout << "Up tick Probability = " << uptick_prob << endl;
	cout << "Down tick Probability = " << down_prob << endl;
	cout << "No tick Probability = " << (1 - uptick_prob - down_prob) << endl;
	cout << "--------------------------------------" << endl;
	float call_price = american_call_option(0, 0, initial_stock_price);
	cout << "trinomial Price of an American Call Option = " << call_price << endl;
	float put_price = american_put_option(0, 0, initial_stock_price);
	cout << "trinomial Price of an American Put Option = " << put_price << endl;
	cout << "--------------------------------------" << endl;
	float euro_call_price = european_call_option(0, 0);
	cout << "trinomial Price of an European Call Option = " << euro_call_price << endl;
	float euro_put_price = european_put_option(0, 0);
	cout << "trinomial Price of an European put Option = " << euro_put_price << endl;
	cout << "--------------------------------------" << endl;
	cout << "Let us verify the Put-Call Parity: S+P-C = Kexp(-r*T) for American Options" << endl;
	cout << initial_stock_price << " + " << put_price << " - " << call_price;
	cout << " = " << strike_price << "exp(-" << risk_free_rate << " * " << expiration_time << ")" << endl;
	cout << initial_stock_price + put_price - call_price << " ?=? " << strike_price*exp(-risk_free_rate*expiration_time) << endl;
	if (abs(initial_stock_price + put_price - call_price - strike_price*exp(-risk_free_rate*expiration_time)) <= 1e-2)
		cout << "Looks like Put-Call Parity holds within three decimal places" << endl;
	else
		cout << "Looks like Put-Call Parity does NOT hold" << endl;
	cout << "--------------------------------------" << endl;

}
