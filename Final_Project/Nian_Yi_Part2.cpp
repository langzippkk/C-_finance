//#include <iostream>
//#include <iomanip>
//#include <cmath>
//#include <fstream>
//#include <cstdlib>
//#include <math.h>
//#include<random>
//#include "normdist.h" 
//using namespace std;
//
//
//double risk_free_rate, strike_price, barrier_price, initial_stock_price, expiration_time, volatility, R;
//int no_of_trials, no_of_divisions;
//
//
//// define the probability that a price path did not touch the barrier
//double dicrete_prob(double current_stock, double period)
//{
//	double mean = initial_stock_price + period / (double)no_of_divisions * (current_stock - initial_stock_price);
//	double variance = period / (double)no_of_divisions*(expiration_time - expiration_time / ((double)no_of_divisions) * period);
//	return   (1 - N((barrier_price - mean) / sqrt(variance)));
//}
//
//double max(double a, double b) {
//	return (b < a) ? a : b;
//}
//int seed = 12345;
//std::default_random_engine generator(seed);
//// i.i.d. uniform distribution generator, search on google.//
//// enabled the program to run trials > 500//
//double get_uniform()
//{
//	std::uniform_real_distribution <double> distribution(0.0, 1.0);
//	double number = distribution(generator);
//	return (number);
//}
//
//
////using the sample code from Professor R.S in simulations
//
//double N(const double& z) {
//	if (z > 6.0) { return 1.0; }; // this guards against overflow 
//	if (z < -6.0) { return 0.0; };
//	double b1 = 0.31938153;
//	double b2 = -0.356563782;
//	double b3 = 1.781477937;
//	double b4 = -1.821255978;
//	double b5 = 1.330274429;
//	double p = 0.2316419;
//	double c2 = 0.3989423;
//	double a = fabs(z);
//	double t = 1.0 / (1.0 + a*p);
//	double b = c2*exp((-z)*(z / 2.0));
//	double n = ((((b5*t + b4)*t + b3)*t + b2)*t + b1)*t;
//	n = 1.0 - b*n;
//	if (z < 0.0) n = 1.0 - n;
//	return n;
//};
////
////double option_price_put_black_scholes(const double& S,      // spot price
////	const double& K,      // Strike (exercise) price,
////	const double& r,      // interest rate
////	const double& sigma,  // volatility
////	const double& time){
////	double time_sqrt = sqrt(time);
////	double d1 = (log(S / K) + r*time) / (sigma*time_sqrt) + 0.5*sigma*time_sqrt;
////	double d2 = d1 - (sigma*time_sqrt);
////	return K*exp(-r*time)*N(-d2) - S*N(-d1);
////};
////
////double option_price_call_black_scholes(const double& S,       // spot (underlying) price
////	const double& K,       // strike (exercise) price,
////	const double& r,       // interest rate
////	const double& sigma,   // volatility 
////	const double& time) {  // time to maturity 
////	double time_sqrt = sqrt(time);
////	double d1 = (log(S / K) + r*time) / (sigma*time_sqrt) + 0.5*sigma*time_sqrt;
////	double d2 = d1 - (sigma*time_sqrt);
////	return S*N(d1) - K*exp(-r*time)*N(d2);
////};
////
//////Use Professor R.S's sample code from lesson7 and change it to double;
////double closed_form_down_and_out_european_call_option()
////{
////	// I took this formula from Wilmott, Howison and Dewynne, "The Mathematics of Financial Derivatives"
////	double K = (2 * risk_free_rate) / (volatility*volatility);
////	double A = option_price_call_black_scholes(initial_stock_price, strike_price,
////		risk_free_rate, volatility, expiration_time);
////	double B = (barrier_price*barrier_price) / initial_stock_price;
////	double C = pow(initial_stock_price / barrier_price, -(K - 1));
////	double D = option_price_call_black_scholes(B, strike_price, risk_free_rate, volatility, expiration_time);
////	return (A - D*C);
////}
////
////
////double closed_form_down_and_in_european_put_option()
////{
////	// just making it easier by renaming the global variables locally
////	double S = initial_stock_price;
////	double r = risk_free_rate;
////	double T = expiration_time;
////	double sigma = volatility;
////	double H = barrier_price;
////	double X = strike_price;
////
////	// Took these formulae from some online reference
////	double lambda = (r + ((sigma*sigma) / 2)) / (sigma*sigma);
////	double temp = 2 * lambda - 2.0;
////	double x1 = (log(S / H) / (sigma*sqrt(T))) + (lambda*sigma*sqrt(T));
////	double y = (log(H*H / (S*X)) / (sigma*sqrt(T))) + (lambda*sigma*sqrt(T));
////	double y1 = (log(H / S) / (sigma*sqrt(T))) + (lambda*sigma*sqrt(T));
////	return (-S*N(-x1) + X*exp(-r*T)*N(-x1 + sigma*sqrt(T)) +
////		S*pow(H / S, 2 * lambda)*(N(y) - N(y1)) -
////		X*exp(-r*T)*pow(H / S, temp)*(N(y - sigma*sqrt(T)) - N(y1 - sigma*sqrt(T))));
////}
////
////double closed_form_down_and_out_european_put_option()
////{
////	double vanilla_put = option_price_put_black_scholes(initial_stock_price, strike_price,
////		risk_free_rate, volatility, expiration_time);
////	double put_down_in = closed_form_down_and_in_european_put_option();
////	return (vanilla_put - put_down_in);
////}
////
//
//int main(int argc, char* argv[])
//{
//	sscanf_s(argv[1], "%lf", &expiration_time);
//	sscanf_s(argv[2], "%lf", &risk_free_rate);
//	sscanf_s(argv[3], "%lf", &volatility);
//	sscanf_s(argv[4], "%lf", &initial_stock_price);
//	sscanf_s(argv[5], "%lf", &strike_price);
//	sscanf_s(argv[6], "%d", &no_of_trials);
//	sscanf_s(argv[7], "%d", &no_of_divisions);
//	sscanf_s(argv[8], "%lf", &barrier_price);
//
//
//	//***********************************************Simulation method**********************************************************//
//	//                         revise professor R.S sample code for simulations of European options(multiperiod)//
//	double t = expiration_time / ((int)no_of_divisions);
//	double r = (risk_free_rate - 0.5*pow(volatility, 2))*t;
//	double sd = volatility*sqrt(t);
//	// I need to cut-off the path in many segment so that I can determine it touches the barrier or not;
//
//	double put_option_price = 0.0;
//	double call_option_price = 0.0;
//
//
//
//	for (int i = 0; i < no_of_trials; i++){
//		double path1 = 1;
//		double path2 = 1;
//		double path3 = 1;
//		double path4 = 1;
//		// initialize S1....S4 for each trial;
//		double S1 = initial_stock_price;
//		double S2 = initial_stock_price;
//		double S3 = initial_stock_price;
//		double S4 = initial_stock_price;
//		// define the direct probability method's variable: if the barrier is above the stock price, then 
//		// the down and out option have no value(1-p=0)'
//
//
//		for (int j = 0; j < no_of_divisions; j++) {
//
//			// generate unit-normals using Box-Muller Transform
//			double x = get_uniform();
//			double y = get_uniform();
//			double a = sqrt(-2.0*log(x)) * cos(6.283185307999998*y);
//			double b = sqrt(-2.0*log(x)) * sin(6.283185307999998*y);
//
//
//			S1 = S1 * exp(r + sd*a);
//			S2 = S2 * exp(r - sd*a);
//			S3 = S3* exp(r + sd*b);
//			S4 = S4 * exp(r - sd*b);
//			if (S1 <= barrier_price)
//				path1 = 0;
//			if (S2 <= barrier_price)
//				path2 = 0;
//			if (S3 <= barrier_price)
//				path3 = 0;
//			if (S4 <= barrier_price)
//				path4 = 0;
//			/*	cout << S1<<endl;
//			cout << path1;*/
//		}
//		//end inside for loop;
//		if (path1 > 0)
//		{
//			put_option_price += max(0.0, strike_price - S1);
//			call_option_price += max(0.0, S1 - strike_price);
//
//		}
//		if (path2 > 0)
//		{
//			put_option_price += max(0.0, strike_price - S2);
//			call_option_price += max(0.0, S2 - strike_price);
//		}
//		if (path3 > 0)
//		{
//			put_option_price += max(0.0, strike_price - S3);
//			call_option_price += max(0.0, S3 - strike_price);
//		}
//		if (path4 > 0)
//		{
//			put_option_price += max(0.0, strike_price - S4);
//			call_option_price += max(0.0, S4 - strike_price);
//		}
//	}
//	call_option_price = exp(-risk_free_rate*expiration_time)*((call_option_price / 4.0) / ((double)no_of_trials));
//	put_option_price = exp(-risk_free_rate*expiration_time)*((put_option_price / 4.0) / ((double)no_of_trials));
//
//	// for simulation method : divided by 4 and discount it back and divide no_of_trials at last//
//	// up to here, the call option price add all the non-zero payoff of 4 paths for no_of_division times;
//	//**************************************************************************************//
//
//	//                                   direct probability method                          //
//	// define the probabilities: if the barrier is above the stock price, then 
//	// the down and out option have no value(1-p=0)'
//	double T = expiration_time;
//	double R = (risk_free_rate - 0.5*pow(volatility, 2))*T;
//	double SD = volatility*sqrt(T);
//
//	double call_probability = 0; double put_probability = 0;
//
//	for (int i = 0; i < no_of_trials; i++){
//		double p1 = 1.0;
//		double p2 = 1.0;
//		double p3 = 1.0;
//		double p4 = 1.0;
//
//		// initialize S1....S4 for each trial;
//		double S1P = initial_stock_price;
//		double S2P = initial_stock_price;
//		double S3P = initial_stock_price;
//		double S4P = initial_stock_price;
//		// define the direct probability method's variable: if the barrier is above the stock price, then 
//		// the down and out option have no value(1-p=0)'
//
//		// generate unit-normals using Box-Muller Transform
//		double x1 = get_uniform();
//		double y1 = get_uniform();
//		double a1 = sqrt(-2.0*log(x1)) * cos(6.283185307999998*y1);
//		double b1 = sqrt(-2.0*log(x1)) * sin(6.283185307999998*y1);
//
//
//		S1P = S1P * exp(R + SD*a1);
//		S2P = S2P * exp(R - SD*a1);
//		S3P = S3P* exp(R + SD*b1);
//		S4P = S4P * exp(R - SD*b1);
//
//		// compute the probability for a path that did not touch the barrier//
//
//		for (int j = 1; j <= no_of_divisions; j++)
//		{
//			if (S1P > barrier_price && initial_stock_price > barrier_price)
//				p1 = p1*dicrete_prob(S1P, j);
//			else
//				p1 = p1 * 0;
//			if (S2P > barrier_price && initial_stock_price > barrier_price)
//				p2 = p2*dicrete_prob(S2P, j);
//			else
//				p2 = p2 * 0;
//			if (S3P > barrier_price && initial_stock_price > barrier_price)
//				p3 = p3*dicrete_prob(S3P, j);
//			else
//				p3 = p3 * 0;
//			if (S4P > barrier_price&& initial_stock_price > barrier_price)
//				p4 = p4*dicrete_prob(S4P, j);
//			else
//				p4 = p4 * 0;
//		}
//
//		/*	for (int j = 1; j <= no_of_divisions; j++){
//		p1 = p1*dicrete_prob(S1P, j);
//		p2 = p2*dicrete_prob(S2P, j);
//		p3 = p3*dicrete_prob(S3P, j);
//		p4 = p4*dicrete_prob(S4P, j);
//
//		}*/
//
//		call_probability += (max(0.0, S1P - strike_price)*(p1)+
//			max(0.0, S2P - strike_price)* (p2)+
//			max(0.0, S3P - strike_price)*(p3)+
//			max(0.0, S4P - strike_price)*(p4)) / 4.0;
//		put_probability += (max(0.0, strike_price - S1P)*(p1)+
//			max(0.0, strike_price - S2P)* (p2)+
//			max(0.0, strike_price - S3P)*(p3)+
//			max(0.0, strike_price - S4P)*(p4)) / 4.0;
//		//	cout << call_probability<<endl;
//	}
//
//	call_probability = exp(-risk_free_rate*expiration_time)*((call_probability) / ((double)no_of_trials));
//	put_probability = exp(-risk_free_rate*expiration_time)*((put_probability) / ((double)no_of_trials));
//
//
//	// for probability method : get the mean of all the trials                                 //
//	//******************************************************************************************************************//
//
//	cout << "European Down-and-Out Discrete Barrier Options Pricing via Monte Carlo Simulation " << endl;
//	cout << "Expiration Time (Years) = " << expiration_time << endl;
//	cout << "Risk Free Interest Rate = " << risk_free_rate << endl;
//	cout << "Volatility (%age of stock value) = " << volatility * 100 << endl;
//	cout << "Initial Stock Price = " << initial_stock_price << endl;
//	cout << "Strike Price = " << strike_price << endl;
//	cout << "Number of trials = " << no_of_trials << endl;
//	cout << "Number of Discrete Barriers = " << no_of_divisions << endl;
//	cout << "barrier price is = " << barrier_price << endl;
//	cout << "--------------------------------------" << endl;
//	cout << "--------------------------------------" << endl;
//	cout << "The average Call Price via explicit simulation of price paths             = " << call_option_price << endl;
//	cout << "The average Call Price with Brownian Bridge correction on the final price = " << call_probability << endl;
//	cout << "The average Put Price via explicit simulation of price paths              = " << put_option_price << endl;
//	cout << "The average Put Price with Brownian Bridge correction on the final price  = " << put_probability << endl;
//	cout << "--------------------------------------" << endl;
//}