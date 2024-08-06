#include "Shooting.h"

#include <iostream> 
#include <fstream> 
#include <vector>
#include <complex>


#define pause 	system("pause");
using std::cout;
using std::endl;




int main()
{
	double Ra_crit = 1e+10, Real_crit = 0, Im_crit = 0, k_crit = 0;
	double Im = 0;
	double Ra = 1000;


	Shooting_Method shoot;
	
	shoot.Le = 130;
	shoot.recalculate_vibro_coefs(90, 90, 0, 0.00); 
	shoot.k_ratio = 1;


	//single computation
	//shoot.k = 3;
	//shoot.solve_with_secant_method(Ra, Im);
	//cout << "result = " << Ra << " " << Im << endl;


	double re, im;
	shoot.Ra = 2000;
	shoot.k = 2.0;

	shoot.Pr = 3;


	
	ofstream w("spectrum.dat"); 
	w << "Ra, re, Im" << endl;
	for (double r = -300; r <= 1000; r = r + 1)
	{
		shoot.Ra = r;

		re = 0;
		im = 0.;
		int error = shoot.solve_with_secant_method(re, im);
		if (!error)
		{
			w << r << " " << -re << " " << im << endl;
		}
		cout << error << " " << r << " " << re << " " << im << endl;
		
	}
	













	auto neutral_curve = [&shoot](double k_start, double k_end, double dk = 0.1)
	{
		cout << "neutral curve" << endl;
		ofstream w("neutral_curve.dat");
		w << "k, Ra, Im" << endl;
		double k = k_start;
		int stop = 0;
		while (!stop)
		{
			if (k <= 0.9) k = k - 0.1 * dk; else 
				k = k - dk;
			
			if (k < k_end + 1e-10) { k = 0.001; stop = 1; }
			shoot.k = k;

			double Im = 0, Ra = 1000;
			int err = shoot.solve_with_secant_method(Ra, Im);

			if (!err)
				w << k << " " << Ra << " " << Im << endl;
			cout << err << ", k = " << k << ", Ra = " << Ra << endl;
		}
	};
	//neutral_curve(4, 1);



	auto tabulate = [&shoot](double R_start, double R_end, double dR)
	{
		cout << "tabulating " << endl;
		complex<double> res;
		for (double R = R_start; R < R_end; R = R + dR)
		{
			shoot.ssystem(res, R, 0);
			cout << R << " " << res << endl;
		}
	};
	//tabulate(0, 2000, 10);






	auto find = [&shoot, &Ra_crit, &Im_crit, &k_crit](double k_start, double k_end, bool show, double dk = 0.1)
	{
		double k = k_start;
		Ra_crit = 1e+10, Im_crit = 0, k_crit = k;
		int stop = 0;
		while (!stop)
		{
			k = k - dk;
			if (k < k_end + 1e-10) { k = 0.001; stop = 1; }
			shoot.k = k;

			double Im = 0, Ra = 1000;
			shoot.solve_with_secant_method(Ra, Im);
			{
				if (abs(Ra) < Ra_crit)
				{
					Ra_crit = Ra;
					k_crit = k;
					Im_crit = Im;
				}
			}
			//cout << "k = " << k << endl;
		}
		if (show) cout << endl << k_crit << " " << Ra_crit << " " << Im_crit << endl;
	};
	//find(4, 0, true);



	//shoot.k_ratio = 0.5;
	//ofstream w("data.dat");
	//w << "a, k, Ra" << endl;
	//for (double a = 90; a > 0; a = a - 2)
	//{
	//	shoot.recalculate_vibro_coefs(a, 90, 0, 0.00);
	//	find(4, 0, false, 0.1);
	//	cout << a << " " << k_crit / 2 << " " << Ra_crit / 16 << " " << Im_crit << endl;
	//	w << a << " " << k_crit / 2 << " " << Ra_crit / 16 << endl;
	//	w << a << " " << k_crit  << " " << Ra_crit  << endl;
	//}

	//cout << endl << k_crit << " " << Ra_crit << " " << Im << endl;

}

