#include "Shooting.h"

#include <iostream> 
#include <fstream> 
#include <vector>
#include <complex>
#include <iomanip>
#include <cmath>

#define pause 	system("pause");
using std::cout;
using std::endl;




int main()
{
	double Pi = acos(-1.0);
	double Ra_crit = 1e+10, Real_crit = 0, Im_crit = 0, k_crit = 0;
	double Im = 0;
	double Ra = 1000;


	Shooting_Method shoot;
	

	shoot.Le = 130;
	shoot.recalculate_vibro_coefs(90, 90, 8000, 0.000); 
	shoot.k = 0;

	////single computation
	auto single = [&shoot](double x, double y)
	{
		double Ra = x;
		double Im = y;
		int err = shoot.solve_with_secant_method(Ra, Im);
		cout << err <<  ", result = " << Ra << " " << Im << endl;
	};


	auto single_kxky = [&shoot](double kx, double ky)
	{
		double Ra = 10;
		double Im = 0.0001;

		shoot.k = sqrt(kx * kx + ky * ky);
		shoot.k_ratio = kx / shoot.k;
		int err = shoot.solve_with_secant_method(Ra, Im);
		cout << err << ", k= " << shoot.k << ", a=" << shoot.k_ratio 
			<< ", kx= " << shoot.k_ratio * shoot.k << ", ky= " << sqrt(1.0 - pow(shoot.k_ratio, 2)) * shoot.k 
			<< ", Ra = " << Ra << ", Im = " << Im << endl;
	};
	auto kxky = [&shoot]()
	{
		double Ra = 100;
		double Im = 0.0001;
		double Ra2 = 100;
		double Im2 = 0.0001;

		double dx = 0.01;
		double dy = 0.01;
		ofstream w("kx_ky.dat");
		w << "k, a, kx, ky, Ra, Im " << endl;
		for (double kx = 0.001; kx < 4; kx+= dx)
			for (double ky = 0.001; ky < 4; ky += dy)
			{
				if (kx < 1 && ky < 0.5 || kx < 0.5) dx = dy = 0.01;
				else dx = dy = 0.1;

				shoot.k = sqrt(kx * kx + ky * ky);
				shoot.k_ratio = kx / shoot.k;
				shoot.recalculate_vibro_coefs(60, 90, 0, 0.001);
				int err = shoot.solve_with_secant_method(Ra, Im);

				//shoot.recalculate_vibro_coefs(60, 90, 0, 0.00);
				int err2 = 0; // shoot.solve_with_secant_method(Ra2, Im2);
				if (err == 0 && err2 == 0)
				{
					w << shoot.k << " " << shoot.k_ratio
						<< " " << shoot.k_ratio * shoot.k << " " << sqrt(1.0 - pow(shoot.k_ratio, 2)) * shoot.k
						<< " " << Ra << " " << Im
						<< " " << Ra2 << " "
						<< " " << ((Ra2 - Ra) / Ra2 > 0.25 ? 1 : 0)
						<< endl;
				}
				cout << "kx = " << kx << ", ky = " << ky << ", err = " << err << endl;
			}


	};

	auto neutral_curve = [&shoot](double k_start, double k_end, double dk = 0.1)
	{
		cout << "neutral curve" << endl;
		ofstream w("neutral_curve.dat");
		w << "k, Ra, Im" << endl;
		double k = k_start + dk;
		int stop = 0;
		while (!stop)
		{
			//if (k <= 0.9) k = k - 0.1 * dk; else 
				k = k - dk;
			
			if (k < 0.001 + 1e-10) { k = 0.001; stop = 1; }
			shoot.k = k;

			double Im = 0.0001, Ra = 1000;
			int err = shoot.solve_with_secant_method(Ra, Im);

			if (!err)
				w << k << " " << Ra << " " << Im << endl;
			cout << err << ", k= " << k << ", kx= " << shoot.k_ratio*k << ", ky= " << sqrt(1.0 - pow(shoot.k_ratio, 2)) * k << ", Ra = " << Ra << ", Im = " << Im << endl;
			
			if (k <= k_end) stop = 1;
		}
	};

	shoot.Le = 130;
	shoot.k_ratio = 0;
	shoot.recalculate_vibro_coefs(20, 90, 0, 0.001);
	
	//shoot.k = 0.01;
	//single(5000, 0.0001);
	
	//neutral_curve(4, 0, 0.1);
	//single_kxky(2 * Pi / 10, 2 * Pi / 5);
	//single_kxky(2 * Pi / 10, 2 * Pi / (5 - 0.05));

	
	//for (double l = 1; l < 3; l = l + 0.1 )		single_kxky(2 * Pi / 10, 2 * Pi / (5 + l));
	//kxky();



	shoot.k = 3.14;
	shoot.recalculate_vibro_coefs(90, 90, 0, 0.00);
	auto tabulate = [&shoot](double R_start, double R_end, double dR)
	{
		cout << "tabulating " << endl;
		complex<double> res;
		for (double R = R_start; R < R_end; R = R + dR)
		{
			shoot.ssystem(res, R, 0);
			cout << fixed << setprecision(2) << R << " " << setprecision(10) << res << endl;
		}
	};
	//tabulate(1000, 3000, 1);



	//shoot.recalculate_vibro_coefs(90, 90, 0, 0.000);
	//shoot.k = 1;


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
				if ((Ra) < Ra_crit)
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
	auto find2 = [&find, &shoot, &Ra_crit, &Im_crit, &k_crit](double k_start, double k_end, bool show, double dk = 0.1)
	{
		double k = k_start;

		find(k_start, k_end, show, dk);
		if (k_crit < 0.01) return;

		while (dk > 0.01)
		{
			k_start = k_crit + dk;
			k_end = k_crit - dk;
			dk = dk * 0.1;
			find(k_start, k_end, show, dk);
			
		}
	};


	



	shoot.k_ratio = 0.95;
	ofstream w("data.dat");
	w << "a, k, Ra" << endl;
	for (double a = 90; a > 89; a = a - 0.01)
	{
		shoot.recalculate_vibro_coefs(a, 90, 0, 0.001);
		find2(4, 0.0, false, 0.1);
		cout << a << " " << k_crit  << " " << Ra_crit  << " " << Im_crit << endl;
		w << a << " " << k_crit  << " " << Ra_crit  << endl;
	}

	



	
	shoot.k = 3.1;
	shoot.recalculate_vibro_coefs(90, 90, 0, -0.1);



	auto input = [&shoot, &single]()
	{
		while (true)
		{
			double x = -1500;
			double y = 0.000;
			cout << "Ra = ";	cin >> x;
			cout << "Im = ";	cin >> y;
			single(x, y);
		}
	};

	shoot.Ra = 1700;
	shoot.spectrum = true;
	//input();
	
	auto spectr = [&shoot]()
	{
		shoot.spectrum = true;
		double re, im;
		ofstream w("spectrum.dat");
		w << "Ra, x, y" << endl;
		re = 1;
		double r = 1500;
		int q = 0;
		for (double r = 2500; r > 1800; r = r - 1)
		{

			shoot.Ra = r;
			re = 0.1;
			im = 0.001;
			int error = shoot.solve_with_secant_method(re, im);
			if (!error)
			{
				if (q % 1 == 0)
					w << r << " " << setprecision(25) << re << " " << im << endl;
			}
			cout << error << " " << setw(4) << r << " " << re << " " << im << endl;
			q++;
		}
	};

	//spectr();












}

