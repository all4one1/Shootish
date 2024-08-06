#pragma once
#include <iostream>
#include <cmath>
#include <complex>
#include <vector>
using namespace std;

//#define PRINTINFO 





template <typename type>
void allocation2(type*** f, int nx, int ny)
{
	*(f) = (type**)malloc(sizeof(type) * nx);
	for (int i = 0; i < nx; i++)
		(*f)[i] = (type*)malloc(sizeof(type) * ny);
}
template <typename type>
void allocation1(type** f, int N)
{
	*f = (type*)malloc(sizeof(type) * N);
	for (int l = 0; l < N; l++)
		(*f)[l] = (type)0.0;
};

struct Shooting_Method
{
#define loop(i, n1, n2) for (int i = n1; i < n2; i++)
	typedef complex<double> cd;
	typedef complex<int> ci;
	double pi = 3.1415926535897932384626433832795;
	double A1, A2, A3, A4, A5, Rav, K, Le, V0_prime;
	double Pr, psi1, psi2, psiS, SC1, SC2;
	double alpha, beta, gamma;
	double cosA, cosB, cosG, sinA, sinB, sinG;
	double Ra = 0, Im = 0, k = 0, k_ratio = 1;
	cd **y, **k1, **k2, **k3, **k4, **k5, **e, f, ld;
	cd **a, **at, *b, *bt, e0, sume, sumy, **dy;
	int n = 0, s = 0;
	int m = 1000, z = 0, *l;
	bool doubleSizeLayer = false;
	long long int u = 0;
	double eps = 1e-6;
	vector<int> bc;
	vector<int> unity;

	struct SecantMethod
	{
		#define NSec 3
		std::complex<double> *d, **W, **wt, *bw, *btw, d4, pp;
		std::complex<double> **a, **at, *b, *bt;
		int n = 3;
		int z = 0, *l, s = 0, NN = 0, s1 = 0, s2 = 0;
		double *x, *y, delx = 0, dely = 0, dd = 0;
		double e = 1e-6, esmall = 1e-10, ebig = 1e+10;
		Shooting_Method* parent = nullptr;
		SecantMethod() {};
		SecantMethod(Shooting_Method* p) : parent(p)
		{
			allocation1<complex<double>>(&d, n);
			allocation2<complex<double>>(&W, n, n);
			allocation2<complex<double>>(&wt, n, n);
			allocation1<complex<double>>(&bw, n);
			allocation1<complex<double>>(&btw, n);

			allocation2<complex<double>>(&a, n, n);
			allocation2<complex<double>>(&at, n, n);
			allocation1<complex<double>>(&b, n);
			allocation1<complex<double>>(&bt, n);

			allocation1<double>(&x, n);
			allocation1<double>(&y, n);
			allocation1<int>(&l, n);
		}

		int solve(double &x4, double &y4) {

			s1 = 1;
			NN = 0;
			x[2] = x4;
			y[2] = y4;

			while (true)
			{
			begin:

				double delta = 0.01;
				double base = 0.1;
				if (s1 != 0)
				{
					if (x[2] == 0.0 || std::abs(x[2]) < esmall)
					{
						x[0] = delta;
					}
					else {
						x[0] = (base + delta) * x[2];
					}
					y[0] = y[2];
					x[1] = x[2];
					if (y[2] == 0.0 || std::abs(y[2]) < esmall)
					{
						y[1] = delta;
					}
					else {
						y[1] = (base + delta) * y[2];
					}
				}

				for (int i = 0; i < n; i++)
				{
					parent->ssystem(d[i], x[i], y[i]);
				}

				s1 = 0;
				std::fill_n(bw, n, 0.0);
				W[0][0] = x[0];
				W[0][1] = x[1];
				W[0][2] = x[2];
				W[1][0] = d[0].real();
				W[1][1] = d[1].real();
				W[1][2] = d[2].real();
				W[2][0] = d[0].imag();
				W[2][1] = d[1].imag();
				W[2][2] = d[2].imag();

				parent->tran(n, W, wt, bw, btw, l, z);
				parent->det(n, wt, pp, z);
				delx = std::real(pp);

				std::fill_n(bw, n, 0.0);
				W[0][0] = y[0];
				W[0][1] = y[1];
				W[0][2] = y[2];
				W[1][0] = d[0].real();
				W[1][1] = d[1].real();
				W[1][2] = d[2].real();
				W[2][0] = d[0].imag();
				W[2][1] = d[1].imag();
				W[2][2] = d[2].imag();
				parent->tran(n, W, wt, bw, btw, l, z);
				parent->det(n, wt, pp, z);
				dely = std::real(pp);

				std::fill_n(bw, n, 0.0);
				W[0][0] = 1.0;
				W[0][1] = 1.0;
				W[0][2] = 1.0;
				W[1][0] = d[0].real();
				W[1][1] = d[1].real();
				W[1][2] = d[2].real();
				W[2][0] = d[0].imag();
				W[2][1] = d[1].imag();
				W[2][2] = d[2].imag();
				parent->tran(n, W, wt, bw, btw, l, z);
				parent->det(n, wt, pp, z);
				dd = std::real(pp);

				if (esmall * std::abs(delx) < ebig * std::abs(dd) && esmall * std::abs(dely) < ebig * std::abs(dd)) {
					x4 = delx / dd;
					y4 = dely / dd;
				}
				else {
					s1 = 1;
					#ifdef PRINTINFO
						PRINTINFO std::cout << "reinitialization" << std::endl;
					#endif 
					goto begin;
				}

				//std::cout << "Ra " << std::abs(x4 - x[2]) << ", " << (std::abs(x4 - x[2])) << std::endl;
				//std::cout << "Im " << std::abs(y4 - y[2]) << ", " << (std::abs(y4 - y[2])) << std::endl;

				if (std::abs(x4 - x[2]) < e && std::abs(y4 - y[2]) < e)
				{
					return 0;
				}

				else
				{
					parent->ssystem(d4, x4, y4);
					NN++;
					
					#ifdef PRINTINFO
						std::cout << "iterations = " << NN << std::endl;
					#endif // PRINTINFO

				
					if (NN > 100  || abs(x4) > 15000 || abs(y4) > 100)
					{
						#ifdef PRINTINFO
							std::cout << "over looped" << std::endl;
							std::cout << x4 << " " << y4 << std::endl;
						#endif // PRINTINFO

						return 1;
					}
					//std::cout << "d4 = " << d4 << std::endl;
					//std::cout << x4 << ", " << y4 << std::endl;
					x[0] = x[1]; x[1] = x[2]; x[2] = x4;
					y[0] = y[1]; y[1] = y[2]; y[2] = y4;
					d[0] = d[1]; d[1] = d[2]; d[2] = d4;
				}

			}
			return 2;
		}
	};
	SecantMethod sec;


	Shooting_Method()
	{

		n = 6; 
		s = n / 2;

		// 0: w, 1: w', 4: T, 6: T'', 8: C, 9: C'
		//bc = { 0, 1, 4, 6, 9};

		bc = { 0, 2, 4 };

		for (int i = 0; i < n; i++)
			if (std::find(bc.begin(), bc.end(), i) == bc.end())
				unity.push_back(i);



		Pr = 7;
		psi1 = 0.2;
		psi2 = -0.4;
		psiS = psi1+psi2;
		//psi2 = psiS - psi1;
		SC2 = 1000;
		SC1 = 100;

		Le = 130;

		k = 0.001;

		alpha = 90;
		gamma = 90;
		beta = 90;
		K = 0.0;
		Rav = 0.0;
		recalculate_vibro_coefs(alpha, beta, Rav, K);


		allocation2<complex<double>>(&y, s, n);
		allocation2<complex<double>>(&k1, s, n);
		allocation2<complex<double>>(&k2, s, n);
		allocation2<complex<double>>(&k3, s, n);
		allocation2<complex<double>>(&k4, s, n);
		allocation2<complex<double>>(&k5, s, n);
		allocation2<complex<double>>(&e, s, n);
		allocation2<complex<double>>(&dy, s, n);
		allocation2<complex<double>>(&a, s, s);
		allocation2<complex<double>>(&at, s, s);
		allocation1<complex<double>>(&b, s);
		allocation1<complex<double>>(&bt, s);
		allocation1<int>(&l, n);

		sec = SecantMethod(this);
	}

	void recalculate_vibro_coefs(double a, double b, double Rav_, double K_)
	{
		alpha = a;
		beta = b;
		gamma = alpha;
		Rav = Rav_;
		K = K_;
		auto trig = [this](double angle, double (*func)(double))
		{	
			return std::floor(func(angle * pi / 180.0) * 1e+10) / 1e+10;	
		};

		cosA = trig(alpha, cos);
		cosB = trig(beta,  cos);
		cosG = trig(gamma, cos);
		sinA = trig(alpha, sin);
		sinB = trig(beta,  sin);
		sinG = trig(gamma, sin);


		A1 = Rav * (1.0 + K) * trig(gamma, cos) * trig(beta, cos);
		A2 = Rav * (1.0 + K) * trig(gamma, cos) * trig(beta, sin);
		A3 = Rav * (1.0 + K) * trig(gamma, sin) * trig(beta, cos);
		A4 = Rav * (1.0 + K) * trig(gamma, sin) * trig(beta, sin);
		V0_prime = (1.0 + K) * trig(beta - gamma, sin);
		A5 = Rav * V0_prime * cosB;
	}


	//void equation(cd ld, double Ra, int i, int n, cd** y, int j, double x, cd& f)
	//{
	//	//y[0] = w, y[1] = w', y[2] = w'', y[3] = w'''
	//	#define TT (y[j][4])
	//	#define TTprime (y[j][5])
	//	#define TTprime2 (y[j][6])

	//	#define CC (y[j][8])
	//	#define CCprime (y[j][9])
	//	#define CCprime2 (y[j][10])

	//	#define FF (TT + K * CC)
	//	#define FFprime (TTprime + K * CCprime)

	//	#define wz (y[j][0])
	//	#define wzprime (y[j][1])
	//	#define wzprime2 (y[j][2])
	//	#define wzprime3 (y[j][3])

	//	#define WZ y[j][8]
	//	#define WZprime y[j][9]
	//	#define Wprime2 (pow(k, 2) * y[j][8] - pow(k, 2) * FF * sinA - I * k * k_ratio * FFprime * cosB)
	//	#define V0(x) (V0_prime * (x - 0.5))
	//	#define VX ((1.0 - pow(k_ratio, 2)) * cosB * FF + I * k_ratio / k * WZprime)
	//	#define VXprime ((1.0 - pow(k_ratio, 2)) * cosB * FFprime + I * k_ratio / k * WZprime2)

	//	#define SOURCE (Ra*cosA*FF + A5 * (x - 0.5) * I * k * k_ratio * FF - cosG * (1.0 + K) * Rav - A1 * VX - A2 * WZ)



	//	std::complex<double> I = { 0, 1 };
	//	f = { 0,0 };

	//	switch (i) {
	//	case 0:			 // w
	//		f = y[j][1]; // w' = w''
	//		break;
	//	case 1:			 // w'
	//		f = y[j][2]; // w''
	//		break;
	//	case 2:		  	 // w''
	//		f = y[j][3]; // w'''
	//		break;
	//	case 3:			// w''''
	//		f =  (ld / Pr) * (wzprime2 - pow(k, 2) * wz)  // ∇^2 * λ / Pr
	//			+ 2 * pow(k, 2) * wzprime2 - pow(k, 4) * wz   // rest of ∇^4
	//			+ pow(k, 2) * Ra * sinA * FF
	//			+ I * k * k_ratio * Ra * cosA * FFprime
	//			- A1 * I * k * k_ratio * (1.0 - pow(k_ratio, 2)) * cosB * FFprime
	//			+ A1 * pow(k_ratio, 2) * Wprime2
	//			- A2 * I * k * k_ratio * y[j][9]
	//			- A3 * pow(k_ratio, 2) * (1 - pow(k_ratio, 2)) * cosB * FF
	//			- A3 * I * k * k_ratio * y[j][9]
	//			- A4 * pow(k, 2) * y[j][8]
	//			- A5 * pow(k, 2) * FF;
	//		break;
	//	case 4:		// T
	//		f = y[j][5];
	//		break;
	//	case 5:	   // T'	
	//		f = y[j][6];
	//		break;
	//	case 6:	   // T''	
	//		f = y[j][7];
	//		break;
	//	case 7:	   // T'''
	//		f = 2 * pow(k, 2) * TTprime2 - pow(k, 4) * TT
	//			+ ld * (1.0 + 1.0 / Pr) * (TTprime2 - pow(k, 2) * TT)
	//			- pow(ld, 2) / Pr * TT
	//			- sinG * (wzprime2 - pow(k, 2) * wz)
	//			- cosG * (I * k_ratio / k * wzprime3 - I * k_ratio * k * wzprime - ((1.0 - pow(k_ratio, 2)) * SOURCE));
	//		break;
	//	
	//	case 8:
	//		f = y[j][9];
	//		break;
	//	case 9: 
	//		f = Le * (TTprime2 - pow(k, 2)*TT + ld * CC - ld * TT) + pow(k, 2) * CC;
	//		break;


	//	//case 8:		
	//	//	f = y[j][9];
	//	//	break;
	//	//case 9:	  	
	//	//	f = y[j][10];
	//	//	break;
	//	//case 10:	   	
	//	//	f = y[j][11];
	//	//	break;
	//	//case 11:	   
	//	//	f = (2 * pow(k, 2) * CCprime2 - pow(k, 4) * CC) / Le
	//	//		+ ld * (1.0 + 1.0 / (Pr * Le)) * (CCprime2 - pow(k, 2) * CC)
	//	//		- pow(ld, 2) / Pr * CC
	//	//		- sinG * (wzprime2 - pow(k, 2) * wz)
	//	//		- cosG * (I * k_ratio / k * wzprime3 - I * k_ratio * k * wzprime - ((1.0 - pow(k_ratio, 2)) * SOURCE));
	//	//	f = f * Le;
	//	//	break;
	//	}
	//}

	void equation(cd ld, double Ra, int i, int n, cd** y, int j, double x, cd& f)
	{
		//y[0] = w, y[1] = w', y[2] = w'', y[3] = w'''
		#define TT y[j][4]
		#define TTprime y[j][5]
		#define FF (TT)
		#define FFprime (y[j][5])


		std::complex<double> I = { 0, 1 };
		f = { 0,0 };

		switch (i) {
		case 0:			 // w
			f = y[j][1]; // w' = w''
			break;
		case 1:			 // w'
			f = y[j][2]; // w''
			break;
		case 2:		  	 // w''
			f = y[j][3]; // w'''
			break;
		case 3:			// w''''
			f = (ld / Pr) * (y[j][2] - pow(k, 2) * y[j][0])  // ∇^2 * λ / Pr
				+ 2 * pow(k, 2) * y[j][2] - pow(k, 4) * y[j][0]   // rest of ∇^4
				+ pow(k, 2) * Ra * sinA * FF
				+ I * k * k_ratio * Ra * cosA * FFprime;
			break;
		case 4: // T
			f = y[j][5]; // T'
			break;
		case 5: //T''
			f = (-I * y[j][1] * cosG / k - y[j][0] * sinG) + pow(k, 2) * y[j][4] + ld * y[j][4];
			break;

		}

	}

	void ssystem(cd &d, double x4, double y4)
	{

		loop(q, 0, s)
			b[q] = 0;
		loop(i, 0, s)
			loop(j, 0, n)
				y[i][j] = 0;




		for (int i = 0; i < s; i++)
			y[i][unity[i]] = 1;

		//y[0][2] = 1;
		//y[1][3] = 1;
		//y[2][4] = 1;
		//y[3][7] = 1;
		//y[4][9] = 1;





		double x = 0.0;
		double xn = 1;
		u = 0;

		double h = 0.001;
		double hh = h;
		ld = cd(x4, y4);

		double ra = Ra;


		while (x < xn)
		{
			if (abs(xn - x) < abs(h))
				h = xn - x;
		continue2:
			u = u + 1;
			hh = h;


			//K1

			loop(j, 0, s) {
				loop(i, 0, n) {
					equation(ld, ra, i, n, y, j, x, f);
					k1[j][i] = h / 3 * f;

				}
			}
			x = x + h / 3.0;


			loop(j, 0, s)
				loop(i, 0, n)
				y[j][i] = y[j][i] + k1[j][i];



			// K2

			loop(j, 0, s) {
				loop(i, 0, n) {
					equation(ld, ra, i, n, y, j, x, f);
					k2[j][i] = h / 3 * f;

				}
			}


			loop(j, 0, s)
				loop(i, 0, n)
				y[j][i] = y[j][i] - 0.5*k1[j][i] + 0.5*k2[j][i];





			// K3

			loop(j, 0, s) {
				loop(i, 0, n) {
					equation(ld, ra, i, n, y, j, x, f);
					k3[j][i] = h / 3 * f;
				}
			}

			x = x + h / 6.0;
			loop(j, 0, s)
				loop(i, 0, n)
				y[j][i] = y[j][i] - 0.125*k1[j][i] - 0.5*k2[j][i] + 1.125*k3[j][i];




			// K4
			loop(j, 0, s) {
				loop(i, 0, n) {
					equation(ld, ra, i, n, y, j, x, f);
					k4[j][i] = h / 3 * f;
				}
			}

			x = x + 0.5*h;
			loop(j, 0, s)
				loop(i, 0, n)
				y[j][i] = y[j][i] + 1.125*k1[j][i] - 5.625*k3[j][i] + 6.0*k4[j][i];

			// K5
			loop(j, 0, s) {
				loop(i, 0, n) {
					equation(ld, ra, i, n, y, j, x, f);
					k5[j][i] = h / 3 * f;
				}
			}

			sume = 0.0;
			sumy = 0.0;

			loop(j, 0, s) {
				loop(i, 0, n) {
					e[j][i] = (k1[j][i] + 4.0*k4[j][i] - 9.0 / 2.0*k3[j][i] - k5[j][i] / 2.0) / 5.0;
					dy[j][i] = 0.5*(k1[j][i] + 4.0*k4[j][i] + k5[j][i]);
					y[j][i] = y[j][i] - 1.5*k1[j][i] + 4.5*k3[j][i] - 6.0*k4[j][i];
					y[j][i] = y[j][i] + dy[j][i];
					sume = sume + pow(e[j][i], 2);
					sumy = sumy + pow(y[j][i], 2);
				}
			}
			//loop(j, 0, s)	loop(i, 0, n)		cout << setprecision(12) << j << " " << i << " " << real(y[j][i]) << endl;	pause
			//cout << sumy << " " << sume << endl;

			e0 = sqrt(sume / sumy);
			h = h*pow((abs(e0) / eps + 0.001), -0.2);

			if (abs(e0) > 5.0*eps) {
				#ifdef PRINTINFO
					printf("otkat: step, x = %llu, %f \n", u, x);
				#endif 

				
				x = x - hh;
				loop(j, 0, s)
					loop(i, 0, n)
					y[j][i] = y[j][i] - dy[j][i];
				goto continue2;
			}


			if (u % m == 0) 			
				ort(y, s, n);

			//	cout << setprecision(12) << real(sumy) << " " << real(sume) << " " << h << endl; 		pause

		}




		//a[0][0] = y[0][0];  a[0][1] = y[1][0];  a[0][2] = y[2][0];  a[0][3] = y[3][0];  a[0][4] = y[4][0]; // w = 0
		//a[1][0] = y[0][1];  a[1][1] = y[1][1];  a[1][2] = y[2][1];  a[1][3] = y[3][1];  a[1][4] = y[4][1]; // w' = 0
		//a[2][0] = y[0][5];  a[2][1] = y[1][5];  a[2][2] = y[2][5];  a[2][3] = y[3][5];  a[2][4] = y[4][5]; // C' = 0
		//a[3][0] = y[0][6];  a[3][1] = y[1][6];  a[3][2] = y[2][6];  a[3][3] = y[3][6];  a[3][4] = y[4][6]; // Т = 0
		//a[4][0] = y[0][8];  a[4][1] = y[1][8];  a[4][2] = y[2][8];  a[4][3] = y[3][8];  a[4][4] = y[4][8]; // W = 0


	
		for (int i = 0; i < s; i++)
		{
			for (int j = 0; j < s; j++)
			{
				//cout << i << "," << j << " = "  << j << "," << bc[i] << "    ";
				a[i][j] = y[j][bc[i]];
			} // cout << endl;
		}
		



		tran(s, a, at, b, bt, l, z);
		det(s, at, d, z);
	}

	void det(int n, cd **a, cd &d, int &z)
	{
		d = a[0][0];
		loop(i, 1, n)
		{
			d = d*a[i][i];
		}
		d = d*pow(-1, z);
	}

	void ort(cd **y, int m, int n)
	{

		cd **z, p, t;
		allocation2<complex<double>>(&z, m, n);
		loop(i, 0, m)
			loop(j, 0, n)
			z[i][j] = 0;

		loop(i, 0, n) z[0][i] = y[0][i];


		loop(k, 1, m) {
			loop(i, 0, k - 1) {
				p = 0;
				t = 0;
				loop(j, 0, n) {
					p = p + y[k][j] * z[i][j];
					t = t + z[i][j] * z[i][j];
				}
				loop(q, 0, n) z[k][q] = z[i][q] * p / t + z[k][q];
			}
			loop(q, 0, n) z[k][q] = y[k][q] - z[k][q];
		}

		loop(k, 0, m) {
			p = 0;
			loop(j, 0, n) {
				p = p + pow(z[k][j], 2);
			}
			p = sqrt(p);
			loop(q, 1, n) y[k][q] = z[k][q] / p;
		}
	}

	void tran(int n, cd **a, cd **at, cd *b, cd *bt, int *l, int &z)
	{
		cd m, *per, h, p;
		allocation1<complex<double>>(&per, n);
		int mj, mi, v;

		loop(i, 0, n) {
			loop(j, 0, n) {
				at[i][j] = a[i][j];
			}
			bt[i] = b[i];
		}

		loop(i, 0, n) l[i] = i;
		z = 0;
		loop(k, 0, n - 1) {
			m = at[k][k];
			mj = k;
			mi = k;


			loop(i, k, n) {
				loop(j, k, n) {
					if (abs(at[i][j]) > abs(m))
					{
						m = at[i][j];
						mi = i;
						mj = j;
					}
				}
			}
			if (abs(m) < 1e-80)
			{
				//printf("!!!BAD INITIAL CONDITIONS %f\n", abs(m));
			}
			v = l[k];
			l[k] = l[mj];
			l[mj] = v;


			if (mj != k)
			{
				loop(q, 0, n) per[q] = at[q][k];
				loop(q, 0, n) at[q][k] = at[q][mj];
				loop(q, 0, n) at[q][mj] = per[q];
				z = z + 1;
			}

			if (mi != k)
			{
				loop(q, 0, n) per[q] = at[k][q];
				loop(q, 0, n) at[k][q] = at[mi][q];
				loop(q, 0, n) at[mi][q] = per[q];
				h = bt[k];
				bt[k] = bt[mi];
				bt[mi] = h;
				z = z + 1;
			}

			// !!! k or k+1
			loop(i, k + 1, n) {
				p = at[i][k] / at[k][k];
				loop(j, k, n) {
					at[i][j] = at[i][j] - p*at[k][j];
				}
				bt[i] = bt[i] - p*bt[k];
			}


		}
	}

	void maxx(int n, int s, cd **e, double e0)
	{
		//complex(8) e(s, n)
		e0 = abs(e[0][0]);
		loop(j, 0, s) {
			loop(i, 0, n) {
				if (abs(e[j][i]) > e0)
					e0 = abs(e[j][i]);
			}
		}
	}

	int solve_with_secant_method(double & in_x, double & in_y)
	{
		return sec.solve(in_x, in_y);
	}


};