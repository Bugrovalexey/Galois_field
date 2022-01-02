#pragma once
#include "BPF.h"
#include <fstream>
#include <vector>

const int p = 2, t = 4;// GF(p^t)

//const int p = 5, t = 1;// GF(p^t)

void wf1()
{
	mt19937 gen;
	gen.seed(time(0));
	//gen.seed(0);

	/*
	Test(5, 2, f55);
	Get_random_quazyfunc(gen, 5, 2, f55);
	Test(5, 2, f55);
	*/

	//Algebraic_degree(p, n, f2);

	/*
	int * pol = BPF(p, 2, Get_B_matrix_for_func_to_pol(p), f);

	cout << Algebraic_degree(p, n, pol) << endl;

	Get_B_matrix_for_func_to_spectr(3);
	*/


	//complex <double> * spectr = BPF(p, n, Get_B_matrix_for_func_to_spectr(p), f);

	//complex <double> * spectr = BPF(2, 3, Get_B_matrix_for_func_to_spectr(2), chi_func(2, 3, f2));

	//cout << Resilient_index(2, 3,spectr);


	//полсчитаем распределение суммы степени функции и ее индекса устойчисвости в случайной выборке
	int p = 5;
	int n = 2;
	int pn = (int)pow(p, n);
	int size_distr_ri_plus_deg = p * n + 2;
	int * distr_ri_plus_deg = new int[size_distr_ri_plus_deg];

	int size_distr_deg = p*n;
	int * distr_deg = new int[size_distr_deg];

	int size_distr_ri = n + 2;
	int * distr_ri = new int[size_distr_ri];


	for (int i = 0; i < size_distr_ri_plus_deg; i++)
		distr_ri_plus_deg[i] = 0;

	for (int i = 0; i < size_distr_ri; i++)
		distr_ri[i] = 0;

	for (int i = 0; i < size_distr_deg; i++)
		distr_deg[i] = 0;



	int deg, ri;
	int MinPolWeight = pn, *MinPol = new int[pn], *MinF = new int[pn];
	int * f = new int[pn];//выделим память для работы
	int * pol = new int[pn];
	complex <double> * spectr = new complex <double>[pn];
	complex <double> * chi_f = new complex <double>[pn];

	complex <double> nonlin;

	complex< double > ** B_spectr = Get_B_matrix_for_func_to_spectr(p);
	int ** B_polinom = Get_B_matrix_for_func_to_pol(p);

	//while (1)
	{


		for (int i = 0; i < 1000; i++)
		{
			Get_random_func(gen, p, n, f); //просто случайная функция

										   //Get_random_quazyfunc(gen, p, n, f);//случайная квазигруппа

			chi_func(p, n, f, chi_f);

			BPF(p, n, B_spectr, chi_f, spectr);

			nonlin = AbsMax(pn, spectr)*pn;

			if (abs(nonlin) < 53)
				cout << endl;

			//	ri = Resilient_index(p, n, spectr);

			ri = Correalation_immunity_index(p, n, spectr);


			BPF(p, n, B_polinom, f, pol);

			deg = Algebraic_degree(p, n, pol);



			if (deg == 6)
				if (MinPolWeight > WeightArray(pn, pol))
				{
					MinPolWeight = WeightArray(pn, pol);
					memcpy(MinPol, pol, sizeof(int) * pn);
					memcpy(MinF, f, sizeof(int) * pn);
				}



			distr_ri_plus_deg[deg + ri]++;
			distr_ri[ri]++;
			distr_deg[deg]++;


			/*
			if (deg + ri > n)
			{
			for (int j = 0; j < pn; j++)
			cout << f[j] << " ";
			cout << endl;
			}
			*/



		}
		cout << "Ci + deg:" << endl;
		for (int i = 0; i < size_distr_ri_plus_deg; i++)
			cout << i << "\t" << distr_ri_plus_deg[i] << endl;
		cout << endl;


		cout << "Correlation immunity:" << endl;
		for (int i = 0; i < size_distr_ri; i++)
			cout << i << "\t" << distr_ri[i] << endl;
		cout << endl;

		cout << "Degree:" << endl;
		for (int i = 0; i < size_distr_deg; i++)
			cout << i << "\t" << distr_deg[i] << endl;
		cout << endl;

		cout << "p = " << p << endl;
		cout << "n = " << n << endl;

		//cout << "Minimal polinom: ";
		//WritePol(p, n, MinPol);
	}


}

void wf2()
{
	//посчитаем распределение нелинейности функции и ее индекса устойчисвости в случайной выборке
	//для простых полей

	mt19937 gen;
	gen.seed(time(0));
	//gen.seed(0);


	int p = 5;
	int n = 4;
	int pn = (int)pow(p, n);

	double *distr = new double[n+1]; // значению устойчивости будем сопоставлять минимальную нелинейность.
	for (int it = 0; it < n+1; it++)
		distr[it] = pn;
	
	int deg, ri;
	int MinPolWeight = pn, *MinPol = new int[pn], *MinF = new int[pn];
	int * f = new int[pn];//выделим память для работы
	int * pol = new int[pn];
	complex <double> * spectr = new complex <double>[pn];
	complex <double> * chi_f = new complex <double>[pn];

	double nonlin, min=pn;

	complex< double > ** B_spectr = Get_B_matrix_for_func_to_spectr(p);
	int ** B_polinom = Get_B_matrix_for_func_to_pol(p);


	int count = 0;
	while (1)
	{
		
		if (count == 0)
		{
			Get_random_func(gen, p, n, f); //просто случайная функция
			count = 1;
		}
			
		else
		{
			Get_random_quazyfunc(gen, p, n, f);//случайная квазигруппа
			count = 0;
		}
			

			chi_func(p, n, f, chi_f);

			BPF(p, n, B_spectr, chi_f, spectr);

			nonlin = AbsMax(pn, spectr)*pn;

			ri = Resilient_index(p, n, spectr);

			if (distr[ri +1] > nonlin)
			{
				distr[ri +1] = nonlin;

				cout << "Resilent index - min C(f)" << endl;
				for (int it = 0; it < n + 1; it++)
					cout << it - 1 << " - " << distr[it] << endl;
				cout << endl;
			}
				
	}
}

//тест функций над полем GF(2^t)
void wf3()
{
	
	InitGF(p, t);
	int q = (int)pow(p, t); 
	int n = 1, // n - количество переменных у функции
		qn = (int)pow(q, n);

	

	element *f = new element[qn],
		*pol = new element[qn],
		*x = new element[n],
		xi(1, p, t),
		one(0, p, t),
		zero(q-1, p, t),
		ebuf(1, p, t);

	for (int i = 0; i < n; i++)
		x[i] = one;

	int w;



	for (int i = 0; i < qn; i++)//зададим функцию многочленом
	{
		w = i;
		for (int j = 0; j < n; j++)
		{
			x[j].exp = w % q;
			w = w / q;
		}
		
		f[i] = x[0] + x[0].ElPower(3) +  x[0].ElPower(12);

	}

		
		/*
		f[i] = x[2]*x[3] + (x[0].ElPower(q-2) + x[1].ElPower((q - 2)/2) + one);
		f[i] = f[i].ElPower(q - 5);
		*/

		/*
		f[i] = (x[0].ElPower(q-2) + x[1].ElPower((q - 2)/2) + one);
		f[i] = f[i].ElPower(q - 5);
		*/
		
		/*
		ebuf = x[0]+one;
		f[i] = ebuf.ElPower(q-2);
		ebuf = xi*(x[0].ElPower(q-2));
		f[i] = f[i] + ebuf;
		*/

	/*
	//зададим функцию поточечно
	for (int i = 0; i < qn; i++)
	{
		f[i] = zero;
	}
	f[qn-1] = one;
	*/

	Test_2t(q, t, n, f);

	//Convert_exp_func_on_nomera_2t(q, n, f, fnomer);
	
}

//тест булевых функций
void wf4()
{
	int n = 4; // количество переменных
	int *x = new int[n]; // для хранения значений переменных
	int q = (int)pow(2, n);
	
	int non, deg, w;

	int *f = new int[q];

	for (int i = 0; i < q; i++)//зададим функцию многочленом
	{
		w = i;
		for (int j = 0; j < n; j++)
		{
			x[j] = w % 2;
			w = w / 2;
		}
		f[i] =  x[0]*x[1]+x[2] + x[2]*x[0]*x[3] ;
		f[i] = f[i] % 2;

	}

	TestBoolFunc(q, n, f, non, deg);

}


//полный тест измельчения функции над полем GF(2^t)
void wf5()
{
	//зададим поле
	InitGF(2, t);
	int q = (int)pow(2, t); 
	int n = 2, // n - количество переменных у функции
		qn = (int)pow(q, n);

	
	//зададим функцию
	element *f = new element[qn],
		*pol = new element[qn],
		*x = new element[n],
		xi(1, p, t),
		one(0, p, t),
		zero(q-1, p, t);

	for (int i = 0; i < n; i++)
		x[i] = one;


	for (int i = 0; i < qn; i++)//зададим функцию многочленом
	{
		int	w = i;
		for (int j = 0; j < n; j++)
		{
			x[j].exp = w % q;
			w = w / q;
		}
		/*
		f[i] = (x[0].ElPower(q-2) + x[1].ElPower((q - 2)/2) + one);
		f[i] = f[i].ElPower(q - 5);
		*/

		f[i] = x[0]*x[0]*x[0]*x[1]*x[1] +  x[0].ElPower(q-2)*x[1] + x[1];  //Падает устойчивость


		

		/*
		f[i] = x[0]*x[0]*x[0] + x[1].ElPower(6) + one;
		f[i] = f[i].ElPower(6);
		*/

	}

	//тест функции над большим полем
	Test_reduce_func_2t(q, t, n, f);





}

//построение одной устойчивой функции из четырех устойчивых нелинейных функций типа как у Таранникова
void wf44()
{
	
	InitGF(p, t);


	int q = (int)pow(p, t);
	int n = 4, w;
	int qn = (int)pow(q, n);

	cout << "q = " << q << ", n = " << n << endl;

	element
		*F = new element[qn * q],
		*f1 = new element[qn],
		*f2 = new element[qn],
		*f3 = new element[qn],
		*f4 = new element[qn],
		*x = new element[n + 1],
		*c = new element[q],
		xi(1, p, t),
		one(0, p, t);

	for (int i = 0; i < n+1; i++)
		x[i] = one;

	for (int i = 0; i < q; i++)
		c[i] = one;

	for (int i = 0; i < q; i++)
		c[i].exp = i;

	for (int i = 0; i < qn; i++)//зададим функцию 1 многочленом
	{
		w = i;
		for (int j = 0; j < n; j++)
		{
			x[j].exp = w % q;
			w = w / q;
		}
		f1[i] = x[0] * x[1] + x[2] + x[3];//вот тут надо написать многочлен от переменных x[0],...,x[n-1]
	}

	for (int i = 0; i < qn; i++)//зададим функцию 2 многочленом
	{
		w = i;
		for (int j = 0; j < n; j++)
		{
			x[j].exp = w % q;
			w = w / q;
		}
		f2[i] = x[0] * x[1] + xi * x[2] + x[3];//вот тут надо написать многочлен от переменных x[0],...,x[n-1]
	}

	for (int i = 0; i < qn; i++)//зададим функцию 3 многочленом
	{
		w = i;
		for (int j = 0; j < n; j++)
		{
			x[j].exp = w % q;
			w = w / q;
		}
		f3[i] = x[0] * x[1] + x[2] + xi * x[3];//вот тут надо написать многочлен от переменных x[0],...,x[n-1]
	}

	for (int i = 0; i < qn; i++)//зададим функцию 4 многочленом
	{
		w = i;
		for (int j = 0; j < n; j++)
		{
			x[j].exp = w % q;
			w = w / q;
		}
		f4[i] = x[0] * x[1] + xi* x[2] + xi * x[3];//вот тут надо написать многочлен от переменных x[0],...,x[n-1]
	}

	for (int i = 0; i < qn*q; i++)//зададим функцию 4 многочленом
	{
		w = i;
		for (int j = 0; j < n + 1; j++)
		{
			x[j].exp = w % q;
			w = w / q;
		}
		F[i] = (one + (x[4] + c[0])*(x[4] + c[0])*(x[4] + c[0])) * f1[i % qn] +
			   (one + (x[4] + c[1])*(x[4] + c[1])*(x[4] + c[1])) * f2[i % qn] +
			   (one + (x[4] + c[2])*(x[4] + c[2])*(x[4] + c[2])) * f3[i % qn] +
			   (one + (x[4] + c[3])*(x[4] + c[3])*(x[4] + c[3])) * f4[i % qn] ;

		;//вот тут надо написать многочлен от переменных x[0],...,x[n-1]
	}
	
	/*
	for (int i = 0; i < qn; i++)
		cout << f1[i].exp << " ";
	cout << endl;
	for (int i = 0; i < qn; i++)
		cout << f2[i].exp << " ";
	cout << endl;
	for (int i = 0; i < qn; i++)
		cout << f3[i].exp << " ";
	cout << endl;
	for (int i = 0; i < qn; i++)
		cout << f4[i].exp << " ";
	cout << endl;

	for (int i = 0; i < qn*q; i++)
		cout << F[i].exp << " ";
	cout << endl;
	*/

	
	Test_2t(q, t, n, f1);
	Test_2t(q, t, n, f2);
	Test_2t(q, t, n, f3);
	Test_2t(q, t, n, f4);
	


	Test_2t(q, t, n + 1, F);

}

//итерационный процесс построения устойчивых функций из x1*x2 от 3к+2 переменных, (2к-1)-устойчивых, (наверное) максимальной нелинейности, почти по Таранникову
void wf55()
{
	

	InitGF(p, t);
	int K = 3;//количество итераций
	int q = (int)pow(p, t);
	int n, w, qn;


	element **f = new element*[K];//здесь будут хранится создаваемые функции
	for (int i = 0; i < K; i++)
		f[i] = new element[(int)pow(q, 3 * i + 2)];

	element ***fbuf = new element**[K];//здесь будут хранится промежуточные создаваемые функции
	for (int i = 0; i < K; i++)
		fbuf[i] = new element*[q];

	for (int i = 0; i < K; i++)
		for (int j = 0; j < q; j++)
			fbuf[i][j] = new element[(int)pow(q, 3 * i + 4)];


	element
		*x = new element[3 * K + 2],
		*c = new element[q],
		xi(1, p, t),
		one(0, p, t),
		zero(q-1, p, t);

	for (int i = 0; i < 3 * K + 2; i++)
		x[i] = one;

	for (int i = 0; i < q; i++)
		c[i] = one;
	for (int i = 0; i < q; i++)
		c[i].exp = i;

	n = 2;
	for (int i = 0; i < q*q; i++)//зададим функцию f_0 многочленом
	{
		w = i;
		for (int j = 0; j < n; j++)
		{
			x[j].exp = w % q;
			w = w / q;
		}
		f[0][i] = x[0] * x[1] ;//вот тут надо написать многочлен от переменных x[0],...,x[n-1]
	}

	Test_2t(q, t, n, f[0]);//тестим ее

	//дальше включаем итеративный процесс
	for (int k = 0; k < K; k++)// k<K - math_style ) // строим функцию f_(k+1)
	{
		//сначала построим q буферных функций

		n = 3 * k + 4;// количество переменных в функции fbuf[k][iq]
		qn = (int)pow(q, n);

		for (int iq = 0; iq < q; iq++)
		{

			for (int i = 0; i < (int)pow(q, n); i++)//зададим функцию fbuf[k][iq] многочленом
			{
				w = i;
				for (int j = 0; j < n; j++)
				{
					x[j].exp = w % q;
					w = w / q;
				}

				fbuf[k][iq][i] =  f[k][i % (qn/(q*q)) ] + c[iq / (q-1)]*x[n-2] + c[iq % (q-1)] * x[n-1];//вот тут надо написать многочлен от переменных x[0],...,x[n-1]
			}
			/*
			for (int i = 0; i < q*q; i++)
				cout << f[k][i].exp << " ";
			cout << endl;

			for (int i = 0; i < (int)pow(q, n); i++)
				cout << fbuf[k][iq][i].exp << " ";
			cout << endl;
			*/
			Test_2t(q, t, n, fbuf[k][iq]);//тестим ее
		}

		//построили q штук вспомогательных функций, теперь склеиваем их
		//инициализируем функцию f[k+1]
		n = 3 * (k+1) + 2;// количество переменных в функции f[k+1]

		qn = (int)pow(q, n);

		for (int i = 0; i < qn; i++)
			f[k + 1][i] = zero;

		for (int i = 0; i < qn; i++)//зададим функцию f[k+1] многочленом
		{
			w = i;
			for (int j = 0; j < n; j++)
			{
				x[j].exp = w % q;
				w = w / q;
			}

			for (int iq = 0; iq < q; iq++)
				f[k + 1][i] = f[k + 1][i] + (one + (x[n-1] + c[iq])*(x[n-1] + c[iq])*(x[n-1] + c[iq])) * fbuf[k][iq][i % (qn/q)];

		}

		//типа задали f[k+1]

		Test_2t(q, t, n, f[k + 1]);//тестим ее

	}
}

//ищем 1-сильно устойчивую функцию типа GF^2(4)->GF(2)
void wf6()
{
	

	InitGF(p, t);
	int q = (int)pow(p, t);

	int n=2, w;

	int qn = (int)pow(q, n);

	bool wb;

	element *f = new element[qn];
	
	int *chi_f = new int[qn], *spectr_f = new int[qn];

	int ** B_spectr = Get_B_matrix_for_func_to_spectr_2t(t);

	element
		*x = new element[n],
		*c = new element[q],
		xi(1, p, t),
		one(0, p, t),
		zero(q - 1, p, t);

	for (int i = 0; i < n; i++)
		x[i] = one;

	for (int i = 0; i < q; i++)
		c[i] = one;
	for (int i = 0; i < q; i++)
		c[i].exp = i;


	int FuncSoul = 32640; // Число в диапозоне от 0 до 2^15-1 = 32767  определяющее функцию f своим двоичным представлением


	for (FuncSoul = 0; FuncSoul < 32767; FuncSoul++)
	{
		if (Weigth1(2, FuncSoul) != qn / 2)
			continue;

	

		w = FuncSoul;
		f[0] = zero;
		for (int i = 1; i < qn; i++)//зададим функцию f 
		{
			if (w % 2 == 0)
				f[i] = zero;
			else
				f[i] = one;
			w = w / 2;
		}

	
		//Test_2t(q, t, n, f);//тестим ее

		//проверим на сильную устойчивость


			wb = Test_strong_resilient_2t(FuncSoul, q, t, n, f);

			if (wb == true)
			{
				for (int i = 0; i < qn; i++)
					cout << f[i].exp << " ";
				cout << endl;
				cout << "Succes!\n";
			}
				


	}


}

void wf7()
{
	//построили 1-сильно устойчивую функцию вида GF^3(4)->GF(2)

	//строим булеву функцию, которая при особом преобразовании сохраняет устойчивость

	int n = 6, w;
	int pn = (int)pow(p, n);
	const int count_basis = 3;//количество базисов
	int res[count_basis];
	bool wb;

	int *h = new int[(int)pow(p, 4)];
	int **f = new int* [count_basis];//память для функций
	for (int i = 0; i < count_basis; i++)
		f[i] = new int[pn];

	int * * mem = new int*[n + 1];//память для БПФ
	for (int i = 0; i < n + 1; i++)
		mem[i] = new int[pn];

	//зададим изначальную функцию числом FuncSoul


	int FuncSoul = 32640; // Число в диапозоне от 0 до 2^15-1 = 32767  определяющее функцию f своим двоичным представлением

	//unsigned __int64 FuncSoul=0, x=3, M64 = (1 << 63)-1; 

	for (FuncSoul = 0; FuncSoul < 32640; FuncSoul++)
	{
		if (Weigth1(2, FuncSoul) != (int)pow(p, 4) / 2)//*функция должна быть сбалансированной
			continue;



		//зададим промежуточную функцию h
		w = FuncSoul;
		h[0] = 0;
		for (int i = 1; i < (int)pow(p, 4); i++)//зададим функцию h
		{
			h[i] = w % 2;
			w = w / 2;
		}
		

		for (int i = 0; i < pn; i++)//зададим функцию f 
			f[0][i] = ( h[i % 16] + (i/16) % 2 + (i / 32) % 2) % 2; // h плюс две линейных переменных




		StatStructBPFBool(mem, n, pn, f[0]);
		res[0] = Resilient_Index(pn, n, mem[n]);

		if (res[0] == -1)
			continue;

		
			
		//построим остальные функции

		int new_basis[2][4] = { {0,1,3,2},{0,3,1,2} };//функция перехода из стандартного базиса в остальные

		for (int ib = 1; ib <= 2; ib++)
		{
			for (int i = 0; i < pn; i++)
				f[ib][new_basis[ib-1][i%4]     + new_basis[ib-1][(i / 4) % 4]*4 + new_basis[ib - 1][(i / 16) % 4] * 16] = f[0][i];

			StatStructBPFBool(mem, n, pn, f[ib]);
			res[ib] = Resilient_Index(pn, n, mem[n]);

		}



		if ((res[0] + res[1] + res[2]) > 8)
		{
			for (int ib = 0; ib <= 2; ib++)
				cout << res[ib] << " ";
			cout << endl;
		}
	


		/*
		for (int i = 0; i < pn; i++)
			cout << f[0][i] << " ";
		cout << endl;

		for (int i = 0; i < pn; i++)
			cout << f[1][i] << " ";
		cout << endl;

		for (int i = 0; i < pn; i++)
			cout << f[2][i] << " ";
		cout<< endl;
		*/

		
	}

}

void wf8()
{
	//2-сильно устойчивая функция

	//строим булеву функцию от 8 переменных, которая при особом преобразовании сохраняет устойчивость

	int n = 8, w;
	int pn = (int)pow(p, n);
	const int count_basis = 3;//количество базисов
	int res[count_basis];


	int *h = new int[(int)pow(p, 4)];
	int **f = new int*[count_basis];//память для функций
	for (int i = 0; i < count_basis; i++)
		f[i] = new int[pn];

	int * * mem = new int*[n + 1];//память для БПФ
	for (int i = 0; i < n + 1; i++)
		mem[i] = new int[pn];

	//зададим изначальную функцию числом FuncSoul


	int FuncSoul = 32640; // Число в диапозоне от 0 до 2^15-1 = 32767  определяющее функцию f своим двоичным представлением

						  //unsigned __int64 FuncSoul=0, x=3, M64 = (1 << 63)-1; 

	for (FuncSoul = 0; FuncSoul < 32640; FuncSoul++)
	{
		if (Weigth1(2, FuncSoul) != (int)pow(p, 4) / 2)//*функция должна быть сбалансированной
			continue;



		//зададим промежуточную функцию h
		w = FuncSoul;
		h[0] = 0;
		for (int i = 1; i < (int)pow(p, 4); i++)//зададим функцию h
		{
			h[i] = w % 2;
			w = w / 2;
		}


		for (int i = 0; i < pn; i++)//зададим функцию f 
			f[0][i] = (h[i % 16] + (i / 16) % 2 + (i / 32) % 2 + (i / 64) % 2 + (i / 128) % 2) % 2; // h плюс четыре линейных переменных




		StatStructBPFBool(mem, n, pn, f[0]);
		res[0] = Resilient_Index(pn, n, mem[n]);

		if (res[0] == -1)
			continue;



		//построим остальные функции

		int new_basis[2][4] = { { 0,1,3,2 },{ 0,3,1,2 } };//функция перехода из стандартного базиса в остальные

		for (int ib = 1; ib <= 2; ib++)
		{
			for (int i = 0; i < pn; i++)
				f[ib][new_basis[ib - 1][i % 4] + new_basis[ib - 1][(i / 4) % 4] * 4 + new_basis[ib - 1][(i / 16) % 4] * 16 + new_basis[ib - 1][(i / 64) % 4] * 64] = f[0][i];

			StatStructBPFBool(mem, n, pn, f[ib]);
			res[ib] = Resilient_Index(pn, n, mem[n]);

		}



		if ((res[0] + res[1] + res[2]) > 11)
		{
			for (int ib = 0; ib <= 2; ib++)
				cout << res[ib] << " ";
			cout << endl;
		}



		/*
		for (int i = 0; i < pn; i++)
		cout << f[0][i] << " ";
		cout << endl;

		for (int i = 0; i < pn; i++)
		cout << f[1][i] << " ";
		cout << endl;

		for (int i = 0; i < pn; i++)
		cout << f[2][i] << " ";
		cout<< endl;
		*/


	}

}

void wf9()
{
	//строим функцию GF(4)^n -> GF(2) (точнее ее булев аналог GF(2)^(2n) -> GF(2)), которая при измельчении переменных (каждой переменной по своему базису) сохраняет относительную устойчивость
	
	ofstream dt("data.txt");

	dt << "asdasd\n";

	srand(time(0));

	const unsigned int timer = 214750;

	int hn=5, n = 8, w;
	int pn = ( int)pow(p, n);
	const int count_basis = 3;//количество базисов
	int *res = new int [81]; //сюда будем сохранять устойчивость измельченных функций
	bool wb;

	int *h = new  int[(int)pow(p, hn)];//промежуточная функция для построения
	 int **f = new int*[count_basis];//память для функций
	for ( int i = 0; i < count_basis; i++)
		f[i] = new  int[pn];

	int * * mem = new int*[n + 1];//память для БПФ
	for (int i = 0; i < n + 1; i++)
		mem[i] = new int[pn];

	//зададим изначальную функцию числом FuncSoul


	unsigned int w64, t, FuncSoul = 32640; // Число в диапозоне от 0 до 2^15-1 = 32767  определяющее функцию f своим двоичным представлением

						  //unsigned __int64 FuncSoul=0, x=3, M64 = (1 << 63)-1; 

	//cin >> t;

	t = 8000000;

	//for (FuncSoul = 65535; FuncSoul != (int)pow(2, 31) ; FuncSoul++)
	for (FuncSoul = t; FuncSoul != (unsigned int) pow(2, 31)-1; FuncSoul++)
	/*
	h[0] = 0;
	for (int i = 1; i < 32; i++)
		if (i < 17)
			h[i] = 1;
		else
			h[i] = 0;
	t = 0;
	while (Get_next_vector_this_weight(h, 32))
	*/
	{
		/*
		t++;
		if (t % 3000 == 0)
			cout << t << endl << endl;
			*/


		if (FuncSoul % timer == 0)
			cout << FuncSoul / timer << endl;
		
		


		if (Weigth1(p, FuncSoul) != (int)pow(2, hn) / 2)//*функция должна быть сбалансированной
			continue;
		
		//cout << "asdasd";
		
		//зададим промежуточную функцию h
		w64 = FuncSoul;

		
		h[0] = 0;
		for (int i = 1; i < (int)pow(p, hn); i++)//зададим функцию h
		{
			h[i] = (int)(w64 % 2);
			w64 = w64 / 2;
		}
		

		/*
		w = 0;
		//зададим функцию h случайно
		for (int i = 1; i < (int)pow(p, 5); i++)//зададим функцию h
		{
			h[i] = rand() % 2;
			w += h[i];
		}
		
		if (w != (int)pow(p, 5) / 2)
			continue;
		*/


		for (int i = 0; i < pn; i++)//зададим функцию f 
			f[0][i] = (h[i % 32]  + (i / 32) % 2 + (i / 64) % 2 + (i / 128) % 2) % 2; // h плюс 3 линейных переменных




		StatStructBPFBool(mem, n, pn, f[0]);
		res[0] = Resilient_Index(pn, n, mem[n]);

		/*
		if (res[0] == -1)
			continue;
		*/
		//cout << "asdasd";

		//построим остальные функции

		int new_basis[3][4] = { {0,1,2,3}, { 0, 1, 3, 2 }, { 0,3,1,2 }};//функция перехода из стандартного базиса в остальные

	

		for (int ib = 1; ib < 81; ib++)//каждая переменная будет измельчаться по своему базису n переменных для каждой count_basis вариантов базисов
		{//f[1] будет буфером для получаемых функций
			//count_basis === 3
			for (int i = 0; i < pn; i++)
				f[1][new_basis[ib % 3][i % 4] + new_basis[(ib/3) % 3][(i / 4) % 4] * 4 + new_basis[(ib / 9) % 3][(i / 16) % 4] * 16 + new_basis[(ib / 27) % 3][(i / 64) % 4] * 64] = f[0][i];

			StatStructBPFBool(mem, n, pn, f[1]);
			res[ib] = Resilient_Index(pn, n, mem[n]);

		}

		w = 0;
		for (int i = 0; i < 81; i++)
		{
			w += res[i];
		}

		if (w >= 81*3)
		{
			for (int ib = 0; ib < 81; ib++)
				cout << res[ib] << " ";
			cout << endl << FuncSoul << endl << endl;
			system("pause");
		}



		/*
		for (int i = 0; i < pn; i++)
		cout << f[0][i] << " ";
		cout << endl;

		for (int i = 0; i < pn; i++)
		cout << f[1][i] << " ";
		cout << endl;

		for (int i = 0; i < pn; i++)
		cout << f[2][i] << " ";
		cout<< endl;
		*/


	}

}

// тест булевой функции
void wf10()
{
	int f[16]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1};
	int non, deg;

	TestBoolFunc(16, 4, f, non, deg);
}

//тест функций в произвольных конечных полях
void wf13()
{
	
	/*
	int f[25] = { 0,1,2,3,4,1,4,3,2,0,2,3,4,0,1,3,0,1,4,2,4,2,0,1,3 };
	Test(5, 2, f);
	*/

	const int q = 3; 
	const int n = 3;
	int f[27], w;

	int x[n];


	for (int i=0; i < 27; ++i)
	{
		w=i;
		for (int j = 0; j < n; j++)
		{
			x[j] = w % 3;
			w = w / 3;
		}
		
		f[i] = 0;

		if (x[0]==0)
			f[i] = 0;

		if (x[1]==0)
			f[i] = 0;

		if (x[2]==0)
			f[i] = 0;

		w = x[0] + x[1] + x[2];

		if (x[0] > 0)
			if (x[1] > 0)
				if (x[2] > 0)
					if (w >= 3)
						f[i] = 1;

		if (x[0] >= 1)
			if (x[1] >= 1)
				if (x[2] >= 1)
					if (w >= 5)
						f[i] = 2;
	}

	Test(3, 3, f);

	//Convert_exp_func_on_nomera_2t(q, n, f, fnomer);

}

//посчитаем распределение нелинейности у квазигрупп
//для простых полей
void wf12()
{
	

	mt19937 gen;
	gen.seed(time(0));
	//gen.seed(0);


	int p = 5;
	int n = 10;
	int pn = (int)pow(p, n);

	cout << "p = " << p << ", n = " << n << endl;

	double *distr = new double[n + 1]; // значению устойчивости будем сопоставлять минимальную нелинейность.
	for (int it = 0; it < n + 1; it++)
		distr[it] = pn;

	int deg, ri;
	int MinPolWeight = pn, *MinPol = new int[pn], *MinF = new int[pn];
	int * f = new int[pn];//выделим память для работы
	int * pol = new int[pn];
	complex <double> * spectr = new complex <double>[pn];
	complex <double> * chi_f = new complex <double>[pn];

	double nonlin, min = pn;

	complex< double > ** B_spectr = Get_B_matrix_for_func_to_spectr(p);
	int ** B_polinom = Get_B_matrix_for_func_to_pol(p);


	while (1)
	{

		Get_random_quazyfunc(gen, p, n, f);//случайная квазигруппа

		chi_func(p, n, f, chi_f);

		BPF(p, n, B_spectr, chi_f, spectr);

		nonlin = AbsMax(pn, spectr)*pn;

		ri = Resilient_index(p, n, spectr);

		if (min > nonlin)
		{
			min = nonlin;

			cout << min << endl;
		
		}

	}
}


// тест булевой функции
void wf14()
{
	

	int n = 6, w;
	int pn = (int)pow(p, n);
	const int count_basis = 3;//количество базисов
	int res[count_basis];
	bool wb;

	int *h = new int[(int)pow(p, 4)];
	int **f = new int*[count_basis];//память для функций
	for (int i = 0; i < count_basis; i++)
		f[i] = new int[pn];

	int * * mem = new int*[n + 1];//память для БПФ
	for (int i = 0; i < n + 1; i++)
		mem[i] = new int[pn];

	//зададим изначальную функцию числом FuncSoul


	int FuncSoul = 32640; // Число в диапозоне от 0 до 2^15-1 = 32767  определяющее функцию f своим двоичным представлением

						  //unsigned __int64 FuncSoul=0, x=3, M64 = (1 << 63)-1; 

	for (FuncSoul = 0; FuncSoul < 32640; FuncSoul++)
	{
		if (Weigth1(2, FuncSoul) != (int)pow(p, 4) / 2)//*функция должна быть сбалансированной
			continue;



		//зададим промежуточную функцию h
		w = FuncSoul;
		h[0] = 0;
		for (int i = 1; i < (int)pow(p, 4); i++)//зададим функцию h
		{
			h[i] = w % 2;
			w = w / 2;
		}


		for (int i = 0; i < pn; i++)//зададим функцию f 
			f[0][i] = (h[i % 16] + (i / 16) % 2 + (i / 32) % 2) % 2; // h плюс две линейных переменных




		StatStructBPFBool(mem, n, pn, f[0]);
		res[0] = Resilient_Index(pn, n, mem[n]);

		if (res[0] == -1)
			continue;



		//построим остальные функции

		int new_basis[2][4] = { { 0,1,3,2 },{ 0,3,1,2 } };//функция перехода из стандартного базиса в остальные

		for (int ib = 1; ib <= 2; ib++)
		{
			for (int i = 0; i < pn; i++)
				f[ib][new_basis[ib - 1][i % 4] + new_basis[ib - 1][(i / 4) % 4] * 4 + new_basis[ib - 1][(i / 16) % 4] * 16] = f[0][i];

			StatStructBPFBool(mem, n, pn, f[ib]);
			res[ib] = Resilient_Index(pn, n, mem[n]);

		}



		if ((res[0] + res[1] + res[2]) > 8)
		{
			for (int ib = 0; ib <= 2; ib++)
				cout << res[ib] << " ";
			cout << endl;
		}



		/*
		for (int i = 0; i < pn; i++)
		cout << f[0][i] << " ";
		cout << endl;

		for (int i = 0; i < pn; i++)
		cout << f[1][i] << " ";
		cout << endl;

		for (int i = 0; i < pn; i++)
		cout << f[2][i] << " ";
		cout<< endl;
		*/


	}

}


void wf_basis(int k)
{

	Build_All_Basis_2t(k);

}

//перебираем функции (x^?+y^?+1)^? в поисках низкой линейности
void wf33()
{
	srand(time(NULL));
	int a,b,c,d,w;

	InitGF(p, t);
	int q = (int)pow(p, t); 
	int n = 2, // n - количество переменных у функции
		qn = (int)pow(q, n);

	

	element *f = new element[qn],
		*x = new element[n],
		xi(1, p, t),
		one(0, p, t),
		zero(q-1, p, t);

	for (int i = 0; i < n; i++)
		x[i] = one;
	/*
	int *fnomer = new int[qn],
		*chi_f2 = new int[qn],
		*spectr_f2 = new int[qn], w;
	*/

	int **chi_f = new int *[q-1], **spectr_f = new int *[q-1];

	for (int i=0; i< q-1; i++)
	{
		chi_f[i] = new int[qn]; // здесь будет храниться столбец значения функции \chi(a*f(x))
		spectr_f[i] = new int[qn]; // здесь будут храниться кросс-корреляционные коэффициенты С(af, <u,x>) функции \chi(a*f(x))
	}
	
	element * pol = new element [qn];

	element ** B_polinom = Get_B_matrix_for_func_to_pol_2t(2, t);

	int ** B_spectr = Get_B_matrix_for_func_to_spectr_2t(t);

	element ela(2,t);

	//int v[]={5,11,13,17,19,23,29,31,37,41,43,47,53};

	//int v[51]={7, 11, 13, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199, 211, 223, 227,229, 233, 239, 241, 251, 254};


	//while (1)
	//for (int i : v)
	{
		/*
		a=(rand()%(q/2-2))*2+1;
		b=(rand()%(q/2-2))*2+1;
		c=(rand()%(q/2-2))*2+1;
		*/


		/*
		a = v[rand()%51];
		b = v[rand()%51];
		c = 254;
		*/

		a = q - 2;
		b = (q - 2) / 2;
		c = q - 5;
		d = (q - 2) / 2;


		for (int i = 0; i < qn; i++)//зададим функцию многочленом
		{
			w = i;
			for (int j = 0; j < n; j++)
			{
				x[j].exp = w % q;
				w = w / q;
			}
			
		
			f[i] = (x[0].ElPower(a) + x[1].ElPower(b) + one);
			f[i] = f[i].ElPower(c);
			//f[i] = f[i] + x[2].ElPower((q - 2) / 2);
			//f[i] = f[i].ElPower(q-2);
			//f[i] = f[i] + x[3].ElPower((q - 2) / 2);
			//f[i] = f[i].ElPower(q - 2);
	
		}

		BPF_2t(q, 2, t, n, B_polinom, f, pol);

		int deg = Algebraic_degree_2t(q, t, n, pol);

		int dl = Nonliniarity_degree_2t(q, t, n, pol);

		for (int i=0; i<q-1; i++)
		{
			cout << i << endl;
			ela.exp = i;
			chi_func_2t(n, ela, f, chi_f[i]);
		
			BPF_2t(q, t, n, B_spectr, chi_f[i], spectr_f[i]);
		}
	

		int res = Resilient_index_2t(q, n, spectr_f);

		int nl = AbsMax(q-1, qn, spectr_f);


		//if (dl + res == n*t - 1)
		{
			cout << a << " "<< b << " " << c << endl;
			cout << "deg = " << deg << " <= " << n*(q-1) << endl;
			cout << "dl = "<< dl <<", res = "<< res << endl;
			cout << "nl : " << pow(qn, 1.0 / 2) << " <= " << nl << " <= " <<  qn << endl;
			cout << endl;
		}

	
	}

}

void wf34()
{
	
	InitGF(p, t);
	int q = (int)pow(p, t); 
	int n = 1, // n - количество переменных у функции
		qn = (int)pow(q, n);

	

	element *f = new element[qn],
		*pol = new element[qn],
		*x = new element[n],
		xi(1, p, t),
		one(0, p, t),
		zero(q-1, p, t),
		we(1, p, t);

	for (int i = 0; i < n; i++)
		x[i] = one;

	int w,w1;

	for (int ex=0; ex < q*q*q*q; ++ex)
		//int ex = 2;
	{

		vector <int> vec;
		w=ex;
		for (int z=1; z < t; ++z)
		{
			vec.push_back(w % p);
			w = w/p;
		}

		for (int i = 0; i < qn; i++)//зададим функцию многочленом
		{
			w = i;
			for (int j = 0; j < n; j++)
			{
				x[j].exp = w % q;
				w = w / q;
			}
		
			f[i] =  x[0].ElPower(q-2);

		
			for (int z=0; z < t-1; ++z)
			{
				we.exp = vec[z];
				f[i] = f[i] + we*x[0].ElPower((int)pow(2, z+1));
			}
			

		}

		if ( ex % q*q == 0)
			cout << ex << endl;


		if (Work2_Test_2t(q, t, n, f))
			cout <<  endl;
	}
	
	

	
}


//Возведем функцию во все степени и посмотрим как изменятся параметры
void wf35()
{
	
	
	/*
	int d[] = {1, 2, 4, 5, 8, 10, 11, 13, 16, 17, 19, 20, 22, 23, 25, 26, 29, 31, 32, 34, 37, 38, 40, 41, 43, 44, 46, 47, 50, 52, 53, 55, 58, 59, 61, 62};
	const int N=36;
	*/

	int d[] = {1, 2, 4, 7, 8, 11, 13, 14, 16, 19, 22, 23, 26, 28, 29, 31, 32, 37, 38, 41, 43, 44, 46, 47, 49, 52, 53, 56, 58, 59, 61, 62, 64, 67, 71, 73,
		74, 76, 77, 79, 82, 83, 86, 88, 89, 91, 92, 94, 97, 98, 101, 103, 104, 106, 107, 109, 112, 113, 116, 118, 121, 122, 124, 127, 128, 131, 133, 134,
		137, 139, 142, 143, 146, 148, 149, 151, 152, 154, 157, 158, 161, 163, 164, 166, 167, 169, 172, 173, 176, 178, 179, 181, 182, 184, 188, 191, 193, 
		194, 196, 197, 199, 202, 203, 206, 208, 209, 211, 212, 214, 217, 218, 223, 224, 226, 227, 229, 232, 233, 236, 239, 241, 242, 244, 247, 248, 251, 253, 254};
	const int N=128;

	vector <int> degs(d,d+N);



	InitGF(p, t);
	int q = (int)pow(p, t); 
	int n = 2, // n - количество переменных у функции
		qn = (int)pow(q, n);

	

	element *f = new element[qn],
		*g = new element[qn],
		*pol = new element[qn],
		*x = new element[n],
		xi(1, p, t),
		one(0, p, t),
		zero(q-1, p, t),
		ebuf(1, p, t);

	for (int i = 0; i < n; i++)
		x[i] = one;

	int w;



	for (int i = 0; i < qn; i++)//зададим функцию многочленом
	{
		w = i;
		for (int j = 0; j < n; j++)
		{
			x[j].exp = w % q;
			w = w / q;
		}
		
		//f[i] = x[0];

		/*
		f[i] = (x[0].ElPower(q-2) + x[1].ElPower((q-2)/2) + one);
		f[i] = f[i].ElPower(q-5);
		*/

		f[i] = (x[0].ElPower(2) + x[1].ElPower(13) + one);
		f[i] = f[i].ElPower(q-2);


	}

		

	for (auto it : degs)	// возведем в степень
	{
		for (int i = 0; i < qn; i++)
		{
			g[i] = f[i].ElPower(it);
		}
		cout << it << endl;
		Work2_Test_2t(q, t, n, g);
		cout << endl;
	}
	
	
}


void main()
{
	//wf13();
	wf33();

}