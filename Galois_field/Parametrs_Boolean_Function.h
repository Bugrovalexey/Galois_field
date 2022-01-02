// Бугров А.Д.
// Хедер для работы с булевыми функциями
// булевы функции задаем столбцом значений через массивы unsigned char
// при лексикографическом упорядочивании переменных


#pragma once

#include <iostream>
//#include <map>
#include<vector>
using namespace std;

// Нахождение модуля элемента массива максимального по модулю
int MaxAbs(int n, int * ar)
{
	int w = 0;
	for (int i = 0; i < n; i++)
	{
		if ( abs(ar[i]) > abs(w) )
			w = abs(ar[i]);
	}
	return w;
}

// БПФ для нахождения коэффициентов Уолша-Адамара (спектральных коэффициентов) булевой функции от n переменных
int * StatStructBPFBool(int **w, int n, int pn, unsigned char * func) 
{
	// n - количество перемнных
	// pn = 2^n
	// **w - указатель на память, выделенную для проведения операций

	int pnn = pn, pnnn, i, j;

	for (i = 0; i < pn; i++)//преобразование булевой в псевдобулеву
	{
		if (func[i] == 0)
			w[0][i] = 1;
		else
			w[0][i] = -1;
	}


	for (i = 1; i < n + 1; i++)
	{
		pnnn = pnn / 2;
		for (j = 0; j < pn; j++)
		{
			if (j%pnn < pnnn)
				w[i][j] = w[i - 1][j] + w[i - 1][j + pnnn];
			else
				w[i][j] = w[i - 1][j] - w[i - 1][j - pnnn];
			//cout<< w[i][j]<<" ";
		}
		//	cout << "\n";
		pnn = pnnn;
	}
	return w[n];
}

// Двоичный вес числа а
int BinWeight(int a)
{
	int b = 0;
	while (a > 0)
	{
		b += a % 2;
		a = a / 2;
	}
	return b;
}

// получение степени многочлена Жигалкина, bw - массив с двоичными весами чисел от 0 до pn-1
int DegPolZheg( int * pol, int pn, int * bw)
{
	int deg = 0;
	for (int i = 0; i < pn; i++)
		if ((pol[i] == 1) && (deg < bw[i]))
			deg = bw[i];

	return deg;
}

// получение многочлена Жегалкина булевой функции через БПФ
void PolZhegBoolFunc(unsigned char **w, int n, int pn, unsigned char* func)
{

	int pnn = pn, pnnn, i, j;

	for (i = 0; i < pn; i++)//заполнение в память столбца значений функции
		w[0][i] = func[i];

	for (i = 1; i < n + 1; i++)
	{
		pnnn = pnn / 2;
		for (j = 0; j < pn; j++)
		{
			if (j%pnn < pnnn)
				w[i][j] = w[i - 1][j];
			else
				w[i][j] = w[i - 1][j] ^ w[i - 1][j - pnnn];
			//	cout<< w[i][j]<<" ";
		}
		//	cout << "\n";
		pnn = pnnn;
	}
}

// вычисление порядка корреляционной иммунности булевой функции по спектру
int Correlation_Immunity_Index(int q, int n, int * spectr)
{
	int w=0;

	vector<bool> cor(n+1);

	for (int i = 0; i < n + 1; i++)
		cor[i] = true;


	for (int i = 1; i < q; i++)
	{
		if (spectr[i] != 0)
			cor[BinWeight(i)] = false;
	}

	for (int i = 1; i < n + 1; i++)
		if (cor[i] == false)
		{
			w = i - 1;
			break;
		}

	return w;
}

int Resilient_Index(int q, int n, int * spectr)
{
	//определить коэффициент устойчивости булевой функции
	int w = -1;

	//bool * cor = new bool[n + 1];
	vector<bool> res(n + 1);

	for (int i = 0; i < n + 1; i++)
		res[i] = true;


	for (int i = 0; i < q; i++)
	{
		if (spectr[i] != 0)
			res[BinWeight(i)] = false;
	}

	for (int i = 0; i < n + 1; i++)
		if (res[i] == false)
		{
			w = i - 1;
			break;
		}

	return w;
}


// Вычисление характеристик булевой функции func
void TestBoolFunc(int n, unsigned char *func, int & nl, int & deg, int& res, bool print_flag=false)
{
	// t - количество переменных
	//q = 2^t

	ui q = ui(pow(2, n));



	int **mem = new int*[n + 1]; 
	for (int i = 0; i < n + 1; i++)
		mem[i] = new int[q];

	int *spectr = new int[q];
	int *pol = new int[q];

	int *binweight;
	binweight = new int[q];
	for (int i = 0; i < q; i++)
		binweight[i] = BinWeight(i);
	 

	StatStructBPFBool(mem, n, q, func);

	for (int i = 0; i < q; i++)
		spectr[i] = mem[n][i];

	nl = MaxAbs(q, spectr);

	
	PolZhegBoolFunc( (unsigned char**)mem, n, q, func);

	for (int i = 0; i < q; i++)
		pol[i] = mem[n][i];

	deg = DegPolZheg(pol, q, binweight);


	res = Resilient_Index(q, n, spectr);

	if (print_flag)
	{
		cout << endl << "n = " << n << endl;
		cout << "nl : " << pow(q, 1.0 / 2) << " <= " << nl << " <= " << q << endl;
		cout << "deg = " << deg << endl;
		cout << "res = " << res << endl;

		cout << endl;
		cout << "Ziegenthaler's inequality for\n";
		cout << "deg: " << deg << " + " << res;
		if (deg + res == n - 1)
			cout << " = ";
		else
			cout << " < ";
		cout << n - 1 << endl;

		
	}


	for (int i = 0; i < n + 1; i++)
		delete[] mem[i];
	delete[] mem;
	delete[] binweight;
}


// Вычисление характеристик подстановки S на двоичных векторах длины length, где 2^length=q
void TestBoolPermut(int length, int *S, int &Nl, int &Deg, int &Def, bool print_flag = false)
{
	int q = (int)pow(2, length);

	//память для хранения нетривиальных линейных комбинаций кординатных функций
	unsigned char* cordfunc = new unsigned char[q];

	int maxdefchar = 0, maxnl = 0, nl, deg, res, mindeg = length;

	// память для вычисления разностной характеристики
	int ** defchar = new int*[q];
	for (int i = 0; i < q; i++)
	{
		defchar[i] = new int[q];
		for (int j = 0; j < q; j++)
			defchar[i][j] = 0;
	}

	// Переберем все нетривиальные линейные комбинации кординатных функций
	for (int i = 1; i < q; i++)
	{
		for (int ii = 0; ii < q; ii++)
			cordfunc[ii] = 0;

		for (int j = 0; j < q; j++) // для задания столбца функции
		{
			for (int k = 0; k < length; k++) // для сложения битиков
				cordfunc[j] = cordfunc[j] ^ (((S[j] & i) >> k) & 1);

			cordfunc[j] = (2 + cordfunc[j] % 2) % 2; //задаем i-ую линейную комбинацию кординатных функций									 
		}

		TestBoolFunc(length, cordfunc, nl, deg, res);

		if (maxnl<nl)
			maxnl = nl;
		if (mindeg>deg)
			mindeg = deg;
	}

	for (int ii = 1; ii < q; ii++)
	{
		for (int j = 0; j < q; j++)
			defchar[ii] [ ((S[((j^ii) % q + q) % q] ^ S[j]) % q + q) % q ]++;
	}

	for (int i = 1; i < q; i++)
	{
		for (int j = 0; j < q; j++)
		{
			//cout<<defchar[i][j]<<" ";
			if (maxdefchar < defchar[i][j])
				maxdefchar = defchar[i][j];
		}
	}

	Nl = maxnl;
	Deg = mindeg;
	Def = maxdefchar;


	if (print_flag)
	{
		cout << "nl=" << maxnl << " deg=" << mindeg << " def=" << maxdefchar << "\n";
	}
	
	delete[] cordfunc;

	for (int i = 0; i < q - 1; i++)
		delete[] defchar[i];
	delete[] defchar;
}


//Проверка на подстановочность
int TestPermutation(int q, int* S)
{
	int* test = new int[q];

	for (int i = 0; i < q; i++)
		test[i] = 0;

	for (int i = 0; i < q; i++)
		if (test[S[i]] == 0)
			test[S[i]] = 1;
		else
		{
			delete[] test;
			return 0;//Не подстановка
		}
	delete[] test;
	return 1;//Подстановка
}


//Расширенный тест на подстановочность
//Если возвращает 0, то подстановка
int Test2Permutation(int q, ui *S)
{
	ui w=0, *test = new ui[q];

	for (int i = 0; i < q; i++)
		test[i] = 0;

	for (int i = 0; i < q; i++)
		test[S[i]] = 1;

	for (int i = 0; i < q; i++)
		if (test[i] == 0)
			w++;

	delete[] test;
	return w;
}

