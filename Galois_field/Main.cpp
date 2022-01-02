// Демонстрация возможностей проведения вычислений в конечных полях



//#include "Galua_Field.h"
//#include "Parametrs_Boolean_Function.h"
#include "BPF.h"
#include <fstream>
#include <vector>


//Вычислим параметры булевой функции
void Info_about_boolean_function()
{
	int n = 4; // количество переменных
	int* x = new int[n]; // для хранения значений переменных
	int q = (int)pow(2, n);

	int nl, deg, res, w;

	unsigned char* f = new unsigned char[q];

	for (int i = 0; i < q; i++) //зададим функцию многочленом
	{
		w = i;
		for (int j = 0; j < n; j++)
		{
			x[j] = w % 2;
			w = w / 2;
		}
		f[i] = x[0] * x[1] + x[2] * x[3];
		f[i] = f[i] % 2;

	}

	TestBoolFunc(n, f, nl, deg, res, true);

}


//Вычислим параметры подстановки x^{-1} в поле GF(256), представив ее как подстановку двоичных векторов длины 8
void Info_about_inverse_permutation()
{	
	const int p = 2, t = 8; // задаем поле GF(p^t)
	GaloisField F(p, t);
	int q = int(pow(p, t));

	element xi(1, F), one(0, F), zero(q - 1, F), x(F);

	ui* f = new ui[q];

	f[q-1] = q-1;
	for (int i = 0; i < q-1; i++) //зададим подстановку на степенях
	{
		f[i] = (i * (q - 2)) % (q-1);
	}

	ui* newf = new ui[q];
	Convert_exp_func_on_nomera_2t(F, 1, f, newf);

	if ( Test2Permutation(q, newf) != 0 )
		cout << "Warning!";

	int Nl, Deg, Def;
	TestBoolPermut(8, (int*) newf, Nl, Deg, Def, true);
}


//Вычислим параметры функции над полем GF(2^t)
void Info_about_bin_function()
{
	const int p = 2, t = 4; // задаем поле GF(p^t)
	GaloisField F(p, t);
	int q = F.q;

	int n = 2, // n - количество переменных у функции
		qn = (int)pow(q, n);



	element* f = new element[qn],
		* pol = new element[qn],
		* x = new element[n],
		xi(1, F),
		one(0, F),
		zero(q - 1, F),
		we(1, F);

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

		f[i] = x[0] * x[1] ;

	}

	Test_2t(F,  n, f);

}

//Вычислим параметры функции над полем GF(p)
void Info_about_p_function()
{
	const int p = 3, t = 1; // задаем поле GF(p^t)
	GaloisField F(p, t);
	int q = F.q;

	int n = 3, // n - количество переменных у функции
		qn = (int)pow(q, n);

	
	int *f = new int[qn], w;

	int *x = new int[n];


	for (int i = 0; i < qn; i++)//зададим функцию многочленом
	{
		w = i;
		for (int j = 0; j < n; j++)
		{
			x[j] = w % q;
			w = w / q;
		}

		f[i] = x[0] * x[1];

	}

	Test(p, n, f);
}



int main()
{

	//Info_about_inverse_permutation();

	//Info_about_boolean_function();

	//Info_about_bin_function();

	Info_about_p_function();

	return 0;
}