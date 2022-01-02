//������ �.�.
//��������� ������ � ��������� � ������ 1) ������� ����� Z_p 2) ������ ����� GF(2^n)
//���������� ������� �������������� �����
//���������� ������ ������ ��� ���
//����������� � ������ ����� �� ����� "������������ � ������������� ���������� �������" �.�. ����������� ���. 39

#pragma once
#define _USE_MATH_DEFINES

#include <math.h>
#include <complex>
#include <ctime>
#include <random>
#include <omp.h>
#include "Galua_Field.h"
#include "Parametrs_Boolean_Function.h"

using namespace std;


const complex< double > I(0, 1);


// ���������� ���������� ������� ��� ��������� �������������� ��������� �������������� ����� 
// ��� ��������� ������������� ���������� �� ������� �������� ������� � ������ �������� ��������� ����
int ** Get_B_matrix_for_func_to_pol(int p)
{
	int xi=get_primitive_root(p);

	int ** B = new int *[p];
	for (int i = 0; i < p; i++)
		B[i] = new int[p];


	B[0][0] = 1; //������ ������
	for (int i = 1; i < p; i++)
		B[0][i] = 0;

	for (int i = 0; i < p; i++)// ��������� ������
		B[p - 1][i] = p - 1;

	for (int i = 1; i < p - 1; i++)//�������� ������� � ������� �������
	{
		B[i][0] = 0;
		B[i][1] = p-1;
	}
		

	for (int i = 1; i < p - 1; i++)//����������� �����
	{
		int wi = int(   pow(xi, ((p-1) - i))) ;//xi ^(-i)
		for (int j = 2; j < p; j++)
			B[i][j] = (B[i][j - 1] * wi) % p;
	}
		

	/*
	for (int i = 0; i < p; i++)
	{ 
		for (int j = 0; j < p; j++)
			cout << B[i][j] << "	";
		cout << endl;
	}
	cout << endl;
	*/


	
	//�������������� ������ � ������� ������� � ������������ �� ��������� �������������
	//����� �������������� ���������� ���� �� ������� � �� �� ������� ������������ ��������
	
	int * spec_permute = new int [p];
	spec_permute[0] = 0;
	spec_permute[1] = 1;
	for (int i = 2; i < p; i++)
		spec_permute[i] = ( (int) pow(xi, i-1)) % p ;


	int ** Buf = new int *[p];
	for (int i = 0; i < p; i++)
		Buf[i] = new int[p];

	for (int i = 0; i < p; i++)
		for (int j = 0; j < p; j++)
			Buf[i][spec_permute[j]] =  B[i][j] ;  //�����

	for (int i = 0; i < p; i++)
		for (int j = 0; j < p; j++)
			B[i][j] =  Buf[i][j] ;  //�����
	
	
	for (int i = 0; i < p; i++)
		delete Buf[i];
	delete Buf;
	

	return B;
}


// ���������� ���������� ������� ��� ��������� �������������� ��������� �������������� ����� ��� ��������� ������������� ���������� �� ������� �������� ������� � ������ �� �������� ���� �������������� 2
element ** Get_B_matrix_for_func_to_pol_2t(GaloisField &F)
{
	int p = F.p;
	int t = F.t;
	int q = F.q;

	element w(F), xi(q-2, F), ixi(1, F), one(0, F), zero(q-1, F);

	element ** B = new element *[q];
	for (int i = 0; i < q; i++)
		B[i] = new element[q];


	B[q-1][q-1] = one; //��������� ������
	for (int i = 0; i < q-1; i++)
		B[q - 1][i] = zero;

	for (int i = 0; i < q; i++) // ������������� ������
		B[q - 2][i] = one;

	for (int i = 0; i < q - 2; i++) //������� ����� ������� � ���������� �������
	{
		B[i][q-1] = zero;
		B[i][0] = one;
	}


	for (int i = 0; i < q - 2; i++)//����������� �����
	{
		w.exp = (q-2 - i);  //xi ^(-i)

		for (int j = 1; j < q-1; j++)
			B[i][j] = B[i][j - 1] * w;
	}

	return B;
}

complex<double> ** Get_B_matrix_for_func_to_spectr(int p)
{
	//����� ������� B ����� �������� �� �������� ���������� ���������� ���� Z_p


	complex< double > ** B = new complex <double> * [p];
	for (int i = 0; i < p; i++)
		B[i] = new complex < double> [p];


	//������ ������
	for (int it = 0; it < p; it++)
		B[0][it] = 1;


	for (int it = 1; it < p; it++)
		for (int j = 0; j < p; j++)
			B[it][j] = pow(M_E, 2 * M_PI * I * (double)it * (double)j / (double)p);

	for (int it = 0; it < p; it++)
		for (int j = 0; j < p; j++)
			B[it][j] = conj(B[it][j]) / (double)p ;

	return B;
}

int ** Get_B_matrix_for_func_to_spectr_2t(GaloisField &F)
{
	//����� ������� B ����� �������� �� �������� ���������� ���������� ���� GF(2^t)
	int t = F.t;
	int q = F.q;

	element a(F), b(F);

	int ** B = new int *[q];
	for (int i = 0; i < q; i++)
		B[i] = new int[q];


	for (int i = 0; i < q; i++)
		for (int j = 0; j < q; j++)
		{
			a.exp = i;
			b.exp = j;
			B[i][j] = chi(a, b);
		}

	//�� ���������
	return B;
}

complex<double> ** Get_B_matrix_for_spectr_to_func(int p)
{
	//����� ������� B ����� �������� �� �������� ���������� ���������� ���� Z_p

	complex< double > I(0, 1);

	complex< double > ** B = new complex <double> *[p];
	for (int i = 0; i < p; i++)
		B[i] = new complex < double>[p];


	//������ ������
	for (int it = 0; it < p; it++)
		B[0][it] = 1;


	for (int it = 1; it < p; it++)
		for (int j = 0; j < p; j++)
			B[it][j] = pow(M_E, 2 * M_PI * I * (double)it * (double)j / (double)p);

	return B;
}

//��� � ���� GF(p)
void BPF(int p, int n, int ** B, int * f, int * ff)
{
	int pn = int(pow(p, n));

	int ** wa = new int*[n + 1];//������� ������
	for (int i = 0; i < n + 1; i++)
		wa[i] = new int[pn];

	for (int i = 0; i < pn; i++)
		wa[0][i] = f[i];//���������� ������� ������� � ������� ������

	for (int i = 1; i < n + 1; i++)// �������������
		for (int j = 0; j < pn; j++)
			wa[i][j] = 0;

	int i, j, k, l;
	for (i = 1; i < n + 1; i++)// n - �������� �� �������� �������� �������
	{
		l = int(pow(p, n - i));

		for (j = 0; j < pn; j++)
		{
			for (k = 0; k < p; k++)//����������� �������� ���� �����
				wa[i][j] = (wa[i][j] + wa[i - 1][ (j / (l*p))*(l*p) + j % l + k * l] * B[(j / l) % p][k] ) % p;

		}
		
	}
	/*
	for (int i = n; i < n + 1; i++)//�����
	{

		for (int j = 0; j < pn; j++)
			cout << wa[i][j] << "  ";
		cout << endl;
	}
	*/
	//int * ff = new int[pn];//���������

	for (j = 0; j < pn; j++)
		ff[j] = wa[n][j];

	for (int i = 0; i < n + 1; i++)
		delete wa[i];

	delete wa;

}

//��� � ���� GF(2^t) ��� ����������
void BPF_2t(GaloisField &F, int n, element ** B, element * f, element * ff)
{
	int p = 2;
	int t = F.t;
	int q = F.q;
	int qn = int(pow(q, n));

	element ** wa = new element*[n + 1];//������� ������
	element zero(q - 1, F);

	for (int i = 0; i < n + 1; i++)
		wa[i] = new element[qn];

	for (int i = 0; i < qn; i++)
		wa[0][i] = f[i];//���������� ������� ������� � ������� ������

	for (int i = 1; i < n + 1; i++)// ������������� �����
		for (int j = 0; j < qn; j++)
			wa[i][j] = zero;

	int i, j, k, l;
	for (i = 1; i < n + 1; i++)// n - �������� �� �������� �������� �������
	{
		l = int(pow(q, n - i));

		for (j = 0; j < qn; j++)
		{
			for (k = 0; k < q; k++)//����������� �������� ���� �����
				wa[i][j] = wa[i][j] + wa[i - 1][(j / (l*q))*(l*q) + j % l + k * l] * B[(j / l) % q][k];
		}
	}
	/*
	for (int i = n; i < n + 1; i++)//�����
	{

	for (int j = 0; j < pn; j++)
	cout << wa[i][j] << "  ";
	cout << endl;
	}
	*/
	//int * ff = new int[pn];//���������

	for (j = 0; j < qn; j++)
		ff[j] = wa[n][j];

	for (int i = 0; i < n + 1; i++)
		delete wa[i];

	delete wa;
}


//��� � ���� GF(2^t) ��� �������
void BPF_2t(GaloisField &F, int n, int ** B, int * f, int * ff)
{
	int p = 2;
	int t = F.t;
	int q = F.q;
	int qn = int(pow(q, n));

	int ** wa = new int*[n + 1];//������� ������
	int zero = 0;

	for (int i = 0; i < n + 1; i++)
		wa[i] = new int[qn];

	for (int i = 0; i < qn; i++)
		wa[0][i] = f[i];//���������� ������� ������� � ������� ������

	for (int i = 1; i < n + 1; i++)// ������������� �����
		for (int j = 0; j < qn; j++)
			wa[i][j] = zero;

	int i, j, k, l;
	for (i = 1; i < n + 1; i++)// n - �������� �� �������� �������� �������
	{
		l = int(pow(q, n - i));

		for (j = 0; j < qn; j++)
		{
			for (k = 0; k < q; k++)//����������� �������� ���� �����
				wa[i][j] = wa[i][j] + wa[i - 1][(j / (l*q))*(l*q) + j % l + k * l] * B[(j / l) % q][k];
		}

	}
	/*
	for (int i = 0; i < n + 1; i++)//�����
	{

		for (int j = 0; j < qn; j++)
		cout << wa[i][j] << "  ";
		cout << endl;
	}
	*/
	//int * ff = new int[pn];//���������

	for (j = 0; j < qn; j++)
		ff[j] = wa[n][j];

	for (int i = 0; i < n + 1; i++)
		delete[] wa[i];

	delete[] wa;

}


//���������������� ��� � ���� GF(2^t) ��� �������
void Parallel_BPF_2t(GaloisField &F, int n, int ** B, int * f, int * ff)
{
	int p = 2;
	int t = F.t;
	int q = F.q;
	int qn = int(pow(q, n));

	int ** wa = new int*[n + 1];//������� ������
	int nol = 0;

	for (int i = 0; i < n + 1; i++)
		wa[i] = new int[qn];

	for (int i = 0; i < qn; i++)
		wa[0][i] = f[i];//���������� ������� ������� � ������� ������

	for (int i = 1; i < n + 1; i++)// ������������� �����
		for (int j = 0; j < qn; j++)
			wa[i][j] = nol;

	int i, j, k, l;
	for (i = 1; i < n + 1; i++)// n - �������� �� �������� �������� �������
	{
		l = int(pow(q, n - i));

		for (j = 0; j < qn; j++)
		{

			for (k = 0; k < q; k++)//����������� �������� ���� �����
				wa[i][j] = wa[i][j] + wa[i - 1][(j / (l*q))*(l*q) + j % l + k * l] * B[(j / l) % q][k];
		}

	}
	/*
	for (int i = 0; i < n + 1; i++)//�����
	{

		for (int j = 0; j < qn; j++)
		cout << wa[i][j] << "  ";
		cout << endl;
	}
	*/
	//int * ff = new int[pn];//���������

	for (j = 0; j < qn; j++)
		ff[j] = wa[n][j];

	for (int i = 0; i < n + 1; i++)
		delete[] wa[i];

	delete[] wa;

}


void * BPF(int p, int n, complex< double > ** B, complex< double > * f, complex< double > * ff)
{
	int pn = int(pow(p, n));

	complex< double > ** wa = new complex< double > * [n + 1];//������� ������
	for (int i = 0; i < n + 1; i++)
		wa[i] = new complex< double > [pn];

	for (int i = 0; i < pn; i++)
		wa[0][i] = f[i];//���������� ������� ������� � ������� ������

	for (int i = 1; i < n + 1; i++)// �������������
		for (int j = 0; j < pn; j++)
			wa[i][j] = 0;

	int i, j, k, l;
	for (i = 1; i < n + 1; i++)// n - �������� �� �������� �������� �������
	{
		l = int(pow(p, n - i));

		for (j = 0; j < pn; j++)
		{


			for (k = 0; k < p; k++)//����������� �������� ���� �����
				wa[i][j] = wa[i][j] + wa[i - 1][(j / (l*p))*(l*p) + j % l + k * l] * B[(j / l) % p][k];

		}

	}

	/*
	for (int i = n; i < n + 1; i++)//�����
	{

		for (int j = 0; j < pn; j++)
		{
			cout.precision(2);
			cout << fixed << wa[i][j] * (double)pn << "  ";
		}
		cout << endl;
	}
	*/

	//complex <double> * ff = new complex <double> [pn];//���������

	for (j = 0; j < pn; j++)
		ff[j] = wa[n][j];

	for (int i = 0; i < n + 1; i++)
		delete wa[i];

	delete wa;


	return (ff);
}

void chi_func(int p, int n, int * f, complex <double> * chi_f)
{
	int pn = (int)pow(p, n);

	//complex <double> * chi_func = new complex <double> [pn];

	for (int j = 0; j < pn; j++)
		chi_f[j] = pow(M_E, 2 * M_PI * I * (double)f[j] / (double)p);
}

//chi_a(f)
void chi_func_2t(int n, element a, element * f, int * chi_f)
{
	int w = (int)pow(a.F->q, n);
	for (int j = 0; j < w; j++)
		chi_f[j] = chi(a, f[j]);
}


void pol_func(int p, int n, int * f, int * ff)
{
	int pn = (int)pow(p, n);
	int xi = get_primitive_root(p);

	//�������������� ������� ������� � ������������ �� ��������� �������������
	//����� �������������� ���������� ���� �� ������� ������������ ��������, � �� �� ������� 
	
	int * spec_permute = new int [p];
	spec_permute[0] = 0;
	spec_permute[1] = 1;
	for (int i = 2; i < p; i++)
		spec_permute[i] = ( (int) pow(xi, i-1)) % p ;

	for (int i = 0; i < pn; i++)
		ff[i] = f[spec_permute[i]];

	delete spec_permute;
}



// ��������� ���(���������� ���������) ����� i � p-����� ������
int Weigth1(unsigned int p, unsigned  int i)
{
	int sum = 0;

	while (i != 0)
	{
		if (i % p)
			sum ++;
		i = i / p;
	}

	return sum;
}

int Weigth1(int p, int i)
{
	int sum = 0;

	while (i != 0)
	{
		if (i % p)
			sum++;
		i = i / p;
	}

	return sum;
}

int Weigth1(int p, unsigned i)
{
	int sum = 0;

	while (i != 0)
	{
		if (i % p)
			sum++;
		i = i / p;
	}

	return sum;
}

// ��������� ���(���������� ��������� ���������) ����� i � p-����� ������ � ������ ����������������� ���� (q-1)
int Weigth1_exp(int q, int n,  int i)
{
	int sum = 0;

	for (int j = 0; j < n; j++)
	{
		if ((i % q) != (q-1))
			sum++;
		i = i / q;
	}

	return sum;
}

// ��������� ���(����� ����) ����� i � p-����� ������
int Weigth2(int p, int i)
{
	int sum = 0;

	while (i != 0)
	{
		sum += i % p;
		i = i / p;
	}

	return sum;
}

int WeightArray(int N, int * array)
{
	int wi = 0;
	for (int j = 0; j < N; j++)
		if (array[j] != 0)
			wi++;

	return wi;
}


void WritePol(int p, int n, int *  pol)
{
	//������� ��������� 
	int pn = int(pow(p, n));
	int wi;

	for (int i = 0; i < pn; i++)
	{
		wi = i;
		if (pol[i] != 0)
		{
			cout << pol[i];
			for (int j = 0; j < n; j++)
			{
				if (pol[i])
					if (wi % p)
					{
						cout << "x_" << j << "^" << wi % p;
					
					}
					
				wi = wi / p;
			}

			if (pol[i])
				if (i < pn - 1)
					cout << " + ";
		}
	}
	cout << endl << endl << endl;
}


int Algebraic_degree(int p, int n, int * polinom)
{
	int max = 0, w;

	for (int i = 0; i < int(pow(p, n)); i++)
	{
		if (polinom[i])
		{
			w = Weigth2(p, i);
			if (w > max)
				max = w;
		}
		
	}
		
	return max;
}

//��������� �������������� ������� ����������
int Algebraic_degree_2t(GaloisField &F, int n, element * polinom)
{
	int p = 2;
	int t = F.t;
	int q = F.q;
	int max = 0, w, wi;

	for (int i = 0; i < int(pow(q, n)); i++)
	{
		if (polinom[i].exp != q-1)//���� ����������� �� ����� ����, ��...
		{			
			w = 0;
			wi = i;
			for (int j = 0; j < n; j++)
			{
				w +=  (wi % q + 1) % q;
				wi = wi / q;
			}

			if (w > max)
				max = w;
		}
	}
	return max;
}

//��������� ������� ������������ dl(f) ��� ����� ����
int Nonliniarity_degree_2t(GaloisField &F, int n, element * polinom)
{
	int p = 2;
	int t = F.t;
	int q = F.q;
	int max = 0, w, wi;

	for (int i = 0; i < int(pow(q, n)); i++)
	{
		if (polinom[i].exp != q - 1)//���� ����������� �� ����� ����, ��...
		{//���������� ������� ������������ ��������� � ������������ � ������������ ����������� (��� ����� ����)
			w = 0;
			wi = i;
			for (int j = 0; j < n; j++)
			{
				w += Weigth2(2, (wi % q + 1) % q);
				wi = wi / q;
			}

			if (w > max)
				max = w;
		}
	}

	return max;
}

int AbsMax(int N, int * spectr)
{
	int max = 0;

	for (int i = 0; i < N; i++)
		if (max < (int) abs(spectr[i]))
			max = (int) abs(spectr[i]);

	return max;
}

int AbsMax(int N1, int N2, int ** spectr)
{
	int max = 0;

	for (int i = 0; i < N1; i++)
		for (int j = 0; j < N2; j++)
		{
			if (max < (int) abs(spectr[i][j]))
			max = (int) abs(spectr[i][j]);
		}
	return max;
}

double AbsMax(int pn, complex <double> * spectr)
{
	double max = 0;

	for (int i = 0; i < pn; i++)
		if (max < abs(spectr[i]))
			max = abs(spectr[i]);

	return max;
}

int Resilient_index(int p, int n, complex <double> * spectr)
{
	//���������� ������ ������������ ������� �� �� �������

	int w=-1, pn = (int)pow(p, n);

	bool * ri = new bool[n+1];
	for (int i = 0; i < n+1; i++)
		ri[i] = true;


	for (int i = 0; i < pn; i++)
	{
		if (abs(spectr[i]) > 0.00001)
			ri[Weigth1(p, i)] = false;
	}

	for (int i = 0; i < n+1; i++)
		if (ri[i] == false)
		{
			w = i - 1;
			break;
		}

	delete ri;
	return w;
}

int Resilient_index_2t(int q, int n, int ** spectr)
{
	//���������� ������ ������������ ������� ��� ����� GF(2^t) �� �� �������

	int w = -1, qn = (int)pow(q, n);

	bool * ri = new bool[n + 1];
	for (int i = 0; i < n + 1; i++)
		ri[i] = true;

	for (int j = 0; j < q-1; j++)
		for (int i = 0; i < qn; i++)
		{
			if (spectr[j][i] != 0)
				ri[Weigth1_exp(q, n, i)] = false;
		}

	for (int i = 0; i < n + 1; i++)
		if (ri[i] == false)
		{
			w = i - 1;
			break;
		}

	delete ri;
	return w;
}

int Correalation_immunity_index_2t(int q, int n, int * spectr)
{
	//���������� ������ �������������� ���������� ������� ��� ����� GF(2^t) �� �� �������

	int w = 0, qn = (int)pow(q, n);

	bool * ri = new bool[n + 1];
	for (int i = 0; i < n + 1; i++)
		ri[i] = true;


	for (int i = 1; i < qn; i++)
	{
		if (spectr[i] != 0)
			ri[Weigth1_exp(q, n, i)] = false;
	}

	for (int i = 1; i < n + 1; i++)
		if (ri[i] == false)
		{
			w = i - 1;
			break;
		}

	delete ri;
	return w;
}

int Correalation_immunity_index(int p, int n, complex <double> * spectr)
{
	//���������� ������ ������������ ������� �� �� �������

	int w = 0, pn = (int)pow(p, n);

	bool * ri = new bool[n + 1];
	for (int i = 0; i < n + 1 ; i++)
		ri[i] = true;


	for (int i = 1; i < pn; i++)
	{
		if (abs(spectr[i]) > 0.00001)
			ri[Weigth1(p, i)] = false;
	}

	for (int i = 1; i < n + 1 ; i++)
		if (ri[i] == false)
			break;
		else
			w ++ ;

	delete ri;
	return w;
}

void Get_random_func(mt19937 &gen, int p, int n, int * f)
{
	int pn = (int)pow(p, n);

	for (int i = 0; i < pn; i++)
		f[i] = gen() % p;

}

void Get_random_permute_Zp(mt19937 &gen, int p, int * permute)
{
	int *w = new int[p];
	int t = p, rand;
	

	for (int i=0; i<p; i++)
		w[i]=i;

	for (int i=0; i<p; i++)
	{
		rand = gen() % t;
		permute[i] = w[rand];
		t--;
		w[rand] = w[t];
	}

}

//������ ������� - ����������� �� ��������� �����������  f(x_1,...,x_n) = g_{n+1}(g_1(x_1)+...+g_n(x_n)), ��� g - ����������� ���� Z_p
void Get_random_quazyfunc(mt19937 &gen, int p, int n, int * f)
{

	int pn = (int)pow(p, n);
	int ws, wi;
	int **g = new int * [n+1];
	for (int i=0; i<n+1; i++)
		g[i] = new int[p];

	for (int i=0; i<n+1; i++)
		Get_random_permute_Zp(gen, p, g[i]);


	for (int i=0; i<pn; i++)
	{
		wi = i;
		ws = 0 ;
		for (int j=0; j<n; j++)
		{
			ws = ( ws + g[j][wi % p]) % p;
			wi = wi / p;
		}

		f[i] = g[n][ws];
	}


	for (int i=0; i<n+1; i++)
		delete g[i];
	delete g;
}

// ������ ��������� ������� ��� ����� GF(p)
void Test(int p, int n, int * f)
{
	//���� ����� �������

	int pn = (int)pow(p, n);

	cout << "Function:\n";
	for (int i = 0; i < pn; i++)
		cout << f[i] << " ";
	cout << endl;

	int deg, ri;

	int * pol = new int[pn];
	int * pol_f = new int[pn];
	complex <double> * spectr = new complex <double>[pn];
	complex <double> * chi_f = new complex <double>[pn];

	complex< double > ** B_spectr = Get_B_matrix_for_func_to_spectr(p);
	int ** B_polinom = Get_B_matrix_for_func_to_pol(p);

	//pol_func(p, n, f, pol_f);
	
	//BPF(p, n, B_polinom, pol_f, pol);

	BPF(p, n, B_polinom, f, pol);

	cout << "Polinom:\n";
	for (int i = 0; i < pn; i++)
		cout << pol[i] << " ";
	cout << endl;

	deg = Algebraic_degree(p, n, pol);

	cout << "Degree = " << deg << endl;

	chi_func(p, n, f, chi_f);

	BPF(p, n, B_spectr, chi_f, spectr);

	cout << "\nSpectr:\n";
	for (int i = 0; i < pn; i++)
		cout << abs(spectr[i]) *pn << endl;
	cout << endl;

	cout << "C(f) = " << AbsMax(pn, spectr)*pn << endl << endl;


	ri = Correalation_immunity_index(p, n, spectr);

	cout << "Correlation immunity = " << ri << endl;


	cout << "Arg\tFunc\tPol\n";
	int w;
	for (int i = 0; i < pn; i++)
	{
		w = i;
		for (int j = 0; j < n; j++)
		{
			cout << w%p;
			w = w / p;
		}
		cout << "\t";
		cout <<f[i] << "\t";
		cout << pol[i] << "\n";

	}


	system("Pause");
}

// ������ ��������� ������� ��� ����� GF(2^t)
void Test_2t(GaloisField &F, int n, element * f)
{
	// q = 2^t - �������� ����
	// n - ���������� ����������
	//��������� ���������� C(f), �������������� ������� deg(f), ������� ������������ dl(f) � ������ ������������ res(f) ��� ������� ��� ����� GF(2^t)

	int p = 2;
	int t = F.t;
	int q = F.q;

	cout << "p = " << 2 << ", t = " << t << ", q = " << q << ", n = " << n << endl << endl;

	int qn = (int)pow(q, n); // �������� ��������� ���������� ������� f

	/*
	//����� ������
	for (int i = 0; i < qn ; i++)
		cout << f[i].exp << " ";
	cout << endl;
	*/

	int **chi_f = new int *[q-1], **spectr_f = new int *[q-1];

	for (int i=0; i< q-1; i++)
	{
		chi_f[i] = new int[qn]; // ����� ����� ��������� ������� �������� ������� \chi(a*f(x))
		spectr_f[i] = new int[qn]; // ����� ����� ��������� �����-�������������� ������������ �(af, <u,x>) ������� \chi(a*f(x))
	}


	element * pol = new element [qn];

	element ** B_polinom = Get_B_matrix_for_func_to_pol_2t(F);

	BPF_2t(F, n, B_polinom, f, pol);

	int deg = Algebraic_degree_2t(F, n, pol);
	cout << "deg = " << deg << " <= " << n*(q-1) << endl;

	int dl = Nonliniarity_degree_2t(F, n, pol);
	cout << "dl = " << dl << " <= " << n*t << endl;

	int ** B_spectr = Get_B_matrix_for_func_to_spectr_2t(F);



	// ���������� ���������� ������������ ������������� ��� ������ ����������
	# pragma omp parallel for shared(chi_f,spectr_f)
	for (int i=0; i<q-1; i++)
	{
		element a(i, F);
		chi_func_2t(n, a, f, chi_f[i]);
		
		BPF_2t(F, n, B_spectr, chi_f[i], spectr_f[i]);
	}
	


	/*
	cout << endl;
	for (int i = 0; i < q; i++)
	{
		for (int j = 0; j < q; j++)
			cout << B_spectr[i][j] << " ";
		cout << endl;
	}
	*/

	int res = Resilient_index_2t(q, n, spectr_f);

	cout << "res = " << res << endl;

	cout << endl;
	cout << "Ziegenthaler's inequality for\n";
	cout << "deg: " << deg << " + " << res;
	if (deg + res == n*(q - 1) - 1)
		cout << " = ";
	else
		cout << " < ";
	cout << n*(q - 1) - 1 << endl;


	cout << "dl: " << dl << " + " << res;
	if (dl + res == n*t - 1)
		cout << " = ";
	else
		cout << " < ";
	cout << n*t - 1 << endl;


	int nl = AbsMax(q-1, qn, spectr_f);
	cout << endl;
	cout << "nl : " << pow(qn, 1.0 / 2) << " <= " << nl << " <= " <<  qn << endl;

	
	/*
	//����� �������
	for (int i=0; i<q-1; i++)
	{
		for (int j=0; j<qn; j++)
			cout << spectr_f[i][j] << " ";
		cout << endl<<endl;
	}
		
	cout <<endl<<endl<<endl;
	*/


	/*
	for (int i = 0; i < qn; i++)
	{
		cout << spectr_f[i] << " ";
		if ((i+1) % (qn / q) == 0)
			cout << endl;
	}
	cout << endl;
	*/

	delete chi_f;
	delete spectr_f;
	delete pol;


}

//���������������� Test_2t
bool Work_Test_2t(GaloisField &F, int n, element * f)
{
	int p = 2;
	int t = F.t;
	int q = F.q;
	int qn = (int)pow(q, n); // �������� ��������� ���������� ������� f

	int **chi_f = new int *[q-1], **spectr_f = new int *[q-1];

	for (int i=0; i< q-1; i++)
	{
		chi_f[i] = new int[qn]; // ����� ����� ��������� ������� �������� ������� \chi(a*f(x))
		spectr_f[i] = new int[qn]; // ����� ����� ��������� �����-�������������� ������������ �(af, <u,x>) ������� \chi(a*f(x))
	}
	
	element * pol = new element [qn];

	element ** B_polinom = Get_B_matrix_for_func_to_pol_2t(F);

	BPF_2t(F, n, B_polinom, f, pol);

	int ** B_spectr = Get_B_matrix_for_func_to_spectr_2t(F);


	int deg = Algebraic_degree_2t(F, n, pol);
	//cout << "deg = " << deg << " <= " << n*(q-1) << endl;

	int dl = Nonliniarity_degree_2t(F, n, pol);
	//cout << "dl = " << dl << " <= " << n*t << endl;

	

	element a(F);
	for (int i=0; i<q-1; i++)
	{
		a.exp = i;
		chi_func_2t(n, a, f, chi_f[i]);
		
		BPF_2t(F, n, B_spectr, chi_f[i], spectr_f[i]);
	}
	
	int res = Resilient_index_2t(q, n, spectr_f);

	int nl = AbsMax(q-1, qn, spectr_f);


	for (int i=0; i< q-1; i++)
	{
		delete[] chi_f[i]; // ����� ����� ��������� ������� �������� ������� \chi(a*f(x))
		delete[] spectr_f[i]; // ����� ����� ��������� �����-�������������� ������������ �(af, <u,x>) ������� \chi(a*f(x))
	}

	delete[] chi_f;
	delete[] spectr_f;
	delete[] pol;


	if (dl + res == n*t - 1)
	{
		cout << "deg = " << deg << " <= " << n*(q-1) << endl;
		cout << "dl = "<< dl <<", res = "<< res << endl;
		cout << "nl : " << pow(qn, 1.0 / 2) << " <= " << nl << " <= " <<  qn << endl;
		return true;
	}

	return false;
}

//���������������� Test_2t
bool Work2_Test_2t(GaloisField &F, int n, element * f)
{
	int p = 2;
	int t = F.t;
	int q = F.q;
	int qn = (int)pow(q, n); // �������� ��������� ���������� ������� f

	int **chi_f = new int *[q-1], **spectr_f = new int *[q-1];

	for (int i=0; i< q-1; i++)
	{
		chi_f[i] = new int[qn]; // ����� ����� ��������� ������� �������� ������� \chi(a*f(x))
		spectr_f[i] = new int[qn]; // ����� ����� ��������� �����-�������������� ������������ �(af, <u,x>) ������� \chi(a*f(x))
	}
	
	element * pol = new element [qn];

	element ** B_polinom = Get_B_matrix_for_func_to_pol_2t(F);

	BPF_2t(F, n, B_polinom, f, pol);

	int ** B_spectr = Get_B_matrix_for_func_to_spectr_2t(F);


	int deg = Algebraic_degree_2t(F, n, pol);
	//cout << "deg = " << deg << " <= " << n*(q-1) << endl;

	/*
	for (int i=0; i< qn; ++i)
		cout << pol[i].exp << " ";
	cout << endl;
	*/

	int dl = Nonliniarity_degree_2t(F, n, pol);
	//cout << "dl = " << dl << " <= " << n*t << endl;

	

	element a(F);
	# pragma omp parallel for shared(chi_f,spectr_f)
	for (int i=0; i<q-1; i++)
	{
		a.exp = i;
		chi_func_2t(n, a, f, chi_f[i]);
		
		BPF_2t(F, n, B_spectr, chi_f[i], spectr_f[i]);
	}
	

	int res = Resilient_index_2t(q, n, spectr_f);

	int nl = AbsMax(q-1, qn, spectr_f);


	for (int i=0; i< q-1; i++)
	{
		delete[] chi_f[i]; // ����� ����� ��������� ������� �������� ������� \chi(a*f(x))
		delete[] spectr_f[i]; // ����� ����� ��������� �����-�������������� ������������ �(af, <u,x>) ������� \chi(a*f(x))
	}

	delete[] chi_f;
	delete[] spectr_f;
	delete[] pol;


	//if (dl + res == n*t - 1)
	{
		cout << "deg = " << deg << " <= " << n*(q-1) << endl;
		cout << "dl = "<< dl << " <= " << n*t <<  endl;
		cout <<"res = "<< res << endl;
		cout << "nl = " <<  nl << endl;
		if (res + dl == n*t-1)
			cout << "\t\t\tOptimal\n";
		return true;
	}

	return false;
}

int Linear_Times(int a, int b) // ��������� ��������� �������� �������� �������������� int'���
{
	int vec = a & b, w=0;
	//������ ����� ������� ��� ������� vec �� ������ 2
	while (vec)
	{
		if (vec % 2)
			w++;
		vec = vec/2;
	}

	return w % 2;
}


//���������� ��������� � ������������������ ������� �������� ������ ������ �� ����
int Get_next_vector_this_weight(int *array, int len)
{
	//���� ����� ����� ������ �������, ������� ����� �������� ������, �������� ��, � ��� ������� ������� ����� �� ��� ���������� � ���

	int w, c;

	for (int i = len - 1; i >=0; i--)
	{
		if(array[i] == 1)
			if (array[i + 1] == 0)
			{
				array[i] = 0;
				array[i+1] = 1;

				c = 0;// ��������� ������� ������ ������
				for (int j = i + 2; j < len; j++)
					c += array[j];

				//������ ��������� ��

				for (int j = i + 2; j < len; j++)
					if (j - i - 1 <= c)
						array[j] = 1;
					else
						array[j] = 0;
				return 1;
			}
	}
	//��� ������, ��� �������� ������ �������� ��������� � ������������������ ������� � ����� ���� ������ ���
	return 0;
}