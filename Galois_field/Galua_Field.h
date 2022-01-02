// ������ �.�.
// ���������� �������� �����
// �������� �������������� � ����������������� (��� ������� ������������ ��������) �����, � �� � ����������(��� ����������)

using ui = unsigned unsigned int;

#include <iostream>
//#include <fstream>
//#include <map>
using namespace std;



ui get_primitive_polinoms(ui p, ui n)//���������� ����������� ��������� ������� n ��� ����� p^n
{
	/*
	if (p==2)//������������ ����� ������ ������� �� �������
	{
		if (n==2) {static ui p2n2[3]={1,1,1};return p2n2;}
		if (n==3) {static ui p2n3[4]={1,1,0,1};return p2n3;}// ��������� x^3+x+1
		if (n==4) {static ui p2n4[5]={1,1,0,0,1};return p2n4;}
		if (n==5) {static ui p2n5[6]={1,0,1,0,0,1};return p2n5;}
		if (n==6) {static ui p2n6[7]={1,1,0,0,0,0,1};return p2n6;}
		if (n==7) {static ui p2n7[8]={1,1,0,0,0,0,0,1};return p2n7;}
		if (n==8) {static ui p2n8[9]={1,0,1,1,1,0,0,0,1};return p2n8;}
		if (n==9) {static ui p2n9[10]={1,0,0,0,1,0,0,0,0,1};return p2n9;}
		if (n==10) {static ui p2n10[11]={1,0,0,1,0,0,0,0,0,0,1};return p2n10;}
	}
	if (p==3)
	{
		if (n==4) {static ui p3n4[5]={2,1,0,0,1}; return p3n4;}
	}
	if (p==5)
	{
		if (n==2) {static ui p5n2[3]={2,1,1}; return p5n2;}
	}
	if (p==7)
	{
		if (n==2) {static ui p7n2[3]={3,1,1}; return p7n2;}
	}

	cout << "Need primitive polynomial!/n";
	system("PAUSE");
	return 0 ;
	*/

	if (p == 2)//������������ ����� ������ ������� �� �������
	{
		if (n == 2) return 7;
		if (n == 3) return 11; // ��������� x^3+x+1
		if (n == 4) return 19;
		if (n == 5) return 37;
		if (n == 6) return 67;
		if (n == 7) return 131;
		if (n == 8) return 285;
		if (n == 9) return 529;
		if (n == 10) return 1033;
	}


	cout << "Need primitive polynomial!/n";
	system("PAUSE");
	return 0;
}

ui get_primitive_root(ui p)
{
	ui xi = 0;

	if (p == 2)
		xi = 1;

	if ((p == 3) || (p == 5) || (p == 11) || (p == 13))
		xi = 2;

	if ((p == 7) || (p == 17))
		xi = 3;

	if (xi == 0)
	{
		cout << "Warning! No xi!\n";
		system("pause");
	}


	return xi;
}

// �������� �������� �� ������ p
ui MagicPlus(ui a, ui b, ui p)
{
	ui c = 0, t = 0;
	while ((a > 0) || (b > 0))
	{
		c += (((a % p) + (b % p)) % p) * (ui)pow(p, t);
		t++;
		a = a / p;
		b = b / p;
	}
	return c;
}

//���������� ���������� w �� ������ ���������� root � ���� GF(p)
ui Reduce(ui w, ui root, ui p, ui q)
{
	// ���������� ������� root � w �� ������ p �� ��� ���, ���� w �� ������ ������ ��� q

	while (w >= q)
		w = MagicPlus(w, root, p);

	return w;
}

class GaloisField {

public:
	ui p, t, q; // GF(p^t=q)

	ui *PolToExp, *ExpToPol; // ����� ����� �������� ��������������� ������� ��������� ����,
	// p-����� ������ ui � ������ ������ ���� ���������, ����������� �������������

	ui** KallyTable;// ����� ����� �������� ������� ����� �� �������� ��������� � ����������������� �����
	
	GaloisField(ui _p, ui _t)
	{
		p = _p;
		t = _t;
		q = (ui)pow(p, t);
		 
		PolToExp = new ui[q]; //��������������� �����
		ExpToPol = new ui[q];

		

		//��������� ����� ����� ���������� ������, ��� p-����� ������ � ����� �����������

		PolToExp[0] = q - 1;//���� ����
		ExpToPol[q - 1] = 0;

		PolToExp[1] = 0;//������� ����
		ExpToPol[0] = 1;

		if (t == 1)// � ������ �������� ���� 
		{
			ui w = 1, root = get_primitive_root(p); //����������� �������

			for (ui i = 1; i < q - 1; i++)
			{
				w = (w * root) % p;
				PolToExp[w] = i;//����������� ������� ����
				ExpToPol[i] = w;
			}
		}
		else // � ������ �� �������� ����
		{ 
			ui w = 1, root = get_primitive_polinoms(p, t); //����������� ������� ���� ���������� ����������� "x" (10 � p-����� ������, p � ����������)

			for (ui i = 1; i < q - 1; i++)
			{
				w = Reduce(w * p, root, p, q);
				PolToExp[w] = i;//����������� ������� ����
				ExpToPol[i] = w;
			}
		}

		//������� ����� �� ��������

		KallyTable = new ui* [q];
		for (ui i = 0; i < q; i++)
			KallyTable[i] = new ui[q];

		for (ui i = 0; i < q; i++)//i � j ����������� �������� ����, ������ �� ����������
			for (ui j = 0; j < q; j++)
				KallyTable[i][j] = PolToExp[MagicPlus(ExpToPol[i], ExpToPol[j], p)];

		cout << "\nGalois field GF("<< q<<") successfully initiated!\n";

	}
};

class element //������� ��������� ����
{
public:
	ui exp; //����������������� ����� �������������, �� ���� �������� �������� �� ��������� ������������ �������� "�", ���� ���� ������������ ����� q-1.
	GaloisField *F; 

	element() {}

	element(ui _exp, GaloisField &_F)
	{
		exp = _exp;
		F = &_F;
	}

	element(ui _exp, GaloisField *_F)
	{
		exp = _exp;
		F = _F;
	}

	element(GaloisField* _F)
	{
		exp = _F->q-1;  // ���� ����
		F = _F;
	}

	element(GaloisField &_F)
	{
		exp = _F.q - 1;  // ���� ����
		F = &_F;
	}

	element Power(ui n) // ���������� �������� � ������� n
	{
		if (this->exp == F->q - 1) 
		{
			element w(F->q - 1, F);
			return w;
		}
		else
		{
			element w( (this->exp * n) % (F->q - 1), F);
			return w;
		}
	}

	friend const element operator *(const element& l, const element& r); // ��������� ���������, ��� �������� �� ����������� �� ������ q-2

	friend const element operator+(const element& left, const element& right);
};

const element operator * (const element& l, const element& r)// ��������� ���������, ��� �������� �� ����������� �� ������ q-2
{
	if ((l.exp == l.F->q - 1) || (r.exp == r.F->q - 1))
		return element(l.F->q - 1, l.F);// ��������� �� ����
	else
		return element((l.exp + r.exp) % (l.F->q - 1), l.F);
}

const element operator + (const element& l, const element& r)
{
	return element (l.F->KallyTable[l.exp][r.exp], l.F);
}

//������� ����������� ����� � ���� GF(2^t)
ui Tr(element x)
{
	element w( x.F->q - 1, x.F);

	for (ui i = 0; i < x.F->t; i++)
	{
		w = w + x;
		x = x * x;
	}

	if (w.exp == 0)
		return 1;//���������� ���������� 0 � 1, � �� ����������������

	if (w.exp == x.F->q - 1)
		return 0;

	cout << "Warning. Trace!";
	system("Pause");

	return(0);
}

//���������� �������� chi_a(b) ���� GF(2^t)
ui chi(element a, element b)
{
	ui w = Tr(a * b); //1 ��� 0

	if (w == 0)
		return 1;

	if (w == 1)
		return -1;

	cout << "Warning. Chi2^t!";
	system("Pause");

	return(0);
}

ui VectorPolToExp(GaloisField &F, ui n, ui x) // ��������� ������ ������� PolToExp
{
	// q-����� ������ ����� n, �������������� �������� ����, �������������� ������������
	// ����������� q-����� ������, �������������� �������� ����, �������������� ��������� ������������ ��������

	ui w = 0, q = F.q;
	for (ui i = 0; i < n; i++)
	{
		w += F.PolToExp[x % q] * (ui)pow(q, i);
		x = x / q;
	}
	return w;
}

ui VectorExpToPol(GaloisField & F, ui n, ui x) // ��������� ������ ������� ExpToPol
{
	ui w = 0, q = F.q;
	for (ui i = 0; i < n; i++)
	{
		w += F.ExpToPol[x % q] * (ui)pow(q, i);
		x = x / q;
	}
	return w;
}


void Convert_exp_func_on_nomera_2t(GaloisField& F, ui n, ui *f, ui *newf)
{
	// ������� f ������ �� ���������-�����������, � ������� newf ������ �� ���������-�����������

	// ���� ��� � ���
	// ��� ������ � ��������� �����(t>1) ������������ ���������������� ������������� ���������
	// ��� ��������� ��������� ������� �� ����������������� ������������� f(xi^{i_1},...,xi^{i_n})=xi^{e * f(i_1,...,i_n)}
	// � ������������� ���������� �� ����������� ��������� ���� - ����������� "��������" f(i_1,...,i_n) = j, ��� �������� ������������� i_s ������������� ����������-�������� ����
	// ��� ����� ������������ ��������������� �������

	ui qn = (ui)pow(F.q, n);
	ui w = 0;

	for (ui i = 0; i < qn; i++)
	{
		//cout << VectorPolToExp(F, n, i) << endl;
		//cout << f[VectorPolToExp(F, n, i)] << endl;
		//cout << F.ExpToPol[f[VectorPolToExp(F, n, i)]] << endl;
		//cout << endl;

		newf[i] = F.ExpToPol[  f[VectorPolToExp(F, n, i)]  ];
	}
}

void Convert_exp_func_on_nomera_2t(GaloisField& F, ui n, element* f, ui* newf)
{
	// ������� f ������ �� ���������-�����������, � ������� newf ������ �� ���������-�����������

	// ���� ��� � ���
	// ��� ������ � ��������� �����(t>1) ������������ ���������������� ������������� ���������
	// ��� ��������� ��������� ������� �� ����������������� ������������� f(xi^{i_1},...,xi^{i_n})=xi^{e * f(i_1,...,i_n)}
	// � ������������� ���������� �� ����������� ��������� ���� - ����������� "��������" f(i_1,...,i_n) = j, ��� �������� ������������� i_s ������������� ����������-�������� ����
	// ��� ����� ������������ ��������������� �������

	ui qn = (ui)pow(F.q, n);
	ui w = 0;

	for (ui i = 0; i < qn; i++)
	{
		//cout << VectorPolToExp(F, n, i) << endl;
		//cout << f[VectorPolToExp(F, n, i)].exp << endl;

		newf[i] = F.ExpToPol[f[VectorPolToExp(F,n, i)].exp];
	}	
}

/*
void Recursive_find_basis_2t(ui k, ui t, ui q, ui* basis, bool* old_ar)
{
	// k - ������� �������� - ���������� �������� ��������� ������
	// t - �������� ������
	element a(2, t), b(2, t);


	bool* ar = new bool[q];
	ui bound;


	if (k == t) // ���� ����� ��������, �� // ������� ��� � ����
	{
		string fn = "All_Basis_2^";
		fn = fn + to_string(t) + ".txt";
		ofstream res(fn, ios_base::app);

		for (ui i = 0; i < t; ++i)
			res << basis[i] << " ";
		res << endl;

		res.close();
	}



	if (k == 0)
		bound = 0;
	else
		bound = basis[k - 1] + 1;


	if (k < t)// ���� ����� �� ������, �� ������� ���� �������
	{
		for (ui i = bound; i < q - t + k; ++i) // ������� ������� ��� ������ �� ��� ������� �� ���������� ��������� ��� ���������
		{
			if (old_ar[i] == false)
				basis[k] = i;
			else
				continue;


			a.exp = i;

			//�������� ������ �����
			for (ui j = 0; j < q; ++j)
				ar[j] = old_ar[j];

			//������� ����� ������������ �������� ������ ����� ������� �� ����� ���������� ���������� ���������� ����
			//�������� ����� �����
			for (ui j = 0; j < q; ++j)
			{
				if (old_ar[j] == true)
				{
					b.exp = j;
					b = a + b;
					ar[b.exp] = true;

				}

			}
			ar[i] = true; // ��������� ������� ���� ��� �������


			//����� ��������
			Recursive_find_basis_2t(k + 1, t, q, basis, ar);
		}
	}
}
*/

/*
// ���������� ���� ������� �� P_{P_0} ��� p=2
void Build_All_Basis_2t(ui t)
{
	//�������� ��� ������ ���������� ������������ P_{P_0}
	// � ��������� �� ������������ ��������
	//������� ����
	InitGF(2, t);
	ui q = (ui)pow(2, t);

	cout << "Start build basis in GF(" << q << ")" << endl;

	ui* basis = new ui[t]; // ��� �������� ������

	for (ui i = 0; i < t; ++i)
		basis[i] = i;

	bool* ar = new bool[q];
	for (ui i = 0; i < q; ++i)
		ar[i] = false;

	Recursive_find_basis_2t(0, t, q, basis, ar);
}
*/



