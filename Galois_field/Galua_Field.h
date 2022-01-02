// Бугров А.Д.
// Реализация конечных полей
// элементы представляются в мультипликативной (как степени примитивного элемента) форме, а не в аддитивной(как многочлены)

using ui = unsigned unsigned int;

#include <iostream>
//#include <fstream>
//#include <map>
using namespace std;



ui get_primitive_polinoms(ui p, ui n)//возвращает примитивный многочлен степени n над полем p^n
{
	/*
	if (p==2)//коэффициенты нужно писать начиная со старших
	{
		if (n==2) {static ui p2n2[3]={1,1,1};return p2n2;}
		if (n==3) {static ui p2n3[4]={1,1,0,1};return p2n3;}// многочлен x^3+x+1
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

	if (p == 2)//коэффициенты нужно писать начиная со старших
	{
		if (n == 2) return 7;
		if (n == 3) return 11; // многочлен x^3+x+1
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

// сложение разрядов по модулю p
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

//приведение многочлена w по модулю многочлена root в поле GF(p)
ui Reduce(ui w, ui root, ui p, ui q)
{
	// прибавляем разряды root к w по модулю p до тех пор, пока w не станет меньше чем q

	while (w >= q)
		w = MagicPlus(w, root, p);

	return w;
}

class GaloisField {

public:
	ui p, t, q; // GF(p^t=q)

	ui *PolToExp, *ExpToPol; // Здесь будет хранится логарифмические таблицы элементов поля,
	// p-ичная запись ui в данном случае есть многочлен, отображения взаимообратны

	ui** KallyTable;// Здесь будет хранится таблица Келли по сложению элементов в мультипликативной форме
	
	GaloisField(ui _p, ui _t)
	{
		p = _p;
		t = _t;
		q = (ui)pow(p, t);
		 
		PolToExp = new ui[q]; //Логарифмические шкалы
		ExpToPol = new ui[q];

		

		//многочлен здесь будем кодировать числом, его p-ичная запись и будет многочленом

		PolToExp[0] = q - 1;//ноль поля
		ExpToPol[q - 1] = 0;

		PolToExp[1] = 0;//единица поля
		ExpToPol[0] = 1;

		if (t == 1)// в случае простого поля 
		{
			ui w = 1, root = get_primitive_root(p); //примитивный элемент

			for (ui i = 1; i < q - 1; i++)
			{
				w = (w * root) % p;
				PolToExp[w] = i;//примитивный элемент поля
				ExpToPol[i] = w;
			}
		}
		else // в случае не простого поля
		{ 
			ui w = 1, root = get_primitive_polinoms(p, t); //примитивный элемент поля задающийся многочленом "x" (10 в p-ичной записи, p в десятичной)

			for (ui i = 1; i < q - 1; i++)
			{
				w = Reduce(w * p, root, p, q);
				PolToExp[w] = i;//примитивный элемент поля
				ExpToPol[i] = w;
			}
		}

		//Таблица Келли по сложению

		KallyTable = new ui* [q];
		for (ui i = 0; i < q; i++)
			KallyTable[i] = new ui[q];

		for (ui i = 0; i < q; i++)//i и j перечисляют элементы поля, точнее их экспоненты
			for (ui j = 0; j < q; j++)
				KallyTable[i][j] = PolToExp[MagicPlus(ExpToPol[i], ExpToPol[j], p)];

		cout << "\nGalois field GF("<< q<<") successfully initiated!\n";

	}
};

class element //Элемент конечного поля
{
public:
	ui exp; //мультипликативная форма представления, то есть логарифм элемента по основанию примитивного элемента "х", ноль поля обозначается через q-1.
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
		exp = _F->q-1;  // Ноль поля
		F = _F;
	}

	element(GaloisField &_F)
	{
		exp = _F.q - 1;  // Ноль поля
		F = &_F;
	}

	element Power(ui n) // возведение элемента в степень n
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

	friend const element operator *(const element& l, const element& r); // умножение элементов, как сложение их показателей по модулю q-2

	friend const element operator+(const element& left, const element& right);
};

const element operator * (const element& l, const element& r)// умножение элементов, как сложение их показателей по модулю q-2
{
	if ((l.exp == l.F->q - 1) || (r.exp == r.F->q - 1))
		return element(l.F->q - 1, l.F);// умножение на ноль
	else
		return element((l.exp + r.exp) % (l.F->q - 1), l.F);
}

const element operator + (const element& l, const element& r)
{
	return element (l.F->KallyTable[l.exp][r.exp], l.F);
}

//Функция абсолютного следа в поле GF(2^t)
ui Tr(element x)
{
	element w( x.F->q - 1, x.F);

	for (ui i = 0; i < x.F->t; i++)
	{
		w = w + x;
		x = x * x;
	}

	if (w.exp == 0)
		return 1;//возвращаем нормальные 0 и 1, а не экспоненциальные

	if (w.exp == x.F->q - 1)
		return 0;

	cout << "Warning. Trace!";
	system("Pause");

	return(0);
}

//Аддитивный характер chi_a(b) поля GF(2^t)
ui chi(element a, element b)
{
	ui w = Tr(a * b); //1 или 0

	if (w == 0)
		return 1;

	if (w == 1)
		return -1;

	cout << "Warning. Chi2^t!";
	system("Pause");

	return(0);
}

ui VectorPolToExp(GaloisField &F, ui n, ui x) // Векторный аналог функции PolToExp
{
	// q-ичный вектор длины n, представляющий элементы поля, закодированные многочленами
	// переводится q-ичный вектор, представляющий элементы поля, закодированные степенями примитивного элемента

	ui w = 0, q = F.q;
	for (ui i = 0; i < n; i++)
	{
		w += F.PolToExp[x % q] * (ui)pow(q, i);
		x = x / q;
	}
	return w;
}

ui VectorExpToPol(GaloisField & F, ui n, ui x) // Векторный аналог функции ExpToPol
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
	// Функция f задана на элементах-экспонентах, а функция newf задана на элементах-многочленах

	// дело вот в чем
	// при работе в непростых полях(t>1) используется экспоненциальное представление элементов
	// эта процедура переводит функцию из экспоненциального представления f(xi^{i_1},...,xi^{i_n})=xi^{e * f(i_1,...,i_n)}
	// в представление основанное на кодировании элементов поля - многочленов "номерами" f(i_1,...,i_n) = j, где двоичное представление i_s соответствует многочлену-элементу поля
	// для этого используется логарифмическая таблица

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
	// Функция f задана на элементах-экспонентах, а функция newf задана на элементах-многочленах

	// дело вот в чем
	// при работе в непростых полях(t>1) используется экспоненциальное представление элементов
	// эта процедура переводит функцию из экспоненциального представления f(xi^{i_1},...,xi^{i_n})=xi^{e * f(i_1,...,i_n)}
	// в представление основанное на кодировании элементов поля - многочленов "номерами" f(i_1,...,i_n) = j, где двоичное представление i_s соответствует многочлену-элементу поля
	// для этого используется логарифмическая таблица

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
	// k - глубина рекурсии - количество заданных элементов базиса
	// t - мощность базиса
	element a(2, t), b(2, t);


	bool* ar = new bool[q];
	ui bound;


	if (k == t) // если базис построен, то // вывести его в файл
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


	if (k < t)// если базис не полный, то добавим туда элемент
	{
		for (ui i = bound; i < q - t + k; ++i) // выберем элемент для базиса из тех которые не образуются сложением уже имеющихся
		{
			if (old_ar[i] == false)
				basis[k] = i;
			else
				continue;


			a.exp = i;

			//выставим старые флаги
			for (ui j = 0; j < q; ++j)
				ar[j] = old_ar[j];

			//допишем новые образующиеся элементы сложив новый элемент со всеми ненулевыми имеющимися элементами поля
			//выставим новые флаги
			for (ui j = 0; j < q; ++j)
			{
				if (old_ar[j] == true)
				{
					b.exp = j;
					b = a + b;
					ar[b.exp] = true;

				}

			}
			ar[i] = true; // выбранный элемент тоже уже имеется


			//виток рекурсии
			Recursive_find_basis_2t(k + 1, t, q, basis, ar);
		}
	}
}
*/

/*
// Построение всех базисов вп P_{P_0} при p=2
void Build_All_Basis_2t(ui t)
{
	//построим все базисы векторного пространства P_{P_0}
	// с точностью до перестановки векторов
	//зададим поле
	InitGF(2, t);
	ui q = (ui)pow(2, t);

	cout << "Start build basis in GF(" << q << ")" << endl;

	ui* basis = new ui[t]; // для хранения базиса

	for (ui i = 0; i < t; ++i)
		basis[i] = i;

	bool* ar = new bool[q];
	for (ui i = 0; i < q; ++i)
		ar[i] = false;

	Recursive_find_basis_2t(0, t, q, basis, ar);
}
*/



