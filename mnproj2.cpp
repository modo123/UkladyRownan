#include <iostream>
#include <fstream>
#include <ctime>

using namespace std;

int const N=913;
double M[N][N];
double NN[N];
double x1[N];
double x2[N];
double **A;
double b[N];
double D[N][N];
double U[N][N];
double L[N][N];
double x[N];


void wypelnijA_A()
{
	int a1=8;
	int a2=-1;
	int a3=-1;
	// Get values of A
	A=new double* [N];
	for (int i = 0; i < N; i++)
	{
		A[i] = new double [N];
	}

	//wypelnianie zerami
	for (int i=0;i<N;i++)
	{
		for (int j=0;j<N;j++)
		{
			A[i][j]=0;
		}

	}
	//wypelnianie a1
	for(int i=0;i<N;i++)
	{
		for (int j=0;j<N;j++)
		{
			if(i==j)
			{
				A[i][j]=a1;
			}
		}
	}
	//wypelnianie a2
	for (int i=1;i<N;i++)
	{
		for (int j=i-1;j<i;j++)
		{
			A[i][j]=a2;
			A[j][i]=a2;
		}

	}
	//wypelnianie a3
	for (int i=2;i<N;i++)
	{
		for (int j=i-2;j<i;j=j+2)
		{
			A[i][j]=a3;
			A[j][i]=a3;
		}
	}

	for (int i=0;i<N;i++)
	{
		b[i]=sin(0.2*i);
	}

}
void Jacobi()
{
	clock_t s,f;
	double czas;
	s=clock();

	// Calculate N = D^-1
	for (int i=0; i<N; i++)
		NN[i] = pow(A[i][i], -1);

	// Calculate M = -D^-1 (L + U)
	for (int i=0; i<N; i++)
		for (int j=0; j<N; j++)
			if (i == j)
				M[i][j] = 0;
			else
				M[i][j] = - (A[i][j] * NN[i]);

	//Initialize x
	for (int i=0; i<N; i++)
		x1[i] = 0;

	//for (k=0; k<iter; k++) {
	for (int i=0; i<N; i++) 
	{
		x2[i] = NN[i]*b[i];
		for (int j=0; j<N; j++)
			x2[i] += M[i][j]*x1[j];

		
	}
	for (int i=0; i<N; i++)
		x1[i] = x2[i];
	//}

	fstream plik ("Jacobi.txt", ios::out);
	f=clock();
	czas = (double)(f - s) / (double)(CLOCKS_PER_SEC);
	plik<<czas;

	if (plik.good() )
	{
		plik<<endl;
		for (int i=0;i<N;i++)
			plik<<"x"<<i+1<<" = "<<x1[i]<<endl;

		plik.close();
	}

}
void GaussS()
{
	clock_t s,f;
	double czas;
	s=clock();
	// Divide A into L + D + U
	for (int i=0; i<N; i++)
		for (int j=0; j<N; j++) 
		{
			if (i < j) 
			{
				U[i][j] = A[i][j];
			}
			else if (i > j) 
			{
				L[i][j] = A[i][j];
			}
			else 
			{
				D[i][j] = A[i][j];
			}
		}

		// Calculate D^-1
		for (int i=0; i<N; i++)
			D[i][i] = 1/D[i][i];

		// Calculate D^-1 * b
		for (int i=0; i<N; i++)
			b[i] *= D[i][i];

		//Calculate D^-1 * L
		for (int i=0; i<N; i++)
			for (int j=0; j<i; j++)
				L[i][j] *= D[i][i];

		//Calculate D^-1 * U
		for (int i=0; i<N; i++)
			for (int j=i+1; j<N; j++)
				U[i][j] *= D[i][i];

		//Initialize x
		for (int i=0; i<N; i++)
			x[i] = 0;

		//printf("Ile iteracji algorytmu wykonac?\n");
		//scanf("%d", &iter);

		//for (k=0; k<iter; k++)
		for (int i=0; i<N; i++) 
		{
			x[i] = b[i];                       // x = D^-1*b -
			for (int j=0; j<i; j++)
				x[i] -= L[i][j]*x[j];    // D^-1*L * x -
			for (int j=i+1; j<N; j++)
				x[i] -= U[i][j]*x[j];    // D^-1*U * x
		}

		fstream plik2 ("GaussS.txt", ios::out);
		f=clock();
		czas = (double)(f - s) / (double)(CLOCKS_PER_SEC);
		plik2<<czas;

		if (plik2.good() )
		{
			plik2<<endl;
			for (int i=0;i<N;i++)
				plik2<<"x"<<i+1<<" = "<<x1[i]<<endl;

			plik2.close();
		}




}
void wypelnijA_C()
{
	int a1=3;
	int a2=-1;
	int a3=-1;
	// Get values of A
	A=new double* [N];
	for (int i = 0; i < N; i++)
	{
		A[i] = new double [N];
	}

	//wypelnianie zerami
	for (int i=0;i<N;i++)
	{
		for (int j=0;j<N;j++)
		{
			A[i][j]=0;
		}

	}
	//wypelnianie a1
	for(int i=0;i<N;i++)
	{
		for (int j=0;j<N;j++)
		{
			if(i==j)
			{
				A[i][j]=a1;
			}
		}
	}
	//wypelnianie a2
	for (int i=1;i<N;i++)
	{
		for (int j=i-1;j<i;j++)
		{
			A[i][j]=a2;
			A[j][i]=a2;
		}

	}
	//wypelnianie a3
	for (int i=2;i<N;i++)
	{
		for (int j=i-2;j<i;j=j+2)
		{
			A[i][j]=a3;
			A[j][i]=a3;
		}
	}

	for (int i=0;i<N;i++)
	{
		b[i]=sin(0.2*i);
	}

}


int main()
{
	wypelnijA_A();
	Jacobi();
	wypelnijA_A();
	GaussS();
	//wypelnijA_C();









	return 0;
}
