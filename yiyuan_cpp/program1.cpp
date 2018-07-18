#include<iomanip>
#include<iostream>
#include<math.h>
#include<stdlib.h>
#include<stdio.h>
#include<fstream>
#include<time.h>
using namespace std;

extern "C" void dsyev_(char *a,char *b,int *n,double *A,int *nn,double *e,double *work,int *lwork,int *info);
void tr(double **h,double **vec,double *e,int dim)
{
	char a='V';
	char b='L';
	int N1;
	int NN;
	int lwork;
	int info;
	lwork=dim*3;
	double *work;
	N1=dim;
	NN=N1;
	info=0;
	work=new double [dim*3];   
	double *h_temp=new double [dim*dim];
	int i,j;
	for(i=0;i<dim;i++)
	{
		for(j=0;j<dim;j++)
		{
			h_temp[i*dim+j]=h[i][j];
		}
	}
	dsyev_(&a,&b,&N1,h_temp,&NN,e,work,&lwork,&info);
	for(i=0;i<dim;i++)
	{
		for(j=0;j<dim;j++)
		{
			vec[i][j]=h_temp[i*dim+j];  
		}
	}
	delete h_temp;
	delete work;
}
// h is a symmetric matrix, vec[i] is the i-th eigenvector, e[i] is the i-th eigenvalue, dim is the dimension of h,vec,e 

int nj;
int npart;
int pair_ex_max;
int dimension;

double fac(int x)
{
	int k;
	double res=1;

	if (x<0) 
	{
		return 0;
	}
	else if(x==0)
	{
		return 1;
	}
	else
	{
		for(k=1;k<=x;k++) 
		{
			res=res*k;
		}

		return res;
	}
}

void ph_config(int *hh, int h_num, int *pp, int p_num, int mode, int ph_num, int **config, int *config_num)
{
	int i,k;
	int kf=npart/2;

	if(mode==0)
	{
		if(h_num>=ph_num)
		{
			ph_config(hh, ph_num, pp, 0, 1, ph_num, config, config_num);
		}
		else
		{
			if(h_num==0)
			{
				for(hh[0]=0;hh[0]<(kf-ph_num+1);hh[0]++)
				{
					ph_config(hh, 1, pp, 0, 0, ph_num, config, config_num);			
				}
			}
			else
			{
				for(hh[h_num]=hh[h_num-1]+1;hh[h_num]<(kf-ph_num+h_num+1);hh[h_num]++)
				{
					ph_config(hh, h_num+1, pp, 0, 0, ph_num, config, config_num);
				}
			}
		}
	}
	else
	{
		if(p_num>=ph_num)
		{
			for(i=0;i<ph_num;i++)
			{
				config[config_num[0]][hh[i]]=0;
			}
			for(i=0;i<ph_num;i++)
			{
				config[config_num[0]][pp[i]]=1;
			}

			cout<<"   pair ex. number = "<<ph_num<<" : ";
			for(k=0;k<nj;k++)
			{
				cout<<config[config_num[0]][k]<<", ";
			}
			cout<<endl;

			config_num[0]=config_num[0]+1;

		}
		else
		{
			if(p_num==0)
			{
				for(pp[0]=kf;pp[0]<(nj-ph_num+1);pp[0]++)
				{
					ph_config(hh, ph_num, pp, 1, 1, ph_num, config, config_num);
				}
			}	
			else
			{
				for(pp[p_num]=pp[p_num-1]+1;pp[p_num]<(nj-ph_num+p_num+1);pp[p_num]++)
				{
					ph_config(hh, ph_num, pp, p_num+1, 1, ph_num, config, config_num);
				}
			}
		}

	}
}

void space(int ***pair_config, int *pair_config_num)
{
	int i,j,k;
	int h,p;

	int kf=npart/2;

	int *hh,*pp;
	hh=new int[npart/2];
	pp=new int[nj-npart/2];

	int *num;
	num= new int [1];
	for(i=1;i<=pair_ex_max;i++)
	{
		num[0]=0;
		ph_config(hh, 0, pp, 0, 0, i, pair_config[i], num);

		if(num[0] != pair_config_num[i])
		{
			cout<<"Error for pair_config_num of pair_ex_num = "<<i<<endl;
		}
	}

	delete hh;
	delete pp;
}

void pair_Hmatrix(int ***pair_config, int *pair_config_num, double *esp, double g, double **Hmatrix)
{
	int i1,i2,j1,j2,k;
	int num1,num2,d_num;

	num1=0;
	for(i1=0;i1<=pair_ex_max;i1++)
	{
		for(j1=0;j1<pair_config_num[i1];j1++)
		{
			num2=0;
			for(i2=0;i2<=pair_ex_max;i2++)
			{
				for(j2=0;j2<pair_config_num[i2];j2++)
				{
					if(num1==num2)
					{
						for(k=0;k<nj;k++)
						{
							Hmatrix[num1][num1]=Hmatrix[num1][num1]+2*esp[k]*pair_config[i1][j1][k];
							if(pair_config[i1][j1][k]==1)
							{
								Hmatrix[num1][num1]=Hmatrix[num1][num1]-g;
							}
						}

					}			
					else
					{
						d_num=0;
						for(k=0;k<nj;k++)
						{
							if(pair_config[i1][j1][k]==1 && pair_config[i2][j2][k]==0)
							{
								d_num=d_num+1;
							}
						}

						if(d_num == 1)
						{
							Hmatrix[num1][num2]=Hmatrix[num1][num2]-g;
						}
					}	
					num2++;
				}
			}		
			num1++;
		}
	}
}

int main()
{
	int i,j,k;

	cout<<"input the number of single-particle states: ";
	cin>>nj;

	cout<<"input the number of particles: ";
	cin>>npart;

	if(npart%2 != 0)
	{
		cout<<"Error! only treat even-number system"<<endl;
		return 0;
	}

	int kf=npart/2;
	if(kf>nj)
	{
		cout<<"Error! need more orbits"<<endl;
		return 0;
	}

	cout<<"input the max number of pairs excited across the Fermi surface: ";
	cin>>pair_ex_max;

	if(pair_ex_max>npart/2)
	{
		cout<<"Error! pair_ex_max is larger than npart/2"<<endl;
		return 0;
	}

	int *pair_config_num;
	pair_config_num= new int [pair_ex_max+1];

	pair_config_num[0]=1;	
	dimension=pair_config_num[0];
	for(i=1;i<=pair_ex_max;i++)
	{
		if(kf<i || (nj-kf)<i)
		{
			pair_config_num[i]=0;
		}
		else
		{
			pair_config_num[i]=fac(kf)*fac(nj-kf)/(fac(i)*fac(kf-i)*fac(i)*fac(nj-kf-i));
			dimension=dimension+pair_config_num[i];

			cout<<"   There are "<<pair_config_num[i]<<" basis states for "<<i<<" pair(s) excited across Fermi surface"<<endl;
		}

	}

	int ***pair_config;
	pair_config= new int **[pair_ex_max+1];
	for(i=0;i<=pair_ex_max;i++)
	{
		pair_config[i]= new int *[1000];
		for(j=0;j<1000;j++)
		{
			pair_config[i][j]= new int [nj];
		}
	}

	for(i=0;i<=pair_ex_max;i++)
	{
		for(j=0;j<1000;j++)
		{
			for(k=0;k<nj;k++)
			{
				if(k<kf)
				{
					pair_config[i][j][k]=1;
				}
				else
				{
					pair_config[i][j][k]=0;
				}
			}
		}
	}

	space(pair_config, pair_config_num);

	double **Hmatrix;
	Hmatrix= new double *[dimension];
	for(i=0;i<dimension;i++)
	{
		Hmatrix[i]= new double [dimension];
	}

	double **vec;
	vec= new double *[dimension];
	for(i=0;i<dimension;i++)
	{
		vec[i]= new double [dimension];
	}

	double *ener;
	ener= new double [dimension];

	double *esp;
	esp= new double [nj];

	double gap=1.0;	
	double g;

	for(i=0;i<nj;i++)
	{
		esp[i]=i*gap;
	}

	int dim=dimension;

	char filename1[90]="output1.dat";
	FILE *fp;

	if((fp=fopen(filename1,"w"))==NULL)
	{
		cout<<"open file error"<<endl;
	}
	else
	{
		for(g=-1.0;g<=1.0;g=g+0.05)
		{
			for(i=0;i<dimension;i++)
			{
				for(j=0;j<dimension;j++)
				{
					Hmatrix[i][j]=0.0;
				}
			}
			pair_Hmatrix(pair_config, pair_config_num, esp, g, Hmatrix);
			/*for(i=0;i<dimension;i++)
			{
				for(j=0;j<dimension;j++)
				{
					cout<<Hmatrix[i][j]<<", ";
				}
				cout<<endl;
			}*/

			tr(Hmatrix, vec, ener, dim);
			fprintf(fp, " %6.3lf   %10.5lf   %10.5lf   %10.5lf\n", g, ener[0], Hmatrix[0][0], ener[0]-Hmatrix[0][0]);
		}
		fprintf(fp, " g-val     tot-ener     ref-ener    corre-ener\n");
		fprintf(fp, " orbit num: %2d\n", nj);
		fprintf(fp, " particle num: %2d\n", npart);
		fprintf(fp, " max number of pair ex.:%2d", pair_ex_max);
		fclose(fp);
	}

	delete pair_config_num;

	for(i=0;i<=pair_ex_max;i++)
	{
		for(j=0;j<1000;j++)
		{
			delete pair_config[i][j];
		}
		delete pair_config[i];
	}
	delete pair_config;

	for(i=0;i<dimension;i++)
	{
		delete Hmatrix[i];
		delete vec[i];
	}
	delete Hmatrix;
	delete vec;

	delete ener;
	delete esp;
}

