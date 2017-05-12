// 航带法解析空中三角测量.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <malloc.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
using namespace std;

#define f 0.153033
#define bx_ 0.200

typedef struct sData
{
	int iCount;
	int *iPointNum; //点号
	double *x1;      
	double *y1;
	double *x2;      
	double *y2;
}sData;

typedef struct ModelPoint   //模型点坐标
{
	double *Xm;
	double *Ym;
	double *Zm;
}ModelPoint;

typedef struct pMPoint   //摄影测量坐标系坐标
{
	double *Xp;
	double *Yp;
	double *Zp;
}pMPoint;

typedef struct tpMPoint   //地面摄影测量坐标系坐标
{
	double *Xtp;
	double *Ytp;
	double *Ztp;
}tpMPoint;

typedef struct tMPoint   //地面坐标系坐标
{
	double *Xt;
	double *Yt;
	double *Zt;
}tMPoint;

typedef struct CPoint   //控制点、检查点坐标
{
	int *iPointNum;
	double *X;
	double *Y;
	double *Z;
	
}CPoint;

double square(double x)              
{
	return x*x;
}

void Scale(sData *Data,CPoint *cPoint,double *m)              //求比例尺系数
{
	int iCount = 0;
	int iComPoint[3];
	double sumM = 0,M,L,photoL;

	for(int i = 0;i < Data->iCount;i++)   //寻找公共点
	{
		for(int j = 0;j < 3;j++)
		{
			if(Data->iPointNum[i] == cPoint->iPointNum[j])
			{
				iComPoint[iCount++] = i;
			}
		}
	}

	int n = 0;
	for(int i = 0;i < 3;i++)
	{
		for(int j = i + 1;j < 3;j++)
		{
			photoL = sqrt(square(Data->x1[ iComPoint[j] ] - Data->x1[ iComPoint[i] ]) + square(Data->y1[ iComPoint[j] ] - Data->y1[ iComPoint[i] ]));
			L = sqrt(square(cPoint->X[j] - cPoint->X[i]) + square(cPoint->Y[j] - cPoint->Y[i]));
			M = L / photoL;
			sumM += M;
			n++;
		}
	}
	*m = sumM/n;
}


//读取数据
void ReadData(char* name,sData *Data)
{
	int temNum;
	double tempData;
	FILE * fp;
	if((fp  = fopen(name,"r")) == NULL)
	{
		printf("ERROR!!!!!!!!!!");
		exit(0);
	}
	else
	{
		do
		{
			fscanf(fp,"%d",&temNum);
			fscanf(fp,"%lf",&tempData);
			fscanf(fp,"%lf",&tempData);
			fscanf(fp,"%lf",&tempData);
			fscanf(fp,"%lf",&tempData);

			Data->iCount++;
		}while(temNum > 0);

		rewind(fp);

		Data->iPointNum = (int*)malloc(sizeof(int) * Data->iCount);
		Data->x1 = (double*)malloc(sizeof(double) * Data->iCount);
		Data->y1 = (double*)malloc(sizeof(double) * Data->iCount);
		Data->x2 = (double*)malloc(sizeof(double) * Data->iCount);
		Data->y2 = (double*)malloc(sizeof(double) * Data->iCount);

		/*读取数据*/
		for(int i = 0;i < Data->iCount;i++)
		{
			fscanf(fp,"%d",&temNum);
			if(temNum < 0)
				Data->iPointNum[i] = -temNum;
			else
				Data->iPointNum[i] = temNum;

			fscanf(fp,"%lf",&tempData);
			Data->x1[i] = tempData;
			fscanf(fp,"%lf",&tempData);
			Data->y1[i] = tempData;
			fscanf(fp,"%lf",&tempData);
			Data->x2[i] = tempData;
			fscanf(fp,"%lf",&tempData);
			Data->y2[i] = tempData;
		}
	}
	fclose(fp);
}

//读取控制点、检查点数据
void ReadCData(char* name,CPoint *cPoint,int x)
{
	int temNum;
	double tempData;
	FILE * fp;
	if((fp  = fopen(name,"r")) == NULL)
	{
		printf("ERROR!!!!!!!!!!");
		exit(0);
	}
	else
	{
		cPoint->iPointNum = (int*)malloc(sizeof(int) * x);
		cPoint->X = (double*)malloc(sizeof(double) * x);
		cPoint->Y = (double*)malloc(sizeof(double) * x);
		cPoint->Z = (double*)malloc(sizeof(double) * x);


		/*读取数据*/
		for(int i = 0;i < x;i++)
		{
			fscanf(fp,"%d",&temNum);
			cPoint->iPointNum[i] = temNum;

			fscanf(fp,"%lf",&tempData);
			cPoint->X[i] = tempData;
			fscanf(fp,"%lf",&tempData);
			cPoint->Y[i] = tempData;
			fscanf(fp,"%lf",&tempData);
			cPoint->Z[i] = tempData;
		}
	}
	fclose(fp);
}


//相对定向
void RelativeOrientation(sData *Data,ModelPoint *mPoint,double *by,double *bz) 
{
	//单位换算，以米为单位 
	for(int i = 0;i < Data->iCount;i++)     
	{
		Data->x1[i] = Data->x1[i] * 0.000001;
		Data->y1[i] = Data->y1[i] * 0.000001;
		Data->x2[i] = Data->x2[i] * 0.000001;
		Data->y2[i] = Data->y2[i] * 0.000001;
	}

	/*确定初始值*/
	int icount = 0;
	double mu,nu,alpha,omega,kappa;
	double *X1,*Y1,*Z1,*X2,*Y2,*Z2;
	double *N1,*N2;

	mu = 0;nu = 0;alpha = 0;omega = 0;kappa = 0;
	X1 = (double*)malloc(sizeof(double) * Data->iCount);
	Y1 = (double*)malloc(sizeof(double) * Data->iCount);
	Z1 = (double*)malloc(sizeof(double) * Data->iCount);
	X2 = (double*)malloc(sizeof(double) * Data->iCount);
	Y2 = (double*)malloc(sizeof(double) * Data->iCount);
	Z2 = (double*)malloc(sizeof(double) * Data->iCount);
	N1 = (double*)malloc(sizeof(double) * Data->iCount);
	N2 = (double*)malloc(sizeof(double) * Data->iCount);


	while(icount < 50)
	{

		/*计算旋转矩阵*/
		double R1[3][3]={{1,0,0},{0,1,0},{0,0,1}},R2[3][3];
		R2[0][0] = cos(alpha)*cos(kappa) - sin(alpha)*sin(omega)*sin(kappa);
		R2[0][1] = -cos(alpha)*sin(kappa) - sin(alpha)*sin(omega)*cos(kappa);
		R2[0][2] = -sin(alpha)*cos(omega);
		R2[1][0] = cos(omega)*sin(kappa);
		R2[1][1] = cos(omega)*cos(kappa);
		R2[1][2] = -sin(omega);
		R2[2][0] = sin(alpha)*cos(kappa) + cos(alpha)*sin(omega)*sin(kappa);
		R2[2][1] = -sin(alpha)*sin(kappa) + cos(alpha)*sin(omega)*cos(kappa);
		R2[2][2] = cos(alpha)*cos(omega); 


		/*组成误差方程*/
		double *Q;
		double **A; 
		double **At;
		double AtA[5][5];
		double AtQ[5];
	
		Q = (double*)malloc(sizeof(double) * Data->iCount);
		A = new double*[Data->iCount];
		for(int i = 0;i < Data->iCount;i++)
			A[i] = new double[5];
		At = new double*[5];
		for(int i = 0;i < 5;i++)
			At[i] = new double[Data->iCount];
	
		for(int i = 0;i < Data->iCount;i++)
		{
			X1[i] = Data->x1[i];
			Y1[i] = Data->y1[i];
			Z1[i] = -f;
			X2[i] = R2[0][0] * Data->x2[i] + R2[0][1] * Data->y2[i] + R2[0][2] * (-f);
			Y2[i] = R2[1][0] * Data->x2[i] + R2[1][1] * Data->y2[i] + R2[1][2] * (-f);
			Z2[i] = R2[2][0] * Data->x2[i] + R2[2][1] * Data->y2[i] + R2[2][2] * (-f);
	
			*by = mu * bx_;
			*bz = nu * bx_;
	
			N1[i] = (bx_*Z2[i] - *bz*X2[i]) / (X1[i]*Z2[i] - X2[i]*Z1[i]);
			N2[i] = (bx_*Z1[i] - *bz*X1[i]) / (X1[i]*Z2[i] - X2[i]*Z1[i]);
	
			Q[i] = N1[i] * Y1[i] - N2[i] * Y2[i] - *by;
			A[i][0] = bx_;
			A[i][1] = -Y2[i]/Z2[i] * bx_;
			A[i][2] = -X2[i]*Y2[i]/Z2[i] * N2[i];
			A[i][3] = -(Z2[i] + Y2[i]*Y2[i]/Z2[i])*N2[i];
			A[i][4] = X2[i]*N2[i];

		}
	
		/*计算法方程*/
		for(int i = 0;i < 5;i++)              //转置
		{
			for(int j = 0;j < Data->iCount;j++)
			{
				At[i][j] = A[j][i];
			}
		}
	
		for(int i = 0;i < 5;i++)                //矩阵相乘  
		{
			for(int j = 0;j < 5;j++) 
			{
				AtA[i][j] = 0;
				AtQ[i] = 0;
				for(int k = 0;k < Data->iCount;k++)
				{
					AtA[i][j] += At[i][k] * A[k][j];
					AtQ[i] += At[i][k] * Q[k];
				}
			}
		}

		/*计算相对定向元素改正数*/
		double N_AtA[5][5];
		double Tem[5][5];

		int k = 0;                               //矩阵求逆
		for(int i = 0;i < 5;i++)
		{
			for(int j = 0;j < 5;j++)
			{
				if(i == k &&j == k)
					Tem[i][j] = 1 / AtA[k][k];
				if(i == k && j!= k)
					Tem[i][j] = AtA[k][j] / AtA[k][k];
				if(i != k && j == k)
					Tem[i][j] = -(AtA[i][k] / AtA[k][k]);
				if(i != k && j != k)
					Tem[i][j] = AtA[i][j] - AtA[i][k] * AtA[k][j] / AtA[k][k];
			}
		}

		for(int k = 1;k < 5;k++)
		{
			for(int i = 0;i < 5;i++)
			{
				for(int j = 0;j < 5;j++)
				{
					if(i == k &&j == k)
						N_AtA[i][j] = 1 /Tem[k][k];
					if(i == k && j!= k)
						N_AtA[i][j] = Tem[k][j] / Tem[k][k];
					if(i != k && j == k)
						N_AtA[i][j] = -(Tem[i][k] / Tem[k][k]);
					if(i != k && j != k)
						N_AtA[i][j] = Tem[i][j] - Tem[i][k] * Tem[k][j] / Tem[k][k];
				}
			}
			for(int i = 0;i < 5;i++)
				for(int j = 0;j < 5;j++)
					Tem[i][j] = N_AtA[i][j];
		}


		double X[5];            //矩阵相乘
		for(int i = 0;i < 5;i++)
		{
			X[i] = 0;
			for(int k = 0;k < 5;k++)
				X[i] += N_AtA[i][k] * AtQ[k];
		}
	
		/*更新相对定向元素*/
		mu += X[0];
		nu += X[1];
		alpha += X[2];
		omega += X[3];
		kappa += X[4];

		icount++;

		//判断限差
		int I;
		for(I = 0;I < 5; I++)
			if(fabs(X[I]) > 3e-5)
				break;
		if(I >= 5)
			break;	
	

		free(Q);
		for(int i = 0;i < Data->iCount;i++)
			delete(A[i]);
		delete(A);
		for(int i = 0;i < 5;i++)
			delete(At[i]);
		delete(At);
	}

	/*模型点坐标*/
	
	mPoint->Xm = (double*)malloc(sizeof(double) * Data->iCount);
	mPoint->Ym = (double*)malloc(sizeof(double) * Data->iCount);
	mPoint->Zm = (double*)malloc(sizeof(double) * Data->iCount);


	for(int i = 0;i < Data->iCount;i++)
	{
		mPoint->Xm[i] = N1[i] * X1[i];
		mPoint->Ym[i] = (N1[i]*Y1[i] + N2[i]*Y2[i] + *by) / 2;
		mPoint->Zm[i] = N1[i] * Z1[i];
	}

	free(X1);free(Y1);free(Z1);
	free(X2);free(Y2);free(Z2);
	free(N1);free(N2);

}



int _tmain(int argc, _TCHAR* argv[])
{
	char name1[] = "像对一.txt";
	char name2[] = "像对二.txt";
	char name3[] = "像对三.txt";
	char name4[] = "控制点数据.txt";
	char name5[] = "检查点数据.txt";
	double temby,tembz;
	double bz[3];
	double by_[3];
	double m;
	double conXtpPoint[4];
	double conYtpPoint[4];
	double conZtpPoint[4];
	sData Data1;
	sData Data2;
	sData Data3;
	ModelPoint mPoint1;  //模型点像空间辅助坐标系
	ModelPoint mPoint2;
	ModelPoint mPoint3;
	pMPoint pPoint1;     //摄影测量坐标系
	pMPoint pPoint2;
	pMPoint pPoint3;
	tpMPoint tpPoint1;   //地面摄影测量坐标系
	tpMPoint tpPoint2;
	tpMPoint tpPoint3;
	tMPoint tPoint1;     //大地坐标系
	tMPoint tPoint2;
	tMPoint tPoint3;
	CPoint conPoint;
	CPoint chPoint;

	Data1.iCount = 0;
	Data2.iCount = 0;
	Data3.iCount = 0;


	ReadData(name1,&Data1);

	RelativeOrientation(&Data1,&mPoint1,&temby,&tembz);  
	by_[0] = temby;
	bz[0] = tembz;


	ReadData(name2,&Data2);

	RelativeOrientation(&Data2,&mPoint2,&temby,&tembz);
	by_[1] = temby;
	bz[1] = tembz;


	ReadData(name3,&Data3);

	RelativeOrientation(&Data3,&mPoint3,&temby,&tembz);
	by_[2] = temby;
	bz[2] = tembz;
	

	ReadCData(name4,&conPoint,4);

	ReadCData(name5,&chPoint,5);


	
	Scale(&Data1,&conPoint,&m);


	/*分配内存空间*/
	pPoint1.Xp = (double*)malloc(sizeof(double) * Data1.iCount);
	pPoint1.Yp = (double*)malloc(sizeof(double) * Data1.iCount);
	pPoint1.Zp = (double*)malloc(sizeof(double) * Data1.iCount);

	pPoint2.Xp = (double*)malloc(sizeof(double) * Data2.iCount);
	pPoint2.Yp = (double*)malloc(sizeof(double) * Data2.iCount);
	pPoint2.Zp = (double*)malloc(sizeof(double) * Data2.iCount);

	pPoint3.Xp = (double*)malloc(sizeof(double) * Data3.iCount);
	pPoint3.Yp = (double*)malloc(sizeof(double) * Data3.iCount);
	pPoint3.Zp = (double*)malloc(sizeof(double) * Data3.iCount);


	
	tpPoint1.Xtp = (double*)malloc(sizeof(double)*Data1.iCount);
	tpPoint1.Ytp = (double*)malloc(sizeof(double)*Data1.iCount);
	tpPoint1.Ztp = (double*)malloc(sizeof(double)*Data1.iCount);

	tpPoint2.Xtp = (double*)malloc(sizeof(double)*Data2.iCount);
	tpPoint2.Ytp = (double*)malloc(sizeof(double)*Data2.iCount);
	tpPoint2.Ztp = (double*)malloc(sizeof(double)*Data2.iCount);

	tpPoint3.Xtp = (double*)malloc(sizeof(double)*Data3.iCount);
	tpPoint3.Ytp = (double*)malloc(sizeof(double)*Data3.iCount);
	tpPoint3.Ztp = (double*)malloc(sizeof(double)*Data3.iCount);


	
	tPoint1.Xt = (double*)malloc(sizeof(double)*Data1.iCount);
	tPoint1.Yt = (double*)malloc(sizeof(double)*Data1.iCount);
	tPoint1.Zt = (double*)malloc(sizeof(double)*Data1.iCount);

	tPoint2.Xt = (double*)malloc(sizeof(double)*Data2.iCount);
	tPoint2.Yt = (double*)malloc(sizeof(double)*Data2.iCount);
	tPoint2.Zt = (double*)malloc(sizeof(double)*Data2.iCount);

	tPoint3.Xt = (double*)malloc(sizeof(double)*Data3.iCount);
	tPoint3.Yt = (double*)malloc(sizeof(double)*Data3.iCount);
	tPoint3.Zt = (double*)malloc(sizeof(double)*Data3.iCount);


	/**************************模型连接**********************************/
	int iCount = 0;
	double k1,k2,k3;
	double k;
	int iCommonPoint1[3] = {0};
	int iCommonPoint2[3] = {0};

	for(int i = 0;i < Data1.iCount;i++)   //寻找公共点
	{
		for(int j = 0;j < Data2.iCount;j++)
		{
			if(Data1.iPointNum[i] == Data2.iPointNum[j])
			{
				iCommonPoint1[iCount] = i;
				iCommonPoint2[iCount] = j;
				iCount++;
			}
		}
	}

	
	for(int i = 0;i < Data1.iCount;i++)  //模型1，转为摄影测量坐标
	{
		pPoint1.Xp[i] = m * mPoint1.Xm[i];
		pPoint1.Yp[i] = m * mPoint1.Ym[i];
		pPoint1.Zp[i] = m * mPoint1.Zm[i] + m*f;
	}
	

	
	k1 = (mPoint1.Zm[ iCommonPoint1[0] ] - bz[0]) / mPoint2.Zm[ iCommonPoint2[0] ];
	k2 = (mPoint1.Zm[ iCommonPoint1[1] ] - bz[0]) / mPoint2.Zm[ iCommonPoint2[1] ];
	k3 = (mPoint1.Zm[ iCommonPoint1[2] ] - bz[0]) / mPoint2.Zm[ iCommonPoint2[2] ];
	k = (k1 + k2 + k3)/3;

	for(int i = 0;i < Data2.iCount;i++)  //模型2，转为摄影测量坐标
	{
		pPoint2.Xp[i] = bx_*m    + m*k*mPoint2.Xm[i];
		pPoint2.Yp[i] = by_[0]*m + m*k*mPoint2.Ym[i];
		pPoint2.Zp[i] = bz[0]*m  + m*k*mPoint2.Zm[i] + m*f;
	}

	for(int i = 0;i < Data3.iCount;i++)   
	{
		pPoint3.Xp[i] = m * k * mPoint3.Xm[i];
		pPoint3.Yp[i] = m * k * mPoint3.Ym[i];
		pPoint3.Zp[i] = m * k * mPoint3.Zm[i];
	}


	iCount = 0;
	for(int i = 0;i < Data2.iCount;i++)   //寻找公共点
	{
		for(int j = 0;j < Data3.iCount;j++)
		{
			if(Data2.iPointNum[i] == Data3.iPointNum[j])
			{
				iCommonPoint1[iCount] = i;
				iCommonPoint2[iCount] = j;
				iCount++;
			}
		}
	}
	
	k1 = (mPoint2.Zm[ iCommonPoint1[0] ] - bz[1]) / mPoint3.Zm[ iCommonPoint2[0] ];
	k2 = (mPoint2.Zm[ iCommonPoint1[1] ] - bz[1]) / mPoint3.Zm[ iCommonPoint2[1] ];
	k3 = (mPoint2.Zm[ iCommonPoint1[2] ] - bz[1]) / mPoint3.Zm[ iCommonPoint2[2] ];
	k  = (k1 + k2 + k3)/3;


	for(int i = 0;i < Data3.iCount;i++)    //模型3，转为摄影测量坐标
	{
		pPoint3.Xp[i] = bx_*m    + bx_*m    + k*pPoint3.Xp[i];
		pPoint3.Yp[i] = by_[1]*m + by_[0]*m + k*pPoint3.Yp[i];
		pPoint3.Zp[i] = bz[1]*m  + bz[0]*m  + k*pPoint3.Zp[i] + m*f;
	}


	/*********************模型的绝对定向**************************/
	int K1 = 0;
	int place[4];
	double a,b,lamda;
	double DeltaXp,DeltaYp;
	double DeltaXt,DeltaYt;
	
	for(int i = 0;i < Data1.iCount;i++)
	{
		for(int j = 0;j < 4;j++)
		{
			if(Data1.iPointNum[i] == conPoint.iPointNum[j])
			{
				place[K1++] = i;
			}
		}
	}

	for(int i = 0;i < Data3.iCount;i++)   //寻找公共点位置
	{
		if(Data3.iPointNum[i] == conPoint.iPointNum[3])
		{
			place[3] = i;
		}
	}

	DeltaXp = pPoint3.Xp[ place[3] ] - pPoint1.Xp[ place[0] ];
	DeltaYp = pPoint3.Yp[ place[3] ] - pPoint1.Yp[ place[0] ];

	DeltaXt = conPoint.X[3] - conPoint.X[0];
	DeltaYt = conPoint.Y[3] - conPoint.Y[0];
	a = (DeltaXt*DeltaYp + DeltaYt*DeltaXp) / (DeltaXt*DeltaXt + DeltaYt*DeltaYt);
	b = (DeltaXt*DeltaXp - DeltaYt*DeltaYp) / (DeltaXt*DeltaXt + DeltaYt*DeltaYt);
	lamda = sqrt(a*a + b*b);

	for(int i = 0;i < 4;i++)
	{
		conXtpPoint[i] = b*(conPoint.X[i] - conPoint.X[0]) + a*(conPoint.Y[i] - conPoint.Y[0]);
		conYtpPoint[i] = a*(conPoint.X[i] - conPoint.X[0]) - b*(conPoint.Y[i] - conPoint.Y[0]);
		conZtpPoint[i] = lamda * (conPoint.Z[i] - conPoint.Z[0]);
	}

	//解求绝对定向元素
	double abL[12];
	double abA[12][7];
	double abAt[7][12];
	double abAtA[7][7];
	double abAtL[7];

	//初始化
	double R2[3][3];
	double DeltaX = 0,DeltaY = 0,DeltaZ = 0,Lamda = 0,Alpha = 0,Omega = 0,Kappa = 0;

	int count = 0;
	while(count < 20)
	{

		/*计算旋转矩阵*/
		R2[0][0] = cos(Alpha)*cos(Kappa) - sin(Alpha)*sin(Omega)*sin(Kappa);
		R2[0][1] = -cos(Alpha)*sin(Kappa) - sin(Alpha)*sin(Omega)*cos(Kappa);
		R2[0][2] = -sin(Alpha)*cos(Omega);
		R2[1][0] = cos(Omega)*sin(Kappa);
		R2[1][1] = cos(Omega)*cos(Kappa);
		R2[1][2] = -sin(Omega);
		R2[2][0] = sin(Alpha)*cos(Kappa) + cos(Alpha)*sin(Omega)*sin(Kappa);
		R2[2][1] = -sin(Alpha)*sin(Kappa) + cos(Alpha)*sin(Omega)*cos(Kappa);
		R2[2][2] = cos(Alpha)*cos(Omega); 

		for(int i = 0;i < 3;i++)
		{
			abA[3*i+0][0] = 1;
			abA[3*i+0][1] = 0;
			abA[3*i+0][2] = 0;
			abA[3*i+0][3] = pPoint1.Xp[ place[i] ];
			abA[3*i+0][4] = -pPoint1.Zp[ place[i] ];
			abA[3*i+0][5] = 0;
			abA[3*i+0][6] = -pPoint1.Yp[ place[i] ];
			abL[3*i+0] = conXtpPoint[i] - Lamda * ( R2[0][0]*pPoint1.Xp[ place[i] ] + R2[0][1]*pPoint1.Yp[ place[i] ] + R2[0][2]*pPoint1.Zp[ place[i] ] ) - DeltaX;

			abA[3*i+1][0] = 0;
			abA[3*i+1][1] = 1;
			abA[3*i+1][2] = 0;
			abA[3*i+1][3] = pPoint1.Yp[ place[i] ];
			abA[3*i+1][4] = 0;
			abA[3*i+1][5] = -pPoint1.Zp[ place[i] ];
			abA[3*i+1][6] = pPoint1.Xp[ place[i] ];
			abL[3*i+1] = conYtpPoint[i] - Lamda *( R2[1][0]*pPoint1.Xp[ place[i] ] + R2[1][1]*pPoint1.Yp[ place[i] ] + R2[1][2]*pPoint1.Zp[ place[i] ]) - DeltaY;

			abA[3*i+2][0] = 0;
			abA[3*i+2][1] = 0;
			abA[3*i+2][2] = 1;
			abA[3*i+2][3] = pPoint1.Zp[ place[i] ];
			abA[3*i+2][4] = pPoint1.Xp[ place[i] ];
			abA[3*i+2][5] = pPoint1.Yp[ place[i] ];
			abA[3*i+2][6] = 0;
			abL[3*i+2] = conZtpPoint[i] - Lamda * ( R2[2][0]*pPoint1.Xp[ place[i] ] + R2[2][1]*pPoint1.Yp[ place[i] ] + R2[2][2]*pPoint1.Zp[ place[i] ]) - DeltaZ;
		}

		abA[9][0] = 1;
		abA[9][1] = 0;
		abA[9][2] = 0;
		abA[9][3] = pPoint3.Xp[ place[3] ];
		abA[9][4] = -pPoint3.Zp[ place[3] ];
		abA[9][5] = 0;
		abA[9][6] = -pPoint3.Yp[ place[3] ];
		abL[9] = conXtpPoint[3] - Lamda * ( R2[0][0]*pPoint3.Xp[ place[3] ] + R2[0][1]*pPoint3.Yp[ place[3] ] + R2[0][2]*pPoint3.Zp[ place[3] ] ) - DeltaX;

		abA[10][0] = 0;
		abA[10][1] = 1;
		abA[10][2] = 0;
		abA[10][3] = pPoint3.Yp[ place[3] ];
		abA[10][4] = 0;
		abA[10][5] = -pPoint3.Zp[ place[3] ];
		abA[10][6] = pPoint3.Xp[ place[3] ];
		abL[10] = conYtpPoint[3] - Lamda *( R2[1][0]*pPoint3.Xp[ place[3] ] + R2[1][1]*pPoint3.Yp[ place[3] ] + R2[1][2]*pPoint3.Zp[ place[3] ]) - DeltaY;

		abA[11][0] = 0;
		abA[11][1] = 0;
		abA[11][2] = 1;
		abA[11][3] = pPoint3.Zp[ place[3] ];
		abA[11][4] = pPoint3.Xp[ place[3] ];
		abA[11][5] = pPoint3.Yp[ place[3] ];
		abA[11][6] = 0;
		abL[11] = conZtpPoint[3] - Lamda * ( R2[2][0]*pPoint3.Xp[ place[3] ] + R2[2][1]*pPoint3.Yp[ place[3] ] + R2[2][2]*pPoint3.Zp[ place[3] ]) - DeltaZ;

		/*计算法方程*/
		for(int i = 0;i < 7;i++)              //转置
		{
			for(int j = 0;j < 12;j++)
			{
				abAt[i][j] = abA[j][i];
			}
		}
		
		for(int i = 0;i < 7;i++)                //矩阵相乘  
		{
			for(int j = 0;j < 7;j++) 
			{
				abAtA[i][j] = 0;
				abAtL[i] = 0;
				for(int k = 0;k < 12;k++)
				{
					abAtA[i][j] += abAt[i][k] * abA[k][j];
					abAtL[i] += abAt[i][k] * abL[k];
				}
			}
		}

		/*计算相对定向元素改正数*/
		double N_abAtA[7][7];
		double abTem[7][7];

		int K = 0;                               //矩阵求逆
		for(int i = 0;i < 7;i++)
		{
			for(int j = 0;j < 7;j++)
			{
				if(i == K &&j == K)
					abTem[i][j] = 1 / abAtA[K][K];
				if(i == K && j!= K)
					abTem[i][j] = abAtA[K][j] / abAtA[K][K];
				if(i != K && j == K)
					abTem[i][j] = -(abAtA[i][K] / abAtA[K][K]);
				if(i != K && j != K)
					abTem[i][j] = abAtA[i][j] - abAtA[i][K] * abAtA[K][j] / abAtA[K][K];
			}
		}

		for( K = 1;K < 7;K++)
		{
			for(int i = 0;i < 7;i++)
			{
				for(int j = 0;j < 7;j++)
				{
					if(i == K &&j == K)
						N_abAtA[i][j] = 1 /abTem[K][K];
					if(i == K && j!= K)
						N_abAtA[i][j] = abTem[K][j] / abTem[K][K];
					if(i != K && j == K)
						N_abAtA[i][j] = -(abTem[i][K] / abTem[K][K]);
					if(i != K && j != K)
						N_abAtA[i][j] = abTem[i][j] - abTem[i][K] * abTem[K][j] / abTem[K][K];
				}
			}
			for(int i = 0;i < 7;i++)
				for(int j = 0;j < 7;j++)
					abTem[i][j] = N_abAtA[i][j];
		}

	
		double abX[7];            //矩阵相乘
		for(int i = 0;i < 7;i++)
		{
			abX[i] = 0;
			for(int k = 0;k < 7;k++)
				abX[i] += N_abAtA[i][k] * abAtL[k];
		}

		/*更新绝对定向元素*/
		DeltaX += abX[0];
		DeltaY += abX[1];
		DeltaZ += abX[2];
		Lamda += abX[3];
		Alpha += abX[4];
		Omega += abX[5];
		Kappa += abX[6];
		


		for(int i = 0;i < 7;i++)
		{
			printf("%lf\n",abX[i]);
		}
		printf("\n");

		count++;
	
		//判断限差
		int I;
		for(I = 0;I < 7; I++)
			if(fabs(abX[I]) > 1e-4)
				break;
		if(I >= 7)
			break;
	}
		
		


	/*求地面摄影测量坐标系*/

	for(int i = 0;i < Data1.iCount;i++)
	{
		tpPoint1.Xtp[i] = DeltaX + Lamda * ( R2[0][0]*pPoint1.Xp[i] + R2[0][1]*pPoint1.Yp[i] + R2[0][2]*pPoint1.Zp[i]);
		tpPoint1.Ytp[i] = DeltaY + Lamda * ( R2[1][0]*pPoint1.Xp[i] + R2[1][1]*pPoint1.Yp[i] + R2[1][2]*pPoint1.Zp[i]);
		tpPoint1.Ztp[i] = DeltaZ + Lamda * ( R2[2][0]*pPoint1.Xp[i] + R2[2][1]*pPoint1.Yp[i] + R2[2][2]*pPoint1.Zp[i]);
	}

	for(int i = 0;i < Data2.iCount;i++)
	{
		tpPoint2.Xtp[i] = DeltaX + Lamda * ( R2[0][0]*pPoint2.Xp[i] + R2[0][1]*pPoint2.Yp[i] + R2[0][2]*pPoint2.Zp[i]);
		tpPoint2.Ytp[i] = DeltaY + Lamda * ( R2[1][0]*pPoint2.Xp[i] + R2[1][1]*pPoint2.Yp[i] + R2[1][2]*pPoint2.Zp[i]);
		tpPoint2.Ztp[i] = DeltaZ + Lamda * ( R2[2][0]*pPoint2.Xp[i] + R2[2][1]*pPoint2.Yp[i] + R2[2][2]*pPoint2.Zp[i]);
	}

	for(int i = 0;i < Data3.iCount;i++)
	{
		tpPoint3.Xtp[i] = DeltaX + Lamda * ( R2[0][0]*pPoint3.Xp[i] + R2[0][1]*pPoint3.Yp[i] + R2[0][2]*pPoint3.Zp[i]);
		tpPoint3.Ytp[i] = DeltaY + Lamda * ( R2[1][0]*pPoint3.Xp[i] + R2[1][1]*pPoint3.Yp[i] + R2[1][2]*pPoint3.Zp[i]);
		tpPoint3.Ztp[i] = DeltaZ + Lamda * ( R2[2][0]*pPoint3.Xp[i] + R2[2][1]*pPoint3.Yp[i] + R2[2][2]*pPoint3.Zp[i]);
	}


	/*求地面坐标系*/

	for(int i = 0;i < Data1.iCount;i++)
	{
		tPoint1.Xt[i] =( b*tpPoint1.Xtp[i] + a*tpPoint1.Ytp[i] ) / (lamda*lamda) + conPoint.X[0];
		tPoint1.Yt[i] =( a*tpPoint1.Xtp[i] - b*tpPoint1.Ytp[i] ) / (lamda*lamda) + conPoint.Y[0];
		tPoint1.Zt[i] = 1/lamda * tpPoint1.Ztp[i] + conPoint.Z[0];
	}

	for(int i = 0;i < Data2.iCount;i++)
	{
		tPoint2.Xt[i] = ( b*tpPoint2.Xtp[i] + a*tpPoint2.Ytp[i] ) / (lamda*lamda) + conPoint.X[0];
		tPoint2.Yt[i] = ( a*tpPoint2.Xtp[i] - b*tpPoint2.Ytp[i] ) / (lamda*lamda) + conPoint.Y[0];
		tPoint2.Zt[i] = 1/lamda * tpPoint2.Ztp[i] + conPoint.Z[0];
	}

	for(int i = 0;i < Data3.iCount;i++)
	{
		tPoint3.Xt[i] = ( b*tpPoint3.Xtp[i] + a*tpPoint3.Ytp[i] ) / (lamda*lamda) + conPoint.X[0];
		tPoint3.Yt[i] = ( a*tpPoint3.Xtp[i] - b*tpPoint3.Ytp[i] ) / (lamda*lamda) + conPoint.Y[0];
		tPoint3.Zt[i] = 1/lamda  * tpPoint3.Ztp[i] + conPoint.Z[0];
	}

	/*计算误差*/
	int checkNum[5];
	double devX[5];
	double devY[5];
	double devZ[5];

	int chK = 0;
	for(int i = 0;i < Data1.iCount;i++)
	{
		for(int j = 0;j < 5;j++)
		{
			if(Data1.iPointNum[i] == chPoint.iPointNum[j])
			{
				checkNum[chK] = Data1.iPointNum[i];
				devX[chK] = tPoint1.Xt[i] - chPoint.X[j];
				devY[chK] = tPoint1.Yt[i] - chPoint.Y[j];
				devZ[chK] = tPoint1.Zt[i] - chPoint.Z[j];
				chK++;
			}
		}
	}
	for(int i = 0;i < Data2.iCount;i++)
	{
		for(int j = 0;j < 5;j++)
		{
			if(Data2.iPointNum[i] == chPoint.iPointNum[j])
			{
				checkNum[chK] = Data2.iPointNum[i];
				devX[chK] = tPoint2.Xt[i] - chPoint.X[j];
				devY[chK] = tPoint2.Yt[i] - chPoint.Y[j];
				devZ[chK] = tPoint2.Zt[i] - chPoint.Z[j];
				chK++;
			}
		}
	}
	for(int i = 0;i < Data3.iCount;i++)
	{
		for(int j = 0;j < 5;j++)
		{
			if(Data3.iPointNum[i] == chPoint.iPointNum[j])
			{
				checkNum[chK] = Data3.iPointNum[i];
				devX[chK] = tPoint3.Xt[i] - chPoint.X[j];
				devY[chK] = tPoint3.Yt[i] - chPoint.Y[j];
				devZ[chK] = tPoint3.Zt[i] - chPoint.Z[j];
				chK++;
			}
		}
	}


	//结果输出
	FILE * fp1;
	if((fp1  = fopen("空中三角测量Result.txt","w")) == NULL)
	{
		printf("ERROR!!!!!!!!!!");
		exit(0);
	}
	else
	{
		fprintf(fp1,"模型点像空间辅助坐标：\n");
		for(int i = 0;i < Data1.iCount;i++)
		{
			fprintf(fp1,"%d         ",Data1.iPointNum[i]);
			fprintf(fp1,"%lf          ",mPoint1.Xm[i]);
			fprintf(fp1,"%lf          ",mPoint1.Ym[i]);
			fprintf(fp1,"%lf        \n",mPoint1.Zm[i]);
		}
		fprintf(fp1,"\n");
		for(int i = 0;i < Data2.iCount;i++)
		{
			fprintf(fp1,"%d          ",Data2.iPointNum[i]);
			fprintf(fp1,"%lf          ",mPoint2.Xm[i]);
			fprintf(fp1,"%lf          ",mPoint2.Ym[i]);
			fprintf(fp1,"%lf        \n",mPoint2.Zm[i]);
		}
		fprintf(fp1,"\n");
		for(int i = 0;i < Data3.iCount;i++)
		{
			fprintf(fp1,"%d          ",Data3.iPointNum[i]);
			fprintf(fp1,"%lf          ",mPoint3.Xm[i]);
			fprintf(fp1,"%lf          ",mPoint3.Ym[i]);
			fprintf(fp1,"%lf        \n",mPoint3.Zm[i]);
		}
		fprintf(fp1,"\n");

		fprintf(fp1,"\n摄影测量坐标：\n");
		for(int i = 0;i < Data1.iCount;i++)
		{
			fprintf(fp1,"%d          ",Data1.iPointNum[i]);
			fprintf(fp1,"%lf          ",pPoint1.Xp[i]);
			fprintf(fp1,"%lf          ",pPoint1.Yp[i]);
			fprintf(fp1,"%lf        \n",pPoint1.Zp[i]);
		}
		fprintf(fp1,"\n");
		for(int i = 0;i < Data2.iCount;i++)
		{
			fprintf(fp1,"%d          ",Data2.iPointNum[i]);
			fprintf(fp1,"%lf          ",pPoint2.Xp[i]);
			fprintf(fp1,"%lf          ",pPoint2.Yp[i]);
			fprintf(fp1,"%lf        \n",pPoint2.Zp[i]);
		}
		fprintf(fp1,"\n");
		for(int i = 0;i < Data3.iCount;i++)
		{
			fprintf(fp1,"%d          ",Data3.iPointNum[i]);
			fprintf(fp1,"%lf          ",pPoint3.Xp[i]);
			fprintf(fp1,"%lf          ",pPoint3.Yp[i]);
			fprintf(fp1,"%lf        \n",pPoint3.Zp[i]);
		}
		fprintf(fp1,"\n");

		fprintf(fp1,"\n地面摄影测量坐标：\n");
		for(int i = 0;i < Data1.iCount;i++)
		{
			fprintf(fp1,"%d          ",Data1.iPointNum[i]);
			fprintf(fp1,"%lf          ",tpPoint1.Xtp[i]);
			fprintf(fp1,"%lf          ",tpPoint1.Ytp[i]);
			fprintf(fp1,"%lf        \n",tpPoint1.Ztp[i]);
		}
		fprintf(fp1,"\n");
		for(int i = 0;i < Data2.iCount;i++)
		{
			fprintf(fp1,"%d          ",Data2.iPointNum[i]);
			fprintf(fp1,"%lf          ",tpPoint2.Xtp[i]);
			fprintf(fp1,"%lf          ",tpPoint2.Ytp[i]);
			fprintf(fp1,"%lf        \n",tpPoint2.Ztp[i]);
		}
		fprintf(fp1,"\n");
		for(int i = 0;i < Data3.iCount;i++)
		{
			fprintf(fp1,"%d          ",Data3.iPointNum[i]);
			fprintf(fp1,"%lf          ",tpPoint3.Xtp[i]);
			fprintf(fp1,"%lf          ",tpPoint3.Ytp[i]);
			fprintf(fp1,"%lf        \n",tpPoint3.Ztp[i]);
		}
		fprintf(fp1,"\n");

		fprintf(fp1,"\n地面坐标：\n");
		for(int i = 0;i < Data1.iCount;i++)
		{
			fprintf(fp1,"%d          ",Data1.iPointNum[i]);
			fprintf(fp1,"%lf          ",tPoint1.Xt[i]);
			fprintf(fp1,"%lf          ",tPoint1.Yt[i]);
			fprintf(fp1,"%lf        \n",tPoint1.Zt[i]);
		}
		fprintf(fp1,"\n");
		for(int i = 0;i < Data2.iCount;i++)
		{
			fprintf(fp1,"%d          ",Data2.iPointNum[i]);
			fprintf(fp1,"%lf          ",tPoint2.Xt[i]);
			fprintf(fp1,"%lf          ",tPoint2.Yt[i]);
			fprintf(fp1,"%lf        \n",tPoint2.Zt[i]);
		}
		fprintf(fp1,"\n");
		for(int i = 0;i < Data3.iCount;i++)
		{
			fprintf(fp1,"%d          ",Data3.iPointNum[i]);
			fprintf(fp1,"%lf          ",tPoint3.Xt[i]);
			fprintf(fp1,"%lf          ",tPoint3.Yt[i]);
			fprintf(fp1,"%lf        \n",tPoint3.Zt[i]);
		}
		fprintf(fp1,"\n");

		//输出误差
		fprintf(fp1,"\n所求大地坐标与检查点坐标的差：\n");
		for(int i = 0;i < 5;i++)
		{
			fprintf(fp1,"%d          ",checkNum[i]);
			fprintf(fp1,"%lf          ",devX[i]);
			fprintf(fp1,"%lf          ",devY[i]);
			fprintf(fp1,"%lf        \n",devZ[i]);
		}

	}
	fclose(fp1);


	//释放内存空间

	free(Data1.iPointNum);
	free(Data1.x1);
	free(Data1.y1);
	free(Data1.x2);
	free(Data1.y2);	
	
	free(Data2.iPointNum);
	free(Data2.x1);
	free(Data2.y1);
	free(Data2.x2);
	free(Data2.y2);	

	free(Data3.iPointNum);
	free(Data3.x1);
	free(Data3.y1);
	free(Data3.x2);
	free(Data3.y2);	

	free(mPoint1.Xm);
	free(mPoint1.Ym);
	free(mPoint1.Zm);

	free(mPoint2.Xm);
	free(mPoint2.Ym);
	free(mPoint2.Zm);
	
	free(mPoint3.Xm);
	free(mPoint3.Ym);
	free(mPoint3.Zm);

	free(pPoint1.Xp);
	free(pPoint1.Yp);
	free(pPoint1.Zp);

	free(pPoint2.Xp);
	free(pPoint2.Yp);
	free(pPoint2.Zp);

	free(pPoint3.Xp);
	free(pPoint3.Yp);
	free(pPoint3.Zp);

	free(tpPoint1.Xtp);
	free(tpPoint1.Ytp);
	free(tpPoint1.Ztp);

	free(tpPoint2.Xtp);
	free(tpPoint2.Ytp);
	free(tpPoint2.Ztp);

	free(tpPoint3.Xtp);
	free(tpPoint3.Ytp);
	free(tpPoint3.Ztp);

	free(tPoint1.Xt);
	free(tPoint1.Yt);
	free(tPoint1.Zt);

	free(tPoint2.Xt);
	free(tPoint2.Yt);
	free(tPoint2.Zt);

	free(tPoint3.Xt);
	free(tPoint3.Yt);
	free(tPoint3.Zt);

	free(conPoint.iPointNum);
	free(conPoint.X);
	free(conPoint.Y);
	free(conPoint.Z);

	free(chPoint.iPointNum);
	free(chPoint.X);
	free(chPoint.Y);
	free(chPoint.Z);
	
	return 0;
}
