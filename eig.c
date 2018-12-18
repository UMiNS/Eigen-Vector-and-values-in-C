//数值计算程序-特征值和特征向量 

#include "stdio.h" 
#include "math.h" 
#define M 4  //方阵的行数 列数
#define e0  0.0000001//ε0为要求的精度
#define N  10//最大迭代次数


//求实对称矩阵的特征值及特征向量的雅格比法 
//利用雅格比(Jacobi)方法求实对称矩阵的全部特征值及特征向量 
//返回值小于0表示超过迭代jt次仍未达到精度要求 
//返回值大于0表示正常返回 
//a-长度为n*n的数组，存放实对称矩阵，返回时对角线存放n个特征值 
//n-矩阵的阶数 
//u-长度为n*n的数组，返回特征向量(按列存储) 
//eps-控制精度要求 
//jt-整型变量，控制最大迭代次数 
int eejcb(double a[],int n,double v[],double eps,int jt) 
{ 
    int i,j,p,q,u,w,t,s,l; 
    double fm,cn,sn,omega,x,y,d; 
    l=1; 
    for (i=0; i<=n-1; i++) 
    { 
        v[i*n+i]=1.0; 
        for (j=0; j<=n-1; j++) 
        { 
            if (i!=j) 
            { 
                v[i*n+j]=0.0; 
            }    
        } 
    } 
    while (1==1) 
    { 
        fm=0.0; 
        for (i=0; i<=n-1; i++) 
        { 
            for (j=0; j<=n-1; j++) 
            { 
                d=fabs(a[i*n+j]); 
                if ((i!=j)&&(d>fm)) 
                { 
                    fm=d; 
                    p=i; 
                    q=j; 
                } 
            } 
        } 
        if (fm<eps) 
        { 
            return(1); 
        } 
        if (l>jt) 
        { 
            return(-1); 
        } 
        l=l+1; 
        u=p*n+q; 
        w=p*n+p; 
        t=q*n+p; 
        s=q*n+q; 
        x=-a[u]; 
        y=(a[s]-a[w])/2.0; 
        omega=x/sqrt(x*x+y*y); 
        if (y<0.0) 
        { 
            omega=-omega; 
        } 
        sn=1.0+sqrt(1.0-omega*omega); 
        sn=omega/sqrt(2.0*sn); 
        cn=sqrt(1.0-sn*sn); 
        fm=a[w]; 
        a[w]=fm*cn*cn+a[s]*sn*sn+a[u]*omega; 
        a[s]=fm*sn*sn+a[s]*cn*cn-a[u]*omega; 
        a[u]=0.0; 
        a[t]=0.0; 
        for (j=0; j<=n-1; j++) 
        { 
            if ((j!=p)&&(j!=q)) 
            { 
                u=p*n+j; 
                w=q*n+j; 
                fm=a[u]; 
                a[u]=fm*cn+a[w]*sn; 
                a[w]=-fm*sn+a[w]*cn; 
            } 
        } 
        for (i=0; i<=n-1; i++) 
        { 
            if ((i!=p)&&(i!=q)) 
            { 
                u=i*n+p; 
                w=i*n+q; 
                fm=a[u]; 
                a[u]=fm*cn+a[w]*sn; 
                a[w]=-fm*sn+a[w]*cn; 
            } 
        } 
        for (i=0; i<=n-1; i++) 
        { 
            u=i*n+p; 
            w=i*n+q; 
            fm=v[u]; 
            v[u]=fm*cn+v[w]*sn; 
            v[w]=-fm*sn+v[w]*cn; 
        } 
    } 
    return(1); 
} 

//矩阵的打印
void printMatrix(double a[M*M], int m, int n)
{
    for (int i = 0; i<m; i++)
    {
        for (int j = 0; j<n; j++)
        {
            printf("%lg15  ", a[i*M+j]);
        }
        printf("\n");
    }
}
//向量的打印
void printVector(double a[], int m)
{
    for (int i = 0; i < m; i++)
    {
        printf("%lf  ", a[i]);
    }
}

//主函数
int main(void)
{
    double a[M*M] = {     9.972485169120526e+000 ,  -2.082338684083218e-003 ,  -8.910478324391178e-001,    4.494330653839068e-006,
    -2.082338684083218e-003,    8.916849167386474e+002,   -3.338046458058518e-002 ,   1.904868847946309e+001,
    -8.910478324391178e-001,   -3.338046458058518e-002,    8.927934869770321e+002 ,   1.204640400767189e+000,
    4.494330653839068e-006 ,   1.904868847946309e+001  ,  1.204640400767189e+000  ,  1.112125299097807e+001 };//待求特征值和特征向量的矩阵
    double v[M*M]={0};
    int ret;
   
    printf("待求特征值和特征向量的矩阵A:\n");
    printMatrix(a, M, M);
    printf("\n");
    //求实对称矩阵的特征值及特征向量的雅格比法 
    //利用雅格比(Jacobi)方法求实对称矩阵的全部特征值及特征向量 
    //返回值小于0表示超过迭代jt次仍未达到精度要求 
    //返回值大于0表示正常返回 
    //a-长度为n*n的数组，存放实对称矩阵，返回时对角线存放n个特征值 
    //n-矩阵的阶数 
    //u-长度为n*n的数组，返回特征向量(按列存储) 
    //eps-控制精度要求 
    //jt-整型变量，控制最大迭代次数 
    ret = eejcb(a,M,v,e0,N) ;
    printf("eigen values:\n");
    printMatrix(a, M, M);
    printf("\neigen vector:\n");
    printMatrix(v, M, M);
}






