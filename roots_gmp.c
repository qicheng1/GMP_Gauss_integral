#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<gmp.h>
#define N 300000
#define PI 3.14159265
#define PREC 1024
mpf_t d[1000],roots[N],ders[N],w[N],erro,p0,p1,p2,dp0,dp1,dp2;
int p[3],q[2];
int sign_p[3],sign_q[2];
int a[N],b[N],c[N],dd[N];
char* flag;
long int r,n,nn;

void taylor2(mpf_t *dux_h,mpf_t h,mpf_t d[],long int m);
void taylor1(mpf_t *ux_h,mpf_t d[],long int m);
double RK(double theta, double value_at_theta,long int emm);
double f(double theta,double x);
void find_first_root(mpf_t xe, mpf_t uxe,long int m);
void m_ders(mpf_t point,mpf_t value_at_point,mpf_t der_at_point,long int m,mpf_t h);
void newton(mpf_t point,mpf_t value_at_prepoint,mpf_t h,mpf_t erro,long int index,long int m);
void iteration_P(mpf_t xk,int m);
void coefs_Legendre();
void coefs_Hermite();
void coefs_Laguerre();

int main(int argc, char **argv)
{
	mpf_set_default_prec(PREC);	
	
        FILE *fp,*fq1,*fq2;
	mpf_t value_of_zero1,value_of_zero2,der_of_zero1,der_of_zero2,temp,zero;

	mpf_init(value_of_zero1);
	mpf_init(value_of_zero2);
	mpf_init(der_of_zero1);
	mpf_init(der_of_zero2);
	mpf_init(temp);
	mpf_init(p0);
	mpf_init(p1);
	mpf_init(p2);
        mpf_init(dp0);
	mpf_init(dp1);
	mpf_init(dp2);
        mpf_init(erro);
        mpf_init(zero);
        mpf_set_ui(zero,0);
	long int i,j,m,t;
        
	for(i=0;i<1000;i++)
	{
		mpf_init(d[i]);
	}
	for(i=0;i<N;i++)
	{
		mpf_init(roots[i]);
		mpf_init(ders[i]);
		mpf_init(w[i]);
	}
  
	mpf_init(erro);
	fp=fopen("data.txt","r");
	fq1=fopen("points.txt","w");
        fq2=fopen("coefs.txt","w");
	fscanf(fp,"%ld",&n);
	gmp_fscanf(fp,"%Ff",erro);
        flag = argv[1];   /*命令行参数*/
        if (strcmp(flag,"Legendre")==0)
           {
               p[0] = 1;  sign_p[0] = 1;
               p[1] = 0;  sign_p[1] = 1;
               p[2] = 1;  sign_p[2] = -1;
               q[0] = 0;  sign_q[0] = 1;
               q[1] = 2;  sign_q[1] = -1;
               r=n*(n+1);
               for(i=0;i<=n;i++)
                   {
                       a[i] = i+1;
                       b[i] = 0;
                       c[i] = 2*i+1;
                       dd[i] = i;
                   }
           }
       else if(strcmp(flag,"Hermite")==0)
           {
               p[0] = 1;  sign_p[0] = 1;
               p[1] = 0;  sign_p[1] = 1;
               p[2] = 0;  sign_p[2] = 1;
               q[0] = 0;  sign_q[0] = 1;
               q[1] = 2;  sign_q[1] = -1;
               r=2*n;
               for(i=0;i<=n;i++)
                   {
                       a[i] = 1;
                       b[i] = 0;
                       c[i] = 2;
                       dd[i] = 2*i;
                   }
               
           }
       else if(strcmp(flag,"Laguerre")==0)
           {
               p[0] = 0; sign_p[0] = 1;
               p[1] = 1; sign_p[1] = 1;
               p[2] = 0; sign_p[2] = 1;
               q[0] = 1; sign_q[0] = 1;
               q[1] = 1; sign_q[1] = -1;
               r=n;
               for(i=0;i<=n;i++)
                   {
                       a[i] = i+1;
                       b[i] = 2*i+1;
                       c[i] = 1; /*符号运算过程中体现*/
                       dd[i] = i;
                   }
           }
	m=30;

    if ((strcmp(flag,"Legendre")==0)||(strcmp(flag,"Hermite")==0))
      {
        mpf_set_ui(temp,0);
        iteration_P(temp,n);
	if(n%2==0)
	{
                mpf_set_ui(roots[0],0);
		mpf_set_ui(ders[0],0);
                /*gmp_printf("%Ff\n",p2);*/ 
		find_first_root(temp,p2,m);
		nn=n/2-1;
	}
	else
	{
		mpf_set_ui(roots[1],0);
		mpf_set(ders[1],dp2);
		nn=n/2;
	}
      }
      if(strcmp(flag,"Laguerre")==0)
      {
          mpf_set_ui(temp,0);
          iteration_P(temp,n);
          if(n%2==0)
	  {
                mpf_set_ui(roots[0],0);
		mpf_set(ders[0],dp2);
		find_first_root(roots[0],p2,m);
		nn=n/2-1;
	   }
	  else
	  {
		mpf_set_ui(roots[0],0);
		mpf_set(ders[0],dp2);
		find_first_root(roots[0],p2,m);
		nn=n/2;
	  }
      }
	for(i=1;i<=nn;i++)
	{
		mpf_set_d(roots[i+1],RK(PI/2,mpf_get_d(roots[i]),0));
		mpf_sub(temp,roots[i+1],roots[i]);
		newton(roots[i+1],zero,temp,erro,i+1,m);
	}
	mpf_set_ui(w[0],0);
        fprintf(fq1,"%ld\n",nn+1);
	for(i=1;i<=nn+1;i++)
	{
		gmp_fprintf(fq1,"%.315Ff\n",roots[i]);
	}

        /*求系数*/
        if(strcmp(flag,"Legendre")==0) coefs_Legendre();
        if(strcmp(flag,"Hermite")==0) coefs_Hermite();
        if(strcmp(flag,"Laguerre")==0) coefs_Laguerre();

	for(i=1;i<=nn+1;i++)
	{
		gmp_fprintf(fq2,"%.315Ff\n",w[i]);
	}
	mpf_clear(value_of_zero1);
	mpf_clear(value_of_zero2);
	mpf_clear(der_of_zero1);
	mpf_clear(der_of_zero2);
	mpf_clear(temp);
	for(i=0;i<100;i++)
	{
		mpf_clear(d[i]);
	}
	for(i=0;i<N;i++)
	{
		mpf_clear(roots[i]);
		mpf_clear(ders[i]);
		mpf_clear(w[i]);
	}
	mpf_clear(erro);
	mpf_clear(p0);
	mpf_clear(p1);
	mpf_clear(p2);
        mpf_clear(dp0);
	mpf_clear(dp1);
	mpf_clear(dp2);
	fclose(fp);
	fclose(fq1);
        fclose(fq2);
	return 0;
}
void iteration_P(mpf_t xk,int m)
{
        int i,j;
        mpf_t temp;
        mpf_init(temp);
        mpf_set_ui(p2,0);
        mpf_set_ui(dp2,0);
        if(strcmp(flag,"Legendre")==0)
          {
            mpf_set_ui(p0,1); mpf_set_ui(dp0,0);
            mpf_set(p1,xk); mpf_set_ui(dp1,1);
          }
        if(strcmp(flag,"Hermite")==0)
          {
            mpf_set_ui(p0,1); mpf_set_ui(dp0,0);
            mpf_mul_ui(p1,xk,2); mpf_set_ui(dp1,2);
          }
        if(strcmp(flag,"Laguerre")==0)
          {
            mpf_set_ui(p0,1); mpf_set_ui(dp0,0);
            mpf_ui_sub(p1,1,xk); mpf_set_si(dp1,-1);
          }
        for(i=1;i<m;i++)
          {
            mpf_mul_ui(temp,xk,c[i]);
            if (strcmp(flag,"Laguerre")==0) mpf_neg(temp,temp);
            mpf_add_ui(temp,temp,b[i]);
            mpf_mul(p2,temp,p1);
            mpf_mul_ui(temp,p0,dd[i]); 
            mpf_sub(p2,p2,temp);
            mpf_div_ui(p2,p2,a[i]);

            mpf_mul_ui(temp,xk,c[i]);
            if (strcmp(flag,"Laguerre")==0) mpf_neg(temp,temp);
            mpf_add_ui(temp,temp,b[i]);
            mpf_mul(dp2,temp,dp1);
            mpf_mul_ui(temp,dp0,dd[i]); 
            mpf_sub(dp2,dp2,temp);
            mpf_mul_ui(temp,p1,c[i]);
            if (strcmp(flag,"Laguerre")==0) mpf_neg(temp,temp);
            mpf_add(dp2,dp2,temp);
            mpf_div_ui(dp2,dp2,a[i]);

            mpf_set(p0,p1);
	    mpf_set(p1,p2);
            mpf_set(dp0,dp1);
            mpf_set(dp1,dp2);
          }
       mpf_clear(temp);
}

void coefs_Legendre()
{ 
    int i;
    for(i=1;i<=nn+1;i++)
       {
          iteration_P(roots[i],n);
          mpf_mul(w[i],roots[i],roots[i]);
	  mpf_ui_sub(w[i],1,w[i]);
          mpf_mul(w[i],w[i],dp2); 
          mpf_mul(w[i],w[i],dp2);
          mpf_ui_div(w[i],2,w[i]);
       }
}
void coefs_Hermite()
{
    int i;
    double t;
    mpf_t temp;
    mpf_init(temp);
    t = 1;
    for (i=1;i<=n;i++);
        t = t*2*i;
    t = t*2*sqrt(PI);
    mpf_set_d(temp,t);
    for(i=1;i<=nn+1;i++)
       {
          iteration_P(roots[i],n);
          mpf_mul(w[i],dp2,dp2); 
          mpf_div(w[i],temp,w[i]); 
       }
    mpf_clear(temp);
}
void coefs_Laguerre()
{
    int i;
    for(i=1;i<=nn+1;i++)
       {
          iteration_P(roots[i],n);
          mpf_mul(w[i],dp2,dp2); 
          mpf_mul(w[i],w[i],roots[i]);
          mpf_ui_div(w[i],1,w[i]); 
       }
}
void taylor2(mpf_t *p,mpf_t h,mpf_t d[],long int m)
{
	long int i;
	mpf_t temp,dux_h;
	mpf_init(temp);
	mpf_init(dux_h);
	mpf_div(dux_h,d[1],h);
	for(i=2;i<=m;i++)
	{
		mpf_mul_ui(temp,d[i],i);
		mpf_div(d[i],temp,h);
		mpf_add(dux_h,dux_h,d[i]);
	}
	mpf_set(*p,dux_h);
	mpf_clear(temp);
	mpf_clear(dux_h);
}

void taylor1(mpf_t *p,mpf_t d[],long int m)
{
	long int i;
	mpf_t ux_h;
	mpf_init(ux_h);
	mpf_set_ui(ux_h,0);
	for(i=0;i<=m;i++)
	{
		mpf_add(ux_h,ux_h,d[i]);
	}
	mpf_set(*p,ux_h);
	mpf_clear(ux_h);
}

double RK(double theta, double value_at_theta,long int emm)
{
	double h,x,k1,k2,k3,k4,nn;
 	long int i;
 	h=-PI/100;	
 	if(emm==0)
 	nn=100; 
 	else
 	nn=50;
 	x=value_at_theta;
 	for(i=1;i<=nn;i++)
 	{
 	 	k1=f(theta,x);
	 	k2=f(theta+h/2,x+h*k1/2);
	 	k3=f(theta+h/2,x+h*k2/2);
	 	k4=f(theta+h,x+h*k3);
	 	x=x+h*(k1+2*k2+2*k3+k4)/6;
	 	theta=theta+h;	
	}
	return x;
}

double f(double theta,double x)
{
	double z,z1;
        z1 = r;
	z=sqrt(z1/(p[0]*sign_p[0]+p[1]*sign_p[1]*x+p[2]*sign_p[2]*x*x))+(-(p[1]*sign_p[1]+2*p[2]*sign_p[2]*x)+2*(q[0]*sign_q[0]+q[1]*sign_q[1]*x))*sin(2*theta)/(4*(p[0]*sign_p[0]+p[1]*sign_p[1]*x+p[2]*sign_p[2]*x*x));
	z=-1.0/z;
	return z;
}

void find_first_root(mpf_t xe, mpf_t uxe,long int m) 
{
	double x,xed;
	mpf_t x1,temp;
	mpf_init(x1);
	mpf_init(temp);
	xed=mpf_get_d(xe);
	x=RK(0,xed,1);
	mpf_set_d(x1,x);
	mpf_sub(temp,x1,xe);
	newton(x1,uxe,temp,erro,1,m);
	mpf_clear(x1);
	mpf_clear(temp);
}

void m_ders(mpf_t point,mpf_t value_at_point,mpf_t der_at_point,long int m,mpf_t h)   
{
	long int k,t;
	mpf_t temp1,temp2;
	mpf_init(temp1);
        mpf_init(temp2);
	mpf_set(d[0],value_at_point);
	mpf_mul(d[1],der_at_point,h);
	if(m<=n) t=m-2; else t=n-2;
	for(k=0;k<=t;k++)
	{
            mpf_mul_ui(temp1,point,2*p[2]);
            if (sign_p[2]==-1)
                   mpf_neg(temp1,temp1);
            if (sign_p[1]==-1)
                   mpf_sub_ui(temp1,temp1,p[1]);
            else
                   mpf_add_ui(temp1,temp1,p[1]);
            mpf_mul_ui(temp1,temp1,k);
            mpf_mul_ui(temp2,point,q[1]);
            if (sign_q[1]==-1)
                   mpf_neg(temp2,temp2);
            if (sign_q[0]==-1)
                   mpf_sub_ui(temp2,temp2,q[0]);
            else
                   mpf_add_ui(temp2,temp2,q[0]);
            mpf_add(temp1,temp1,temp2);
            mpf_div_ui(temp1,temp1,(k+2));
            mpf_neg(temp1,temp1);
            mpf_mul(temp1,temp1,h);
            mpf_mul(temp1,temp1,d[k+1]);
            mpf_set_d(temp2,(k*(k-1)*p[2]*sign_p[2]+k*q[1]*sign_q[1]+r)*1.0/((k+2)*(k+1)));
            mpf_mul(temp2,h,temp2);
            mpf_mul(temp2,temp2,h);
            mpf_mul(temp2,temp2,d[k]);
            mpf_sub(d[k+2],temp1,temp2);
            mpf_mul(temp1,point,point);
            mpf_mul_ui(temp1,temp1,p[2]);
            if (sign_p[2]==-1)
                mpf_neg(temp1,temp1);
            mpf_mul_ui(temp2,point,p[1]);
            if (sign_p[1]==-1)
                mpf_neg(temp2,temp2);
            mpf_add(temp1,temp1,temp2);
            if (sign_p[0]==-1)
                mpf_sub_ui(temp1,temp1,p[0]);
            else
                mpf_add_ui(temp1,temp1,p[0]);
            mpf_div(d[k+2],d[k+2],temp1);
    }
    mpf_clear(temp1);
    mpf_clear(temp2);
}

void newton(mpf_t point,mpf_t value_at_prepoint,mpf_t h,mpf_t erro,long int index,long int m)
{
	mpf_t t,t1,t2,temp;
	mpf_init(t);
	mpf_init(t1);
	mpf_init(t2);
	mpf_init(temp);
	m_ders(roots[index-1],value_at_prepoint,ders[index-1],m,h);
	mpf_set(t,point);
	taylor1(&t1,d,m);
	taylor2(&t2,h,d,m);
	mpf_div(temp,t1,t2);
	mpf_sub(point,point,temp);
	mpf_sub(temp,point,t);
	mpf_abs(temp,temp);

	while(mpf_cmp(temp,erro)>0)
	{
	    mpf_sub(temp,point,t);
	    mpf_add(h,h,temp);
	    m_ders(roots[index-1],value_at_prepoint,ders[index-1],m,h);
	    mpf_set(t,point);
	    taylor1(&t1,d,m);
	    taylor2(&t2,h,d,m);
	    mpf_div(temp,t1,t2);
	    mpf_sub(point,point,temp);
	    mpf_sub(temp,point,t);
	    mpf_abs(temp,temp);
        }
        mpf_set(roots[index],point);
        mpf_sub(temp,point,t);
	mpf_add(h,h,temp);
	m_ders(roots[index-1],value_at_prepoint,ders[index-1],m,h);
	taylor2(&t2,h,d,m);
	mpf_set(ders[index],t2);
	mpf_clear(t);
	mpf_clear(t1);
	mpf_clear(t2);
	mpf_clear(temp);

}
