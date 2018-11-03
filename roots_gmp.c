#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<gmp.h>
#define N 300000
#define PI 3.14159265
#define PREC 1024
mpf_t d[1000],roots[N],ders[N],erro;
int p[3]={1,0,-1},q[2]={0,-2};
long int r,n;

void taylor2(mpf_t *dux_h,mpf_t h,mpf_t d[],long int m);
void taylor1(mpf_t *ux_h,mpf_t d[],long int m);
double RK(double theta, double value_at_theta,long int emm);
double f(double theta,double x);
void find_first_root(mpf_t xe, mpf_t uxe,long int m);
void m_ders(mpf_t point,mpf_t value_at_point,mpf_t der_at_point,long int m,mpf_t h);
void newton(mpf_t point,mpf_t value_at_prepoint,mpf_t h,mpf_t erro,long int index,long int m);

int main(int argc, char *argv)
{
	mpf_set_default_prec(PREC);	
	
        FILE *fp,*fq1,*fq2;
	mpf_t value_of_zero1,value_of_zero2,der_of_zero1,der_of_zero2,temp,p0,p1,p2,w[N],zero;

	mpf_init(value_of_zero1);
	mpf_init(value_of_zero2);
	mpf_init(der_of_zero1);
	mpf_init(der_of_zero2);
	mpf_init(temp);
	mpf_init(p0);
	mpf_init(p1);
	mpf_init(p2);
        mpf_init(erro);
        mpf_init(zero);
        mpf_set_ui(zero,0);
	long int i,j,m,t,nn;
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
	r=n*(n+1);
	m=120;
	mpf_set_ui(value_of_zero1,1);
	mpf_set_ui(value_of_zero2,0);
	mpf_set_ui(der_of_zero1,0);
	mpf_set_ui(der_of_zero2,1);
	if (n%2==0)
	{
		for(i=0;i<n-1;i=i+2)
		{
			mpf_mul_ui(value_of_zero1,value_of_zero1,i+1);
	    	        mpf_div_ui(value_of_zero1,value_of_zero1,i+2);
			mpf_neg(value_of_zero1,value_of_zero1); 
		}
	}
	else
	{
		for(i=1;i<n-1;i=i+2)
		{
			mpf_mul_ui(value_of_zero1,value_of_zero1,i);
			mpf_div_ui(value_of_zero1,value_of_zero1,i+1);
			mpf_neg(value_of_zero1,value_of_zero1);
			mpf_set(der_of_zero2,value_of_zero1);
			mpf_mul_ui(der_of_zero2,der_of_zero2,2*i+3);
			mpf_div_ui(der_of_zero2,der_of_zero2,i+2);
			mpf_set(temp,der_of_zero2);
			mpf_mul_ui(temp,temp,i+1);
			mpf_div_ui(temp,temp,i+2);
			mpf_sub(der_of_zero2,der_of_zero2,temp);
		}
	}
       
	
	if(n%2==0)
	{
                mpf_set_ui(roots[0],0);
		mpf_set_ui(ders[0],0);
		mpf_set_ui(temp,0);
		find_first_root(temp,value_of_zero1,m);
		nn=n/2-1;
	}
	else
	{
		mpf_set_ui(roots[1],0);
		mpf_set(ders[1],der_of_zero2);
		nn=n/2;
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
		mpf_set_ui(p0,1);
		mpf_set(p1,roots[i]);
		for(j=2;j<n;j++)
		{
			mpf_mul_ui(p2,roots[i],2*j-1);
			mpf_mul(p2,p2,p1);
			mpf_mul_ui(temp,p0,j-1);
			mpf_sub(p2,p2,temp);
			mpf_div_ui(p2,p2,j);
			mpf_set(p0,p1);
			mpf_set(p1,p2);
		}
		mpf_mul_ui(p2,p2,n);
		mpf_div_ui(p2,p2,n+1);
		mpf_mul(w[i],roots[i],roots[i]);
		mpf_ui_sub(w[i],1,w[i]);
		mpf_mul_ui(w[i],w[i],2);
		mpf_div_ui(w[i],w[i],n+1);
		mpf_div_ui(w[i],w[i],n+1);
		mpf_div(w[i],w[i],p2);
		mpf_div(w[i],w[i],p2);
	}
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
	fclose(fp);
	fclose(fq1);
        fclose(fq2);
	return 0;
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
	double z,z1,z2;
	z1=n;
	z=sqrt(z1*(z1+1)/(1-x*x))-x*sin(2*theta)/(2*(1-x*x));
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
	mpf_t temp;
	mpf_init(temp);
	mpf_set(d[0],value_at_point);
	mpf_mul(d[1],der_at_point,h);
	if(m<=n) t=m-2; else t=n-2;
	for(k=0;k<=t;k++)
	{
	    mpf_set_ui(d[k+2],2);
	    mpf_mul(d[k+2],d[k+2],h);
	    mpf_mul_ui(d[k+2],d[k+2],k+1);
	    mpf_mul(d[k+2],d[k+2],point);
	    mpf_div_ui(d[k+2],d[k+2],k+2);
	    mpf_mul(d[k+2],d[k+2],d[k+1]);
	    mpf_mul(temp,h,h);
	    mpf_mul_ui(temp,temp,n-k);
	    mpf_neg(temp,temp);
	    mpf_mul_ui(temp,temp,k+n+1);
	    mpf_div_ui(temp,temp,k+1);
	    mpf_div_ui(temp,temp,k+2);
	    mpf_mul(temp,temp,d[k]);
	    mpf_add(d[k+2],d[k+2],temp);
	    
	    mpf_mul(temp,point,point);
	    mpf_ui_sub(temp,1,temp);
	    mpf_div(d[k+2],d[k+2],temp);
    }
    mpf_clear(temp);
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
