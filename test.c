#include<stdio.h>
#include<stdlib.h>
#include<gmp.h>
#define N 30000
int main(int argc,char* argv[])
{
    FILE *fp1,*fp2,*fq;
    int i,j,k,l,m,n,nn;
    mpf_t roots[N],coefs[N],temp,sum,erro,intl;
    mpf_set_default_prec(1024);
    for(i=0;i<N;i++)               //初始化各变量
       {
            mpf_init(roots[i]);
            mpf_init(coefs[i]);
       }
    mpf_init(erro);
    mpf_init(temp);
    mpf_init(sum);
    mpf_init(intl);
    fp1=fopen("points.txt","r");
    fp2=fopen("coefs.txt","r");
    fq=fopen("erros.txt","w");
    fscanf(fp1,"%d",&n);
    printf("hahaha\n");
    for(i=0;i<n;i++)
       {
            gmp_fscanf(fp1,"%Ff",roots[i]);
            gmp_fscanf(fp2,"%Ff",coefs[i]);
       }
    if(mpf_cmp_ui(roots[0],0)==0)
        nn=2*n-1;
    else
        nn=2*n;
    //计算不超过2*nn+1次多项式的数值积分结果（只计算偶次）
    for(i=2;i<=2*nn+1;i=i+2)
        {
            mpf_set_ui(sum,0);
            for(j=0;j<n;j++)
                {
                    if(mpf_cmp_ui(roots[j],0)==0)
                         {
                             mpf_pow_ui(temp,roots[j],i);
                             mpf_mul(temp,temp,coefs[j]);
                             mpf_add(sum,sum,temp);
                         }
                    else
                         {
                             mpf_pow_ui(temp,roots[j],i);
                             mpf_mul(temp,temp,coefs[j]);
                             mpf_mul_ui(temp,temp,2);
                             mpf_add(sum,sum,temp);
                         }
                }
            mpf_set_ui(intl,2);
            mpf_div_ui(intl,intl,i+1);
            mpf_sub(erro,sum,intl);
            mpf_abs(erro,erro);
            gmp_fprintf(fq,"n = %d, erro = %.16Fe\n",i,erro);
        }

    fclose(fp1);
    fclose(fp2);
    fclose(fq);
} 
