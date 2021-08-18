#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include "ran2.c"

#define boxl 100
#define L boxl
#define boxl2 (boxl/2.0)
#define N 16
#define DIR 2
#define pi (22.0/7.0)
#define nmax 17
#define vmax 1.0
#define r_r 1.0
#define Knum 4
#define eps 0.000001
#define delt 0.01
#define tstar 20
#define Ntrain 16000

long int seed=7705832;
int nrep,no,na,nn,listrep[N],listo[N],lista[N],listnn[N];
double r[N][DIR],vel[N][DIR],rinit[N][DIR],vinit[N][DIR],d[N][DIR],A[Knum],B[Knum],Q0,Q1,Q2,r_a,r_o,testf,MAX[Ntrain],MIN[Ntrain];
char outfile[100];

/** generates a Gaussian random variable with mean 'u' and standard deviation 'd'*/

FILE *fpr;
double gaussian(double u, double d)
{
    static double t = 0.0;
    double x, w1, w2, r;
    if( t == 0 )
    {
        do{
            w1 = 2.0 * ran2(&seed) - 1.0;
            w2 = 2.0 * ran2(&seed) - 1.0;
            r = w1 * w1 + w2 * w2;
        }
        while( r >= 1.0 );
        r = sqrt( -2.0*log(r) / r );
        t = w2 * r;
        return( u + w1 * r * d );
    }
    else
    {
        x = t;
        t = 0.0;
        return( u + x * d);
    }
}

double anint(double x)
{
    double d;
    if(x>=0.50)
        d=1.00;
    else
    {
        if(x<=(-0.50))
            d=-1.00;
        else
            d=0.00;
    }
    return d;
}

double correct_pbc(double rr)
{
    double rr_n;
    if(rr > boxl2)
        rr_n=rr-boxl;
    else if(rr < -boxl2)
        rr_n=rr+boxl;
    else
        rr_n=rr;
    return rr_n;
}

void print_configuration(int time)
{
    int i,j;
    FILE *fpw;
    sprintf(outfile,"./output/initial_positions_t%d.dat", time);
    fpw=fopen(outfile,"w");
    for(i=0;i<N;i++)
        fprintf(fpw,"%lf\t%lf\n",r[i][0],r[i][1]);
    fclose(fpw);
}

void get_field_amps(int time)
{
    int i,c0;
    double c1,c2,c3,c4,c5,c6,c7,c8;
    FILE *fpr;
    fpr=fopen("field_amp01.dat","r");
    fscanf(fpr,"%d%lf%lf%lf%lf%lf%lf%lf%lf\n",&c0,&c1,&c2,&c3,&c4,&c5,&c6,&c7,&c8);
    A[0]=c1; A[1]=c2; A[2]=c3; A[3]=c4;
    B[0]=c5; B[1]=c6; B[2]=c7; B[3]=c8;

    /*double Bmean,Bsd;
    Bmean=0.4; Bsd=0.1;
    for(i=0;i<Knum;i++)
    {
        B[i]=gaussian(Bmean,Bsd);
        A[i]=gaussian(Bmean,Bsd);
    }*/

    /*A[0]=gaussian(0,Bsd);
     A[1]=gaussian(0,Bsd);
     A[2]=gaussian(0,Bsd);
     A[3]=gaussian(Bmean,Bsd);*/
    /*fprintf(fpw,"%d\t",time);
    for(i=0;i<Knum;i++)
        fprintf(fpw,"%.8lf\t",A[i]);
    for(i=0;i<Knum-1;i++)
        fprintf(fpw,"%.8lf\t",B[i]);
    fprintf(fpw,"%.8lf\n",B[3]);
    fclose(fpw);*/
}

/*double sigmoid(double x)
{
    double g,sl,sh;
    sl=2.0; sh=0.2;
    g=1.0/(1+exp(-sl*(x-sh)));
    return g;
}*/

double evaluate_field(int s, int rs)
{
    int j;
    int kx,ky,kmax,kgap;
    double F, kdotx,Ffin;
    kmax=nmax; kgap=1;
    j=0;
    F=0;
    for(kx=0;kx<=kmax;kx=kx+kgap)
    {
        for(ky=0;ky<=kmax;ky=ky+kgap)
        {
            kdotx=2.0*pi*(kx*r[s][0]+ky*r[s][1])/L;
            F=F+A[j]*cos(kdotx)+B[j]*sin(kdotx);
            j++;
        }
    }

    /*F=A[0] + A[1]*cos(2.0*pi*r[s][1]/L) + B[1]*sin(2.0*pi*r[s][1]/L) +
    A[2]*cos(2.0*pi*r[s][0]/L) + B[2]*sin(2.0*pi*r[s][0]/L) +
    A[3]*cos(2.0*pi*(r[s][1]+r[s][0])/L) + B[3]*sin(2.0*pi*(r[s][1]+r[s][0])/L);*/

    //Ffin=sigmoid(F);
    Ffin=(F-MIN[rs])/MAX[rs];
    //printf("rs=%d F=%lf Ffin=%lf r=(%lf %lf) A=%lf %lf %lf %lf, B=%lf %lf %lf %lf\n",rs,F, Ffin,r[s][0],r[s][1],A[0],A[1],A[2],A[3],B[0],B[1],B[2],B[3]);
    // if(fabs(F) > 1.0)
    //   printf("A=%lf %lf %lf %lf, B=%lf %lf %lf %lf, r=(%lf %lf), F=%lf\n",A[0],A[1],A[2],A[3],B[0],B[1],B[2],B[3],r[s][0],r[s][1],F);
    //printf("A11=%lf B11=%lf %lf %lf F1=%lf F=%lf\n",A[3],B[3],A[3]*cos(2*pi*(r[s][1]+r[s][0])/L),B[3]*sin(2*pi*(r[s][1]+r[s][0])/L),A[3]*cos(2*pi*(r[s][1]+r[s][0])/L)+B[3]*sin(2*pi*(r[s][1]+r[s][0])/L)+A[0]+A[1]*cos(2*pi*r[s][1]/L)+B[1]*sin(2*pi*r[s][1]/L)+A[2]*cos(2*pi*r[s][0]/L)+B[2]*sin(2*pi*r[s][0]/L),F);
   // if(fabs(Ffin) < eps)
     //   printf("s=%d kdotx=%lf F=%lf Ffin=%lf A=(%lf %lf %lf %lf) B=(%lf %lf %lf %lf)\n",s,kdotx,F,Ffin,A[0],A[1],A[2],A[3],B[0],B[1],B[2],B[3]);
    return Ffin;
}

double * evaluate_grad(int s)
{
    int j;
    double kx,ky,kmax,kgap,Fx,Fy,fac,mod,kdotx;
    static double grad[DIR];
    kmax=2*pi*nmax/L; kgap=2*pi/L;
    j=0;
    Fx=0; Fy=0;
    for(kx=0;kx<(kmax+eps);kx=kx+kgap)
    {
        for(ky=0;ky<(kmax+eps);ky=ky+kgap)
        {
            kdotx=kx*r[s][0]+ky*r[s][1];
            fac=-A[j]*sin(kdotx)+B[j]*cos(kdotx);
            Fx=Fx+kx*fac;
            Fy=Fy+ky*fac;
            j++;
        }
    }
    mod=sqrt(Fx*Fx+Fy*Fy);
    grad[0]=Fx/mod; grad[1]=Fy/mod;

    return grad;
}

/*void initialize_positions()
 {
 int i,j;
 double center[DIR],theta,d,rad;
 rad=(sqrt(N/pi));

 for(j=0;j<DIR;j++)
 center[j]=(boxl-rad)*(ran2(&seed)-0.50);

 for(i=0;i<N;i++)
 {
 d=rad*ran2(&seed);
 theta=2*pi*ran2(&seed);
 r[i][0]=correct_pbc(center[0]+d*cos(theta));
 r[i][1]=correct_pbc(center[1]+d*sin(theta));

 for(j=0;j<DIR;j++)
 rinit[i][j]=r[i][j];
 }
 }*/

void initialize_positions()
{
    FILE *fpw;
    fpw=fopen("initial_pos.dat","w");
    int i,j,k,overlap;
    double center[DIR],theta,d,rad,xt[DIR],sepc[DIR],thrs,sepcsq,sepcsqrt;
    thrs=0.2;
    rad=(sqrt(N/pi));

    for(j=0;j<DIR;j++)
        center[j]=(boxl-rad)*(ran2(&seed)-0.50)+60*(1-j);
    // center[0]=(boxl-rad)*(ran2(&seed)-0.50)+75;
    //center[1]=(boxl-rad)*(ran2(&seed)-0.50)+33;
    j=0;
    while(j<N)
    {

        d=rad*ran2(&seed);
        theta=2*pi*ran2(&seed);
        xt[0]=center[0]+d*cos(theta);
        xt[1]=center[1]+d*sin(theta);

        overlap=0;
        for(i=0;i<j;i++)
        {
            sepcsq=0.0;
            for(k=0;k<DIR;k++)
            {
                sepc[k]=xt[k]-r[i][k];
                sepc[k]=sepc[k]-(boxl)*anint(sepc[k]/(boxl));
                sepcsq=sepcsq+sepc[k]*sepc[k];
            }
            sepcsqrt=sqrt(sepcsq);
            if(sepcsqrt < thrs)
                overlap++;
            if(overlap > 0)
                break;
        }
        if(overlap == 0)
        {
            for(k=0;k<DIR;k++)
            {
                r[j][k]=xt[k];
                rinit[j][k]=r[j][k];
            }
            fprintf(fpw,"%lf %lf\n",r[j][0],r[j][1]);
            j++;

        }
    }
    fclose(fpw);
}

void initialize_velocities(int rs)
{
    int i,j;
    double field,*p,graddir[DIR];

    for(i=0;i<N;i++)
    {
        field=evaluate_field(i,rs);
        p=evaluate_grad(i);
        for(j=0;j<DIR;j++)
        {
            graddir[j]=*(p+j);
            vel[i][j]=-vmax*field*graddir[j];
        }
        for(j=0;j<DIR;j++)
            vinit[i][j]=vel[i][j];
    }
}

double * normalize(double V[DIR])
{
    int k;
    double rdiffsqrt,rdiffsq;
    static double NV[DIR];

    rdiffsq=0;
    for(k=0;k<DIR;k++)
        rdiffsq=rdiffsq+V[k]*V[k];
    rdiffsqrt=sqrt(rdiffsq);
    if(rdiffsqrt > eps)
        for(k=0;k<DIR;k++)
            NV[k]=V[k]/rdiffsqrt;
    else
        for(k=0;k<DIR;k++)
            NV[k]=V[k];
    return NV;
}

void find_direction(int p, int n1, int n2, int n3)
{
    int j,k;
    double rdiff[DIR],rdiffsqrt,rdiffsq,*tmp,tmpd[DIR],tmpv[DIR];

    for(k=0;k<DIR;k++)
        d[p][k]=0;

    if(n1 > 0)
    {
        for(j=0; j<n1; j++)
        {
            rdiffsq=0.0;
            for(k=0;k<DIR;k++)
            {
                rdiff[k]=r[listrep[j]][k]-r[p][k];
                rdiff[k]=rdiff[k]-(boxl)*anint(rdiff[k]/(boxl));
                rdiffsq=rdiffsq+rdiff[k]*rdiff[k];
            }
            rdiffsqrt=sqrt(rdiffsq);
            if(rdiffsqrt > eps)
                for(k=0;k<DIR;k++)
                    rdiff[k]=rdiff[k]/rdiffsqrt;

            for(k=0;k<DIR;k++)
                d[p][k]=d[p][k]-rdiff[k];
        }
    }
    else if((n2+n3) > 0)
    {
        if(n2 > 0)
        {
            for(j=0; j<n2; j++)
            {
                rdiffsq=0.0;
                for(k=0;k<DIR;k++)
                    rdiffsq=rdiffsq+vel[listo[j]][k]*vel[listo[j]][k];
                rdiffsqrt=sqrt(rdiffsq);
                if(rdiffsqrt > eps)
                {
                    //printf("Orien: p=%d v[11]=(%lf %lf) v[15]=%lf %lf rdiffsqrt=%lf listo=%d\n",p,vel[11][0],vel[11][1],vel[15][0],vel[15][1],rdiffsqrt,listo[j]);
                    for(k=0;k<DIR;k++)
                        tmpv[k]=vel[listo[j]][k]/rdiffsqrt;
                }
                else
                    for(k=0;k<DIR;k++)
                        tmpv[k]=vel[listo[j]][k];
                for(k=0;k<DIR;k++)
                    d[p][k]=d[p][k]+tmpv[k];
            }
        }
        if(n3 > 0)
        {
            for(j=0; j<n3; j++)
            {
                rdiffsq=0.0;
                for(k=0;k<DIR;k++)
                {
                    rdiff[k]=r[lista[j]][k]-r[p][k];
                    rdiff[k]=rdiff[k]-(boxl)*anint(rdiff[k]/(boxl));
                    rdiffsq=rdiffsq+rdiff[k]*rdiff[k];
                }
                rdiffsqrt=sqrt(rdiffsq);
                if(rdiffsqrt > eps)
                    for(k=0;k<DIR;k++)
                        rdiff[k]=rdiff[k]/rdiffsqrt;
                for(k=0;k<DIR;k++)
                    d[p][k]=d[p][k]+rdiff[k];
            }

        }
    }
    else
    {
        for(k=0;k<DIR;k++)
            d[p][k]=vel[p][k];
    }
    for(k=0;k<DIR;k++)
        tmpd[k]=d[p][k];
    tmp=normalize(tmpd);
    for(k=0;k<DIR;k++)
        d[p][k]=*(tmp+k);
}

void check_neighbour(int s, double ro, double ra)
{
    int i,j,k;
    double rdiff[DIR],rdiffsqrt,rdiffsq;
    nrep=0; no=0; na=0; nn=0;
    for(j=0;j<N;j++)
    {
        if(j!=s)
        {
            rdiffsq=0.0;
            for(k=0;k<DIR;k++)
            {
                rdiff[k]=r[j][k]-r[s][k];
                rdiff[k]=rdiff[k]-(boxl)*anint(rdiff[k]/(boxl));
                rdiffsq=rdiffsq+rdiff[k]*rdiff[k];
            }
            rdiffsqrt=sqrt(rdiffsq);
            if(rdiffsqrt < r_r)
            {
                listrep[nrep]=j;
                nrep++;
            }
            else if(rdiffsqrt < ro)
            {
                listo[no]=j;
                no++;
            }
            else if(rdiffsqrt < ra)
            {
                lista[na]=j;
                na++;
            }
            else
            {
                listnn[nn]=j;
                nn++;
            }
        }
    }
    find_direction(s,nrep,no,na);
}

double generate_noise()
{
    double eta,thetan;
    eta=0.1*sqrt(delt);
    thetan=2.0 * eta * (ran2(&seed) - 0.5);
    return thetan;
}

void update_posvel(double ro, double ra, int rs)
{
    int i,k;
    double fl,thetan,oldvel[N][DIR],tmpr;

    //thetan=generate_noise();

    for(i=0;i<N;i++)
    {
        thetan=generate_noise();
        check_neighbour(i, ro, ra);
        fl=evaluate_field(i,rs);

        for(k=0;k<DIR;k++)
            oldvel[i][k]=vel[i][k];

        vel[i][0]=vmax*fl*(cos(thetan)*d[i][0]-sin(thetan)*d[i][1]);
        vel[i][1]=vmax*fl*(sin(thetan)*d[i][0]+cos(thetan)*d[i][1]);

        //if(fabs(vel[i][0]) < eps || fabs(vel[i][1]) < eps)
        //   printf("time=%d vel[%d]=(%lf %lf) testf=%lf fl=%lf vmax=%lf d=(%lf %lf) thetan=%lf\n",time,i,vel[i][0],vel[i][1],testf,fl,vmax,d[i][0],d[i][1],thetan);
        for(k=0;k<DIR;k++)
        {
            tmpr=r[i][k]+0.5*(oldvel[i][k]+vel[i][k]);
            r[i][k]=correct_pbc(tmpr);

        }
    }
}

double max(double a, double b)
{
    if(a >= b)
        return a;
    else
        return b;
}

double train(double ro, double ra, int rs)
{
    int i,k;
    double tot_fl0,tot_fl1,Q;
    tot_fl0=0;
    for(i=0;i<N;i++)
        tot_fl0=tot_fl0+evaluate_field(i,rs);
    tot_fl0=tot_fl0/N;

    for(i=0;i<tstar;i++)
        update_posvel(ro, ra, rs);

    tot_fl1=0;

    for(i=0;i<N;i++)
        tot_fl1=tot_fl1+evaluate_field(i,rs);

    tot_fl1=tot_fl1/N;

    Q=max(0,tot_fl0-tot_fl1);
    return Q;
}

void restore_initialize()
{
    int i,j;
    for(i=0;i<N;i++)
    {
        for(j=0;j<DIR;j++)
        {
            r[i][j]=rinit[i][j];
            vel[i][j]=vinit[i][j];
        }
    }

}

void max_min()
{
    int y,j;
    double a,b;
    FILE *fpw;
    fpw=fopen("max_min02.dat","r");
    j=0;
    while(fscanf(fpw,"%d%lf%lf",&y,&a,&b)!=EOF)
    {
        MIN[j]=a; MAX[j]=b;
        j++;
    }
}

int main()
{
    int i,j,k,row,col;
    char outfile1[100],outfile2[100],outfile3[100];
    FILE *fpw;
    sprintf(outfile3,"ttavg_reward_vmax%1.2lft%d.dat",vmax,tstar);
    double delta,gamma,reward;
    double Mx,Mn,Gp;
    Mx=3; Mn=r_r; Gp=0.05;
    int dim;
    dim=(Mx-Mn)/Gp+1;
    double avg_rw[dim+1][dim+1],count[dim+1][dim+1];

    for(i=0;i<=dim;i++)
    {
        for(j=0;j<=dim;j++)
        {
            avg_rw[i][j]=0;
            count[i][j]=0;
        }
    }
    max_min();
    fpr=fopen("field_amp01.dat","r");
    for(i=0;i<1600;i++)
    {
        get_field_amps(i);
        initialize_positions();
        initialize_velocities(i);
        printf("i=%d\n",i);
        for(r_o=Mn;r_o<(Mx+eps+Gp);r_o=r_o+Gp)
        {
            for(r_a=Mn;r_a<(Mx+eps+Gp);r_a=r_a+Gp)
            {
                row=(r_a-Mn)/Gp; col=(r_o-Mn)/Gp;
                restore_initialize();
                reward=train(r_o,r_a,i);
                avg_rw[col][row]=avg_rw[col][row]+reward;
                count[col][row] = count[col][row] + 1;
            }
        }
        fpw=fopen(outfile3,"w");
        fprintf(fpw,"#i=%d r_a r_o\n",i);
        for(k=0;k<dim;k++)
        {
            for(j=0;j<dim;j++)
            {
                if(count[j][k] > 0)
                    fprintf(fpw,"%lf %lf %lf %lf\n",j*Gp+Mn,k*Gp+Mn,avg_rw[j][k]/count[j][k],count[j][k]);
            }
            fprintf(fpw,"\n");
        }
        fclose(fpw);
    }
}
