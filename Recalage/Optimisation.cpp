//---------------------------------------------------------------------------
// Fichier Optimisation.cpp
// Optimisation locale par méthode quasi-Newton
// et optimisation globale par clustering stochastique
//---------------------------------------------------------------------------

#include "stdafx.h"

#include "Optimisation.h"

#define MAXNC 12


#define MAX_nparm 52  
#define MAX_120 1378  


void InitParameterWindow(const CImg<double>& t, const double angle, double transl_ecart, double rot_ecart, int dim, double *p, double* mmax, double* mmin)
{
    p[1] = t(0); p[2] = t(1); p[3] = angle;
	mmin[1] = p[1] - transl_ecart; mmin[2] = p[2] - transl_ecart; mmin[3] = p[3] - rot_ecart;
    mmax[1] = p[1] + transl_ecart; mmax[2] = p[2] + transl_ecart; mmax[3] = p[3] + rot_ecart;

	for (int i=1;i<=dim;i++)
    {
      mmax[i]= (mmax[i]-mmin[i]) / 2;
      mmin[i]= mmin[i]+mmax[i];
    }
}

CImg<double> QuasiNewton( CostFunction& CF, CImg<double> &t, double &angle )
{
    int dim = 3;

    double iter_Max= 1000;
    int iter= 0;
    double fret= 0.0;

	// Détermination des intervalles d'initialisation
	double * p = new double[dim];
	double *mmin = new double[dim+1];
	double *mmax = new double[dim+1];
	InitParameterWindow( t, angle, 10., 10., dim, p, mmax,  mmin );

    double* r = new double[dim+1];
	//for (int i=0;i<dim;i++) r[i] = 0.5;
	memset(r+1,0,dim*sizeof(double));

	CF.set_Start(1);
    double coutref = CF.Evaluate(p);
    
    local1(dim,1E-2,iter_Max,r,&fret,&iter,mmin,mmax,CF);

    for (int i=0;i<dim;i++) p[i] = mmin[i+1] + r[i+1] * mmax[i+1];

    t(0) = p[0];
	t(1) = p[1];
	angle = p[2];

    CImg<double> resume(2);
	resume(0) = fret;
	resume(1) = coutref;

	return resume;
 }

CImg<double> SimulatedAnnealing( CostFunction& CF, CImg<double> &t, double &angle )
{
    int dim = 3;

    double iter_Max= 400;
    int iter= 0;
    double fret= 0.0;

	// Détermination des intervalles d'initialisation
	double * p = new double[dim];
	p[0] = t(0);
	p[1] = t(1);
	p[2] = angle;

	CF.set_Start(0);
    double coutref = CF.Evaluate(p);
    
	SimAnneal sa(&CF,3);

	if ( !sa )
    {
       	cout << "problem initializing SimAnneal object\n";
        exit(1);
    }

	sa.initial(p);	// set initial condition

	sa.melt();		// melt the system

	// make it a bit warmer than melting temperature
	double temp0 = 1.2 * sa.temperature();

	sa.temperature( temp0 );
	
	double temp = sa.anneal( iter_Max );
	sa.optimum( p );

    t(0) = p[0];
	t(1) = p[1];
	angle = p[2];

    CImg<double> resume(2);
	resume(0) = fret;
	resume(1) = coutref;

	delete[] p;

	return resume;
}

CImg<double> StochasticClustering( CostFunction& CF, CImg<double> &t, double &angle, double delta_translation, double delta_rotation )
{
    int dim = 3;

    double iter_Max= 1000;
    int iter= 0;
    double fret= 0.0;

	// Détermination des intervalles d'initialisation
	double * p = new double[dim];
	double *mmin = new double[dim+1];
	double *mmax = new double[dim+1];
	InitParameterWindow( t, angle, delta_translation, delta_translation, dim, p, mmax,  mmin );

    double* r = new double[dim+1];
	//for (int i=0;i<dim;i++) r[i] = 0.5;
	memset(r+1,0,dim*sizeof(double));

	CF.set_Start(1);
    double coutref = CF.Evaluate(p);
    
    double f0[21];
	for (int it= 0;it<21;it++) f0[it]=1;
	int nsampl= 500, nsel= 10, nsig=6, nc, fe;
    TOMB_nx21 x0 = (TOMB_nx21)new double[21*(dim+1)];
	FILE* out = fopen("adat.dat", "w+");
	global(mmin,mmax,dim,nsampl,nsel,out,nsig,x0,&nc,f0,&fe,CF);
	fclose(out);

    for (int i=0;i<dim;i++) p[i] = mmin[i+1] + r[i+1] * mmax[i+1];

    t(0) = (*x0)[1][1];
	t(1) = (*x0)[2][1];
	angle = (*x0)[3][1];

    CImg<double> resume(2);
	resume(0) = fret;
	resume(1) = coutref;

	return resume;
 }

static double   zero=0.0,one=1.0,two=2.0,ten=10.0;
static  double  ax=0.1,op1=0.1,op001=0.001,half=0.5,reps=1.1921e-07,four=4.0,five=5.0,seven=7.0,twelve=12.0;

int global(double *amin,double *amax,int nparm,int n100,int ng0,FILE *fw,int nsig,TOMB_nx21 x0,int *nc,double *f0,int *nfe,CostFunction& CF)
{
             int toto =0;
	     int     i,ng,ns,ncp,n0,n1,im,ig,n,inum,inum1,inum2,nn100,nfe1;
             int     ng10,maxfn,nm,i1,j,jj,icc,in1,iii,ii,l1,icj,iv,it;
             int     ic[101],ic1[21];

	     double    b1,alfa,fm,bb,fc,a,relcon,ff,b;
             double    f[101],f1[21],fcl[101];

             TOMB_nx101 x = (TOMB_nx101)new double[101*(nparm+1)];
             TOMB_nx21 x1 = (TOMB_nx21 )new double[21*(nparm+1)];
             TOMB_nx101 xcl = (TOMB_nx101)new double[101*(nparm+1)];
             TOMB_nx101 r = (TOMB_nx101)new double[101*(nparm+1)];

             double *w = (double *)new double[nparm+1];
             double *y = (double *)new double[nparm+1];
             double *minx = (double *)new double[nparm+1];
             double *maxx = (double *)new double[nparm+1];

             if (nparm <= 0) goto label_460;
             if (nparm  > MAX_nparm) goto label_455;
             for(i=1;i<=nparm;i++) {
                 minx[i] = amin[i];
                 maxx[i] = amax[i];
                 if (minx[i] == maxx[i]) goto label_460;
             }
             b1 = one/(double)nparm;
             if (ng0  < 1) ng0 = 1;
             if (ng0  > 20) ng0 = 20;
             if (n100  < 20) n100 = 20;
             if (n100  > 10000) n100 = 10000;
             if (n100 >= 100) goto label_10;
             nn100 = n100;
             n = 1;
             goto label_MAX_nparm;
  label_10:  nn100 = 100;
             n = n100/100;
             n100 = n*100;
  label_MAX_nparm:  ng10 = 100;
             for(i=1;i<=ng10;i++) {
                 f[i] =9.9e10;
                 ic[i] = 0;
             }
             for(i=1;i<=nparm;i++) {
                 maxx[i]= (maxx[i]-minx[i])/two;
                 minx[i]= minx[i]+maxx[i];
             }
             alfa = .01;
             *nfe = 0;
             ng = 0;
             ns = 0;
             *nc = 0;
             ncp = 1;
             n0 = 0;
             n1 = 0;
             im = 1;
             ig = 0;
             fm = 9.9e10;
             maxfn = 500*nparm;
             relcon = pow(ten,-nsig);

// SAMPLING
  label_20:  n0 = n0+n100;
             nm = n0-1;
             ng = ng+ng0;
             ns = ns+1;
             if (ns*ng0 > 100) goto label_465;
             b = pow(one-pow(alfa,one/(double)nm),b1);
             bb = 0.1*b;
             for(i1=1;i1<=n;i1++) {
                 urdmn(r,nparm);
                 for(j=1;j<=nn100;j++) {
                     for(i=1;i<=nparm;i++) {
                         y[i] = two* (*r)[i][j]-one;
                     }
                     fun(y,&fc,nparm,minx,maxx,CF);
                     if (fc >= fm) continue;
                     f[im] = fc;
                     for(i=1;i<=nparm;i++) {
                         (*x)[i][im] = y[i];
                     }
                     if (im <= ng &&  ic[im]  > 0) ig = ig-1;
                     ic[im] = 0;
                     im = 1;
                     fm = f[1];
                     for(i=2;i<=ng10;i++) {
                         if (f[i]  < fm) continue;
                         im = i;
                         fm = f[i];
                     }
                 }
             }
             *nfe = *nfe+n100;
             if (fw != NULL) //!!!!!
             fprintf(fw,"\n %5d  function evaluations used for sampling",n100);
             //printf("\n %5d  function evaluations used for sampling\n",n100);

// SORTING
             inum = ng10-1;
             for(i=1;i<=inum;i++) {
                 im = i;
                 fm = f[i];
                 inum1 = i+1;
                 for(j=inum1;j<=ng10;j++) {
                     if (f[j] >= fm) continue;
                     im = j;
                     fm = f[j];
                 }
                 if (im <= i) continue;
                 a = fm;
                 for(j=1;j<=nparm;j++) {
                     y[j] = (*x)[j][im];
                 }
                 if (i  > ng || im <= ng) goto label_55;
                 if (ic[ng] == 0 &&  ic[im]  > 0) ig = ig+1;
                 if (ic[ng]  > 0 &&  ic[im] == 0) ig = ig-1;
   label_55:     icc = ic[im];
                 inum1 = im-i;
                 for(j=1;j<=inum1;j++) {
                     inum2 = im-j;
                     f[inum2+1] = f[inum2];
                     ic[inum2+1]= ic[inum2];
                     for(jj=1;jj<=nparm;jj++) {
                         (*x)[jj][inum2+1] = (*x)[jj][inum2];
                     }
                 }
                 f[i] =a;
                 for(j=1;j<=nparm;j++) {
                     (*x)[j][i] = y[j];
                 }
                 ic[i] = icc;
             }
             if (*nc <= 0) goto label_200;

// CLUSTERING TOX*
             for(iii=1;iii<=*nc;iii++) {
                 i = 1;
                 in1 = i;
                 fcl[i]= f0[iii];
                 for(j=1;j<=nparm;j++) {
                     (*xcl)[j][i] = (*x0)[j][iii];
                 }
                 for(j=1;j<=ng;j++) {
                     if (ic[j] != iii) continue;
                     in1 = in1+1;
                     for(ii=1;ii<=nparm;ii++) {
                         (*xcl)[ii][in1]= (*x)[ii][j];
                     }
                 }
  label_95:      for(j=1;j<=ng;j++) {
                     if (ic[j] != 0)  continue;
                     if (fcl[i] >= f[j]) continue;
                     for(l1=1;l1<=nparm;l1++) {
                         w[l1] = fabs((*xcl)[l1][i]-(*x)[l1][j]);
                     }
                     a =zero;
                     for(l1=1;l1<=nparm;l1++) {
                         if (w[l1]  > a) a=w[l1];
                     }
                     if (a >= b) continue;
                     if (fw != NULL)//!!!!!
                     fprintf(fw,"\n sample point added to the cluster no. %2d",iii);
                     for(ii=1;ii<=nparm;ii++) {
                         w[ii]=(*x)[ii][j]*maxx[ii]+minx[ii];
                     }
                     if (fw != NULL)//!!!!!
                     {
                     fprintf(fw,"\n %14.8g\n  ",f[j]);;
                     for(ii=1;ii<=nparm;ii++) {
                        fprintf(fw," %14.8g",w[ii]);;
                        if(ii==5 || ii==10 || ii==MAX_nparm) fprintf(fw,"\n  ");
                     }
                     }
                     ig = ig+1;
                     if (ig >= ng) goto label_395;
                     in1= in1+1;
                     fcl[in1] = f[j];
                     for(ii=1;ii<=nparm;ii++) {
                         (*xcl)[ii][in1]=(*x)[ii][j];
                     }
                     ic[j] = iii;
                 }
                 i = i+1;
                 if (i <= in1) goto label_95;
             }
             if (n1 <= 0) goto label_200;

// CLUSTERING TOX1
             for(iii=1;iii<=n1;iii++) {
                 i = 1;
                 in1 = i;
                 fcl[i]= f1[iii];
                 for(j=1;j<=nparm;j++) {
                     (*xcl)[j][i] = (*x1)[j][iii];
                 }
  label_155:     for(j=1;j<=ng;j++) {
                     if (ic[j] != 0) continue;
                     if (fcl[i] >= f[j]) continue;
                     for(l1=1;l1<=nparm;l1++) {
                         w[l1] = fabs((*xcl)[l1][i]-(*x)[l1][j]);
                     }
                     a =zero;
                     for(l1=1;l1<=nparm;l1++) {
                         if (w[l1]  > a) a=w[l1];
                     }
                     if (a >= b) continue;
                     if (fw != NULL)//!!!!!
                     fprintf(fw,"\n sample point added to the cluster no. %d",ic1[iii]);
		     //printf("\n sample point added to the cluster no. %d\n",ic1[iii]);
                     for(ii=1;ii<=nparm;ii++) {
                         w[ii] = (*x)[ii][j]*maxx[ii]+minx[ii];
                     }
                     if (fw != NULL)//!!!!!
                     {
                     fprintf(fw,"\n %14.8g\n  ",f[j]);;
                     for(ii=1;ii<=nparm;ii++) {
                         fprintf(fw," %14.8g",w[ii]);
                         if(ii==5 || ii==10 || ii==MAX_nparm) fprintf(fw,"\n  ");
                     }
                     }
                     ig = ig+1;
                     if (ig >= ng) goto label_395;
                     in1= in1+1;
                     fcl[in1] = f[j];
                     for(ii=1;ii<=nparm;ii++) {
                         (*xcl)[ii][in1] = (*x)[ii][j];
                     }
                     ic[j]=ic1[iii];
                 }
                 i = i+1;
                 if (i <= in1) goto label_155;
             }

// LOCALSEARCH
  label_200: it = 0;
             for(i1=1;i1<=ng;i1++) {
                 if (ic[i1] != 0) continue;
                 for(i=1;i<=nparm;i++) {
                     y[i] = (*x)[i][i1];
                 }
                 ff = f[i1];
		 //printf("   >> Local Search\n");
                 local1(nparm,relcon,maxfn,y,&ff,&nfe1,minx,maxx,CF);
                 //local2(nparm,&relcon,&maxfn,y,&ff,&nfe1,minx,maxx,func,param);
                 if (*nc <= 0) goto label_290;
                 for(iv=1;iv<=*nc;iv++) {
                     for(l1=1;l1<=nparm;l1++) {
                         w[l1] = fabs((*x0)[l1][iv]-y[l1]);
                     }
                     a =zero;
                     for(l1=1;l1<=nparm;l1++) {
                         if (w[l1]  > a) a=w[l1];
                     }
                     if (a  < bb) goto label_255;
                 }
                 goto label_290;

// NEWSEED-POINT
  label_255:     n1 = n1+1;
                 if (fw != NULL)//!!!!!
                 fprintf(fw,"\n new seed point added to the cluster no %d nfev= %d",iv,nfe1);
                 for(ii=1;ii<=nparm;ii++) {
                     w[ii] = (*x)[ii][i1]*maxx[ii]+minx[ii];
                 }
                 if (fw != NULL)//!!!!!
                 {
                 fprintf(fw,"\n %14.8g\n  ",ff);
                 for(ii=1;ii<=nparm;ii++) {
                      fprintf(fw," %14.8g",w[ii]);
                      if(ii==5 || ii==10 || ii==MAX_nparm) fprintf(fw,"\n  ");
                 }
                 }
                 if (ff >= f0[iv]) goto label_280;
                 if (fw != NULL)//!!!!!
                 fprintf(fw," *** improvement on the local minimum no. %d : %14.8f for %14.8f",iv,f0[iv],ff);
                 for(ii=1;ii<=nparm;ii++) {
                     w[ii] = y[ii]*maxx[ii]+minx[ii];
                 }
                 if (fw != NULL)//!!!!!
                 {
                 fprintf(fw,"\n %14.8g\n  ",ff);
                 for(ii=1;ii<=nparm;ii++){;
                         if (fw != NULL)//!!!!!
                         fprintf(fw," %14.8g",w[ii]);
                     if(ii==5 || ii==10 || ii==MAX_nparm) fprintf(fw,"\n  ");
                 }
                 }
                 f0[iv]= ff;
                 for(ii=1;ii<=nparm;ii++) {
                     (*x0)[ii][iv] =y[ii];
                 }
  label_280:     if (n1  > 20) goto label_470;
                 for(ii=1;ii<=nparm;ii++) {
                     (*x1)[ii][n1] =(*x)[ii][i1];
                     (*xcl)[ii][1] =(*x)[ii][i1];
                 }
                 f1[n1]= f[i1];
                 fcl[1]= f[i1];
                 ic1[n1] = iv;
                 icj = iv;
                 goto label_305;

// NEW LOCAL MINIMUM
  label_290:     *nc = *nc+1;
                 ncp = ncp+1;
                 if (fw != NULL)//!!!!!
                 fprintf(fw, "\n *** the local minimum no. %d : %14.8g   nfev=%d",*nc,ff,nfe1);
                 //printf("\n *** the local minimum no. %d : %14.8g   nfev=%d",*nc,ff,nfe1);
                 for(ii=1;ii<=nparm;ii++) {
                     w[ii] = y[ii]*maxx[ii]+minx[ii];
                 }
                 if (fw != NULL)//!!!!!
                 {
                 fprintf(fw,"\n %14.8g\n  ",ff);
                 //printf("\n  X = ");
                 for(ii=1;ii<=nparm;ii++) {
                         if (fw != NULL)//!!!!!
                         fprintf(fw," %14.8g",w[ii]);
                         //printf(" %14.8g",w[ii]);
                         if(ii==5 || ii==10 || ii==MAX_nparm) fprintf(fw,"\n  ");
				 }
                  }

				 //printf("\n");
                 for(ii=1;ii<=nparm;ii++) {
                     (*x0)[ii][*nc] =y[ii];
                     (*xcl)[ii][1] =y[ii];
                 }
                 fcl[1]= ff;
                 f0[*nc]= ff;
				 if (*nc >= MAXNC) goto label_475;
                 it = 1;
                 icj = *nc;

// CLUSTERING TO THE NEW POINT
  label_305:     *nfe = *nfe+nfe1;
                 ic[i1]= icj;
                 ig = ig+1;
                 if (ig >= ng) goto label_390;
                 i = 1;
                 in1 = i;
  label_310:     for(j=1;j<=ng;j++) {
                     if (ic[j] != 0) continue;
                     if (fcl[i] >= f[j]) continue;
                     for(l1=1;l1<=nparm;l1++) {
                         w[l1] = fabs((*xcl)[l1][i]-(*x)[l1][j]);
                     }
                     a = zero;
                     for(l1=1;l1<=nparm;l1++) {
                         if (w[l1]  > a) a=w[l1];
                     }
                     if (a >= b) continue;
                     in1 = in1+1;
                     for(ii=1;ii<=nparm;ii++) {
                         (*xcl)[ii][in1]= (*x)[ii][j];
                     }
                     fcl[in1] = f[j];
                     ic[j] = icj;
                     if (fw != NULL)//!!!!!
                     fprintf(fw,"\n sample point added to the cluster no. %d",icj);
                     for(ii=1;ii<=nparm;ii++) {
                         w[ii] = (*x)[ii][j]*maxx[ii]+minx[ii];
                     }
                         if (fw != NULL)//!!!!!
                     {
                         fprintf(fw,"\n %14.8g\n  ",f[j]);
                     for(ii=1;ii<=nparm;ii++) {
                         if (fw != NULL)//!!!!!
                         fprintf(fw," %14.8g",w[ii]);
                         if(ii==5 || ii==10 || ii==MAX_nparm) fprintf(fw,"\n  ");
                     }
                     }
                     ig = ig+1;
                     if (ig >= ng) goto label_390;
                 }
                 i = i+1;
                 if (i  < in1) goto label_310;
             }
  label_390: if (it != 0) goto label_20;

// PRINT RESULTS
  label_395: if (fw != NULL)//!!!!!
             fprintf(fw, "\n\n\n\n local minima found: \n\n");
             if (*nc <= 1) goto label_430;
             inum = *nc-1;
             for(i=1;i<=inum;i++) {
                 im = i;
                 fm = f0[i];
                 inum1 = i+1;
                 for(j=inum1;j<=*nc;j++) {
                     if (f0[j] >= fm) continue;
                     im = j;
                     fm = f0[j];
                 }
                 if (im <= i) continue;
                 a = fm;
                 for(j=1;j<=nparm;j++) {
                     y[j] = (*x0)[j][im];
                 }
                 inum1 = im-i;
                 for(j=1;j<=inum1;j++) {
                     inum2 = im-j;
                     f0[inum2+1]= f0[inum2];
                     for(jj=1;jj<=nparm;jj++) {
                         (*x0)[jj][inum2+1] =(*x0)[jj][inum2];
                     }
                 }
                 f0[i] = a;
                 for(j=1;j<=nparm;j++) {
                     (*x0)[j][i] = y[j];
                 }
             }
  label_430: if (*nc <= 0) goto label_445;
             for(i=1;i<=*nc;i++) {
                 for(ii=1;ii<=nparm;ii++) {
                     (*x0)[ii][i] = (*x0)[ii][i]*maxx[ii]+minx[ii];
                 }
                 if (fw != NULL)//!!!!!
                 {
                 fprintf(fw,"\n F(x)= %14.8g\n  ",f0[i]);
                 fprintf(fw,"\n   X =");
                 for(ii=1;ii<=nparm;ii++) {
                         fprintf(fw," %14.8g",(*x0)[ii][i]);
                     if(ii==5 || ii==10 || ii==MAX_nparm) fprintf(fw,"\n  ");
                 }
                 }
             }
  label_445: delete[] maxx;
             delete[] minx;
             delete[] y;
             delete[] w;
             delete[] r;
             delete[] xcl;
             delete[] x1;
             delete[] x;
             if (fw != NULL)//!!!!!
             fprintf(fw,"\n\n normal termination after %d function evaluations \n\n",*nfe);
			 return 0;
  label_455: if (fw != NULL)//!!!!!
             fprintf(fw," ***   too many parameters, abnormal termination");
             exit(3);
  label_460: if (fw != NULL)//!!!!!
             fprintf(fw,"  ***   data error");
             exit(3);
  label_465: if (fw != NULL)//!!!!!
             fprintf(fw," ***   too many sampling");
             goto label_395;
  label_470:  if (fw != NULL)//!!!!!
             fprintf(fw," ***   too many new seed points");
             goto label_395;
  label_475: if (fw != NULL)//!!!!!
             fprintf(fw," ***   too many clusters");
             goto label_395;
}


void fun(double *r,double *f,int nparm,double *minx,double *maxx,CostFunction &CF)
{
             double x[MAX_nparm];
             //memset(x,0,MAX_nparm*sizeof(double));

             for(int i=1; i<=nparm; i++) x[i] = maxx[i] * r[i] + minx[i];
			 *f = CF.Evaluate(x);
}


void urdmn(TOMB_nx101 t,int nparm)
{
        register int i,j;
        int            g;
        double         k;

        RaNdOmIzE;
        for(i=1;i<=100;i++)
           for(j=1;j<=nparm;j++){
              g=rand();
              k=(double)g/(double)RAND_MAX;
              (*t)[j][i]=k;
           }
}


void local1(int n,double eps,int maxfn,double *x,double *f,int *nfev,double *minx,double *maxx,CostFunction &CF)
{
// SPECIFICATIONS FOR LOCAL VARIABLES

               int      ig,igg,is,idiff,ir,ij,i,j,nm1,jj,jp1,l,kj,
                          k,link,itn,ii,im1,jnt,np1,jb,nj,ier;

               double     hh,hjj,v,df,relx,gs0,diff,aeps,alpha,ff,
                          tot,f1,f2,z,gys,dgs,sig,zz,hhh,ghh,iopt;

	       static double h[MAX_120];

               double *g, *w;
               g = (double *)new double[n+1];
               w = (double *)new double[3*(n+1)];

//                        INITIALIZATION
//                        FIRST EXECUTABLE STATEMENT
                iopt = 0;
//            IOPT     - OPTIONS SELECTOR. (INPUT)
//                  IOPT = 0 CAUSES LOCAL     TO INITIALIZE THE
//                    HESSIAN MATRIX H TO     THE IDENTITY MATRIX.
//                  IOPT = 1 INDICATES THAT H HAS     BEEN INITIALIZED
//                    BY THE USER     TO A POSITIVE DEFINITE MATRIX.
//                  IOPT = 2 CAUSES LOCAL     TO COMPUTE THE DIAGONAL
//                    VALUES OF THE HESSIAN MATRIX AND SET H TO
//                    A DIAGONAL MATRIX CONTAINING THESE VALUES.
//                  IOPT = 3 CAUSES LOCAL     TO COMPUTE AN ESTIMATE
//                    OF THE HESSIAN IN H.
                ier = 0;
                hh = sqrt(reps);
                ig = n;
                igg = n+n;
                is = igg;
                idiff = 1;
                ir = n;
                w[1] = -one;
                w[2] = zero;
                w[3] = zero;

// EVALUATE FUNCTION AT     STARTING POINT
                for(i=1;i<=n;i++)
                     g[i] = x[i];
                fun(g,f,n,minx,maxx,CF);
                *nfev = 1;
                if (iopt == 1) goto label_45;

//SET OFF-DIAGONAL ELEMENTS OF H TO 0.0
                if (n == 1) goto label_25;
                ij = 2;
                for(i=2;i<=n;i++) {
                     for(j=2;j<=i;j++) {
                          h[ij] = zero;
                          ij = ij+1;
                     }
                     ij = ij+1;
                }
                if (iopt != 0) goto label_25;

// SET DIAGONAL     ELEMENTS OF H TO one
                ij = 0;
                for(i=1;i<=n;i++) {
                     ij = ij+i;
                     h[ij] = one;
                }
                goto label_80;

// GET DIAGONAL     ELEMENTS OF HESSIAN
  label_25:     im1 = 1;
                nm1 = 1;
                np1 = n+1;
                for(i=2;i<=np1;i++) {
                     hhh = hh*MaX(fabs(x[im1]),ax);
                     g[im1] = x[im1]+hhh;
                     fun (g,&f2,n,minx,maxx,CF);
                     g[im1] = g[im1]+hhh;
                     fun (g,&ff,n,minx,maxx,CF);
                     h[nm1] = (ff-f2+*f-f2)/(hhh*hhh);
                     g[im1] = x[im1];
                     im1 = i;
                     nm1 = i+nm1;
                }
                *nfev = *nfev+n+n;
                if ((iopt != 3) || (n == 1)) goto label_45;

// GET THE REST     OF THE HESSIAN
                jj = 1;
                ii = 2;
                for(i=2;i<=n;i++) {
                     ghh = hh*MaX(fabs(x[i]),ax);
                     g[i] = x[i]+ghh;
                     fun (g,&f2,n,minx,maxx,CF);
                     for(j=1;j<=jj;j++) {
                          hhh = hh*MaX(fabs(x[j]),ax);
                          g[j] = x[j]+hhh;
                          fun (g,&ff,n,minx,maxx,CF);
                          g[i] = x[i];
                          fun (g,&f1,n,minx,maxx,CF);

// H(II) = (FF-F1-F2+F)*SQREPS
                          h[ii] = (ff-f1-f2+*f)/(hhh*ghh);
                          ii = ii+1;
                          g[j] = x[j];
                     }
                     jj = jj+1;
                     ii = ii+1;
                }
                *nfev = *nfev+((n*n-n)/2);

// FACTOR H TO L*D*L-TRANSPOSE
  label_45:     ir = n;
                if (n  > 1) goto label_50;
                if (h[1]  > zero)     goto label_80;
                h[1] = zero;
                ir = 0;
                goto label_75;
  label_50:     nm1 = n-1;
                jj = 0;
                for(j=1;j<=n;j++) {
                     jp1 = j+1;
                     jj = jj+j;
                     hjj = h[jj];
                     if (hjj  > zero) goto label_55;
                     h[jj] = zero;
                     ir = ir-1;
                     continue;
  label_55:          if (j == n) continue;
                     ij = jj;
                     l = 0;
                     for(i=jp1;i<=n;i++) {
                          l = l+1;
                          ij = ij+i-1;
                          v = h[ij]/hjj;
                          kj = ij;
                          for(k=i;k<=n;k++) {
                               h[kj+l] = h[kj+l]-h[kj]*v;
                               kj = kj+k;
                          }
                          h[ij] = v;
                     }
                }
  label_75:     if (ir == n) goto label_80;
                ier = 129;
                goto label_9000;
  label_80:     itn = 0;
                df = -one;

// EVALUATE GRADIENT W(IG+I),I=1,...,N
  label_85:     link = 1;
                goto label_260;
  label_90:     if (*nfev >= maxfn) goto label_235;

// BEGIN ITERATION LOOP
                itn = itn+1;
                for(i=1;i<=n;i++)
                     w[i] = -w[ig+i];

// DETERMINE SEARCH DIRECTION W
// BY     SOLVING     H*W = -G WHERE
// H = L*D*L-TRANSPOSE
                if (ir  < n) goto label_125;

// N .EQ. 1
                g[1] = w[1];
                if (n  > 1) goto label_100;
                w[1] = w[1]/h[1];
                goto label_125;

// N .GT. 1
  label_100:    ii = 1;

// SOLVE L*W = -G
                for(i=2;i<=n;i++) {
                     ij = ii;
                     ii = ii+i;
                     v = w[i];
                     im1 = i-1;
                     for(j=1;j<=im1;j++) {
                          ij = ij+1;
                          v = v-h[ij]*w[j];
                     }
                     g[i] = v;
                     w[i] = v;
                }

// SOLVE (D*LT)*Z = W WHERE
// LT = L-TRANSPOSE
                w[n] = w[n]/h[ii];
                jj = ii;
                nm1 = n-1;
                for(nj=1;nj<=nm1;nj++) {

// J = N-1,N-2,...,1
                     j = n-nj;
                     jp1 = j+1;
                     jj = jj-jp1;
                     v = w[j]/h[jj];
                     ij = jj;
                     for(i=jp1;i<=n;i++) {
                          ij = ij+i-1;
                          v = v-h[ij]*w[i];
                     }
                     w[j] = v;
                }

// DETERMINE STEP LENGTH ALPHA
  label_125:    relx = zero;
                gs0 = zero;
                for(i=1;i<=n;i++) {
                     w[is+i] = w[i];
                     diff = fabs(w[i])/MaX(fabs(x[i]),ax);
                     relx = MaX(relx,diff);
                     gs0 = gs0+w[ig+i]*w[i];
                }
                if (relx == zero)     goto label_230;
                aeps = eps/relx;
                ier = 130;
                if (gs0 >= zero) goto label_230;
                if (df == zero) goto label_230;
                ier = 0;
                alpha = (-df-df)/gs0;
                if (alpha <= zero) alpha = one;
                alpha = MiN(alpha,one);
                if (idiff == 2) alpha = MaX(op1,alpha);
                ff = *f;
                tot = zero;
                jnt = 0;

// SEARCH ALONG     X+ALPHA*W
  label_135:    if (*nfev >= maxfn) goto label_235;
                for(i=1;i<=n;i++)
                     w[i] = x[i]+alpha*w[is+i];
                fun (w,&f1,n,minx,maxx,CF);
                *nfev = *nfev+1;
                if (f1 >= *f) goto label_165;
                f2 = *f;
                tot = tot+alpha;
  label_145:    ier = 0;
                *f = f1;
                for(i=1;i<=n;i++)
                    x[i] = w[i];
                if (jnt-1 <0)      goto label_155;
                else if (jnt-1==0) goto label_185;
                else               goto label_190;
  label_155:    if (*nfev >= maxfn) goto label_235;
                for(i=1;i<=n;i++)
                     w[i] = x[i]+alpha*w[is+i];
                fun (w,&f1,n,minx,maxx,CF);
                *nfev = *nfev+1;
                if (f1 >= *f) goto label_190;
                if ((f1+f2 >= *f+*f) && (seven*f1+five*f2  > twelve*(*f))) jnt = 2;
                tot = tot+alpha;
                alpha = alpha+alpha;
                goto label_145;
  label_165:    if ((*f == ff) && (idiff == 2) && (relx  > eps)) ier = 130;
                if (alpha  < aeps) goto label_230;
                if (*nfev >= maxfn) goto label_235;
                alpha = half*alpha;
                for(i=1;i<=n;i++)
                     w[i] = x[i]+alpha*w[is+i];
                fun (w,&f2,n,minx,maxx,CF);
                *nfev = *nfev+1;
                if (f2 >= *f) goto label_180;
                tot = tot+alpha;
                ier = 0;
                *f = f2;
                for(i=1;i<=n;i++)
                     x[i] = w[i];
                goto label_185;
  label_180:    z = op1;
                if (f1+*f  > f2+f2) z = one+half*(*f-f1)/(*f+f1-f2-f2);
                z = MaX(op1,z);
                alpha = z*alpha;
                jnt = 1;
                goto label_135;
  label_185:    if (tot  < aeps) goto label_230;
  label_190:    alpha = tot;

// SAVE     OLD GRADIENT
                for(i=1;i<=n;i++)
                     w[i] = w[ig+i];

// EVALUATE GRADIENT W(IG+I), I=1,...,N
                link = 2;
                goto label_260;
  label_200:    if(*nfev >= maxfn) goto label_235;
                gys = zero;
                for(i=1;i<=n;i++) {
                     gys = gys+w[ig+i]*w[is+i];
                     w[igg+i] = w[i];
                }
                df = ff-*f;
                dgs = gys-gs0;
                if (dgs <= zero) goto label_90;
                if (dgs+alpha*gs0  > zero) goto label_215;

// UPDATE HESSIAN H USING
// COMPLEMENTARY DFP FORMULA
                sig = one/gs0;
                ir = -ir;
                update (h,n,w,&sig,g,&ir,0,zero);
                for(i= 1;i<=n;i++)
                     g[i] = w[ig+i]-w[igg+i];
                sig = one/(alpha*dgs);
                ir = -ir;
                                update (h,n,g,&sig,w,&ir,0,reps);
                goto label_90;

// UPDATE HESSIAN USING
// DFP FORMULA
  label_215:    zz = alpha/(dgs-alpha*gs0);
                sig = -zz;
                                update (h,n,w,&sig,g,&ir,0,reps);
                z = dgs*zz-one;
                for(i=1;i<=n;i++)
                     g[i] = w[ig+i]+z*w[igg+i];
                sig = one/(zz*dgs*dgs);
                                update (h,n,g,&sig,w,&ir,0,zero);
                goto label_90;

// MaxFN FUNCTION EVALUATIONS
  label_230:    if (idiff == 2) goto label_235;

// CHANGE TO CENTRAL DIFFERENCES
                idiff = 2;
                goto label_85;
  label_235:    if ((relx  > eps) && (ier == 0)) goto label_85;

// COMPUTE H = L*D*L-TRANSPOSE AND
// OUTPUT
                if (n == 1) goto label_9000;
                np1 = n+1;
                nm1 = n-1;
                jj = (n*(np1))/2;
                for(jb=1;jb<=nm1;jb++) {
                     jp1 = np1-jb;
                     jj = jj-jp1;
                     hjj = h[jj];
                     ij = jj;
                     l = 0;
                     for(i=jp1;i<=n;i++) {
                          l = l+1;
                          ij = ij+i-1;
                          v = h[ij]*hjj;
                          kj = ij;
                          for(k=i;k<=n;k++) {
                               h[kj+l] = h[kj+l]+h[kj]*v;
                               kj = kj+k;
                          }
                          h[ij] = v;
                     }
                     hjj = h[jj];
                }
                goto label_9000;

// EVALUATE GRADIENT
  label_260:    if (idiff == 2) goto label_270;

// FORWARD DIFFERENCES
// GRADIENT =     W(IG+I), I=1,...,N
                for(i=1;i<=n;i++) {
                     z = hh*MaX(fabs(x[i]),ax);
                     zz = x[i];
                     x[i] = zz+z;
                     fun(x,&f1,n,minx,maxx,CF);
                     w[ig+i] = (f1-*f)/z;
                     x[i] = zz;
                }
                *nfev = *nfev+n;
                if (link==1) goto label_90;
                 else if (link==2) goto label_200;

// CENTRAL DIFFERENCES
// GRADIENT =     W(IG+I), I=1,...,N
  label_270:    for(i=1;i<=n;i++) {
                     z = hh*MaX(fabs(x[i]),ax);
                     zz = x[i];
                     x[i] = zz+z;
                     fun (x,&f1,n,minx,maxx,CF);
                     x[i] = zz-z;
                     fun(x,&f2,n,minx,maxx,CF);
                     w[ig+i] = (f1-f2)/(z+z);
                     x[i] = zz;
                }
                *nfev = *nfev+n;
                if (link==1) goto label_90;
                 else if (link==2) goto label_200;
// RETURN
label_9000:     delete[] g;
                delete[] w;
                return;
}



void  update(double *a,int n,double *z,double *sig,double *w,int *irp,int mk,double eps)
{
                int       j,jj,ij,jp1,i,ii,mm;
                double      ti,v,tim,al,r,b,gm,y;

// UPDATE FACTORS GIVEN     IN A
// SIG*Z*Z-TRANSPOSE IS ADDED
// FIRST EXECUTABLE STATEMENT
                if (n  > 1) goto label_5;

// N .EQ. 1
                a[1] = a[1]+*sig*z[1]*z[1];
                                *irp = 1;
                if (a[1]  > zero)     goto label_9005;
                a[1] = zero;
                                *irp = 0;
                goto label_9005;

// N .GT. 1
  label_5:      if (*sig  > zero) goto label_65;
                if ((*sig == zero) || (*irp == 0)) goto label_9005;
                                ti = one / *sig;
                jj = 0;
                                if (mk == 0) goto label_MAX_nparm;

// L*W = Z ON INPUT
                                for (j=1; j <= n ;j++) {
                     jj = jj+j;
                     if (a[jj] != zero) ti = ti+(w[j]*w[j])/a[jj];
                                }
                goto label_40;

// SOLVE L*W = Z
  label_MAX_nparm:     for(j=1;j<=(n);j++)
                     w[j] = z[j];
                                for(j=1;j<=(n);j++) {
                     jj = jj+j;
                     v = w[j];
                     if (a[jj]  > zero) goto label_25;
                     w[j] = zero;
                     goto label_40;
  label_25:          ti = ti+(v*v)/a[jj];
                     if (j == n) goto label_40;
                     ij = jj;
                     jp1 = j+1;
                     for(i=jp1;i<=n;i++) {
                          ij = ij+i-1;
                          w[i] = w[i]-v*a[ij];
                     }
                }

// SET     TI, TIM     AND W
  label_40:     if (*irp <= 0) goto label_45;
                if (ti  > zero) goto label_50;
                if (mk-1 <= 0) goto label_65;
                   else        goto label_55;
  label_45:     ti = zero;
                                *irp = -*irp-1;
                goto label_55;
  label_50:     ti = eps / *sig;
                                if (eps == zero) *irp = *irp-1;
  label_55:     tim = ti;
                ii = jj;
                i = n;
                for(j=1;j<=n;j++) {
                     if (a[ii] != zero) tim = ti-(w[i]*w[i])/a[ii];
                     w[i] = ti;
                     ti = tim;
                     ii = ii-i;
                     i = i-1;
                }
                mm = 1;
                goto label_70;
  label_65:     mm = 0;
                tim = one / *sig;
  label_70:     jj = 0;

// UPDATE A
                for(j=1;j<=(n);j++) {
                     jj = jj+j;
                     ij = jj;
                     jp1 = j+1;

// UPDATE A(J,J)
                     v = z[j];
                     if (a[jj]  > zero) goto label_85;

// A(J,J) .EQ. zero
                         if ((*irp  > 0) || (*sig  < zero) || (v == zero)) goto label_80;
                                         *irp = 1-*irp;
                     a[jj] = (v*v)/tim;
                     if (j == n) goto label_9005;
                     for(i=jp1;i<=(n);i++) {
                          ij = ij+i-1;
                          a[ij] = z[i]/v;
                     }
                     goto label_9005;
  label_80:          ti = tim;
                     goto label_115;

// A(J,J) .GT. zero
  label_85:          al = v/a[jj];
                     ti = w[j];
                     if (mm == 0) ti = tim+v*al;
                     r = ti/tim;
                     a[jj] = r*a[jj];
                     if (r == zero)     goto label_115;
                     if (j == n) goto label_115;

// UPDATE REMAINDER OF COLUMN J
                     b = al/ti;
                     if (r  > four)     goto label_95;
                     for(i=jp1;i<=(n);i++) {
                          ij = ij+i-1;
                          z[i] = z[i]-v*a[ij];
                          a[ij] = a[ij]+b*z[i];
                     }
                     goto label_105;
  label_95:          gm = tim/ti;
                     for(i=jp1;i<=(n);i++) {
                          ij = ij+i-1;
                          y = a[ij];
                          a[ij] = b*z[i]+y*gm;
                          z[i] = z[i]-v*y;
                     }
  label_105:         tim = ti;
                }
  label_115:    if (*irp  < 0) *irp = -*irp;
label_9005:     return;
}

