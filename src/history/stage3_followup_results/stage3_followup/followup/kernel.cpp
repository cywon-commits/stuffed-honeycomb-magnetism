#include <cstdint>
#include <cmath>
#include <complex>
#include <omp.h>
static int C[32][32];
extern "C" {
void init(){for(int n=0;n<32;n++){C[n][0]=1;for(int k=1;k<=n;k++)C[n][k]=(k==n?1:C[n-1][k-1]+C[n-1][k]);}}
int rankm(uint32_t m){int r=0,k=1;while(m){int p=__builtin_ctz(m);r+=C[p][k++];m&=m-1;}return r;}
int orbits(int n,int up,int dim,int ng,const uint32_t* lo,const uint32_t* hi,int split,int* rid,uint8_t* shift,uint32_t* reps,uint8_t* lens){
 uint32_t m=(1u<<up)-1;int nr=0;uint32_t lowmask=(1u<<split)-1;int nl=1<<split,nh=1<<(n-split);
 for(int r=0;r<dim;r++){
 if(rid[r]<0){reps[nr]=m;int len=0;for(int g=0;g<ng;g++){uint32_t t=lo[g*nl+(m&lowmask)]|hi[g*nh+(m>>split)];int rt=rankm(t);if(rid[rt]<0){rid[rt]=nr;shift[rt]=g;len++;}}lens[nr++]=len;}
 if(r+1<dim){uint32_t c=m&-m,z=m+c;m=(((z^m)>>2)/c)|z;}}
 return nr;}
void mv(int nb,int nbond,const uint32_t* masks,const uint8_t* lens,const double* diag,const uint32_t* flips,const double* js,const int* rid,const uint8_t* shift,const int* rawbasis,const uint8_t* pop,const int* rl,const int* rh,int split,int nhi,const std::complex<double>* phase,const std::complex<double>* x,std::complex<double>* y,int threads){
 uint32_t lm=(1u<<split)-1;
 #pragma omp parallel for num_threads(threads) schedule(static)
 for(int i=0;i<nb;i++){uint32_t m=masks[i];std::complex<double> v=diag[i]*x[i];for(int b=0;b<nbond;b++){if(js[b]==0.0)continue;uint32_t bits=m&flips[b];if(bits && bits!=flips[b]){uint32_t t=m^flips[b];int l=t&lm;int r=rl[l]+rh[pop[l]*nhi+(t>>split)];int j=rawbasis[rid[r]];if(j>=0)v+=0.5*js[b]*(lens[i]==lens[j]?1.0:sqrt(double(lens[i])/lens[j]))*phase[shift[r]]*x[j];}}y[i]=v;}}
}
