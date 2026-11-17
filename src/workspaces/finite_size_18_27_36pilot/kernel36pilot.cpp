#include <cstdint>
#include <complex>
#include <cmath>
#include <omp.h>
using Z=std::complex<double>;static uint64_t C[40][40];
uint64_t rankm(uint64_t x){uint64_t r=0;int k=1;while(x){int p=__builtin_ctzll(x);r+=C[p][k++];x&=x-1;}return r;}
extern "C" {
void init36(){for(int n=0;n<40;n++){C[n][0]=1;for(int k=1;k<=n;k++)C[n][k]=k==n?1:C[n-1][k-1]+C[n-1][k];}}
int orbits36(int n,int up,uint64_t dim,int ng,const int*perms,int*rid,uint8_t*shift,uint64_t*reps,uint8_t*lens){uint64_t mask=(1ULL<<up)-1;int nr=0;for(uint64_t r=0;r<dim;r++){
 if(rid[r]<0){reps[nr]=mask;int len=0;for(int g=0;g<ng;g++){uint64_t t=0,m=mask;while(m){int p=__builtin_ctzll(m);m&=m-1;t|=1ULL<<perms[g*n+p];}uint64_t rt=rankm(t);if(rid[rt]<0){rid[rt]=nr;shift[rt]=g;len++;}}lens[nr++]=len;}
 if(r+1<dim){uint64_t c=mask&-mask,z=mask+c;mask=(((z^mask)>>2)/c)|z;}}
 return nr;}
void mv36(int nb,int ne,const uint64_t*masks,const uint8_t*lens,const uint64_t*flips,const double*js,const double*diag,const int*rid,const uint8_t*shift,const Z*phase,const Z*x,Z*y,int threads){
#pragma omp parallel for num_threads(threads) schedule(static)
 for(int i=0;i<nb;i++){uint64_t m=masks[i];Z v=diag[i]*x[i];for(int b=0;b<ne;b++){uint64_t bits=m&flips[b];if(bits && bits!=flips[b]){uint64_t r=rankm(m^flips[b]);int j=rid[r];v+=.5*js[b]*phase[shift[r]]*sqrt(double(lens[i])/lens[j])*x[j];}}y[i]=v;}}
}
