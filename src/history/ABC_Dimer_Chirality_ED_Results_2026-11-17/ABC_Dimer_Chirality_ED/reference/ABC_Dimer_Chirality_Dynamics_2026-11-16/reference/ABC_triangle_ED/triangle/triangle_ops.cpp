#include <cstdint>
#include <complex>
#include <vector>
#include <cmath>
#include <omp.h>
using Z=std::complex<double>;
static inline int ranker(uint32_t t,int split,int nh,const uint8_t*pop,const int*rl,const int*rh){int l=t&((1u<<split)-1);return rl[l]+rh[pop[l]*nh+(t>>split)];}
extern "C" void triangle_seeds(const int*c,const uint32_t*masks,const uint8_t*tlens,const int*rid,const uint8_t*shift,const int*rb,const uint8_t*slens,const uint8_t*pop,const int*rl,const int*rh,const Z*phase,const Z*psi,const int*faces,const Z*weights,const Z*ops,Z*out,Z*rdm){
 int nb=c[0],nc=c[1],split=c[2],nh=c[3],threads=c[4],fast=c[5],do_rdm=c[6];
 #pragma omp parallel num_threads(threads)
 {
 std::vector<Z> local(6*64,Z(0));
 #pragma omp for schedule(static)
 for(int i=0;i<nb;i++){
  uint32_t mask=masks[i];Z acc[18]={};
  for(int R=0;R<nc;R++)for(int mu=0;mu<6;mu++){
   int f=R*6+mu,a=faces[3*f],b=faces[3*f+1],d=faces[3*f+2];uint32_t bits=(1u<<a)|(1u<<b)|(1u<<d);int row=((mask>>a)&1)|(((mask>>b)&1)<<1)|(((mask>>d)&1)<<2);int number=__builtin_popcount((unsigned)row);Z val[3]={};
   for(int col=0;col<8;col++)if(__builtin_popcount((unsigned)col)==number){
    uint32_t pre=(mask&~bits)|((col&1)?1u<<a:0)|((col&2)?1u<<b:0)|((col&4)?1u<<d:0);Z amplitude;
    if(fast && pre==mask)amplitude=psi[i];else{int rank=ranker(pre,split,nh,pop,rl,rh),j=rb[rid[rank]];amplitude=(j<0)?Z(0):psi[j]*phase[shift[rank]]*(fast?1.:std::sqrt(double(tlens[i])/slens[j]));}
    for(int k=0;k<3;k++)val[k]+=ops[((mu*3+k)*8+row)*8+col]*amplitude;
    if(do_rdm)local[(mu*8+col)*8+row]+=std::conj(psi[i])*amplitude/double(nc);
   }
   for(int k=0;k<3;k++)acc[k*6+mu]+=weights[R]*val[k];
  }
  for(int k=0;k<18;k++)out[k*nb+i]=acc[k];
 }
 if(do_rdm){
 #pragma omp critical
  {for(int k=0;k<6*64;k++)rdm[k]+=local[k];}
 }
}

}
