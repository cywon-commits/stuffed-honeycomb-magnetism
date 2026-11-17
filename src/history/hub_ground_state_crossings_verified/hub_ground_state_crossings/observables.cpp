#include <cstdint>
#include <complex>
#include <cmath>
#include <vector>
#include <omp.h>
using Z=std::complex<double>;
static inline int rank(uint32_t t,int split,int nh,const uint8_t*pop,const int*rl,const int*rh){int l=t&((1u<<split)-1);return rl[l]+rh[pop[l]*nh+(t>>split)];}
extern "C" {
void pair_stats(const int*c,const uint32_t*masks,const uint8_t*lens,const int*rid,const uint8_t*shift,const int*rb,const uint8_t*pop,const int*rl,const int*rh,const int*aa,const int*bb,const int*cls,const Z*phase,const Z*psi,double*dot,double*zz,double*mag){
 int nb=c[0],np=c[1],nc=c[2],ns=c[3],ncell=c[4],split=c[5],nh=c[6],threads=c[7];
 #pragma omp parallel num_threads(threads)
 {
 std::vector<double>d(nc,0.),z(nc,0.);double mgn[3]={0,0,0};
 #pragma omp for schedule(static)
 for(int i=0;i<nb;i++){
  uint32_t m=masks[i];double w=std::norm(psi[i]);Z vc=std::conj(psi[i]);
  for(int p=0;p<ns;p++)mgn[p%3]+=w*(((m>>p)&1)?0.5:-0.5)/ncell;
  for(int b=0;b<np;b++){
   uint32_t flip=(1u<<aa[b])|(1u<<bb[b]),bits=m&flip;bool anti=(bits && bits!=flip);int k=cls[b];double diag=w*(anti?-0.25:0.25);z[k]+=diag;d[k]+=diag;
   if(anti){int r=rank(m^flip,split,nh,pop,rl,rh);int j=rb[rid[r]];if(j>=0)d[k]+=0.5*(lens[i]==lens[j]?1.0:std::sqrt(double(lens[i])/lens[j]))*std::real(vc*phase[shift[r]]*psi[j]);}
  }
 }
 #pragma omp critical
 {for(int j=0;j<nc;j++){dot[j]+=d[j];zz[j]+=z[j];}for(int a=0;a<3;a++)mag[a]+=mgn[a];}
 }
}
void spin_seed(const int*c,const uint32_t*target_masks,const uint8_t*target_lens,const int*rid,const uint8_t*shift,const int*rb,const uint8_t*source_lens,const uint8_t*pop,const int*rl,const int*rh,const Z*phase,const Z*weights,const Z*psi,Z*out){
 int nb=c[0],ns=c[1],split=c[2],nh=c[3],threads=c[4],kind=c[5];
 #pragma omp parallel for num_threads(threads) schedule(static)
 for(int i=0;i<nb;i++){
  uint32_t m=target_masks[i];Z val=0;
  if(kind==0){
   int r=rank(m,split,nh,pop,rl,rh),j=rb[rid[r]];
   if(j>=0){Z w=0;for(int p=0;p<ns;p++)w+=weights[p]*(((m>>p)&1)?0.5:-0.5);val=w*psi[j]*phase[shift[r]]*std::sqrt(double(target_lens[i])/source_lens[j]);}
  }else{
   for(int p=0;p<ns;p++)if(((m>>p)&1)==(kind==1?1:0)){
    uint32_t t=m^(1u<<p);int r=rank(t,split,nh,pop,rl,rh),j=rb[rid[r]];
    if(j>=0)val+=weights[p]*psi[j]*phase[shift[r]]*std::sqrt(double(target_lens[i])/source_lens[j]);
   }
  }
  out[i]=val;
 }
}
}
