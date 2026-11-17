#include <cstdint>
#include <complex>
#include <cmath>
#include <omp.h>
using Z=std::complex<double>;
extern "C" void action(int nb,int nedge,const uint32_t* masks,const uint8_t* lens,const uint32_t* flips,const double* js,const double* diag,const int* reps,const uint8_t* shift,const int* rb,const Z* phase,const Z* x,Z* y,int threads){
#pragma omp parallel for num_threads(threads) schedule(static)
for(int i=0;i<nb;i++){uint32_t m=masks[i];Z val=diag[i]*x[i];for(int b=0;b<nedge;b++){uint32_t bit=m&flips[b];if(bit && bit!=flips[b]){uint32_t t=m^flips[b];int j=rb[reps[t]];if(j>=0)val+=.5*js[b]*phase[shift[t]]*std::sqrt(double(lens[i])/lens[j])*x[j];}}y[i]=val;}}
