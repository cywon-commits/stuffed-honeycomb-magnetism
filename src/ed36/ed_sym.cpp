// ed_sym: symmetry-adapted matrix-free Lanczos ED for the spin-1/2
// stuffed-honeycomb Heisenberg model  H = J sum_rim S.S + J' sum_hub-rim S.T  (J = 1, J' = r).
//
// Sector specifications are produced by symmetry.py.  Bits 0..NR-1 are rim
// sites, bits NR..N-1 are hub sites.  Every symmetry element maps hub sites to
// hub sites, which allows a fast representative search:
//   minhub[h]  = min over elements of the image of hub pattern h
//   cand[h]    = elements achieving that minimum
// so only cand[h] (usually 1..few elements) must be applied to the rim bits.
//
// Basis state |r~> = P|r>/||P|r>||,  P = sum_g chi(g)^* g,  n_r = sum_{g r = r} chi(g)^*.
// Matrix element (pull form):  y[r] += h * chi(g) * sqrt(n_r'/n_r) * x[r'],  g s = r'.
//
// Modes:  dim | bench | lanczos | dense (small sectors only).
#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <complex>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <random>
#include <sstream>
#include <string>
#include <vector>
#include <omp.h>
#include <sys/resource.h>

using u64 = uint64_t;
using cplx = std::complex<double>;

static double now() { return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count(); }
static double peak_rss_mib() { struct rusage ru; getrusage(RUSAGE_SELF, &ru); return ru.ru_maxrss / 1024.0; }

static inline double conj_(double x) { return x; }
static inline cplx conj_(cplx x) { return std::conj(x); }
static inline double re_(double x) { return x; }
static inline double re_(cplx x) { return x.real(); }
template <class T> static T from_c(cplx c);
template <> double from_c<double>(cplx c) { return c.real(); }
template <> cplx from_c<cplx>(cplx c) { return c; }

// ------------------------------------------------------------------ spec
struct Spec {
  int N = 0, NH = 0, NR = 0, NUP = 0;
  std::vector<int> bi, bj, bt;              // bonds (bit positions, type)
  std::vector<int> flip;                    // per element
  std::vector<cplx> chi;                    // per element
  std::vector<std::vector<int>> perm;       // per element: bit b -> perm[b]
  bool real = true;
  std::string header;
};

static Spec read_spec(const char* path) {
  std::ifstream f(path);
  if (!f) { fprintf(stderr, "cannot open %s\n", path); exit(1); }
  Spec s; std::string line, tok;
  while (std::getline(f, line)) {
    if (line.empty()) continue;
    if (line[0] == '#') { s.header += line + "\n"; continue; }
    std::istringstream is(line); is >> tok;
    if (tok == "N") { is >> s.N >> tok >> s.NH >> tok >> s.NUP; s.NR = s.N - s.NH; }
    else if (tok == "BONDS") {
      int nb; is >> nb;
      for (int k = 0; k < nb; k++) { int a, b, t; f >> a >> b >> t; s.bi.push_back(a); s.bj.push_back(b); s.bt.push_back(t); }
    } else if (tok == "ELEMENTS") {
      int ne; is >> ne;
      for (int e = 0; e < ne; e++) {
        int fl; double re, im; f >> fl >> re >> im;
        std::vector<int> p(s.N); for (int b = 0; b < s.N; b++) f >> p[b];
        s.flip.push_back(fl); s.chi.push_back(cplx(re, im)); s.perm.push_back(p);
        if (std::abs(im) > 1e-12) s.real = false;
      }
    }
  }
  if (s.N > 64 || s.NH > 16 || s.NR > 48) { fprintf(stderr, "unsupported sizes\n"); exit(1); }
  return s;
}

// ------------------------------------------------------------------ symmetry tables
struct Sym {
  int N, NH, NR, NE, nchunk;
  u64 rimmask, hubmask;
  std::vector<uint16_t> hubimg;     // [e][h]
  std::vector<u64> rimtab;          // [e][chunk][256]
  std::vector<uint16_t> minhub;     // [h]
  std::vector<uint32_t> coff;       // CSR offsets into cand
  std::vector<uint16_t> cand;       // element indices
  std::vector<uint8_t> flip;
  std::vector<cplx> chi;

  void build(const Spec& s) {
    N = s.N; NH = s.NH; NR = s.NR; NE = (int)s.perm.size();
    nchunk = (NR + 7) / 8;
    rimmask = (NR == 64) ? ~0ULL : ((1ULL << NR) - 1);
    hubmask = (1ULL << NH) - 1;
    flip.assign(s.flip.begin(), s.flip.end()); chi = s.chi;
    const u64 H = 1ULL << NH;
    hubimg.assign((size_t)NE * H, 0);
    rimtab.assign((size_t)NE * nchunk * 256, 0);
    for (int e = 0; e < NE; e++) {
      const auto& p = s.perm[e];
      for (int b = 0; b < NR; b++) if (p[b] >= NR) { fprintf(stderr, "element maps rim to hub\n"); exit(1); }
      for (u64 h = 0; h < H; h++) {
        u64 o = 0;
        for (int b = 0; b < NH; b++) if (h >> b & 1) o |= 1ULL << (p[NR + b] - NR);
        if (flip[e]) o ^= hubmask;
        hubimg[(size_t)e * H + h] = (uint16_t)o;
      }
      for (int c = 0; c < nchunk; c++)
        for (int v = 0; v < 256; v++) {
          u64 o = 0;
          for (int k = 0; k < 8; k++) { int b = 8 * c + k; if (b < NR && (v >> k & 1)) o |= 1ULL << p[b]; }
          rimtab[((size_t)e * nchunk + c) * 256 + v] = o;
        }
    }
    minhub.assign(H, 0); coff.assign(H + 1, 0); cand.clear();
    for (u64 h = 0; h < H; h++) {
      uint16_t m = 0xFFFF;
      for (int e = 0; e < NE; e++) m = std::min(m, hubimg[(size_t)e * H + h]);
      minhub[h] = m;
      for (int e = 0; e < NE; e++) if (hubimg[(size_t)e * H + h] == m) cand.push_back((uint16_t)e);
      coff[h + 1] = (uint32_t)cand.size();
    }
  }
  inline u64 apply_rim(int e, u64 r) const {
    const u64* t = &rimtab[(size_t)e * nchunk * 256];
    u64 o = 0;
    for (int c = 0; c < nchunk; c++) o |= t[c * 256 + ((r >> (8 * c)) & 255)];
    return flip[e] ? (o ^ rimmask) : o;
  }
  // representative of s and one element g with g s = rep
  inline u64 rep(u64 s, int& g) const {
    u64 h = s >> NR, r = s & rimmask, best = ~0ULL;
    for (uint32_t k = coff[h]; k < coff[h + 1]; k++) {
      int e = cand[k]; u64 o = apply_rim(e, r);
      if (o < best) { best = o; g = e; }
    }
    return ((u64)minhub[h] << NR) | best;
  }
};

// ------------------------------------------------------------------ basis
static u64 binom(int n, int k) {
  if (k < 0 || k > n) return 0;
  u64 r = 1; for (int i = 1; i <= k; i++) r = r * (n - k + i) / i; return r;
}
// m-th (0-based) k-subset of n bits in increasing integer order (colex)
static u64 unrank(int n, int k, u64 m) {
  u64 s = 0;
  for (int b = n - 1; b >= 0 && k > 0; b--) {
    u64 c = binom(b, k);       // subsets not using bit b (all below)
    if (m >= c) { s |= 1ULL << b; m -= c; k--; }
  }
  return s;
}
static inline u64 next_comb(u64 v) { u64 t = v | (v - 1); return (t + 1) | (((~t & -~t) - 1) >> (__builtin_ctzll(v) + 1)); }

struct Basis {
  std::vector<u64> reps;
  std::vector<double> norm;          // n_r (real, > 0)
  std::vector<int32_t> hubidx;       // minimal hub -> block id (-1 otherwise)
  std::vector<u64> off;              // bucket offsets
  int rb = 0;                        // rim key bits
  int NR = 0;

  void build(const Sym& S, int NUP) {
    NR = S.NR;
    const u64 H = 1ULL << S.NH;
    std::vector<u64> hubs;
    for (u64 h = 0; h < H; h++) if (S.minhub[h] == h) {
      int kr = NUP - __builtin_popcountll(h);
      if (kr >= 0 && kr <= S.NR) hubs.push_back(h);
    }
    // work items: (hub, chunk of rim combinations)
    struct Item { u64 h; int kr; u64 m0, m1; };
    std::vector<Item> items;
    for (u64 h : hubs) {
      int kr = NUP - __builtin_popcountll(h); u64 C = binom(S.NR, kr);
      u64 nch = std::max<u64>(1, C / 200000);
      for (u64 c = 0; c < nch; c++) items.push_back({h, kr, C * c / nch, C * (c + 1) / nch});
    }
    std::vector<std::vector<u64>> R(items.size());
    std::vector<std::vector<double>> Nn(items.size());
#pragma omp parallel for schedule(dynamic, 1)
    for (size_t it = 0; it < items.size(); it++) {
      const Item& I = items[it];
      u64 r = unrank(S.NR, I.kr, I.m0);
      for (u64 m = I.m0; m < I.m1; m++) {
        bool ok = true; cplx nrm = 0;
        for (uint32_t k = S.coff[I.h]; k < S.coff[I.h + 1]; k++) {
          int e = S.cand[k]; u64 o = S.apply_rim(e, r);
          if (o < r) { ok = false; break; }
          if (o == r) nrm += std::conj(S.chi[e]);
        }
        if (ok && nrm.real() > 1e-8) { R[it].push_back((I.h << S.NR) | r); Nn[it].push_back(nrm.real()); }
        if (I.kr > 0 && m + 1 < I.m1) r = next_comb(r);
      }
    }
    size_t tot = 0; for (auto& v : R) tot += v.size();
    reps.reserve(tot); norm.reserve(tot);
    for (size_t it = 0; it < items.size(); it++) {
      reps.insert(reps.end(), R[it].begin(), R[it].end()); norm.insert(norm.end(), Nn[it].begin(), Nn[it].end());
      std::vector<u64>().swap(R[it]); std::vector<double>().swap(Nn[it]);
    }
    // lookup buckets: (hub block, top rb rim bits)
    hubidx.assign(H, -1); int nb = 0;
    for (u64 h : hubs) hubidx[h] = nb++;
    double per = (double)tot / std::max(1, nb);
    rb = 0; while (rb < S.NR && (double)(1ULL << rb) * 64 < per) rb++;
    size_t nbk = (size_t)nb << rb;
    off.assign(nbk + 1, 0);
    for (u64 s : reps) off[key(s) + 1]++;
    for (size_t k = 0; k < nbk; k++) off[k + 1] += off[k];
  }
  inline size_t key(u64 s) const {
    return ((size_t)hubidx[s >> NR] << rb) | (size_t)((s & ((1ULL << NR) - 1)) >> (NR - rb));
  }
  inline int64_t find(u64 s) const {
    int32_t hb = hubidx[s >> NR]; if (hb < 0) return -1;
    size_t k = key(s);
    const u64* lo = reps.data() + off[k]; const u64* hi = reps.data() + off[k + 1];
    const u64* p = std::lower_bound(lo, hi, s);
    return (p != hi && *p == s) ? (int64_t)(p - reps.data()) : -1;
  }
};

// ------------------------------------------------------------------ Hamiltonian
struct Ham {
  const Sym* S; const Basis* B;
  std::vector<u64> bmask; std::vector<double> bJ;
  std::vector<double> sqn;     // sqrt(norm)
  void init(const Spec& sp, const Sym& s, const Basis& b, double r) {
    S = &s; B = &b;
    for (size_t k = 0; k < sp.bi.size(); k++) { bmask.push_back((1ULL << sp.bi[k]) | (1ULL << sp.bj[k])); bJ.push_back(sp.bt[k] == 0 ? 1.0 : r); }
    sqn.resize(b.norm.size());
    for (size_t i = 0; i < sqn.size(); i++) sqn[i] = std::sqrt(b.norm[i]);
  }
  template <class T> void apply(const T* x, T* y) const {
    const size_t D = B->reps.size(); const int nb = (int)bmask.size();
#pragma omp parallel for schedule(dynamic, 4096)
    for (size_t i = 0; i < D; i++) {
      const u64 s = B->reps[i];
      double diag = 0; T acc = 0;
      for (int k = 0; k < nb; k++) {
        const u64 m = s & bmask[k];
        if (m == 0 || m == bmask[k]) { diag += 0.25 * bJ[k]; continue; }
        diag -= 0.25 * bJ[k];
        int g = 0; const u64 rp = S->rep(s ^ bmask[k], g);
        const int64_t j = B->find(rp); if (j < 0) continue;
        acc += (0.5 * bJ[k] * sqn[j]) * from_c<T>(S->chi[g]) * x[j];
      }
      y[i] = diag * x[i] + acc / sqn[i];
    }
  }
};

// ------------------------------------------------------------------ linear algebra helpers
template <class T> static double dotr(const std::vector<T>& a, const std::vector<T>& b) {  // Re <a,b>
  double s = 0; const size_t n = a.size();
#pragma omp parallel for reduction(+ : s)
  for (size_t i = 0; i < n; i++) s += re_(conj_(a[i]) * b[i]);
  return s;
}
template <class T> static cplx dotc(const std::vector<T>& a, const std::vector<T>& b) {
  double sr = 0, si = 0; const size_t n = a.size();
#pragma omp parallel for reduction(+ : sr, si)
  for (size_t i = 0; i < n; i++) { cplx v = cplx(conj_(a[i])) * cplx(b[i]); sr += v.real(); si += v.imag(); }
  return cplx(sr, si);
}
template <class T> static void random_vec(std::vector<T>& v, uint64_t seed) {
  const size_t n = v.size();
#pragma omp parallel
  {
    int t = omp_get_thread_num(), nt = omp_get_num_threads();
    std::mt19937_64 rng(seed * 1000003ULL + t);
    std::normal_distribution<double> nd;
    size_t a = n * t / nt, b = n * (t + 1) / nt;
    for (size_t i = a; i < b; i++) { if constexpr (std::is_same<T, double>::value) v[i] = nd(rng); else v[i] = cplx(nd(rng), nd(rng)); }
  }
  double nr = std::sqrt(dotr(v, v));
#pragma omp parallel for
  for (size_t i = 0; i < n; i++) v[i] /= nr;
}

// symmetric tridiagonal eigen-decomposition (implicit QL, from standard tqli), returns sorted evals and evecs (column-major n x n)
static void tqli(std::vector<double> d, std::vector<double> e, std::vector<double>& evals, std::vector<double>& Z, bool wantz = true) {
  int n = (int)d.size(); Z.assign(wantz ? (size_t)n * n : 0, 0); if (wantz) for (int i = 0; i < n; i++) Z[(size_t)i * n + i] = 1;
  e.push_back(0.0);  // e[i] couples i,i+1
  for (int l = 0; l < n; l++) {
    int iter = 0, m;
    do {
      for (m = l; m < n - 1; m++) { double dd = std::fabs(d[m]) + std::fabs(d[m + 1]); if (std::fabs(e[m]) <= 1e-15 * dd) break; }
      if (m != l) {
        if (iter++ == 300) break;
        double g = (d[l + 1] - d[l]) / (2.0 * e[l]); double r = std::hypot(g, 1.0);
        g = d[m] - d[l] + e[l] / (g + (g >= 0 ? std::fabs(r) : -std::fabs(r)));
        double s = 1, c = 1, p = 0; int i;
        for (i = m - 1; i >= l; i--) {
          double f = s * e[i], b = c * e[i];
          e[i + 1] = (r = std::hypot(f, g));
          if (r == 0.0) { d[i + 1] -= p; e[m] = 0.0; break; }
          s = f / r; c = g / r; g = d[i + 1] - p; r = (d[i] - g) * s + 2.0 * c * b; d[i + 1] = g + (p = s * r); g = c * r - b;
          if (wantz) for (int k = 0; k < n; k++) { double fz = Z[(size_t)(i + 1) * n + k]; Z[(size_t)(i + 1) * n + k] = s * Z[(size_t)i * n + k] + c * fz; Z[(size_t)i * n + k] = c * Z[(size_t)i * n + k] - s * fz; }
        }
        if (r == 0.0 && i >= l) continue;
        d[l] -= p; e[l] = g; e[m] = 0.0;
      }
    } while (m != l);
  }
  std::vector<int> idx(n); for (int i = 0; i < n; i++) idx[i] = i;
  std::sort(idx.begin(), idx.end(), [&](int a, int b) { return d[a] < d[b]; });
  evals.resize(n); std::vector<double> Z2(wantz ? (size_t)n * n : 0);
  for (int i = 0; i < n; i++) { evals[i] = d[idx[i]]; if (wantz) for (int k = 0; k < n; k++) Z2[(size_t)i * n + k] = Z[(size_t)idx[i] * n + k]; }
  Z.swap(Z2);  // Z[i*n+k] = component k of eigenvector i
}

// ------------------------------------------------------------------ checkpoint I/O
template <class T> static void save_vec(FILE* f, const std::vector<T>& v) { fwrite(v.data(), sizeof(T), v.size(), f); }
template <class T> static bool load_vec(FILE* f, std::vector<T>& v) { return fread(v.data(), sizeof(T), v.size(), f) == v.size(); }

struct Opts {
  std::string mode = "lanczos", spec, ckpt, out;
  double r = 0.2, tol = 1e-10;
  int maxit = 400, nbench = 3, ckevery = 20, nev = 3; uint64_t seed = 12345; bool vec = true;
};

template <class T>
static void lanczos(const Opts& o, const Ham& H, size_t D, FILE* log) {
  std::vector<T> v0(D), v1(D), w(D);
  std::vector<double> al, be;  // be[j] couples v_j, v_{j+1}
  int start = 0; double t0 = now();
  std::string ck = o.ckpt;
  bool resumed = false;
  if (!ck.empty()) {
    FILE* f = fopen((ck + ".lz").c_str(), "rb");
    if (f) {
      uint64_t d, it; fread(&d, 8, 1, f); fread(&it, 8, 1, f);
      if (d == D) {
        al.resize(it); be.resize(it); fread(al.data(), 8, it, f); fread(be.data(), 8, it, f);
        if (load_vec(f, v0) && load_vec(f, v1)) { start = (int)it; resumed = true; }
      }
      fclose(f);
      if (resumed) fprintf(log, "resumed pass-1 checkpoint at iteration %d\n", start);
    }
  }
  if (!resumed) { random_vec(v1, o.seed); std::fill(v0.begin(), v0.end(), T(0)); }
  std::vector<double> ev, Z; double Eold = 1e300; int conv_it = -1; double resid_est = 1;
  for (int j = start; j < o.maxit; j++) {
    double tm = now();
    H.apply(v1.data(), w.data());
    double a = dotr(v1, w); double bprev = j > 0 ? be[j - 1] : 0.0;
#pragma omp parallel for
    for (size_t i = 0; i < D; i++) w[i] -= a * v1[i] + bprev * v0[i];
    double b = std::sqrt(dotr(w, w));
    al.push_back(a); be.push_back(b);
    // eigenvalues every iteration (O(n^2)); Ritz vectors only when E0 has settled or periodically (O(n^3))
    std::vector<double> sub(be.begin(), be.end() - 1);
    tqli(al, sub, ev, Z, false);
    int n = (int)al.size();
    if (std::fabs(ev[0] - Eold) < 1e-12 * std::max(1.0, std::fabs(ev[0])) || (j + 1) % 25 == 0 || j + 1 == o.maxit) {
      std::vector<double> ev2; tqli(al, sub, ev2, Z, true);
      resid_est = std::fabs(b * Z[(size_t)0 * n + (n - 1)]);
    }
    fprintf(log, "it %4d  E0 % .14f  E1 % .14f  ritz_resid %.3e  dt %.2fs  t %.1fs\n", j + 1, ev[0], n > 1 ? ev[1] : 0.0, resid_est, now() - tm, now() - t0);
    fflush(log);
    if (b < 1e-14) { conv_it = j + 1; break; }
#pragma omp parallel for
    for (size_t i = 0; i < D; i++) { v0[i] = v1[i]; v1[i] = w[i] / b; }
    if (!ck.empty() && (j + 1) % o.ckevery == 0) {
      std::string tmp = ck + ".lz.tmp"; FILE* f = fopen(tmp.c_str(), "wb");
      uint64_t d = D, it = al.size(); fwrite(&d, 8, 1, f); fwrite(&it, 8, 1, f);
      fwrite(al.data(), 8, it, f); fwrite(be.data(), 8, it, f); save_vec(f, v0); save_vec(f, v1); fclose(f);
      rename(tmp.c_str(), (ck + ".lz").c_str());
    }
    if (resid_est < o.tol) { conv_it = j + 1; break; }
    Eold = ev[0];
  }
  int n = (int)al.size();
  { std::vector<double> ev2; tqli(al, std::vector<double>(be.begin(), be.end() - 1), ev2, Z, true); }
  fprintf(log, "PASS1 iterations %d converged %s E0 %.15f ritz_resid %.3e\n", n, conv_it > 0 ? "yes" : "no", ev[0], resid_est);
  fprintf(log, "RITZ");
  for (int k = 0; k < std::min(o.nev * 2, n); k++) fprintf(log, " %.14f", ev[k]);
  fprintf(log, "\n");
  if (!o.vec) return;
  // pass 2: rebuild Krylov vectors and accumulate ground-state vector
  std::vector<T> psi(D, T(0));
  random_vec(v1, o.seed); std::fill(v0.begin(), v0.end(), T(0));
  for (int j = 0; j < n; j++) {
    double c = Z[(size_t)0 * n + j];
#pragma omp parallel for
    for (size_t i = 0; i < D; i++) psi[i] += c * v1[i];
    if (j == n - 1) break;
    H.apply(v1.data(), w.data());
    double bprev = j > 0 ? be[j - 1] : 0.0;
#pragma omp parallel for
    for (size_t i = 0; i < D; i++) { T t = (w[i] - al[j] * v1[i] - bprev * v0[i]) / be[j]; v0[i] = v1[i]; v1[i] = t; }
  }
  double nrm = std::sqrt(dotr(psi, psi));
#pragma omp parallel for
  for (size_t i = 0; i < D; i++) psi[i] /= nrm;
  H.apply(psi.data(), w.data());
  double E = dotr(psi, w);
#pragma omp parallel for
  for (size_t i = 0; i < D; i++) w[i] -= E * psi[i];
  double res = std::sqrt(dotr(w, w));
  fprintf(log, "GROUND E0 %.15f residual %.3e norm_before %.12f (pass-2 total %.1fs)\n", E, res, nrm, now() - t0);
  if (!o.out.empty()) {
    FILE* f = fopen(o.out.c_str(), "wb"); uint64_t d = D; fwrite(&d, 8, 1, f); fwrite(&E, 8, 1, f); save_vec(f, psi); fclose(f);
    fprintf(log, "wrote ground-state vector to %s\n", o.out.c_str());
  }
  if (!ck.empty()) remove((ck + ".lz").c_str());
}

template <class T>
static void run(const Opts& o, const Spec& sp, const Sym& S, const Basis& B) {
  Ham H; H.init(sp, S, B, o.r);
  const size_t D = B.reps.size();
  if (o.mode == "bench") {
    std::vector<T> x(D), y(D), x2(D), y2(D);
    random_vec(x, 1); random_vec(x2, 2);
    for (int k = 0; k < o.nbench; k++) { double t = now(); H.apply(x.data(), y.data()); printf("BENCH H_action %.3f s\n", now() - t); fflush(stdout); }
    H.apply(x2.data(), y2.data());
    cplx a = dotc(x2, y), b = dotc(y2, x);
    printf("HERMITICITY |<x2,Hx>-<Hx2,x>| = %.3e (|<x2,Hx>| = %.3e)\n", std::abs(a - b), std::abs(a));
  } else if (o.mode == "dense") {
    if (D > 6000) { fprintf(stderr, "dense mode limited to D<=6000\n"); exit(1); }
    std::vector<T> e(D), col(D);
    printf("DENSE %zu\n", D);
    for (size_t j = 0; j < D; j++) {
      std::fill(e.begin(), e.end(), T(0)); e[j] = 1; H.apply(e.data(), col.data());
      for (size_t i = 0; i < D; i++) { cplx c = cplx(col[i]); printf("%zu %zu %.17g %.17g\n", i, j, c.real(), c.imag()); }
    }
  } else {
    lanczos<T>(o, H, D, stdout);
  }
}

int main(int argc, char** argv) {
  Opts o;
  for (int i = 1; i < argc; i++) {
    std::string a = argv[i];
    auto nx = [&]() { return std::string(argv[++i]); };
    if (a == "--mode") o.mode = nx();
    else if (a == "--r") o.r = std::stod(nx());
    else if (a == "--tol") o.tol = std::stod(nx());
    else if (a == "--maxit") o.maxit = std::stoi(nx());
    else if (a == "--nbench") o.nbench = std::stoi(nx());
    else if (a == "--ckpt") o.ckpt = nx();
    else if (a == "--ckevery") o.ckevery = std::stoi(nx());
    else if (a == "--out") o.out = nx();
    else if (a == "--seed") o.seed = std::stoull(nx());
    else if (a == "--novec") o.vec = false;
    else o.spec = a;
  }
  if (o.spec.empty()) { fprintf(stderr, "usage: ed_sym SPEC [--mode dim|bench|lanczos|dense] [--r R] [--maxit N] [--tol T] [--ckpt PREFIX] [--out FILE]\n"); return 1; }
  double t0 = now();
  Spec sp = read_spec(o.spec.c_str());
  printf("%sSPEC N=%d NH=%d NUP=%d bonds=%zu elements=%zu real=%d threads=%d r=%.10g\n", sp.header.c_str(), sp.N, sp.NH, sp.NUP, sp.bi.size(), sp.perm.size(), (int)sp.real, omp_get_max_threads(), o.r);
  Sym S; S.build(sp);
  double avgc = (double)S.cand.size() / (1ULL << S.NH);
  printf("TABLES built %.2fs  mean candidates per hub pattern %.2f\n", now() - t0, avgc);
  double t1 = now();
  Basis B; B.build(S, sp.NUP);
  size_t mem = B.reps.size() * 16 + B.off.size() * 8;
  printf("BASIS dim %zu  build %.2fs  basis+lookup %.3f GiB  rim_key_bits %d  peak_RSS %.1f MiB\n", B.reps.size(), now() - t1, mem / 1073741824.0, B.rb, peak_rss_mib());
  fflush(stdout);
  if (o.mode == "dim") return 0;
  if (sp.real) run<double>(o, sp, S, B); else run<cplx>(o, sp, S, B);
  printf("DONE total %.1fs peak_RSS %.1f MiB\n", now() - t0, peak_rss_mib());
  return 0;
}
