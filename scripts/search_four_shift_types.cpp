// Exact bounded center search for Theorem 14, types 4, 5, 15 and 16.
// g++ -O3 -std=c++17 -fopenmp scripts/search_four_shift_types.cpp -o build/search_four_shift_types
// build/search_four_shift_types HEIGHT THREADS [DENOMINATOR_FIRST DENOMINATOR_LAST]
// Searches c=(n/d)^2 with 0<=n<=HEIGHT, gcd(n,d)=1. Finite evidence only.
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>
#include <omp.h>

using I = __int128_t;
struct Family {
    int type, scale;
    // Scaled squared centers X=D*d^2*center^2 have monic polynomials
    // X^2-(s0*n^2+s1*d^2)X+(t0*n^4+t1*n^2*d^2+t2*d^4).
    std::array<int, 2> s, u;
    std::array<int, 3> t, v;
};
const std::array<Family, 4> families{{
    {4, 12, {40,50}, {44,91}, {720,120,501}, {1024,388,1864}},
    {5, 16, {20,-85}, {12,-83}, {80,-840,589}, {16,-232,441}},
    {15,24, {60,45}, {-12,207}, {960,-120,294}, {96,-1848,8556}},
    {16,12, {-20,195}, {-16,144}, {640,-5340,9216}, {384,-3168,5184}}
}};

I square_root(I a) {
    if (a < 0) return -1;
    I r = static_cast<I>(sqrtl(static_cast<long double>(a)));
    // Floating point only supplies a guess; the decision is entirely integral.
    while (r*r > a) --r;
    while ((r+1)*(r+1) <= a) ++r;
    return r*r == a ? r : -1;
}

bool split_into_squares(I S, I T, int D) {
    I r = square_root(S*S-4*T);
    if (r < 0 || (S+r)%2 != 0) return false;
    return square_root(D*((S+r)/2)) >= 0 &&
           square_root(D*((S-r)/2)) >= 0;
}

bool exact(const Family& f, int n, int d) {
    I c=I(n)*n, z=I(d)*d;
    auto lin = [&](auto a) { return a[0]*c+a[1]*z; };
    auto quad = [&](auto a) { return a[0]*c*c+a[1]*c*z+a[2]*z*z; };
    return split_into_squares(lin(f.s),quad(f.t),f.scale) &&
           split_into_squares(lin(f.u),quad(f.v),f.scale);
}

int mod(long long a, int p) { int r=static_cast<int>(a%p); return r<0?r+p:r; }
int inverse(int a, int p) {
    for (int b=1;b<p;++b) if (a*b%p==1) return b;
    throw std::runtime_error("nonunit inverse");
}
struct Filter {
    int p;
    std::vector<uint8_t> allowed;
    // Row d, column n; primitive projective points only.
    bool permits(int n,int d) const {return allowed[(d%p)*p+n%p];}
};

Filter make_filter(const Family& f,int p) {
    std::vector<uint8_t> square(p), splits(p*p), ratio(p);
    for(int r=0;r<p;++r) square[f.scale*r*r%p]=1;
    for(int a=0;a<p;++a) if(square[a])
        for(int b=0;b<p;++b) if(square[b]) splits[((a+b)%p)*p+a*b%p]=1;
    auto at=[&](long long c,long long z) {
        auto lin=[&](auto a){return mod(a[0]*c+a[1]*z,p);};
        auto quad=[&](auto a){return mod(a[0]*c*c+a[1]*c*z+a[2]*z*z,p);};
        return splits[lin(f.s)*p+quad(f.t)] && splits[lin(f.u)*p+quad(f.v)];
    };
    for(int r=0;r<p;++r) ratio[r]=at(r*r%p,1);
    bool infinity=at(1,0);
    Filter result{p,std::vector<uint8_t>(p*p)};
    for(int d=0;d<p;++d) for(int n=0;n<p;++n)
        result.allowed[d*p+n]=d?ratio[n*inverse(d,p)%p]:(n && infinity);
    return result;
}

constexpr int wheel_modulus=7*11*13*17;
struct Sieve {
    std::vector<Filter> filters;
    std::vector<std::vector<uint16_t>> wheel;
};

Sieve make_sieve(const Family& f) {
    Sieve result;
    for(int p:{7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,
               71,73,79,83,89,97,101,103,107,109,113,127,131,137,139,149})
        result.filters.push_back(make_filter(f,p));
    result.wheel.resize(wheel_modulus);
    for(int d=0;d<wheel_modulus;++d) {
        std::vector<int> residues{0};
        int modulus=1;
        for(int j=0;j<4;++j) {
            const auto& filter=result.filters[j]; int p=filter.p;
            int inv=inverse(modulus%p,p);
            std::vector<int> next;
            for(int r:residues) for(int n=0;n<p;++n) if(filter.permits(n,d))
                next.push_back(r+mod(n-r,p)*inv%p*modulus);
            modulus*=p; residues.swap(next);
        }
        for(int r:residues) result.wheel[d].push_back(static_cast<uint16_t>(r));
        std::sort(result.wheel[d].begin(),result.wheel[d].end());
    }
    return result;
}

bool passes(const Sieve& sieve,int n,int d) {
    for(size_t j=4;j<sieve.filters.size();++j)
        if(!sieve.filters[j].permits(n,d)) return false;
    return true;
}

int main(int argc,char** argv) {
    try {
        if(argc<3 || argc>5) throw std::runtime_error(
            "usage: search_four_shift_types HEIGHT THREADS [DENOMINATOR_FIRST DENOMINATOR_LAST]");
        int H=std::stoi(argv[1]), threads=std::stoi(argv[2]);
        int first=argc>=4?std::stoi(argv[3]):1, last=argc>=5?std::stoi(argv[4]):H;
        // This cap bounds every exact intermediate below 2^127 (coefficients <10000).
        if(H<1 || H>10000000 || threads<1 || first<1 || last<first || last>H)
            throw std::runtime_error("require 1<=first<=last<=HEIGHT<=10000000 and THREADS>=1");
        omp_set_num_threads(threads);
        std::array<Sieve,4> sieves;
        for(int k=0;k<4;++k) sieves[k]=make_sieve(families[k]);
        // Small independent loop: no wheel or finite-field rejection may lose a hit.
        for(int k=0;k<4;++k) for(int d=1;d<=100;++d) for(int n=0;n<=100;++n)
            if(std::gcd(n,d)==1 && exact(families[k],n,d)) {
                const auto& w=sieves[k].wheel[d%wheel_modulus];
                if(!std::binary_search(w.begin(),w.end(),n%wheel_modulus) || !passes(sieves[k],n,d))
                    throw std::runtime_error("sieve dropped an exact small control");
            }
        std::cout<<"BOUND "<<H<<" DENOMINATORS "<<first<<" "<<last<<"\n"<<std::flush;
        unsigned long long candidates=0, survivors=0, hits=0;
        double start=omp_get_wtime();
        #pragma omp parallel for schedule(dynamic,16) reduction(+:candidates,survivors,hits)
        for(int d=first;d<=last;++d) for(int k=0;k<4;++k) {
            const auto& sieve=sieves[k];
            for(int r:sieve.wheel[d%wheel_modulus]) for(int n=r;n<=H;n+=wheel_modulus) {
                ++candidates;
                if(!passes(sieve,n,d) || std::gcd(n,d)!=1) continue;
                ++survivors;
                if(!exact(families[k],n,d)) continue;
                ++hits;
                #pragma omp critical
                std::cout<<"HIT "<<families[k].type<<" "<<n<<" "<<d<<"\n";
            }
        }
        std::cout<<"COMPLETE candidates="<<candidates<<" modular_survivors="<<survivors
                 <<" exact_hits="<<hits<<" seconds="<<omp_get_wtime()-start<<"\n";
        std::cout<<"Finite height box only; not a rational-point classification.\n";
    } catch(const std::exception& e) {std::cerr<<e.what()<<"\n";return 1;}
}
