// Exact k-mer counts of reads against genomes, for checking search-fmin's output.
//
// build: g++ -O2 -std=c++17 exact_counts.cpp -o exact_counts
// usage: exact_counts refs.list queries.fa k result.txt t [found]
//
// result.txt is search-fmin's output for queries.fa at threshold t (without --print-denominator).
// Let n be the number of a read's k-mers made of A, C, G, T only, and n_pos the number of
// those that occur in some reference (either strand). A color "should be reported" when it
// holds at least max(1, floor(t * n)) of the read's k-mers (default mode), or, with the
// word "found" as the last argument, when count / n_pos >= t and count >= 1.
// Prints, and returns non-zero on, any violation of:
//   - every color that should be reported is reported (no false negatives)
//   - every reported score is at least the exact count
#include <bits/stdc++.h>
using namespace std;
static int code(char c){ switch(c){case 'A':case 'a':return 0;case 'C':case 'c':return 1;case 'G':case 'g':return 2;case 'T':case 't':return 3;} return -1; }
static vector<string> read_fasta(const string& p){ vector<string> v; ifstream in(p); string l; bool open=false; while(getline(in,l)){ if(!l.empty()&&l.back()=='\r')l.pop_back(); if(l.empty())continue; if(l[0]=='>'){v.emplace_back();open=true;} else if(open) v.back()+=l; } return v; }
template<class F> static void kmers(const string& s,int k,F f){ uint64_t mask=(k==32)?~0ULL:((1ULL<<(2*k))-1),fw=0,rc=0; int run=0; for(size_t i=0;i<s.size();i++){ int c=code(s[i]); if(c<0){run=0;fw=rc=0;continue;} fw=((fw<<2)|c)&mask; rc=(rc>>2)|((uint64_t)(3-c)<<(2*(k-1))); if(++run>=k) f(i+1-k,min(fw,rc)); } }
int main(int argc,char**argv){ string refs=argv[1],q=argv[2]; int k=atoi(argv[3]); string res=argv[4]; float t=atof(argv[5]); double td=atof(argv[5]); bool found_rule=(argc>6 && string(argv[6])=="found"); long same_as_exact_set=0;
  vector<string> reads=read_fasta(q); unordered_map<uint64_t,uint32_t> id; vector<vector<uint32_t>> rk(reads.size());
  for(size_t r=0;r<reads.size();r++) kmers(reads[r],k,[&](size_t,uint64_t x){ auto it=id.emplace(x,id.size()).first; rk[r].push_back(it->second); });
  vector<string> paths; { ifstream in(refs); string l; while(getline(in,l)) if(!l.empty()) paths.push_back(l); }
  size_t nc=paths.size(); vector<vector<bool>> has(nc, vector<bool>(id.size(),false));
  for(size_t c=0;c<nc;c++) for(const string& s: read_fasta(paths[c])) kmers(s,k,[&](size_t,uint64_t x){ auto it=id.find(x); if(it!=id.end()) has[c][it->second]=true; });
  ifstream in(res); string line; size_t r=0; long fn=0,below=0,above=0,equal=0,fp=0,tp=0,lines=0;
  while(getline(in,line)){ istringstream ss(line); uint64_t rid; ss>>rid; if(rid!=r){ printf("line %zu has read number %llu\n",r,(unsigned long long)rid); return 2;} map<int,int> rep; string tok; while(ss>>tok){ size_t p=tok.find(':'); rep[stoi(tok.substr(0,p))]=stoi(tok.substr(p+1)); }
    if(r>=reads.size()){ printf("more result lines than reads\n"); return 2; }
    long n=rk[r].size(); long thr=max<long>(1,(long)(n*t));
    long npos=0; if(found_rule){ for(uint32_t x: rk[r]){ bool any=false; for(size_t c=0;c<nc&&!any;c++) any=has[c][x]; npos+=any; } }
    for(size_t c=0;c<nc;c++){ long e=0; for(uint32_t x: rk[r]) e+=has[c][x]; bool want= found_rule ? (npos>0 && e>=1 && (double)e/(double)npos>=td) : (n>0&&e>=thr); auto it=rep.find(c);
      if(want&&it==rep.end()){ fn++; if(fn<=5) printf("FALSE NEGATIVE read %zu color %zu exact %ld of %ld\n",r,c,e,n); }
      if(it!=rep.end()){ if(it->second<e){ below++; if(below<=5) printf("SCORE BELOW EXACT read %zu color %zu score %d exact %ld\n",r,c,it->second,e);} else if(it->second>e) above++; else equal++; if(want) tp++; else fp++; if(it->second<=0){ printf("ZERO SCORE reported read %zu\n",r); return 2;} } }
    r++; lines++; }
  if(lines!=(long)reads.size()){ printf("result has %ld lines for %zu reads\n",lines,reads.size()); return 2; }
  printf("reads %zu colors %zu | tp %ld fp %ld fn %ld | scores equal %ld above %ld below %ld\n",reads.size(),nc,tp,fp,fn,equal,above,below);
  return (fn||below)?1:0; }
