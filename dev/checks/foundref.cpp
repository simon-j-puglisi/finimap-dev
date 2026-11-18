// Reference for the found-in-index mode, independent of Finimap's code.
// build: g++ -O2 -std=c++17 foundref.cpp -o foundref
// usage: foundref refs.list queries.fa k kmers_out.fna npos_out.txt [min_run]
// Writes every k-mer of every read that occurs in the references (either strand) as its own
// record, in read order, with the read's number as its name; and the number of such k-mers
// per read (one line per read). With min_run > 1, only k-mers that lie in a run of at least
// min_run consecutive such k-mers of the read are kept; a k-mer with a character other than
// A, C, G, T, or the end of the read, ends a run.
#include <bits/stdc++.h>
using namespace std;
static int code(char c){ switch(c){case 'A':case 'a':return 0;case 'C':case 'c':return 1;case 'G':case 'g':return 2;case 'T':case 't':return 3;} return -1; }
static vector<string> read_seqs(const string& p){ vector<string> v; ifstream in(p); string l; bool fq=false, first=true; int state=0; bool open=false;
  while(getline(in,l)){ if(!l.empty()&&l.back()=='\r')l.pop_back(); if(first){ if(l.empty()) continue; fq=(l[0]=='@'); first=false; }
    if(fq){ if(state==0){ if(l.empty()) continue; v.emplace_back(); state=1;} else if(state==1){v.back()=l;state=2;} else if(state==2){state=3;} else {state=0;} }
    else { if(!l.empty()&&l[0]=='>'){v.emplace_back();open=true;} else if(open) v.back()+=l; } }
  return v; }
template<class F> static void kmers(const string& s,int k,F f){ uint64_t mask=((1ULL<<(2*k))-1),fw=0,rc=0; int run=0; for(size_t i=0;i<s.size();i++){ int c=code(s[i]); if(c<0){run=0;continue;} fw=((fw<<2)|c)&mask; rc=(rc>>2)|((uint64_t)(3-c)<<(2*(k-1))); if(++run>=k) f(i+1-k,min(fw,rc)); } }
int main(int argc,char**argv){ string refs=argv[1],q=argv[2]; int k=atoi(argv[3]);
  vector<string> reads=read_seqs(q); unordered_map<uint64_t,bool> present;
  for(auto& r: reads) kmers(r,k,[&](size_t,uint64_t x){ present.emplace(x,false); });
  ifstream in(refs); string path; while(getline(in,path)){ if(path.empty()) continue; for(const string& s: read_seqs(path)) kmers(s,k,[&](size_t,uint64_t x){ auto it=present.find(x); if(it!=present.end()) it->second=true; }); }
  long min_run = argc > 6 ? atol(argv[6]) : 1;
  ofstream ko(argv[4]), no(argv[5]);
  for(size_t r=0;r<reads.size();r++){
    vector<pair<size_t,bool>> v; kmers(reads[r],k,[&](size_t pos,uint64_t x){ v.push_back({pos,present[x]}); });
    vector<bool> keep(v.size(),false);
    for(size_t i=0;i<v.size();){ if(!v[i].second){ i++; continue; } size_t e=i+1; while(e<v.size() && v[e].second && v[e].first==v[e-1].first+1) e++; if((long)(e-i)>=min_run) for(size_t j=i;j<e;j++) keep[j]=true; i=e; }
    long n=0; for(size_t i=0;i<v.size();i++) if(keep[i]){ n++; string km=reads[r].substr(v[i].first,k); for(char& c: km) c=toupper(c); ko<<'>'<<r<<'\n'<<km<<'\n'; }
    no<<n<<'\n'; }
  return 0; }
