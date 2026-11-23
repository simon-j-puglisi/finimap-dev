#ifndef _PREF_TAB_H_
#define _PREF_TAB_H_

#include <algorithm>
#include "mapped_array.hh"
#include "packed_ints.hh"
#include <bit>
#include <chrono>
#include <cstdint>
#include <ctime>
#include <fstream>
#include <iostream>
#include <random>
#include <ratio>
#include <set>
#include <stdio.h>
#include <stdlib.h>
#include <unordered_map>
#include <unordered_set>
#include <vector>

using namespace std;
using namespace std::chrono;

#include "Util.hh"


//#define BIN_SEARCH_THRESHOLD 65536
#define BIN_SEARCH_THRESHOLD 128

class PrefTab {
   public:
    PrefTab() : _preflen(0), _nprefs(0), _nfins(0) {
    }

    ~PrefTab() = default;

    PrefTab(const vector<std::string_view> &finstrings, uint64_t preflen) {

       vector<uint64_t> data;
       vector<uint8_t> lengths;
       for(const string_view& s : finstrings){ 
          uint64_t fini = 0;
          for(int i=0;i<s.length();i++){
             //pack a binary representation of the finimizer in the lower bits of fini
             fini |= ((uint64_t)char2bits(s[i])) << ((uint64_t)(2*(31-i)));
          }
          data.push_back(fini);
          lengths.push_back(s.length());
       }

       _nfins = data.size();
       _nprefs = ((uint64_t)1)<<(preflen<<1);
       uint64_t shift = 64-(preflen<<1);
       _preflen = preflen;

       _X.resize(_nprefs+1); //+1 so we can eventually store _n in _X[_nprefs] and save a branch 
       for(uint64_t i=0;i<_nprefs+1;i++) _X[i] = 0;
       //Count number of items in each bucket
       for(uint64_t i=0;i<data.size();i++){
          _X[data[i]>>shift]++;
       }

       //Collect some statistics
       uint64_t maxcount = 0, emptycount = 0, nbigs = 0, nvbigs = 0;
       for(uint64_t i=0;i<_nprefs;i++){
          if(_X[i] > maxcount){ maxcount = _X[i]; }
          if(_X[i] > 128){ nbigs++; }
          if(_X[i] > 1024){ nvbigs++; }
          if(_X[i] == 0){ emptycount++; }
       }
       cerr << "_nprefs: "<<_nprefs<< '\n';
       cerr << "maxcount emptycount: "<< maxcount << ' ' << emptycount << '\n';
       cerr << "avg. items per bucket: "<<data.size()/(_nprefs-emptycount)<<'\n';
       cerr << "nbigs nvbigs: "<< nbigs << ' ' << nvbigs << '\n';

       //Prefix sum the counts in _X, creating pointers into _F.
       //Empty buckets are treated in a special way.
       uint32_t sum = 0;
       uint32_t lastnonempty = 0;
       uint32_t lastnonemptysz = 0;
       for(uint64_t i=0;i<_nprefs+1;i++){
          uint32_t v = _X[i];
          if(v){
             _X[i] = sum; //point to the first item in _X that's in this bucket
             lastnonempty = sum;
             lastnonemptysz = v;
             _X[i] = _X[i] << 1; //0 in LSB indicates non-empty
          }else{
             //empty bucket
             _X[i] = sum + v - 1; //point to the predecessor (in _F) of anything that lands in the bucket
             _X[i] = (_X[i] << 1) | 1; //1 in LSB indicates empty
          }
          sum += v;
       }
       _X[_nprefs] = (_nfins << 1) | 1;

       //Prepare the array of finimizers in _F, removing there prefixes and tagging them with their lengths in the least significant 5 bits
       _F.resize(_nfins);
       for(uint64_t i=0;i<_nfins;i++){
          //shift off prefix and tag with len
          _F[i] = (data[i] << (preflen<<1)) | lengths[i];
       }

       //run some tests
       //first for present elements
       for(uint64_t i=0;i<data.size();i++){
          pair<int64_t,bool> x = getPred(data[i]);
          if(i != x.first && x.second != true){
             cerr << i << ' ' << x.first << ' ' << x.second << '\n';
          }
       }
       cerr << "Test 1 passed.\n";
/*
       //now test for elements not in the set
       for(uint64_t i=1;i<data.size();i++){
          uint64_t missing = data[i-1]+1;
          while(missing < data[i]){
             pair<int64_t,bool> x = getPred(missing);
             if(i-1 != x.first && x.second != false){
                cerr << i << ' ' << x.first << ' ' << x.second << '\n';
             }
             missing++;
          }
       }
       cerr << "Test 2 passed.\n";
*/
    }

    //  p is the index of the predecessor in the set
    //  bool is 1 if the value at index p is equal to key, else 0
    pair<int64_t, bool> inline getPred(uint64_t key) const{
        uint64_t preflenbits = _preflen<<1;
        uint64_t shift = 64-preflenbits; 
        uint64_t p = key>>shift; //p has the prefix of the key of length preflenbits in its least significant bits
        key = key << preflenbits; //shift the prefix off the key
        uint32_t v = _X[p];
        if(v&1){
           //cerr << "Empty bucket\n";
           //empty bucket
           if(p!=0){
              return {v>>1, false};
           }
           return {-1, false};
        }
        //non-empty bucket
        v = v >> 1;
        //cerr << "Non-empty bucket\n";
        uint32_t blen = (_X[p+1]>>1) + (_X[p+1]&1) - v; //compute # items in bucket
        //cerr << "blen: " << blen << '\n';
        //cerr << "v: " << v << '\n';
        uint32_t i = v;
        for(;i<v+blen;i++){
           //cerr << "key: " << key << " curr: " << ((F(i)>>preflenbits)<<preflenbits) <<'\n';
           if(((F(i)>>preflenbits)<<preflenbits)>=key){
              break;
           }
        }
        return {i-(((F(i)>>preflenbits)<<preflenbits)==key), (((F(i)>>preflenbits)<<preflenbits)==key)};
    }

    pair<int64_t, uint64_t> finiLookup(uint64_t key, uint64_t keylen) const{
       uint64_t preflenbits = _preflen<<1;
       //cerr << "preflenbits: " << preflenbits << '\n';
       uint64_t shift = 64-preflenbits;
       uint64_t p = key>>shift; //p has the prefix of the key of length preflenbits in its least significant bits
       //cerr << "p: " << p << '\n';
       key = key << preflenbits; //shift the prefix off the key
       uint32_t v = _X[p];
       if(v&1){
          //cerr << "Empty bucket\n";
          //empty bucket
          return {-1,0};
       }
       //non-empty bucket
       //cerr << "Non-empty bucket\n";
       v = v >> 1;
       //check if we are smaller than everything in our bucket
       if(((F(v)>>preflenbits)<<preflenbits)>key){
          //cerr << "No predecessor in bucket\n";
          return {-1,0};
       }
       uint32_t blen = (_X[p+1]>>1) + (_X[p+1]&1) - v; //compute # items in bucket
       //cerr << "blen: " << blen << '\n';
       //cerr << "v: " << v << '\n';
       uint32_t i = v;
       //cerr << "v v+blen : " << v << ' ' << v+blen << '\n';
       
       if(blen < BIN_SEARCH_THRESHOLD){
          //cerr << "Scanning search: " << blen << " items \n";
          for(;i<v+blen;i++){
             //cerr << "Item: " << (i-v) << ", key: " << key << " curr: " << ((F(i)>>preflenbits)<<preflenbits) <<'\n';
             if(((F(i)>>preflenbits)<<preflenbits)>=key){
                //cerr << "breaking\n";
                break;
             }
          }
          //cerr << "Normal i: "<<i<<'\n';
          //cerr << "Equal to key? " << (((F(i)>>preflenbits)<<preflenbits)==key)<<'\n';
       }
       else{
          //cerr << "Binary search...\n";
          uint32_t low = v;
          uint32_t high = v+blen-1;
          while(low <= high){
             uint32_t mid = (low + high) >> 1;
             uint64_t midVal = ((F(mid)>>preflenbits)<<preflenbits);
             if (midVal < key){
                low = mid + 1;
                i = low; //key not found, low is the index of the first element greater than key
             }else if (midVal > key){
                high = mid - 1;
                i = low; //key not found, low is the index of the first element greater than key
             }else{
                i = mid; break; // key found at mid
             }
          }
          //cerr << "Binary i: "<<i<<'\n';
       }
       //cerr << "i: " << i << '\n';
       if(i == (v + blen)){
          i--; //i now points to the predecessor
       }else{
          i -= (1 - (((F(i)>>preflenbits)<<preflenbits)==key)); //i now points to the predecessor
       }
       //cerr << "Adjust i: " << i << '\n';
       //compute the LCP between the query and the predecessor
       uint64_t el = (F(i) & (((uint64_t)1)<<preflenbits)-((uint64_t)1)); //the length of the finimizer
       uint64_t pred = ((F(i)>>preflenbits)<<preflenbits);
       //uint64_t lcp = (__builtin_clzll((F(i)>>preflenbits) ^ (~(key >> preflenbits))) - preflenbits)/2;
       uint64_t lcp = countl_zero(~((pred) ^ (~key)));
       lcp = _preflen + (lcp >> 1); //convert from bit matches to symbol matches
       //cerr << "pred el lcp: " << i << ' ' << el << ' ' << lcp << ' ' << ((lcp >= el) ? '*' : ' ') << '\n';
       if(lcp >= el && el <= keylen){
          return {i,el};
       }
       return {-1,0};
       //{i-(((F(i)>>preflenbits)<<preflenbits)==key), (((F(i)>>preflenbits)<<preflenbits)==key)};
    }

    // Finimizer i: its characters after the prefix, from the top bit, and its length in the
    // lowest bits. In index format 4 these are stored packed: the characters after the prefix
    // (62 - 2 * preflen bits, as the last two bits are always zero) and the length (5 bits).
    uint64_t F(uint64_t i) const
    {
       if(!_packed) return _F[i];
       const uint64_t v = _Fp[i];
       return ((v >> 5) << _fshift) | (v & 31);
    }

    // Whether index format 4 stores the finimizers packed with this prefix length
    static bool packs(uint64_t preflen) { return preflen >= 5 && preflen <= 30; }

    uint64_t sizeInBytes() const{
       uint64_t sz = 0;
       sz += 3*sizeof(uint64_t);
       sz += _nprefs*sizeof(uint32_t);
       sz += _packed ? _Fp.size_in_bytes() : _nfins*sizeof(uint64_t);
       return sz;
    }

    int64_t serialize(std::ostream& os, const uint64_t format) const{
       cerr << "PrefTab.serialize()...";
       int64_t written = 0;
       os.write((char *)&_preflen, sizeof(uint64_t));
       os.write((char *)&_nprefs, sizeof(uint64_t));
       os.write((char *)&_nfins, sizeof(uint64_t));
       os.write((char *)_X.data(), (_nprefs+1)*sizeof(uint32_t));
       written += 3*sizeof(uint64_t);
       written += _nprefs*sizeof(uint32_t);
       if(format >= 4 && packs(_preflen)){
          const uint64_t shift = 2*_preflen + 2;
          vector<uint64_t> packed(_nfins);
          for(uint64_t i=0;i<_nfins;i++){
             const uint64_t f = F(i);
             if((f & 31) != (f & ((((uint64_t)1)<<(2*_preflen))-1)) || ((f >> 5) << 5 & ((((uint64_t)1)<<shift)-1)) != 0)
                throw std::logic_error("PrefTab: a finimizer does not have the expected layout");
             packed[i] = ((f >> shift) << 5) | (f & 31);
          }
          PackedInts p;
          p.build(packed.data(), packed.size(), 67 - 2*_preflen);
          p.serialize(os);
          written += p.size_in_bytes() + 24;
       }else{
          os.write((char *)_F.data(), _nfins*sizeof(uint64_t));
          written += _nfins*sizeof(uint64_t);
       }
       cerr << '\n';
       return written;
    }

    void load(std::istream& is, const uint64_t format){
       cerr << "PrefTab.load()...";
       is.read((char *)&_preflen, sizeof(uint64_t));
       is.read((char *)&_nprefs, sizeof(uint64_t));
       is.read((char *)&_nfins, sizeof(uint64_t));
       cerr << "preflen nprefs nfins: "<<_preflen<<' '<<_nprefs<<' '<<_nfins<<'\n';
       if(_preflen < 1 || _preflen > 31 || _nprefs != (((uint64_t)1) << (2*_preflen)))
          throw std::runtime_error("The index is inconsistent (long finimizers)");
       _X.load(is, _nprefs+1);
       _packed = format >= 4 && packs(_preflen);
       _fshift = 2*_preflen + 2;
       if(_packed){
          _Fp.load(is);
          if(is && (_Fp.size() != _nfins || _Fp.bits() != 67 - 2*_preflen))
             throw std::runtime_error("The index is inconsistent (long finimizers)");
       }else{
          _F.load(is, _nfins);
       }
    }
   

   
    //instrumenting statistics
 
   private:
    uint64_t _preflen; //prefix length used for lookup table
    uint64_t _nprefs; //the number of prefixes == 1<<(_preflen<<1);
    uint64_t _nfins; //the number of finimizers == |_F|
    Array<uint32_t> _X;
    Array<uint64_t> _F;   // index formats 1 to 3, and while building
    PackedInts _Fp;       // index format 4
    bool _packed = false;
    uint64_t _fshift = 0;
};

#endif
