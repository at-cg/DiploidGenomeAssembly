#include<bits/stdc++.h>
#include <divsufsort.h>
using namespace std;
#define ll long long int
#define pb push_back
#define rb pop_back
#define ti tuple<int, int, int>
#define pii pair<int, int>
#define pli pair<ll, int>
#define pll pair<ll, ll>
#define mp make_pair
#define mt make_tuple
#define F first
#define S second

using namespace std;
const int maxlogdepth = 30;

struct node{
    int l, r, startIndex, depth; // depth is before starting index of this node
    ll *colocated = new ll[4];
    vector<node*> p;
    vector<node*> child;
    node(int left, int right, int start, int dist, node *parent):l(left), r(right), startIndex(start), depth(dist){
        p.resize(maxlogdepth);
        p[0] = parent;
        for(int i = 0; i < 4; i++){
            colocated[i] = 0;
        } 
    }
};

class SuffixTree{// diploid
public:
    node *root;
    const char *text;
    int *SA;
    int *lcp;
    int *closestloci;
    int len, l0, r0, l1, r1, H;
    bool in_range(int x, vector<int> v){
        for(int i = 0; i < v.size() - 1; i += 2){
            if(x >= v[i] && x <= v[i + 1])return true;
        }
        return false;
    }
    SuffixTree(const char *input_text, int *suffix_array, int *lcp_array):text(input_text), SA(suffix_array), lcp(lcp_array){
        root = new node(-1, -1, -1, -1, NULL);
        len = strlen(input_text);
        int double_length = (len - 2) / 2;
        assert(double_length % 2 == 0);
        int single_length = double_length / 2;

        // only suffixes with starting point in these ranges must be added
        l0 = 0; r0 = double_length / 2 - 1; H = r0 + 1;
        l1 = double_length + 1; r1 = l1 + r0;

        closestloci = new int[single_length];
        for(int i = 0; i < single_length; i++)closestloci[i] = -1;
        int pst = -1;
        for(int i = double_length - 1; i >= 0; i--){
            if(text[i] != text[i + 2 * H + 1]){
                pst = i;
                if(pst >= H)pst -= H;
                continue;
            }
            int loc = i;
            if(loc >= H)loc -= H;
            if(closestloci[loc] != -1)continue;
            closestloci[loc] = pst;
        }

        // finding the first suffix to be added
        int start = 0;
        while(!in_range(SA[start], {l0, r0, l1, r1})){++start;}
        // cout << "start :" << start << endl;
        node *first = new node((int)SA[start], SA[start] + H - 1, SA[start], 0, root);
        root->child.pb(first);

        // Adding other valid suffixes
        node *itr = first; 
        ++start; 
        int LCP = start < len ? min(H, lcp[start]) : -1; 
        
        while(start < len){ // iterating over the SA
            while(!in_range(SA[start], {l0, r0, l1, r1})){
                ++start; 
                if(start >= len)break;
                LCP = min(LCP, lcp[start]);
            };
            if(start >= len)break;
            if(LCP == H){ // haploid
                if(itr->startIndex >= 2 * H + 1 && SA[start] < 2 * H + 1)itr->startIndex = SA[start];
            }else if(LCP == 0){
                node *child = new node(SA[start] + LCP, SA[start] + H - 1, SA[start], 0, root);
                root->child.pb(child);
                itr = child;
            }else{
                int moveUp = H - LCP;
                while(itr->r - itr->l + 1 <= moveUp){// finding the correct parent by traversing up from the previous leaf
                    moveUp -= (itr->r - itr->l + 1);
                    itr = itr->p[0];
                }
                if(moveUp == 0){
                    node *child = new node(SA[start] + LCP, SA[start] + H - 1, SA[start], itr->depth + itr->r - itr->l + 1, itr);
                    itr->child.pb(child);
                    itr = child;  
                }else{
                    node *child1 = new node(itr->r - moveUp + 1, itr->r, itr->startIndex, itr->depth + itr->r - itr->l + 1 - moveUp, itr);
                    child1->child = itr->child;
                    node *child2 = new node(SA[start] + LCP, SA[start] + H - 1, SA[start], itr->depth + itr->r - itr->l + 1 - moveUp, itr);
                    itr->r = itr->r - moveUp;
                    itr->child.clear();
                    itr->child.pb(child1); itr->child.pb(child2);                
                    itr = child2;
                }
            }
            ++start; 
            if(start >= len)break;
            LCP = min(H, lcp[start]);
        }
    }
};

map<char, int> char_map = {{'A', 0}, {'C', 1}, {'G', 2}, {'T', 3}};

ll* doublerepeat(SuffixTree *ST, node *u, map<int, ll> &m, bool *reduce, int depth){
    // precomputation for finding the parent
    for(int i = 1; i < maxlogdepth; i++){
        if(u->p[i - 1] == ST->root)u->p[i] = ST->root;
        else u->p[i] = u->p[i - 1]->p[i - 1];
    }

    ll *freqU = new ll[4]; 
    for(int i = 0; i < 4; i++){
        freqU[i] = 0;
    } 
    if(u->child.empty()){// leaf node
        int baseindex = u->startIndex;
        if(baseindex >= 2 * ST->H + 1)baseindex -= 2 * ST->H + 1;
        int beforeindex = (baseindex - 1 + ST->H) % ST->H;
        if(u->startIndex >= 2 * ST->H + 1)beforeindex += 2 * ST->H + 1;
        int base = char_map[ST->text[beforeindex]];
        freqU[base] = 1;
        if(ST->closestloci[baseindex] != -1 && reduce[baseindex]){ // contributes to colocated repeats
            int searchdepth = (ST->closestloci[baseindex] - baseindex + ST->H) % ST->H;
            node *ptr = u;
            for(int i = maxlogdepth - 1; i >= 0; i--){
                if(ptr->p[i] == ST->root)continue;
                if(ptr->p[i]->depth >= searchdepth){ // There will be a node with exactly that value of depth
                    ptr = ptr->p[i];
                }
            }
            ptr->colocated[base]++;
        }
        if(!reduce[baseindex])reduce[baseindex] = true;

        return freqU;
    }

    ll tsumU = 0;  
    ll contribution = 0;
    for(node *v : u->child){
        ll *freqV = doublerepeat(ST, v, m, reduce, depth + v->r - v->l + 1); // unordered
        
        // removing colocated double repeats
        for(int i = 0; i < 4; i++){
            freqV[i] -= v->colocated[i];  
        }
        
        ll tsumV = 0;
        for(int i = 0; i < 4; i++)tsumV += freqV[i];

        // find contribution
        for(int i = 0; i < 4; i++){
            // contribution of 1 from v
            contribution += (tsumU - freqU[i]) * freqV[i];
        }

        // update freqU
        for(int i = 0; i < 4; i++){
            freqU[i] += freqV[i];
            tsumU += freqV[i];
        }
    }
    
    if(depth != 0 && contribution != 0){
        if(m.find(depth) == m.end())m[depth] = contribution;
        else m[depth] += contribution;
    }
    
    u->p.clear();

    return freqU;
}

void doublerepeat_bruteHaploid(string& s, map<int, ll> &repeatstats_brute){
    int l = s.length();
    s += s;
    for(int i = 0; i < l - 1; i++){
        for(int j = i + 1; j < l; j++){
            if(s[(i - 1 + l) % l] == s[j - 1])continue; // can't be maximal
            for(int sl = l; sl >= 1; sl--){
                if(s[i + sl] == s[j + sl])continue; // can't be maximal
                string s1 = s.substr(i, sl), s2 = s.substr(j, sl);
                if(s1 == s2){
                    if(repeatstats_brute.find(sl) == repeatstats_brute.end())repeatstats_brute[sl] = 1;
                    else repeatstats_brute[sl] = repeatstats_brute[sl] + 1;
                    break; // other smaller lengths won't be maximal
                }
            }
        }
    }
}

void doublerepeat_bruteDiploid(string H0, string H1, map<int, ll> &repeatstats_brute){
    int l = H0.length();
    doublerepeat_bruteHaploid(H0, repeatstats_brute);
    doublerepeat_bruteHaploid(H1, repeatstats_brute);
    assert(H0.length() == 2 * l);

    for(int i = 0; i < l; i++){
        for(int j = 0; j < l; j++){ // unordered
            if(i == j)continue;
            if(H0[(i - 1 + l) % l] == H1[j - 1])continue; // can't be maximal
            for(int sl = l; sl >= 1; sl--){
                if(H0[i + sl] == H1[j + sl])continue; // can't be maximal
                string s1 = H0.substr(i, sl), s2 = H1.substr(j, sl);
                if(s1 == s2){
                    if(repeatstats_brute.find(sl) == repeatstats_brute.end())repeatstats_brute[sl] = 1;
                    else repeatstats_brute[sl] = repeatstats_brute[sl] + 1;
                    break; // other smaller lengths won't be maximal
                }
            }
        }
    }
}

fstream f;

int main(int argc, char* argv[])
{   
    if(argc != 4){
        return 1;
    }

    const char* mat_fastafile = argv[1];
    const char* pat_fastafile = argv[2];
    const char* outputfile = argv[3];

    freopen(outputfile, "w", stdout);

    f.open(mat_fastafile, ios::in);
    string s0;
    if(f.is_open()){
        while(getline(f, s0)){
            if(s0[0] == '>')continue;
            break;
        }
        f.close();
    }

    f.open(pat_fastafile, ios::in);
    string s1;
    if(f.is_open()){
        while(getline(f, s1)){
            if(s1[0] == '>')continue;
            break;
        }
        f.close();
    }
    
    // cerr << "STAGE 1\n";

    // s = s.substr(0, 200000000); s1 = s1.substr(0, 200000000);
    vector<string> v; v.pb(s0); v.pb(s1); // diploid
    string query = v[0] + v[0] + '0' + v[1] + v[1] + '1';
    const char *text = query.c_str();
    int len = strlen(text);

    int *SA = (int *)malloc(len * sizeof(int));
    int *iSA = (int *)malloc(len * sizeof(int));
    int *lcp = (int *)malloc((len - 1) * sizeof(int));

    // cerr << "STAGE 2\n";

    // SA
    divsufsort((unsigned char *)text, SA, len);
    
    // iSA
    for(int i = 0; i < len; i++){
        iSA[SA[i]] = i;
    }

    // cerr << "STAGE 3\n";

    // LCP
    int val = 0; 
    for(int i = 0; i < len; i++){
        if(iSA[i] == 0)continue;
        int start = max(val - 1, 0);
        while(text[i + start] == text[SA[iSA[i] - 1] + start]){
            ++start;
        }
        val = start;
        lcp[iSA[i]] = val;
    }
    free(iSA);

    // cerr << "STAGE 4\n";

    // ST
    SuffixTree *ST = new SuffixTree(text, SA, lcp);
    free(SA);
    
    // cerr << "STAGE 5\n";

    map<int, ll> repeatstats; 
    bool *reduce = new bool[ST->H];
    for(int i = 0; i < ST->H; i++)reduce[i] = false;
    for(node *v : ST->root->child){
        doublerepeat(ST, v, repeatstats, reduce, v->r - v->l + 1);
    }
    free(ST);

    // cerr << "STAGE 6\n";

    for(auto stats : repeatstats){
        cout << stats.first << " " << stats.second << endl;
    }

    map<int, ll> repeatstats_brute;
    doublerepeat_bruteDiploid(s0, s1, repeatstats_brute);

    assert(repeatstats == repeatstats_brute);
} 
