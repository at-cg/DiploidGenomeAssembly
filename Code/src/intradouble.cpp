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

struct node{
    int l, r; // depth is before starting index of this node
    vector<node*> child;
    node* p;
    node(int left, int right, node* parent):l(left), r(right), p(parent){
    }
};

class SuffixTree{// haploid
public:
    node *root;
    int *SA;
    int *lcp;
    int len, l0, r0, H;
    bool in_range(int x, vector<int> v){
        for(int i = 0; i < v.size() - 1; i += 2){
            if(x >= v[i] && x <= v[i + 1])return true;
        }
        return false;
    }
    SuffixTree(int *suffix_array, int *lcp_array):SA(suffix_array), lcp(lcp_array){
        root = new node(-1, -1, NULL);
        len = sizeof(suffix_array) / sizeof(suffix_array[0]);
        int single_length = (len - 1) / 2;

        // only suffixes with starting point in these ranges must be added
        l0 = 0; r0 = single_length - 1; H = r0 + 1;

        // finding the first suffix to be added
        int start = 0;
        while(!in_range(SA[start], {l0, r0})){++start;}
        // cout << "start :" << start << endl;
        node *first = new node((int)SA[start], SA[start] + H - 1, root);
        root->child.pb(first);

        // Adding other valid suffixes
        node *itr = first; 
        ++start; 
        int LCP = start < len ? min(H, lcp[start]) : -1; 
        
        while(start < len){ // iterating over the SA
            while(!in_range(SA[start], {l0, r0})){
                ++start; 
                if(start >= len)break;
                LCP = min(LCP, lcp[start]);
            };
            if(start >= len)break;
            assert(LCP < H);
            if(LCP == 0){
                node *child = new node(SA[start], SA[start] + H - 1, root);
                root->child.pb(child);
                itr = child;
            }else{
                int moveUp = H - LCP;
                while(itr->r - itr->l + 1 <= moveUp){// finding the correct parent by traversing up from the previous leaf
                    moveUp -= (itr->r - itr->l + 1);
                    itr = itr->p;
                }
                if(moveUp == 0){
                    node *child = new node(SA[start] + LCP, SA[start] + H - 1, itr);
                    itr->child.pb(child);
                    itr = child;  
                }else{
                    node *child1 = new node(itr->r - moveUp + 1, itr->r, itr);
                    child1->child = itr->child;
                    node *child2 = new node(SA[start] + LCP, SA[start] + H - 1, itr);
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

ll doublerepeat(SuffixTree *ST, node *u, map<int, ll> &m, int depth){
    if(u->child.empty()){// leaf node
        return 1;
    }

    int cntU = 0;  
    ll contribution = 0;
    for(node *v : u->child){
        int cntV = doublerepeat(ST, v, m, depth + v->r - v->l + 1); // unordered
        contribution += (ll)cntU * cntV;
        cntU += cntV;
    }
    
    if(depth != 0 && contribution != 0){
        if(m.find(depth) == m.end())m[depth] = contribution;
        else m[depth] += contribution;
    }

    return cntU;
}

void doublerepeat_bruteHaploid(string s, map<int, ll> &repeatstats_brute){
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

fstream f;

int main(int argc, char* argv[])
{   
    if(argc != 3){
        return 1;
    }

    const char* fastafile = argv[1];
    const char* outputfile = argv[2];

    freopen(outputfile, "w", stdout);

    f.open(fastafile, ios::in);
    string s;
    if(f.is_open()){
        while(getline(f, s)){
            if(s[0] == '>')continue;
            break;
        }
        f.close();
    } // haploid
    string query = s + s + '1'; // Sentinal symbol for building a generalised suffix tree
    const char *text = query.c_str();
    int len = strlen(text);

    int *SA = (int *)malloc(len * sizeof(int));
    int *iSA = (int *)malloc(len * sizeof(int));
    int *lcp = (int *)malloc((len - 1) * sizeof(int));

    // SA
    divsufsort((unsigned char *)text, SA, len);
    
    // iSA
    for(int i = 0; i < len; i++){
        iSA[SA[i]] = i;
    }

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

    // ST
    SuffixTree *ST = new SuffixTree(SA, lcp);
    free(SA);
    
    map<int, ll> repeatstats; 
    for(node *v : ST->root->child){
        doublerepeat(ST, v, repeatstats, v->r - v->l + 1);
    }
    free(ST);

    map<int, ll> repeatstats_brute;
    doublerepeat_bruteHaploid(s, repeatstats_brute);

    // assert(repeatstats == repeatstats_brute);

    for(auto stats : repeatstats){
        cout << stats.first << " " << stats.second << endl;
    }
} 
