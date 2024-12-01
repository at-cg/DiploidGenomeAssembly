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
const int threshold = 0;

struct node{
    int l, r, startIndex, depth; // depth is before starting index of this node
    node* parent;
    vector<node*> child;
    node(int left, int right, int start, int dist, node* p):l(left), r(right), startIndex(start), depth(dist), parent(p){}
};

class SuffixTree{// diploid
public:
    node *root;
    const char *text;
    int *SA;
    int *lcp;
    int *closestloci;
    int len, l0, r0, H;
    bool in_range(int x, vector<int> v){
        for(int i = 0; i < v.size() - 1; i += 2){
            if(x >= v[i] && x <= v[i + 1])return true;
        }
        return false;
    }
    SuffixTree(const char *input_text, int *suffix_array, int *lcp_array):text(input_text), SA(suffix_array), lcp(lcp_array){
        root = new node(-1, -1, -1, -1, NULL);
        len = strlen(input_text);
        assert((len - 1) % 2 == 0);
        
        // only suffixes with starting point in these ranges must be added
        l0 = 0; r0 = (len - 1) / 2 - 1; H = r0 + 1;

        // finding the first suffix to be added
        int start = 0;
        while(!in_range(SA[start], {l0, r0})){++start;}
        // cout << "start :" << start << endl;
        node *first = new node((int)SA[start], SA[start] + H - 1, SA[start], 0, root);
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
            if(LCP == H){ // haploid
            }else if(LCP == 0){
                node *child = new node(SA[start] + LCP, SA[start] + H - 1, SA[start], 0, root);
                root->child.pb(child);
                itr = child;
            }else{
                int moveUp = H - LCP;
                while(itr->r - itr->l + 1 <= moveUp){// finding the correct parent by traversing up from the previous leaf
                    moveUp -= (itr->r - itr->l + 1);
                    itr = itr->parent;
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

vector<int>* interleavedrepeat(SuffixTree *ST, node *u, map<int, vector<pii>*> &m, int depth){
    if(u->child.empty()){// leaf node
        vector<int> *v = new vector<int>; v->pb(u->startIndex);
        return v;
    }

    vector<int> *contU = new vector<int>;
    for(node *v : u->child){
        vector<int> *contV = interleavedrepeat(ST, v, m, depth + v->r - v->l + 1); // unordered
        if(depth >= threshold){
            for(int p1 : *contV){
                for(int p2 : *contU){
                    if(ST->text[(p1 - 1 + ST->H) % ST->H] == ST->text[(p2 - 1 + ST->H) % ST->H])continue;
                    if(m.find(depth) == m.end()){
                        vector<pii> *cont = new vector<pii>;
                        m[depth] = cont;
                    }
                    m[depth]->pb(mp(min(p1, p2), max(p1, p2)));
                }
            }
        }
        for(int p : *contV){
            contU->pb(p);
        }
    }
    return contU;
}


fstream f;

int main()
{   
    // f.open("/home/daanish/data/HG002/chr1_mat.fasta", ios::in);
    
    // string s;
    // if(f.is_open()){
    //     while(getline(f, s)){
    //         if(s[0] == '>')continue;
    //         break;
    //     }
    //     f.close();
    // }
    string s = "AAGTTGAATT";
    vector<string> v; v.pb(s); // haploid
    string query = v[0] + v[0] + '0';
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
    SuffixTree *ST = new SuffixTree(text, SA, lcp);
    free(SA);
    
    map<int, vector<pii>*> repeatstats; 
    for(node *v : ST->root->child){
        interleavedrepeat(ST, v, repeatstats, v->r - v->l + 1);
    }
    free(ST);

    int ans = -1;

    for(auto it1 = repeatstats.begin(); it1 != repeatstats.end(); it1++){
        vector<pii> *v1 = it1->second;
        // for(pii p : *v1){
        //     cout << "(" << p.first << "," << p.second << ")" << endl;
        // }
        for(auto it2 = repeatstats.begin(); it2 != repeatstats.end(); it2++){
            if(it2->first < it1->first)continue;
            vector<pii> *v2 = it2->second;
            for(pii p1 : *v1){
                for(pii p2 : *v2){
                    if(p1.first < p2.first && p2.first < p1.second && p1.second < p2.second){
                        ans = max(ans, it1->first);
                    }
                    if(p2.first < p1.first && p1.first < p2.second && p2.second < p1.second){
                        ans = max(ans, it1->first);
                    }
                }
            }
        }
    }

    cout << ans << endl;
} 
