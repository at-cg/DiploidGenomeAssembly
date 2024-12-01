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

fstream f;
int *closestloci_forward, *closestloci_backward, H;

bool in_range(int x, vector<int> v){
    for(int i = 0; i < v.size() - 1; i += 2){
        if(x >= v[i] && x <= v[i + 1])return true;
    }
    return false;
}

int distanceforward(int x){
    if(x >= 2 * H + 1)x -= 2 * H + 1;
    assert(x < H);
    return (closestloci_forward[x] - x + H) % H;
}

bool colocated(int x1, int x2){
    if(abs(x2 - x1) == 2 * H + 1)return true;
    return false;
}

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
    string s0, line;
    if(f.is_open()){
        while(getline(f, line)){
            if(line[0] == '>')continue;
            s0 += line;
        }
        f.close();
    }

    f.open(pat_fastafile, ios::in);
    string s1;
    if(f.is_open()){
        while(getline(f, line)){
            if(line[0] == '>')continue;
            s1 += line;
        }
        f.close();
    }
    assert(s0.length() == s1.length());
    cout << "Length of each haplotype: " << s0.length() << endl;
    
    string query = s0 + s0 + '0' + s1 + s1 + '1';
    const char *text = query.c_str();
    const int len = strlen(text);
    const int double_length = (len - 2) / 2;
    assert(double_length % 2 == 0);
    H = double_length / 2;
    const int l0 = 0, r0 = double_length / 2 - 1, l1 = double_length + 1, r1 = l1 + r0;

    closestloci_forward = new int[H];
    for(int i = 0; i < H; i++)closestloci_forward[i] = -1;
    int pst = -1, maxgap = 0;
    for(int i = double_length - 1; i >= 0; i--){
        int loc = i;
        if(loc >= H)loc -= H;
        if(text[i] != text[i + 2 * H + 1]){
            if(pst != -1){
                int ppst = i >= H ? i - H : i;
                maxgap = max(maxgap, abs(pst - ppst) - 1);
            }
            pst = i;
            if(pst >= H)pst -= H;
        }
        if(closestloci_forward[loc] != -1)continue;
        closestloci_forward[loc] = pst;
    }
    cout << "Maximum gap between two adjacent heterozygous loci: " << maxgap << endl;
   
    closestloci_backward = new int[H];
    for(int i = 0; i < H; i++)closestloci_backward[i] = -1;
    pst = -1;
    for(int i = 0; i <= double_length - 1; i++){
        int loc = i;
        if(loc >= H)loc -= H;
        if(text[i] != text[i + 2 * H + 1]){
            pst = i;
            if(pst >= H)pst -= H;
        }
        if(closestloci_backward[loc] != -1)continue;
        closestloci_backward[loc] = pst;
    }

    for(int i = 0; i < H; i++)cout << closestloci_forward[i] << " ";
    cout << endl;
    for(int i = 0; i < H; i++)cout << closestloci_backward[i] << " ";
    cout << endl;

    int *SA = (int *)malloc(len * sizeof(int));
    int *iSA = (int *)malloc(len * sizeof(int));
    int *lcp = (int *)malloc(len * sizeof(int));

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

    int offset = 50;
    
    // max length double repeat
    int maxrepeat = 0; 
    int maxrepeat_mat = 0, maxrepeat_pat = 0; 
    int maxtriplerepeat_mat = 0, maxtriplerepeat_pat = 0;
    int minlength_forwellbridging = H;
    // stats of double repeats for greedy
    int bothcover = 0; 
    int onecover = 0; 
    int nonecover = 0; 
    for(int i = 0; i < len - 1; i++){
        bool mat = false;
        if(!in_range(SA[i], {l0, r0, l1, r1}))continue;
        if(in_range(SA[i], {l0, r0}))mat = true;
        else mat = false;
        int val = H;
        int dist1 = distanceforward(SA[i]);
        int cntmat = 0, cntpat = 0;
        for(int j = i + 1; j < min(len, i + offset); j++){
            val = min(val, lcp[j]);
            
            if(!in_range(SA[j], {l0, r0, l1, r1}))continue;
            if(colocated(SA[i], SA[j]))continue;
            
            maxrepeat = max(maxrepeat, val);
            
            if(mat && in_range(SA[j], {l0, r0})){
                maxrepeat_mat = max(maxrepeat_mat, val);
                cntmat += 1;
                if(cntmat >= 2){
                    maxtriplerepeat_mat = max(maxtriplerepeat_mat, val);
                }
            }
            if(!mat && in_range(SA[j], {l1, r1})){
                maxrepeat_pat = max(maxrepeat_pat, val);
                cntpat += 1;
                if(cntpat >= 2){
                    maxtriplerepeat_pat = max(maxtriplerepeat_pat, val);
                }
            }
            
            int dist2 = distanceforward(SA[j]);
            if(dist1 > val && dist2 > val){// both copies do not cover the locus
                nonecover = max(nonecover, val);

                // int gap1 = findgap(SA[i] + val, closestloci_backward[SA[i]]);
            }else if(dist1 <= val && dist2 <= val){// both copies cover the locus
                bothcover = max(bothcover, val);
            }else{// exactly one copy covers the locus
                onecover = max(onecover, val);
            }
        }
    } 
    cout << "Maximum length of double repeat: " << maxrepeat << endl;
    cout << "Maximum length of double repeat for maternal haplotype: " << maxrepeat_mat << endl;
    cout << "Maximum length of double repeat for paternal haplotype: " << maxrepeat_pat << endl;
    cout << "Maximum length of triple repeat for maternal haplotype: " << maxtriplerepeat_mat << endl;
    cout << "Maximum length of triple repeat for paternal haplotype: " << maxtriplerepeat_pat << endl;
    cout << "Maximum length of double repeat where both copies cover the het loci: " << bothcover << endl;
    cout << "Maximum length of double repeat where exactly one copy covers the het loci: " << onecover << endl;
    cout << "Maximum length of double repeat where both copies do not cover the het loci: " << nonecover << endl;
} 
