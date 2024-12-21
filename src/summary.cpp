#include<bits/stdc++.h>
#include<filesystem>
#include <divsufsort.h>
#define ll long long int
#define pb push_back
#define rb pop_back
#define ti tuple<int, int, int>
#define pii pair<int, int>
#define piii pair<int, pii>
#define pli pair<ll, int>
#define pll pair<ll, ll>
#define mp make_pair
#define mt make_tuple
#define F first
#define S second

using namespace std;
namespace fs = std::filesystem;

fstream f;
int *closestloci_forward, *closestloci_backward, H, l0, r0, l1, r1;
const int offset = 1000;
    
bool in_range(int& x, vector<int> v){
    for(int i = 0; i < v.size() - 1; i += 2){
        if(x >= v[i] && x <= v[i + 1])return true;
    }
    return false;
}

int distance(int x, int y){// end points inclusive, order of x, y matters (distance from x, x + 1, x + 2 -> y)
    return (y - x + H) % H + 1;
}

int distanceloci(int x){
    if(x >= 2 * H + 1)x -= 2 * H + 1;
    return distance(closestloci_backward[x], closestloci_forward[x]) - 2;
}

int distanceforward(int x){
    if(x >= 2 * H + 1)x -= 2 * H + 1;
    assert(x < H);
    return distance(x, closestloci_forward[x]);
}

bool colocated(int& x1, int& x2){
    if(abs(x2 - x1) == 2 * H + 1)return true;
    return false;
}

char find(int& x, string& s0, string& s1){
    if(in_range(x, {l1, r1})){
        x -= 2 * H + 1;
        return s1[(x - 1 + H) % H];
    }else{
        return s0[(x - 1 + H) % H];
    }
}
bool left_maximal(int x1, int x2, string& s0, string& s1){
    return find(x1, s0, s1) != find(x2, s0, s1);
}

int gap(int x, int& val){
    if(x >= 2 * H + 1)x -= 2 * H + 1;
    return min(distance(closestloci_backward[x], x + val - 1), distance(x, closestloci_forward[x])) + 1;
}

piii find_interleaved(vector<piii>& v, string s){ // intra-repeat calculation only
    s += s;
    sort(v.begin(), v.end(), [](const auto& a, const auto& b){return a.first > b.first;});
    int val1 = 0, val2 = 0; int freq = 0;
    for(int i = 0; i < v.size() - 1; i++){
        pii& p1 = v[i].second;
        for(int j = i + 1; j < v.size() && v[j].first >= val1; j++){
            pii& p2 = v[j].second;
            if(s.substr(p2.first, v[j].first) == s.substr(p1.first, v[j].first))continue; // will account for triple repeat 
            if((p1.first < p2.first && p2.first < p1.second && p1.second < p2.second) || (p2.first < p1.first && p1.first < p2.second && p2.second < p1.second)){
                if(v[j].first > val1){
                    freq = 1;
                    val1 = v[j].first; val2 = v[i].first;
                }else if(v[j].first == val1){
                    freq++;
                    val2 = min(val2, v[i].first);
                }
            }
        }
    }
    return {freq, {val1, val2}};
}

// len, (sum, freq)
vector<pii> find_interleaved_stats(vector<piii>& v, string s, pii p){ // intra-repeat calculation only
    s += s;
    vector<pii> o;
    int sum = p.first + p.second;
    for(int i = 0; i < v.size() - 1 && v[i].first >= p.second; i++){
        pii& p1 = v[i].second;
        for(int j = i + 1; j < v.size() && v[j].first + v[i].first >= sum; j++){
            pii& p2 = v[j].second;
            if(s.substr(p2.first, v[j].first) == s.substr(p1.first, v[j].first))continue; // will account for triple repeat 
            if((p1.first < p2.first && p2.first < p1.second && p1.second < p2.second) || (p2.first < p1.first && p1.first < p2.second && p2.second < p1.second)){
                o.pb({v[j].first, v[i].first});
            }
        }
    }
    return o;
}

piii find_specialinterleaved(vector<piii>& v){
    sort(v.begin(), v.end(), [](const auto& a, const auto& b){return a.first > b.first;});
    int val1 = 0, val2 = 0, freq = 0;
    for(int i = 0; i < v.size() - 1; i++){
        pii& p1 = v[i].second;
        for(int j = i + 1; j < v.size() && v[j].first >= val1; j++){
            pii& p2 = v[j].second;
            if(p2.first - p1.first == p2.second - p1.second){
                if(v[j].first > val1){
                    val1 = v[j].first; val2 = v[i].first;
                    freq = 1;
                }else if(v[j].first == val1){
                    val2 = min(val2, v[i].first);
                    freq++;
                }
            }
        }
    }
    return {freq, {val1, val2}};
}

vector<pii> find_specialinterleaved_stats(vector<piii>& v, pii p){
    sort(v.begin(), v.end(), [](const auto& a, const auto& b){return a.first > b.first;});
    vector<pii> o;
    int sum = p.first + p.second;
    for(int i = 0; i < v.size() - 1 && v[i].first >= p.second; i++){
        pii& p1 = v[i].second;
        for(int j = i + 1; j < v.size() && v[j].first + v[i].first >= sum; j++){
            pii& p2 = v[j].second;
            if(p2.first - p1.first == p2.second - p1.second){
                o.pb({v[j].first, v[i].first});
            }
        }
    }
    return o;
}

void printgap(string outputpath, int chk){
    freopen(outputpath.c_str(), "w", stdout);
    for(int i = 0; i < H; i++){
        if(closestloci_forward[i] == i){
            int val = distance(i, closestloci_forward[(i + 1) % H]) - 2;
            assert(val <= chk);
            cout << val << endl;
        }
    }   
}

void print(vector<pii> v, string outputpath){
    freopen(outputpath.c_str(), "w", stdout);
    for(pii p : v){
        cout << p.first << " " << p.second << endl;
    }
}

int main(int argc, char* argv[])
{   
    if(argc != 6){
        return 1;
    }

    const char* mat_fastafile = argv[1];
    const char* pat_fastafile = argv[2];
    string outputpath = argv[3];
    int threshold = stoi(argv[4]);
    string il = argv[5];
    if(!fs::exists(outputpath))fs::create_directories(outputpath);
    string summarypath = outputpath + "/summary.txt";
    // freopen(summarypath.c_str(), "w", stdout);

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
    l0 = 0, r0 = double_length / 2 - 1, l1 = double_length + 1, r1 = l1 + r0;

    closestloci_forward = new int[H];
    for(int i = 0; i < H; i++)closestloci_forward[i] = -1;
    int pst = -1, maxgap = 0, cntloci = 0;
    for(int i = double_length - 1; i >= 0; i--){
        if(text[i] != text[i + 2 * H + 1]){
            cntloci += 1;
            if(pst != -1){
                int val = pst - i - 1;
                maxgap = max(maxgap, val);
            }
            pst = i;
        }
        int loc = i;
        if(loc >= H)loc -= H;
        if(closestloci_forward[loc] != -1)continue;
        closestloci_forward[loc] = pst;
        if(closestloci_forward[loc] >= H)closestloci_forward[loc] -= H;
    }
    cout << "Number of heterozygous loci: " << (cntloci / 2) << endl;
    cout << "Maximum gap between two adjacent heterozygous loci: " << maxgap << endl;

    closestloci_backward = new int[H];
    for(int i = 0; i < H; i++)closestloci_backward[i] = -1;
    pst = -1;
    for(int i = 0; i <= double_length - 1; i++){
        if(text[i] != text[i + 2 * H + 1]){
            pst = i;
            if(pst >= H)pst -= H;
        }
        int loc = i;
        if(loc >= H)loc -= H;
        if(closestloci_backward[loc] != -1)continue;
        closestloci_backward[loc] = pst;
    }

    // for(int i = 0; i < H; i++)cout << closestloci_forward[i] << " ";
    // cout << endl;
    // for(int i = 0; i < H; i++)cout << closestloci_backward[i] << " ";
    // cout << endl;

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

    // max length double repeat
    int maxrepeat = 0; int maxrepeat_F = 0; 
    int maxrepeat_mat = 0, maxrepeat_pat = 0; int maxrepeat_mat_F = 0, maxrepeat_pat_F = 0; 
    int maxtriple_mat = 0, maxtriple_pat = 0; int maxtriple_mat_F = 0, maxtriple_pat_F = 0;
    int minlength_wellbridging = 0; 
    // stats of double repeats for greedy
    int bothcover = 0; int bothcover_F = 0;
    int onecover = 0; int onecover_F = 0;
    int nonecover = 0; int nonecover_F = 0;
    vector<piii> double_mat, double_pat, double_both;
    for(int i = 0; i < len - 1; i++){
        bool mati = false;
        if(!in_range(SA[i], {l0, r0, l1, r1}))continue;
        if(in_range(SA[i], {l0, r0}))mati = true;
        int val = H;
        int dist1 = distanceforward(SA[i]);
        int repeatcnt_mat = 0, repeatcnt_pat = 0; 
        bool repeatcnt_mat_lm = false, repeatcnt_pat_lm = false;
        for(int j = i + 1; j < min(len, i + offset); j++){
            int cntmat = mati; 
        
            val = min(val, lcp[j]);
            if(val < threshold)break;
            
            if(!in_range(SA[j], {l0, r0, l1, r1}))continue;
            if(in_range(SA[j], {l0, r0}))cntmat++;
            
            if(colocated(SA[i], SA[j]))continue;
            bool lm = left_maximal(SA[i], SA[j], s0, s1);
            
            // Precomputation for interleaved repeats
            if(lm){
                int id1 = SA[i], id2 = SA[j];
                if(id1 >= 2 * H + 1)id1 -= 2 * H + 1;
                if(id2 >= 2 * H + 1)id2 -= 2 * H + 1;
                pii id = {min(id1, id2), max(id1, id2)};
                if(cntmat == 2)double_mat.pb({val, id});
                else if(cntmat == 1){// order of indices matter here
                    if(!mati)swap(id1, id2); // want i to be from mat
                    double_both.pb({val, {id1, id2}}); 
                }else double_pat.pb({val, id});
            }

            // Maximum double repeat
            if(lm){
                if(val > maxrepeat){
                    maxrepeat = val;
                    maxrepeat_F = 1;
                }else if(val == maxrepeat)maxrepeat_F++;
            }

            // Maximum intra-double + intra-triple repeat
            if(cntmat == 2){
                if(lm){
                    if(val > maxrepeat_mat){
                        maxrepeat_mat = val;
                        maxrepeat_mat_F = 1;
                    }else if(val == maxrepeat_mat)maxrepeat_mat_F++;
                }
                if(lm || !repeatcnt_mat_lm){
                    repeatcnt_mat += 1;
                    if(repeatcnt_mat >= 2){
                        if(val > maxtriple_mat){
                            maxtriple_mat = val;
                            maxtriple_mat_F = 1;
                        }else if(val == maxtriple_mat)maxtriple_mat_F++;
                    }
                    if(!lm)repeatcnt_mat_lm = true;
                }
            }else if(cntmat == 0){
                if(lm){
                    if(val > maxrepeat_pat){
                    maxrepeat_pat = val;
                    maxrepeat_pat_F = 1;
                    }else if(val == maxrepeat_pat)maxrepeat_pat_F++;
                }
                if(lm || !repeatcnt_pat_lm){
                    repeatcnt_pat += 1;
                    if(repeatcnt_pat >= 2){
                        if(val > maxtriple_pat){
                            maxtriple_pat = val;
                            maxtriple_pat_F = 1;
                        }else if(val == maxtriple_pat)maxtriple_pat_F++;
                    }
                    if(!lm)repeatcnt_pat_lm = true;
                }
            }
            
            // Stats for Well-bridging
            if(lm){
                int dist2 = distanceforward(SA[j]);
                if(dist1 > val && dist2 > val){// both copies do not cover the locus
                    int upd = 2 * val + distanceloci(SA[i]) + distanceloci(SA[j]); 
                    if(upd > nonecover){
                        nonecover = upd;
                        nonecover_F = 1;
                    }else if(upd == nonecover)nonecover_F++;
                    int len1 = gap(SA[i], val);
                    int len2 = gap(SA[j], val);
                    minlength_wellbridging = max(minlength_wellbridging, min(len1, len2));
                }else if(dist1 <= val && dist2 <= val){// both copies cover the locus
                    int upd = val; 
                    if(upd > bothcover){
                        bothcover = upd;
                        bothcover_F = 1;
                    }else if(upd == bothcover)bothcover_F++;
                    minlength_wellbridging = max(minlength_wellbridging, val + 2);
                }else{// exactly one copy covers the locus
                    int upd = 3 * val + 2 * (dist1 > val ? distanceloci(SA[i]) : distanceloci(SA[j]));
                    if(upd > onecover){
                        onecover = upd;
                        onecover_F = 1;
                    }else if(upd == onecover)onecover_F++;
                    int lenc = dist1 > val ? gap(SA[i], val) : gap(SA[j], val);
                    minlength_wellbridging = max(minlength_wellbridging, min(lenc, val + 2));
                }   
            }
        }
    } 
    
    cout << "Maximum length of double repeat: " << maxrepeat << " with freq: " << maxrepeat_F << endl;
    cout << "Minimum read length required for well-bridging of double repeats: " << minlength_wellbridging << endl;
    cout << "Maximum length of double repeat for maternal haplotype: " << maxrepeat_mat << " with freq: " << maxrepeat_mat_F << endl;
    cout << "Maximum length of double repeat for paternal haplotype: " << maxrepeat_pat << " with freq: " << maxrepeat_pat_F << endl;
    cout << "Maximum length of triple repeat for maternal haplotype: " << maxtriple_mat << " with freq: " << maxtriple_mat_F << endl;
    cout << "Maximum length of triple repeat for paternal haplotype: " << maxtriple_pat << " with freq: " << maxtriple_pat_F << endl;
    cout << "Maximum stat of double repeat where both copies cover the het loci: " << bothcover << " with freq: " << bothcover_F << endl;
    cout << "Maximum stat of double repeat where exactly one copy covers the het loci: " << onecover << " with freq: " << onecover_F << endl;
    cout << "Maximum stat of double repeat where both copies do not cover the het loci: " << nonecover << " with freq: " << nonecover_F << endl;
    if(il == "1"){
        piii maxinterleaved_mat = find_interleaved(double_mat, s0);
        vector<pii> interleaved_mat_stats = find_interleaved_stats(double_mat, s0, maxinterleaved_mat.second);
        piii maxinterleaved_pat = find_interleaved(double_pat, s1);
        vector<pii> interleaved_pat_stats = find_interleaved_stats(double_pat, s1, maxinterleaved_pat.second);
        piii specialinterleaved = find_specialinterleaved(double_both);
        vector<pii> specialinterleaved_stats = find_specialinterleaved_stats(double_both, specialinterleaved.second);

        cout << "Maximum length of interleaved repeat for maternal haplotype: " << maxinterleaved_mat.second.first << " " << maxinterleaved_mat.second.second << " with freq: " << maxinterleaved_mat.first << endl;
        cout << "Maximum length of interleaved repeat for paternal haplotype: " << maxinterleaved_pat.second.first << " " << maxinterleaved_pat.second.second << " with freq: " << maxinterleaved_pat.first << endl;
        cout << "Maximum length of repeat of type I2: " << specialinterleaved.second.first << " " << specialinterleaved.second.second << " with freq: " << specialinterleaved.first << endl;
    
        string out = outputpath + "/interleaved_mat_stats.txt";
        print(interleaved_mat_stats, out);
        
        out = outputpath + "/interleaved_pat_stats.txt";
        print(interleaved_pat_stats, out);
        
        out = outputpath + "/specialinterleaved_stats.txt";
        print(specialinterleaved_stats, out);
    }
    string gappath = outputpath + "/gapstats.txt";
    printgap(gappath, maxgap);
} 
