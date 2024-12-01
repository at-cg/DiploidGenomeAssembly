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

int main()
{   
    string query = "ACTGAT";
    const char *text = query.c_str();
    int len = strlen(text);
    
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
    // free(iSA);

    for(int i = 0; i < len; i++)cout << SA[i] << " ";
    cout << endl;

    for(int i = 0; i < len; i++)cout << iSA[i] << " ";
    cout << endl;

    for(int i = 0; i < len; i++)cout << lcp[i] << " ";
    cout << endl;
} 
