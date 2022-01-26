#include <bits/stdc++.h>
using namespace std;
double x[]={0, 15, 30, 45, 60, 75, 90};
double y[]={0, 0.258819, 0.5, 0.707107, 0.866025, 0.965926, 1}; //example 382
double delta(int idx, int n)
{
   if(n==0) return y[idx];
   return delta(idx+1, n-1)-delta(idx, n-1);
}
int main()
{
  cout<<delta(1, 5)<<endl;
  return 0;
}
