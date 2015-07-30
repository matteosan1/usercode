/*
 *  binerr.cpp
 *  
 *
 *  Created by Avi Yagil on 7/23/07.
 *
 */

#include <iostream>
#include <iomanip>
using namespace std;

float binerr(float denom, float num)
{
  float e = num/denom;
//   cout <<  setiosflags(ios::fixed) ;
//   cout << setw(5) << setprecision(4) << e << "+/- " ;
  float err = sqrt(e*(1.-e)/denom);
//   cout << setw(5) << setprecision(4) << err << endl;
  return err;
}                                                            
