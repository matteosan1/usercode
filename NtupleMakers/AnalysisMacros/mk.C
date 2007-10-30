/*
 *  mk.C
 *  
 *
 *  Created by Avi Yagil on 7/23/07.
 *
 */
 #include "TH1F.h"
                                                
 void mk(TH1F &fr, TH1F &num, TH1F &denom)
{
  float numf, denomf;
  for (int i1 = 1; i1 <= fr.GetNbinsX(); ++i1 ) {    
    numf = num.GetBinContent(i1);
    denomf = denom.GetBinContent(i1);
    if ( denomf != 0 ) {
      fr.SetBinContent(i1,numf/denomf );
      float err = binerr(denomf, numf);
      fr.SetBinError(i1, err);
    }
  }
  return;
}                                  
