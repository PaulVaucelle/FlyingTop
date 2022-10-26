#include <TMath.h>
double Deltar(double eta1, double phi1, double eta2, double phi2) {
  double DeltaPhi = TMath::Abs(phi2 - phi1);
  if (DeltaPhi > 3.141593 ) DeltaPhi = 2.*3.141593 - DeltaPhi;
  return TMath::Sqrt( (eta2-eta1)*(eta2-eta1) + DeltaPhi*DeltaPhi );
}
double Deltaphi(double phi1, double phi2) {
  double DeltaPhi = phi1 - phi2;
  if (abs(DeltaPhi) > 3.141593 ) {
    DeltaPhi = 2.*3.141593 - abs(DeltaPhi);
    DeltaPhi = -DeltaPhi * (phi1 - phi2) / abs(phi1 - phi2);
  }
  return DeltaPhi;
}