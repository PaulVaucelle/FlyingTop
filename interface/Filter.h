#include "TLorentzVector.h"
std::vector<float> DiLeptonMass(bool SSLep, int nLep, int nAllLepton, std::vector<bool> ID, std::vector<bool> ISO , std::vector<bool> isPrompt ,std::vector<float> PT,std::vector<float> ETA ,std::vector<float> PHI, std::vector<int> CHARGE, int* index, float* Lepmass) 
    {
        TLorentzVector v1, v2, v;
        std::vector<float> OutputData ;
        float Mmumu = 0, MmumuSameSign = 0;
        int imu1 = -1, imu2 = -1;
        int imu1_SS = -1, imu2_SS = -1;
        int Q1 = 0, Q2 = 0;
        float mupt1, mueta1, muphi1, mupt2, mueta2, muphi2;
        if (nLep >= 2 ) 
            {
                for (int mu = 0; mu < nAllLepton - 1; mu++) {
                    if (!isPrompt[index[mu]]) continue;
                    if (!ID[index[mu]]) continue;
                    if (!ISO[index[mu]]) continue;
                    if (!isPrompt[index[mu]]) continue;
                    double mupt1 = PT[index[mu]];
                    if (mupt1 < 25.) continue;
                    double mueta1 = ETA[index[mu]];
                    double muphi1 = PHI[index[mu]];
                    v1.SetPtEtaPhiM(mupt1, mueta1, muphi1, Lepmass[0]);
                    Q1 = CHARGE[index[mu]];

                    for (int mu2 = mu + 1; mu2 < nAllLepton; mu2++) {
                        if (!isPrompt[index[mu2]]) continue;
                        if (!ID[index[mu2]]) continue;
                        if (!ISO[index[mu2]]) continue;
                        double mupt2 = PT[index[mu2]];
                        if (mupt2 < 10.) continue;
                        Q2 = CHARGE[index[mu2]];
                        if (Q1 == Q2 && !SSLep) continue;
                        mueta2 = ETA[index[mu2]];
                        muphi2 = PHI[index[mu2]];
                        v2.SetPtEtaPhiM(mupt2, mueta2, muphi2, Lepmass[1]);
                        v = v1 + v2;
                        if (v.Mag() > Mmumu && Q1 != Q2) {
                            Mmumu = v.Mag();
                            imu1 = index[mu];
                            imu2 = index[mu2];
                        }
                        if (v.Mag() > MmumuSameSign && SSLep && Q1 == Q2) {
                            MmumuSameSign = v.Mag();
                            imu1_SS = index[mu];
                            imu2_SS = index[mu2];
                        }
                    }
                }
                
                OutputData.push_back(Mmumu);
                OutputData.push_back(MmumuSameSign);
                OutputData.push_back(imu1);
                OutputData.push_back(imu2);
                OutputData.push_back(imu1_SS);
                OutputData.push_back(imu2_SS);
                return OutputData;
            }// if at least 2 leptons
        else    
            {
                return {0.};
            }
    }

    std::vector<float> EMuMass(bool SSLep, int nLep, int nAllLepton, std::vector<bool> ID, std::vector<bool> ISO , std::vector<bool> isPrompt ,std::vector<float> PT,std::vector<float> ETA ,std::vector<float> PHI, std::vector<int> CHARGE, int* index, 
    int nLep2, int nAllLepton2, std::vector<bool> ID2, std::vector<bool> ISO2 , std::vector<bool> isPrompt2 ,std::vector<float> PT2,std::vector<float> ETA2 ,std::vector<float> PHI2, std::vector<int> CHARGE2, int* index2, float* Lepmass) 
    {
        TLorentzVector v1, v2, v;
        std::vector<float> OutputData ;
        float Mmumu = 0, MmumuSameSign = 0;
        int imu1 = -1, imu2 = -1;
        int imu1_SS = -1, imu2_SS = -1;
        int Q1 = 0, Q2 = 0;
        float mupt1, mueta1, muphi1, mupt2, mueta2, muphi2;
        if (nLep >= 1 && nLep2 >= 1 ) 
            {
                for (int mu = 0; mu < nAllLepton - 1; mu++) {
                    if (!isPrompt[index[mu]]) continue;
                    if (!ID[index[mu]]) continue;
                    if (!ISO[index[mu]]) continue;
                    if (!isPrompt[index[mu]]) continue;
                    double mupt1 = PT[index[mu]];
                    if (mupt1 < 25.) continue;
                    double mueta1 = ETA[index[mu]];
                    double muphi1 = PHI[index[mu]];
                    v1.SetPtEtaPhiM(mupt1, mueta1, muphi1, Lepmass[0]);
                    Q1 = CHARGE[index[mu]];

                    for (int mu2 = 0; mu2 < nAllLepton2; mu2++) {
                        if (!isPrompt2[index2[mu2]]) continue;
                        if (!ID2[index2[mu2]]) continue;
                        if (!ISO2[index2[mu2]]) continue;
                        double mupt2 = PT2[index2[mu2]];
                        if (mupt2 < 10.) continue;
                        Q2 = CHARGE2[index2[mu2]];
                        if (Q1 == Q2 && !SSLep) continue;
                        mueta2 = ETA2[index2[mu2]];
                        muphi2 = PHI2[index2[mu2]];
                        v2.SetPtEtaPhiM(mupt2, mueta2, muphi2, Lepmass[1]);
                        v = v1 + v2;
                        if (v.Mag() > Mmumu && Q1 != Q2) {
                            Mmumu = v.Mag();
                            imu1 = index[mu];
                            imu2 = index2[mu2];
                        }
                        if (v.Mag() > MmumuSameSign && SSLep && Q1 == Q2) {
                            MmumuSameSign = v.Mag();
                            imu1_SS = index[mu];
                            imu2_SS = index2[mu2];
                        }
                    }
                }
                
                OutputData.push_back(Mmumu);
                OutputData.push_back(MmumuSameSign);
                OutputData.push_back(imu1);
                OutputData.push_back(imu2);
                OutputData.push_back(imu1_SS);
                OutputData.push_back(imu2_SS);
                return OutputData;
            }// if at least 2 leptons
        else    
            {
                return {0.};
            }
    }



