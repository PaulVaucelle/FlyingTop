class PDFVar { // Becareful, give the vector of only the PDF varation weights, so avoid the first 9 weights in the ntuples
public:
    PDFVar(std::vector<float> lhew) {
        // Constructor code here
        LHEWeight = lhew;
    }

   // Add member functions and variables here
    float GetPDFOriginalFast()
        {
            return LHEWeight[0];
        }
    float GetPDFOriginal()
        {
            float sum =0;
            float n = LHEWeight.size();
            for (unsigned int i = 1 ; i < LHEWeight.size(); i++)
                {
                    sum += abs(LHEWeight[i]);
                    std::cout<<"LHEWeight["<<i<<"] = "<<LHEWeight[i]<<std::endl;
                }
            return sum/n;
        }
    float GetPDFVar()
        {
            float sum = 0;
            float n = LHEWeight.size();
            float mean = GetPDFOriginal();
            for (unsigned int i = 1 ; i < LHEWeight.size(); i++)
                {
                    sum += (abs(LHEWeight[i])-mean)*(abs(LHEWeight[i])-mean);
                }
            sum = sum/(n-1);
            return sqrt(sum);
        }
    float GetPDFVarUp()
        {
            return GetPDFOriginal()+GetPDFVar();
        }
    float GetPDFVarDown()
        {
            return GetPDFOriginal()-GetPDFVar();
        }

    ~PDFVar() {
        // Destructor code here
    }

 

private:
    // Add private member variables here
    std::vector<float> LHEWeight;

};



class ScaleVars {
public:
    ScaleVars(std::vector<float> lhew) {
        // Constructor code here
        LHEWeight = lhew;
    }

    ~ScaleVars() {
        // Destructor code here
    }

    // Add member functions and variables here
    float GetScaleOriginal()
        {
            return LHEWeight[0];
        }

    float GetScaleVarUp()
        {
            float meanup=0;
            float count;
            for (unsigned int i = 1 ; i < LHEWeight.size(); i++)
                {
                    if (LHEWeight[i] >= LHEWeight[0])
                        {
                            meanup += LHEWeight[i];
                            count++;
                        }
                }
            return (meanup/count)/LHEWeight[0];
        }
    float GetScaleVarDown()
        {
            float meandown=0;
            float count;
            for (unsigned int i = 1 ; i < LHEWeight.size(); i++)
                {
                    if (LHEWeight[i] <= LHEWeight[0])
                        {
                            meandown += LHEWeight[i];
                            count++;
                        }
                }
            return (meandown/count)/LHEWeight[0];
        }
private:
    // Add private member variables here
     std::vector<float> LHEWeight;



};
