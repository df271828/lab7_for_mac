

#include <iostream>
#include "matrixfunctions1.h"
using namespace std;
int main(int argc, const char * argv[] )
{
 
    //double r = 4.0;
    //cout << CircleArea(r);

    cout.precision(4);
    vector<vector<double>> data ={{198,92,-1,48,48,45,420,115,-1,98,-1,100},{184,84,-1,44,33,33,350,102,-1,92,-1,130},{183,83,-1,44,37,34,320,98,-1,91,-1,127},{182,80,-1,42,35,30,398,65,-1,85,-1,140},{180,80,-1,43,36,30,388,63,-1,84,-1,129},{183,81,-1,42,37,35,345,45,-1,90,-1,105},{180,82,-1,44,43,37,355,82,-1,88,-1,109},{180,81,-1,44,46,42,362,90,-1,86,-1,113},{185,82,-1,45,26,16,295,180,-1,92,1,109},{187,84,-1,46,27,16.5,299,178,-1,95,1,119},{177,65,-1,41,26,18,209,160,-1,86,1,120},{180,72,-1,43,33,19,236,175,-1,85,1,115},{181,75,-1,43,42,31,198,161,-1,83,1,105},{176,68,-1,42,50,36,195,177,-1,82,1,96},{175,67,1,42,55,38,185,187,-1,80,1,105},{178,75,-1,42,30,24,203,208,-1,81,1,118},{166,47,-1,36,32,28,270,78,1,75,-1,112},{170,60,1,38,23,20,312,99,1,81,-1,110},{172,64,1,39,24,22,308,91,1,82,-1,102},{169,51,1,36,24,23,250,89,1,78,-1,98},{168,52,1,37,27,23.5,260,86,1,78,-1,100},{157,47,1,36,32,32,235,92,1,70,-1,127},{164,50,1,38,41,34,255,134,1,76,-1,101},{162,49,1,37,40,34,265,124,1,75,-1,108},{168,50,1,37,49,34,170,162,1,76,1,135},{166,49,1,36,21,14,150,245,1,75,1,123},{158,46,1,34,30,18,120,120,1,70,1,119},{163,50,1,36,18,11,143,136,1,75,1,102},{162,50,1,36,20,11.5,133,146,1,74,1,132},{165,51,1,36,36,26,121,129,1,76,1,126},{161,48,1,35,41,31.5,116,196,1,75,1,120},{160,48,1,35,40,31,118,198,1,74,1,129}};
    
    Matrix a(data);
    PCA pca(a);
    auto e = pca.ResMatrix(4);

    auto res = e.GetMatrix();
    /* for (int i = 0; i < res.size(); i++){
          for (int j = 0; j < res[i].size(); j++){
              cout << res[i][j] << " ";
          }
         cout << endl;
    }*/
    cout << pca.ERVP(4);
    return 0;
    return 0;
 return 0;
 
}
