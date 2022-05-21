#pragma once
#include <iostream>
#include <vector>
#include <map>
#include <set>
#include <cmath>
#include <initializer_list>
#include <bitset>
#include <exception>
#include <algorithm>
#include <fstream>
#include <iomanip>
#include <string>
#include <complex>
using namespace std;

#define PI 3.1415;
 
double CircleArea (double radius);
double CircleCircum (double radius);


const double EPS = 1E-8;
class Matrix{
    friend Matrix InvertMatrix(Matrix t);
protected:
    vector<vector<double>> matrix;
    int n;
    int m;
public:
    Matrix (){
    }

    vector < double > v ;
    Matrix(initializer_list < double > l )  : v ( l )  {
        matrix.resize(1, vector<double> (1, 0));
        matrix[0] = v;
    }


    Matrix (int a, int b){
        if (a < 0 || b < 0)
            throw invalid_argument("Invalid matrix size");
        n = a;
        m = b;
        matrix.resize(a, vector<double> (b, 0));
    }

    Matrix(const vector<vector<double>>vec){
        set<int> st;
        for (int i = 0; i < vec.size(); i++){
            st.insert(vec[i].size());
        }
        matrix = vec;
        n = vec.size();
        for (auto x: st)
            m = x;
    }

    Matrix (const Matrix& t){
        n = t.GetN();
        m = t.GetM();
        vector<vector<double>> res = t.GetMatrix();
        matrix.resize(n, vector<double> (m, 0));
        for (int i = 0; i < n; i++){
            for (int j = 0; j < m; j++){
                matrix[i][j] = res[i][j];
            }
        }
    }

    int GetN() const{
        return n;
    }

    int GetM() const{
        return m;
    }

    vector<vector<double>> GetMatrix() const{
        return matrix;
    }

    void SetCoef(const vector<vector<double>> &coef){
        matrix = coef;
        n = coef.size();
        set<int> st;
        for (int i = 0; i < coef.size(); i++){
            st.insert(coef[i].size());
        }

        for (auto x: st)
            m = x;
    }

    Matrix operator+ (Matrix& t){
        if (n != t.GetN() || m != t.GetM())
            throw invalid_argument("Can't sum this matrix");
        vector<vector<double>> a = GetMatrix();
        vector<vector<double>> b = t.GetMatrix();
        vector<vector<double>> res;
        res.resize(n, vector<double> (m, 0));
        for (int i = 0; i < n; i++){
            for (int j = 0; j < m; j++){
                res[i][j] = a[i][j] + b[i][j];
            }
        }

        return Matrix(res);
    }

    Matrix operator- (Matrix& t){
        if (n != t.GetN() || m != t.GetM())
            throw invalid_argument("Can't sum this matrix");
        vector<vector<double>> a = GetMatrix();
        vector<vector<double>> b = t.GetMatrix();
        vector<vector<double>>res;
        res.resize(n, vector<double> ());
        for (int i = 0; i < n; i++){
            for (int j = 0; j < m; j++){
                res[i].push_back(a[i][j] - b[i][j]);
            }
        }
        return Matrix(res);
    }

    Matrix operator* (Matrix& t){
        if (GetM() != t.GetN())
            throw invalid_argument("Can't sum this matrix");
        int k = t.GetM();
        vector<vector<double>> a = GetMatrix();
        vector<vector<double>> b = t.GetMatrix();
        vector<vector<double>>res;
        res.resize(n, vector<double> (k, 0));

        for (int i = 0; i < n; i++){
            for (int j = 0; j < k; j++){
                for (int z = 0; z < m; z++)
                    res[i][j] += a[i][z] * b[z][j];
            }
        }
        return Matrix(res);
    }

    Matrix operator/ (Matrix& t){
        Matrix t_inv = InvertMatrix(t);
        return (*this) * t_inv;
    }

    void operator/ (double t){
        if (t == 0)
            throw invalid_argument("Lol kek cheburek 0");
        NumberMultiply(1./t);
    }

    vector<vector<double>> GetLineByNumber(int c){
        if (c > GetN() || c < 0)
            throw invalid_argument("invalid line number");
        vector<double> result;
        for (int i = 0; i < GetM(); i++)
            result.push_back(matrix[c][i]);
        vector<vector<double>> answer;
        answer.push_back(result);
        return answer;
    }

    vector<vector<double>> GetColumnByNumber(int c){
        if (c > GetM() || c < 0)
            throw invalid_argument("invalid column number");
        vector<double> result;
        for (int i = 0; i < GetN(); i++)
            result.push_back(matrix[i][c]);
            vector<vector<double>> answer;
        for (int i = 0; i < result.size(); i++){
            vector<double> res;
            res.push_back(result[i]);
            answer.push_back(res);
        }
        return answer;
    }

    void AddLine(vector<double> vec){
        matrix.push_back(vec);
        n++;
    }

    void AddColumn(vector<double> vec){
        for(int i = 0; i < GetN(); i++){
            matrix[i].push_back(vec[i]);
        }
        m++;
    }

    void AddMatrixColumn(Matrix &t){
        vector<vector<double>> vec = t.GetMatrix();
        if (vec[0].size() != 1)
            AddColumn(vec[0]);
        else{
            vector<double> help;
            for (int i = 0; i < vec.size(); i++){
                help.push_back(vec[i][0]);
            }
            AddColumn(help);
        }
    }




    void NumberMultiply(double d){
        for (int i = 0; i < n; i++){
            for (int j = 0; j < m; j++){
                matrix[i][j] *= d;
            }
        }
    }

    Matrix AdamarMultiply (Matrix& t){
        if (n != t.GetN() || m != t.GetM())
            throw invalid_argument("Can't sum this matrix");
        vector<vector<double>> a = GetMatrix();
        vector<vector<double>> b = t.GetMatrix();
        vector<vector<double>> res;
        res.resize(n, vector<double> ());
        for (int i = 0; i < n; i++){
            for (int j = 0; j < m; j++){
                res[i].push_back(a[i][j] * b[i][j]);
            }
        }
        return Matrix(res);
    }

    double MatrixTrace(){
        double answer = 0;
        if (n > m){
            for (int i = 0; i < m; i++)
                answer += matrix[i][i];
        }
        else {
            for (int i = 0; i < n; i++)
                answer += matrix[i][i];
        }
        return answer;
    }

    double ScalarMultiply(Matrix& t){
        if (n != 1 && m != 1)
            throw invalid_argument("Invalid operation");
         if (t.GetN() != 1 && t.GetM() != 1)
            throw invalid_argument("Invalid operation");
        double answer = 0;
        Matrix c = (*this).AdamarMultiply(t);
        vector<vector<double>> res = c.GetMatrix();
        for (int i = 0; i < res.size(); i++){
            for (int j = 0; j < res[i].size(); j++){
                answer += res[i][j];
            }
        }
        return answer;
    }

    double VectorNorma(){
        if (n != 1 && m != 1)
            throw invalid_argument("Invalid operation");
        double answer = 0;
        for (int i = 0; i < n; i++){
            for (int j = 0; j < m; j++){
                answer += matrix[i][j] * matrix[i][j];
            }
        }
        return sqrt(answer);
    }
    double MaxNorma(){
        if (n != 1 && m != 1)
            throw invalid_argument("Invalid operation");
        double answer = 0;
        for (int i = 0; i < n; i++){
            for (int j = 0; j < m; j++){
                answer = max(answer, abs(matrix[i][j]));
            }
        }
        return answer;
    }
    double FrobeniusNorma(){
        double answer = 0;
        for (int i = 0; i < n; i++){
            for (int j = 0; j < m; j++){
                answer += matrix[i][j] * matrix[i][j];
            }
        }
        return sqrt(answer);
    }

    double Gauss(){
        double det = 1;
        for (int i = 0; i < n; i++) {
            int k = i;
            for (int j = i + 1; j < m; j++)
                if (abs(matrix[j][i]) > abs(matrix[k][i]))
                    k = j;
            if (abs(matrix[k][i]) < EPS) {
                det = 0;
                break;
            }
            swap(matrix[i], matrix[k]);
            if (i != k)
                det = -det;
            det *= matrix[i][i];
            for (int j = i + 1; j < m; j++)
                matrix[i][j] /= matrix[i][i];
            for (int j = 0; j < m; j++)
                if (j != i && abs(matrix[j][i]) > EPS)
                    for (int k = i + 1; k < n; k++)
                        matrix[j][k] -= matrix[i][k] * matrix[j][i];
        }

        return det;
    }

    double Rank(){
        int rank = max(n,m);
        vector<char> used (n);
        vector<vector<double>> res = matrix;
        for (int i = 0; i < n; ++i) {
            int j;
            for (j = 0; j < m; ++j)
                if (!used[j] && abs(res[j][i]) > EPS)
                    break;
            if (j == n)
                rank--;
            else {
                used[j] = true;
                for (int p = i + 1; p < m; ++p)
                    res[j][p] /= res[j][i];
                for (int k = 0; k < n; ++k)
                    if (k != j && abs(res[k][i]) > EPS)
                        for (int p = i + 1; p < m; ++p)
                            res[k][p] -= res[j][p] * res[k][i];
            }
        }
        return rank;
    }

    double GetNumber(int i, int j){
        if (i >= n || i < 0)
            throw invalid_argument("Invalid line at GetNumber");
        if (j >= m || j < 0)
            throw invalid_argument("Invalid column at GetNumber");
        return matrix[i][j];
    }

};




complex<double> VectorsAngle(Matrix &s, Matrix &t);

Matrix Transposition (Matrix &t);


ostream& operator<< (ostream &out, const Matrix &t);

Matrix& operator >> (istream& in, Matrix& t);
    

class PCA{
private:
    Matrix matrix;

public:

    PCA (Matrix a){
        vector<vector<double>> c = a.GetMatrix();
        matrix.SetCoef(c);
    }

    Matrix GetMatrix(){
        return matrix;
    }
    void center(){
        vector<vector<double>> c = matrix.GetMatrix();
        int n = matrix.GetN();
        int m = matrix.GetM();
        for (int j = 0; j < m; j++){
            double mid = 0;
            for (int i = 0; i < n; i++){
                mid += c[i][j];
            }
            mid /= (n);
            for (int i = 0; i < n; i++){
                c[i][j] -= mid;
            }
        }
        matrix = Matrix(c);
    }

    void normalize(){
        vector<vector<double>> c = matrix.GetMatrix();
        int n = matrix.GetN();
        int m = matrix.GetM();
        for (int j = 0; j < m; j++){
            double mid = 0;
            for (int i = 0; i < n; i++){
                mid += c[i][j];
            }
            mid /= (n-1);
            double sum = 0;
            for (int i = 0; i < n; i++){
                sum += (c[i][j] - mid) * (c[i][j] - mid);
            }
            sum /= (n-1);
            sum = sqrt(sum);
             for (int i = 0; i < n; i++){
                c[i][j] /= sum;
            }
        }
        matrix = Matrix(c);
    }

    pair<Matrix, Matrix> Algo(int PC){
        center();
        normalize();
        Matrix E(matrix.GetMatrix());
        Matrix P;
        Matrix T;
        Matrix d;
        Matrix p_trans;
        Matrix p;
        double diff;
        for (int h = 0; h < PC; h++){
            vector<vector<double>> q = E.GetColumnByNumber(h);
            Matrix t(q);
            int cnt = 0;
            do{
                Matrix r = Transposition(t);
                Matrix res1 = r * E;
                Matrix res2 = r * t;
                //cout << h <<  endl;
                double res0 = res2.GetNumber(0,0);
               // cout << res0 << endl;
                res1 / res0;
                p = Matrix(Transposition(res1));
                p / p.VectorNorma();
                Matrix t_old(t.GetMatrix());
                res1 = E * p;
                p_trans = Matrix(Transposition(p));
                res2 = p_trans * p;
                res0 = res2.GetNumber(0,0);
                res1 / res0;
                t = res1;
                Matrix res3;
                res3 = t_old - t;
                d = Matrix(res3.GetMatrix());
                diff = d.VectorNorma();
                //cout << h << " " << diff << endl;
            } while (diff > EPS);
            Matrix res1 = t * p_trans;
            E = E - res1;

            if (h == 0){
                P = p;
                T = t;
            }
            else{
                P.AddMatrixColumn(p);
                T.AddMatrixColumn(t);
            }
        }
        return {P, T};
    }

    vector<double> Eigenvalues(int h){
        auto pr = Algo(h).second;
        Matrix trans = Transposition(pr);
        Matrix res = trans * pr;
        vector<vector<double>> vec = res.GetMatrix();
        vector<double> ans;
        for (int i = 0; i < vec.size(); i++){
            ans.push_back(vec[i][i]);
        }
        return ans;
    }

    vector<double> Leverage(int h){
        auto t = Algo(h).second;
        vector<double> vec = Eigenvalues(h);
        vector<vector<double>> matrix1 = t.GetMatrix();
        vector<double> result;
        for (int i = 0; i < t.GetN(); i++){
            double sum = 0;
            for(int j = 0; j < t.GetM(); j++){
                sum += matrix1[i][j] * matrix1[i][j] * 1. / vec[j];
            }
            result.push_back(sum);
        }
        return result;
    }

    Matrix ResMatrix(int h){
        auto p = Algo(h).first;
        auto t = Algo(h).second;
        auto p_trans = Transposition(p);
        auto res1 = t * p_trans;
        return Matrix(matrix) - res1;
    }

    vector<double> Sigma(int h){
        auto e = ResMatrix(h);

        vector<vector<double>> matrix1 = e.GetMatrix();
        vector<double> result;
        for (int i = 0; i < e.GetN(); i++){
            double sum = 0;
            for(int j = 0; j < e.GetM(); j++){
                sum += matrix1[i][j] * matrix1[i][j];
            }
            result.push_back(sum);
        }
        return result;
    }

    double TRVP(int h){
        auto e = ResMatrix(h);

        vector<vector<double>> matrix1 = e.GetMatrix();
        double sum = 0;
        for (int i = 0; i < e.GetN(); i++){
            for(int j = 0; j < e.GetM(); j++){
                sum += matrix1[i][j] * matrix1[i][j];
            }
        }

        return sum / (e.GetN() * e.GetM());
    }

    double ERVP(int h){
        auto e = ResMatrix(h);
        vector<vector<double>> matrix1 = matrix.GetMatrix();
        double sum1 = 0;
        for (int i = 0; i < e.GetN(); i++){
            for(int j = 0; j < e.GetM(); j++){
                sum1 += matrix1[i][j] * matrix1[i][j];
            }
        }

        vector<vector<double>> matrix2 = e.GetMatrix();
        double sum2 = 0;
        for (int i = 0; i < e.GetN(); i++){
            for(int j = 0; j < e.GetM(); j++){
                sum2 += matrix2[i][j] * matrix2[i][j];
            }
        }

        return 1 - sum2 * 1. / sum1;
    }

};
