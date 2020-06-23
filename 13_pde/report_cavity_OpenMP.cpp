//
// Created by 张骁尧 on 2020/06/24.
//

#include <iostream>
#include <vector>
#include <fstream>
#include <math.h>
#include <omp.h>

using namespace std;

typedef vector<vector<double> > vector_2d;
typedef vector<double> vector_1d;
void copy_2d(vector_2d &from, vector_2d &to){
    #pragma omp parallel for
    for(int i=0; i<from.size(); i++){
        for(int j=0; j<from.at(0).size(); j++){
            to[i][j] = from[i][j];
        }
    }

}
vector_2d build_up_b(vector_2d &b, double rho, double dt, vector_2d u, vector_2d v, double dx, double dy){
    #pragma omp parallel for
    for(int i=1;i<b.size()-1; i++){
        for(int j=1; j<b.at(0).size()-1; j++){
            b[i][j] = (rho * (1 / dt *
                              ((u[i][j+1] - u[i][j-1]) /
                               (2 * dx) + (v[i+1][j] - v[i-1][j]) / (2 * dy)) -
                              pow(((u[i][j+1] - u[i][j-1]) / (2 * dx)),2)-
                              2 * ((u[i+1][j] - u[i-1][j])/( 2 * dy) *
                                   (v[i][j+1] - v[i][j-1]) / (2 * dx)) -
                              pow(((v[i+1][j] - v[i-1][j]) / (2 * dy)),2)));
        }
    }
    return b;
};

vector_2d pressure_poisson(vector_2d p, double dx, double dy, vector_2d b, int nit){
    vector_2d pn(p.size(),vector_1d(p.at(0).size(),0));
    copy_2d(p, pn);
    for(int q=0; q<nit; q++){
        copy_2d(p, pn);
        #pragma omp parallel for
        for(int i=1; i<p.size()-1; i++){
            for(int j=0; j<p.at(0).size()-1; j++){
                p[i][j] = ((pn[i][j+1] + pn[i][j-1]) * pow(dy, 2) +
                           (pn[i+1][j] + pn[i-1][j]) * pow(dx, 2)) /
                          (2 * (pow(dx, 2) + pow(dy, 2))) -
                          pow(dx, 2) * pow(dy, 2) / (2 * (pow(dx, 2) + pow(dy, 2))) *
                          b[i][j];
            }
        }
        #pragma omp parallel for
        for(int i=0; i<p.size(); i++){
            p[i][p.at(0).size()-1] = p[i][p.at(0).size()-2];
            p[i][0] = p[i][1];
        }
        #pragma omp parallel for
        for(int j=0; j<p.at(0).size(); j++){
            p[0][j] = p[1][j];
            p[p.size()-1][j] = 0;
        }
    }
    return p;
}

void cavity_flow(int nt, vector_2d &u, vector_2d &v, double dt, double dx, double dy, vector_2d &p, double rho, double nu, int nit){
    vector_2d un(u.size(),vector_1d(u.at(0).size(),0));
    vector_2d vn(v.size(),vector_1d(v.at(0).size(),0));
    vector_2d b(u.size(),vector_1d(u.at(0).size(),0));

    for(int i=0; i<nt; i++) {
        copy_2d(u, un);
        copy_2d(v, vn);

        b = build_up_b(b, rho, dt, u, v, dx, dy);
        p = pressure_poisson(p, dx, dy, b, nit);

        #pragma omp parallel for
        for (int i = 1; i < u.size()-1; i++) {
            for (int j = 1; j < u.at(0).size()-1; j++) {
                u[i][j] = (un[i][j] -
                           un[i][j] * dt / dx *
                           (un[i][j] - un[i][j - 1]) -
                           vn[i][j] * dt / dy *
                           (un[i][j] - un[i - 1][j]) -
                           dt / (2 * rho * dx) * (p[i][j + 1] - p[i][j - 1]) +
                           nu * (dt / pow(dx, 2) *
                                 (un[i][j + 1] - 2 * un[i][j] + un[i][j - 1]) +
                                 dt / pow(dy, 2) *
                                 (un[i + 1][j] - 2 * un[i][j] + un[i - 1][j])));

                v[i][j] = (vn[i][j] -
                           un[i][j] * dt / dx *
                           (vn[i][j] - vn[i][j - 1]) -
                           vn[i][j] * dt / dy *
                           (vn[i][j] - vn[i - 1][j]) -
                           dt / (2 * rho * dy) * (p[i + 1][j] - p[i - 1][j]) +
                           nu * (dt / pow(dx, 2) *
                                 (vn[i][j + 1] - 2 * vn[i][j] + vn[i][j - 1]) +
                                 dt / pow(dy, 2) *
                                 (vn[i + 1][j] - 2 * vn[i][j] + vn[i - 1][j])));

            }
        }
        #pragma omp parallel for
        for (int i = 0; i < u.size(); i++) {
            u[i][0] = 0;
            u[i][u.at(0).size() - 1] = 0;
            v[i][0] = 0;
            v[i][v.at(0).size() - 1] = 0;
        }
        #pragma omp parallel for
        for (int j = 0; j < u.at(0).size(); j++) {
            u[0][j] = 0;
            u[u.size() - 1][j] = 1;
            v[0][j] = 0;
            v[v.size() - 1][j] = 0;
        }
    }
}

int main(int argc, char** argv){
    int nx = 41;
    int ny = 41;
    int nt = 500;
    int nit = 50;
    int c = 1;
    double dx = (double)2 / (nx - 1);
    double dy = (double)2 / (ny - 1);

    double rho = 1;
    double nu = 0.1;
    double dt = 0.001;
    vector_2d u(ny,vector_1d(nx,0));
    vector_2d v(ny,vector_1d(nx,0));
    vector_2d p(ny,vector_1d(nx,0));
    vector_2d b(ny,vector_1d(nx,0));
    cavity_flow(nt, u, v, dt, dx, dy, p, rho, nu, nit);
    ofstream u_out("u_openmp.txt");
    ofstream v_out("v_openmp.txt");
    ofstream p_out("p_openmp.txt");
    for(int i=0; i<ny; i++){
        for(int j=0; j<nx; j++){
            u_out<<u[i][j]<<" ";
            v_out<<v[i][j]<<" ";
            p_out<<p[i][j]<<" ";
        }
        u_out<<endl;
        v_out<<endl;
        p_out<<endl;
    }
    u_out.close();
    v_out.close();
    p_out.close();

    return 0;


}
