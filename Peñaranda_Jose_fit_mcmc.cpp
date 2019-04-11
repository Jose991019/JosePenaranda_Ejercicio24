#include <fstream>
#include <iostream>
#include <string.h>
#include <stdlib.h>
#include <cmath>
#include <ctime>
#include <algorithm>

using namespace std; 

float *read_file(string filename, int *n_points);
void model(float *y, float*x, int n_points, float *c, int poly_degree);
float loglikelihood(float *x_obs, float *y_obs, int n_points, float *c, int poly_degree);
float logprior(float *c, int poly_degree);
void MCMC_polynomial(float *x_obs, float *y_obs, int n_points, int n_steps, int poly_degree);

int main(){
  float *x=NULL;
  float *y=NULL;
  int n_x=0;
  int n_y=0;

  x = read_file("valores_x.txt", &n_x);
  y = read_file("valores_y.txt", &n_y);
  
  MCMC_polynomial(x, y, n_x, 1000000, 3);
  
  return 0;
}

float * read_file(string filename, int *n_points){
  int n_lines=0;
  ifstream infile; 
  string line;
  int i;
  float *array;

  /*cuenta lineas*/
  infile.open(filename); 
  getline(infile, line);
  while(infile){
    n_lines++;
    getline(infile, line);
  }
  infile.close();
  *n_points = n_lines;

  /*reserva la memoria necesaria*/
  array = new float[n_lines];

  /*carga los valores*/
  i=0;
  infile.open(filename); 
  getline(infile, line);  
  while(infile){
    array[i] = atof(line.c_str());
    i++;
    getline(infile, line);
  }
  infile.close();

  return array;
}

void model(float *y, float*x, int n_points, float *c, int poly_degree){
    for(int i = 0; i < n_points; i ++){
        for(int j = 0; i < poly_degree; i ++){
            y[i] += c[j]*pow(x[i],j);
        }
    }
}
float loglikelihood(float *x_obs, float *y_obs, int n_points, float *c, int poly_degree){
    float suma = 0;
    float sigma = 0.001;
    float *y = new float[n_points];
    model(y,x_obs,n_points,c,poly_degree);
    for(int i = 0; i < n_points; i ++){
        suma += pow((y_obs[i] - y[i])/sigma,2);
    }
    return(-0.5*suma);
}
float logprior(float *c, int poly_degree){
    return 0;
}
void MCMC_polynomial(float *x_obs, float *y_obs, int n_points, int n_steps, int poly_degree){
    float* c = new float[poly_degree];
    for(int i = 0; i < poly_degree; i ++){
            c[i] = 0.0;
        }
    cout << 0 << " " << 0 << " " << 0 << " " << 0 << endl;
    float sigma = 0.1;
    for(int i = 0; i < n_steps; i ++ ){
        float *propuestas = new float[poly_degree];
        for(int j = 0; j <= poly_degree; j ++){
            srand48(time(0));
            propuestas[j] = c[j] + (0.2*drand48()-0.1);
        }
        float log_viejo = loglikelihood(x_obs,y_obs,n_points,c,poly_degree);
        float log_nuevo = loglikelihood(x_obs,y_obs,n_points,propuestas,poly_degree);
        float aceptacion = fmin(1.0,exp(log_nuevo-log_viejo));
        srand48(time(0));
        float random = drand48();
        if(random < aceptacion){
            for(int j = 0; j < 4; j ++){
                c[j] = propuestas[j];               
            }
        }
        cout << c[0] << " " << c[1] << " " << c[2] << " " << c[3] << endl;
    }
}

