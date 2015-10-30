#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <stack>
#include <list>
#include "stdarg.h"
#include <math.h>
#include "utils.h"

using namespace std;

double without_constraint_relative(int,int,int* &,int* &,int* &,double* &,double* &,double* &,double&,double,double);

double without_constraint_rate_relative(double,int,int,int* &,int* &,int* &,double* &,double* &,double* &,double);

double starting_point_rate_relative(double,int,int,int* &,int* &,int* &,double* &,double* &,double* &,list<int>&,bool* &,double);

list<double> without_constraint_active_set_rate_relative(double,int,int,bool* &,int* &,int* &,int* &,double* &,double* &,double* &,list<int>&,double&,double);
  
double with_constraint_active_set_rate_relative(double,int,int,int* &,int* &,int* &,double* &,double* &,double* &,double);

double starting_point_relative(int,int,int* &,int* &,int* &,double* &,double* &,double* &,double&,double,list<int>&,bool* &,double);

list<double> without_constraint_active_set_relative(int,int,bool* &,int* &,int* &,int* &,double* &,double* &,double* &,double&,double,list<int>&,double&,double);

double with_constraint_active_set_relative(int,int,int* &,int* &,int* &,double* &,double* &,double* &,double&,double,double);
