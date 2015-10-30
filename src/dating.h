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

double without_constraint(int,int,int* &,int* &,int* &,double* &,double* &,double* &,double&,double);

double without_constraint_rate(double,int,int,int* &,int* &,int* &,double* &,double* &,double* &);

double starting_point_rate(double,int,int,int* &,int* &,int* &,double* &,double* &,double* &,list<int>&);

list<double> without_constraint_active_set_rate(double,int,int,bool* &,int* &,int* &,int* &,double* &,double* &,double* &,list<int>&,double&);
  
int remove_ne_lambda(list<double> &,list<int>&);

bool conditions(int,list<double>&,int* &,double* &);

double with_constraint_active_set_rate(double,int,int,int* &,int* &,int* &,double* &,double* &,double* &);

double starting_point(int,int,int* &,int* &,int* &,double* &,double* &,double* &,double&,double,list<int>&);

list<double> without_constraint_active_set(int,int,bool* &,int* &,int* &,int* &,double* &,double* &,double* &,double&,double,list<int>&,double&);

double with_constraint_active_set(int,int,int* &,int* &,int* &,double* &,double* &,double* &,double&,double);
