#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <stack>
#include <list>
#include "stdarg.h"
#include <math.h>
#include "utils.h"
#include "dating.h"
using namespace std;
//////LD////////////////
double without_constraint_lambda(int,int,double,int* &,int* &,int* &,double* &,double* &,double* &,double&,double,double&);
double estimate_root_without_constraint_local_rooted(int,int,int*&,int*&,int* &,double* &,double* &,double* &,int&,double&,double&,double);
double estimate_root_without_constraint_rooted(int,int,int,int,int* &,double* &,double* &,double* &,int&,double&,double&,double);
/////LD given rate///////
double without_constraint_lambda_rate(double,int,int,double,int* &,int* &,int* &,double* &,double* &,double* &,double&);
double estimate_root_without_constraint_local_rate_rooted(double,int,int,int*&,int*&,int* &,double* &,double* &,double* &,int&,double&);
double estimate_root_without_constraint_rate_rooted(double,int,int,int,int,int* &,double* &,double* &,double* &,int&,double&);
//////QPD/////////////    
double starting_point_lambda(int,int,double,int* &,int* &,int* &,double* &,double* &,double* &,double,double,list<int>&);
list<double> without_constraint_active_set_lambda(int,int,double,int* &,bool* &,int* &,int* &,double* &,double* &,double* &,double,double,list<int>&,double);
bool conditions_lambda(int,list<double>&,int* &,double* &,double);
double with_constraint_active_set_lambda(int,int,double,int* &,int* &,int* &,double* &,double* &,double* &,double&,double,double&);
double estimate_root_with_constraint_local_rooted(int,int,int*&,int*&,int* &,double* &,double* &,double* &,int&,double&,double&,double);
double estimate_root_with_constraint_fast_rooted(int,int,int,int,int* &,double* &,double* &,double* &,int&,double&,double&,double);
double estimate_root_with_constraint_rooted(int,int,int,int,int* &,double* &,double* &,double* &,int&,double&,double&,double);
/////QPD given rate///////////
double starting_point_lambda_rate(double,int,int,double,int* &,int* &,int* &,double* &,double* &,double* &,double&,list<int>&);
list<double> without_constraint_active_set_lambda_rate(double,int,int,double,int* &,bool* &,int* &,int* &,double* &,double* &,double* &,double&,list<int>&,double&);
double with_constraint_active_set_lambda_rate(double,int,int,double,int* &,int* &,int* &,double* &,double* &,double* &,double&);
double estimate_root_with_constraint_local_rate_rooted(double,int,int,int*&,int*&,int* &,double* &,double* &,double* &,int&,double&);
double estimate_root_with_constraint_fast_rate_rooted(double,int,int,int,int,int*&,double* &,double* &,double* &,int&,double&);
double estimate_root_with_constraint_rate_rooted(double,int,int,int,int,int* &,double* &,double* &,double* &,int&,double&);
