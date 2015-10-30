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
#include "estimate_root_rooted.h"
#include "estimate_root_unrooted.h"
#include "dating_relative.h"
using namespace std;
/////////LD ////////////////////
double without_constraint_lambda_relative(int,int,double,int* &,int* &,int* &,double* &,double* &,double* &,double&,double,double&,double);
double estimate_root_without_constraint_local_rooted_relative(int,int,int*&,int*&,int* &,double* &,double* &,double* &,int&,double&,double&,double,double);
double estimate_root_without_constraint_rooted_relative(int,int,int,int,int* &,double* &,double* &,double* &,int&,double&,double&,double,double);
////////LD given rate///////////
double without_constraint_lambda_rate_relative(double,int,int,double,int* &,int* &,int* &,double* &,double* &,double* &,double&,double);
double estimate_root_without_constraint_local_rate_rooted_relative(double,int,int,int*&,int*&,int* &,double* &,double* &,double* &,int&,double&,double);
double estimate_root_without_constraint_rate_rooted_relative(double,int,int,int,int,int* &,double* &,double* &,double* &,int &,double &,double);
//////QPD///////////////////////     
double starting_point_lambda_relative(int,int,double,int* &,int* &,int* &,double* &,double* &,double* &,double,double,list<int>&,bool* &,double mrca);
list<double> without_constraint_active_set_lambda_relative(int,int,double,int* &,bool* &,int* &,int* &,double* &,double* &,double* &,double,double,list<int>&,double,double);
double with_constraint_active_set_lambda_relative(int,int,double,int* &,int* &,int* &,double* &,double* &,double* &,double&,double,double&,double);
double estimate_root_with_constraint_local_rooted_relative(int,int,int*&,int*&,int* &,double* &,double* &,double* &,int&,double&,double&,double,double);
double estimate_root_with_constraint_fast_rooted_relative(int,int,int,int,int* &,double* &,double* &,double* &,int&,double&,double&,double,double);
double estimate_root_with_constraint_rooted_relative(int,int,int,int,int* &,double* &,double* &,double* &,int&,double&,double&,double,double);
//////QPD given rate///////////////////////     
double starting_point_lambda_rate_relative(double,int,int,double,int* &,int* &,int* &,double* &,double* &,double* &,double&,list<int>&,double mrca);
list<double> without_constraint_active_set_lambda_rate_relative(double,int,int,double,int* &,bool* &,int* &,int* &,double* &,double* &,double* &,double&,list<int>&,double&,double);
double with_constraint_active_set_lambda_rate_relative(double,int,int,double,int* &,int* &,int* &,double* &,double* &,double* &,double&,double);
double estimate_root_with_constraint_local_rate_rooted_relative(double,int,int,int*&,int*&,int* &,double* &,double* &,double* &,int&,double&,double);
double estimate_root_with_constraint_fast_rate_rooted_relative(double,int,int,int,int,int* &,double* &,double* &,double* &,int &,double &,double);
double estimate_root_with_constraint_rate_rooted_relative(double,int,int,int,int,int* &,double* &,double* &,double* &,int&,double&,double);
