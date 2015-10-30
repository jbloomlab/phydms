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
using namespace std;
    
double estimate_root_without_constraint_unrooted(int,int,int* &,double* &,double* &,double* &,int&,double&,double&,double);

double estimate_root_without_constraint_rate_unrooted(double,int,int,int* &,double* &,double* &,double* &,int&,double&);

double estimate_root_with_constraint_unrooted(int,int,int* &,double* &,double* &,double* &,int&,double& ,double&,double);

double estimate_root_with_constraint_rate_unrooted(double,int,int,int* &,double* &,double* &,double* &,int&,double&);
