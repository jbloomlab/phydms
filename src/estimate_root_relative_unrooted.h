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
#include "estimate_root_relative_rooted.h"
#include "dating_relative.h"
using namespace std;
     
double estimate_root_without_constraint_unrooted_relative(int n,int,int* &,double* &,double* &,double* &,int&,double&,double&,double,double);

double estimate_root_without_constraint_rate_unrooted_relative(double,int,int,int* &,double* &,double* &,double* &,int&,double&,double);

double estimate_root_with_constraint_unrooted_relative(int,int,int* &,double* &,double* &,double* &,int&,double& ,double&,double,double);

double estimate_root_with_constraint_rate_unrooted_relative(double,int,int,int* &,double* &,double* &,double* &,int&,double&,double);

