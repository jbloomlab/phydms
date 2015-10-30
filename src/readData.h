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

bool tree2data(FILE *,string,int,int,int* &,double* &,string* &,double* &,string* &);
void tree2dataS(FILE *,int,int,int* &,double* &,string* &,string* &);
void extrait_outgroup(string,string,list<string>&,int);

