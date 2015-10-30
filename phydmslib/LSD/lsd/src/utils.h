#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <stack>
#include <list>
#include <cstdlib>
#include <string.h>
#include "stdarg.h"
#include <math.h>
#ifndef maxIter
#define maxIter 1000
#endif
using namespace std;
  
string readWord(FILE *,string);
char readChar(FILE *,string);
double readdouble(FILE *,string);
int readInt(FILE *,string);
string readTaxon(FILE *);
double str2Double(char*);
double readResult(char*,char*,char*,int);
string readSupport(FILE *,string);
void concat(list<int>&,list<int>);
void concatPos(list<int>,list<int>&);
void concat(stack<int>&,list<int>);
list<int> arrayToList(int* &,list<int>,int);
void listToArray(list<int>,int* &);
void listToArray(list<string>,string* &);
void listToArray(list<double>,double* &);
void listToArray(list<string>&,list<string>&,string* &);
int getPosition(string* &,string,int,int);
bool comparecs(char*,char*);
bool comparec(char*,char*,int);
double minV(double* &,int);
double Min(double,double);
double Min(double,double,double);
double Max(double,double);
double Maxi(double* &,int);
double myabs(double);
bool isAncestor(int* &,int,int);
int lca(list<int>, list<int>);
int mrca(int* &,int,int);
int index(list<int>,int);
int index(string,string* &,int);
bool contain(string,list<string>);
string readLabel(FILE *,FILE *);
string readLabel(char c,FILE *);
char readBracket(FILE *,string);
char read2P(FILE *,string);
void counting(string,int&,int&,int);
bool finish(bool* &,int);
int choose(bool* &,int* &,int* &,int);
int choose(bool* &,list<int>*,int);
void computeSuc(int* &,int* &,int* &,int,int);
int computeSuc_unrooted(int* &,int* &,int* &,int,int);
void computeSuc(int* &,list<int>*,list<int>*,int,int);
list<int> pos(int,int* &,int* &,int* &,int);
list<int> pos_NB(int,int* &,list<int>*,int);
list<int> postorder(int* &,int* &,int* &,int);
list<int> postorder_NB(int* &,list<int>*,int);
list<int> postorder_unrooted(int,int* &,int* &,int* &,int);
list<int> pre(int,int* &,int* &,int* &,int);
list<int> preorder(int* &,int* &,int* &,int);
list<int> preorder_unrooted(int r,int* &,int* &,int* &,int);
void myExit(string, ... );
void myErrorMsg(string, ... );
bool isReal( const char*);
bool isInteger( const char*);
void newicktree(int,int* &,int* &,int* &,string* &,double* &,string* &,FILE*);
void nexustree(int,int* ,int*,int* ,string* ,double*,string* ,double* ,FILE*);
void nexustreeIC(int,int*,int*,int*,string*,double*,string*,double*,double*,double*,FILE *);
bool contain(int,list<int>);
list<int> sub(int,int,int,int* &,int* &,int* &);
void sort(int* &,int);
int index(int* &,int,int);
void subTree(int,int,int,int* &,int* &,int* &,double* &,string* &,string* &,list<int>&,list<double>&,list<string>&,list<string>&);
int mrca(int* &,list<int>);
void leaves(int,int* &,int* &,int,int,int);
double phi(int,int,int* &,double* &,double* &,double* &,double);
void computeSuc(int,int,int* &,int* &,bool* &,bool* &,list<int>&,list<int>&);
list<int> suc(int,int,int* &,int* &,int* &,bool* &,bool* &,double* &);
list<int> suc(int,int,int,int* ,int* ,int* ,bool* ,bool* ,int* &,list<int>&,list<int>&);
list<int> suc(int,int,int,int* ,double* &,int* ,int* ,bool* ,bool* &,int* &,list<int>&,list<int>&);
void computeFeuilles(int,int,int* &,int* &,int* &,double* &,bool* &,bool* &,stack<int>* &,list<int>);
void reduceTree(int,int,int* &,int* &,int* &,int* &,list<int>* &,list<int>* &,bool* ,bool* ,list<int>*&,list<int>);
void reduceTree(int,int,int* &,double* &,int* &,int* &,int* &,list<int>* &,list<int>* &,bool* ,bool* &,list<int>*&,list<int>);
double root_mesure(int,int,int* &,double* &,char**,int* &,double* &,char**);
void nexus(FILE *,char**,int,int);
double reroot_unrootedtree(int,int,int* &,double* &,double* &,int* &,double* &,double* &);
void reroot_unrootedtree(int,int,double lambda,int* &,double* &,double* &,int* &,double* &,double* &);
void reroot_unrootedtree(int,int,double lambda,int* &,double* &,string* &,double* &,int* &,double* &,string* &,double* &);
double reroot_rootedtree(int,int,int* &,int,int,double* &,int* &,double* &);
double reroot_rootedtree(int,int,int* &,int,int,double* &,int* &,double* &,int* &,int* &);
void reroot_rootedtree(int,int,double,int* &,int,int,double* &,int* &,double* &);
void reroot_rootedtree(int,int,double,int* &,int,int,double* &,string* &,int* &,double* &,string* &);
void rooted2unrooted(int,int,int* &,int* &,int* &,double* &,double* &,int* &,double* &,double* &);
void rooted2unrooted(int,int,int* &,int* &,int* &,double* &,string* &);
void reroot(int,int,int,int* &,int&,int* &,int* &,double* &,string* &,int* &,int* &,int* &,double* &,string* &);
void unrooted2rooted(int,int,int,double,int* &,int* &,int* &,double* &,double* &,int* &,int* &,int* &,double* &,double* &);
double unrooted2rooted(int,int,int,int* &,int* &,int* &,double* &,double* &,int* &,int* &,int* &,double* &,double* &);
double phi(int,int,int* &,double* &,double* &,double* &,double);
double phi(int,int,double &,int* &,double* &,double* &,double* &);
void output(bool,bool,FILE*,FILE*,FILE*,FILE*,int,int,double,int*,int*,int*,double*,double*,string*,string*,double);
double* variance(bool,int,double*,int,int);
list<string> getOutgroup(string);
