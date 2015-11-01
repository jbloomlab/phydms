//    LSD - Least Square Dating for etimating substitution rate and divergence dates
//    Copyright (C) <2015> <Thu-Hien To>

//    This program is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.

//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.

//    You should have received a copy of the GNU General Public License
//    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#include "utils.h"
  
string readWord(FILE *f,string fn){
       string s;
       char c=readChar(f,fn);
       int i=0;
       while (c<=33 || c>=126) c=readChar(f,fn);
       while (c>=33 && c<=126){
             s=s+c;
             i++;
             c=getc(f);
       }
       return s;
}

char readChar(FILE *f,string fn){
    char r;
    if (fscanf(f,"%c",&r)==1) return r;
    else {
        cout<<"Error in file "<<fn<<": file ended unexpected"<<endl;
        exit(EXIT_FAILURE);
    }
}

double readdouble(FILE *f,string fn){
    double r;
    if (fscanf(f,"%lf",&r)==1) return r;
    else {
        cout<<"Error in the file "<<fn<<": real expected"<<endl;
        exit(EXIT_FAILURE);
    }
}

int readInt(FILE *f,string msg){
    int d;
    if (fscanf(f,"%d",&d)==1) return d;
    else {
        cout<<msg<<endl;
        exit(EXIT_FAILURE);
    }
}

double readResult(string result,string rate,string date,int nb_tree){
       FILE *f = fopen(result.c_str(),"rt");
       if (f!=NULL){
          FILE *w = fopen(rate.c_str(),"wt");  
          FILE *wd = fopen(date.c_str(),"wt");
          string s= readWord(f,result);
          while (s.compare("Exit)")!=0) s = readWord(f,result);          
          for (int y=0;y<nb_tree;y++){
              s = readWord(f,result);
              while (s.compare("rate")!=0) s = readWord(f,result);
              double r=readdouble(f,result);
              while (s.compare("tMRCA")!=0) s = readWord(f,result);
              double t=readdouble(f,result);
              fprintf(w,"%f\n",r);
              fprintf(wd,"%f\n",t);
          }
          while (s.compare("time:")!=0) s = readWord(f,result);       
          double time=readdouble(f,result);
          fclose(f);
          fclose(w);
          fclose(wd);
          return time;
       }
       else{
            cout<<"Impossible to open file "<<result<<endl;
            exit (EXIT_FAILURE);
       }
}

string readSupport(FILE *f,string fn){
    string s="";
    char c=readChar(f,fn);
    while (c!=':' && c!=';') {
        s=s+c;
        c=readChar(f,fn);
    }
    return s;
}

void concat(list<int> & l1,list<int> l2){
     for (list<int>::iterator iter = l2.begin();iter != l2.end();iter++){
         l1.push_back(*iter);
     }
}

void concatPos(list<int> l1,list<int> &l2){
     for (list<int>::iterator iter = l1.begin();iter != l1.end();iter++){
         l2.push_back(*iter);
     }
}

void concat(stack<int> & l1,list<int> l2){
     for (list<int>::iterator iter = l2.begin();iter != l2.end();iter++){
         l1.push(*iter);
     }
}



list<int> arrayToList(int* & arr,list<int> li,int n){
     for (int i=0;i<n;i++){
         li.push_back(arr[i]);
     }
     return li;
}
void listToArray(list<int> li,int* & arr){
     int count=0;
     for (list<int>::iterator i=li.begin();i!=li.end();i++){
         arr[count]=(*i);
         count++;
     }
}
void listToArray(list<string> li,string* & arr){
     int count=0;
     for (list<string>::iterator i=li.begin();i!=li.end();i++){
         arr[count]=(*i);
         count++;
     }
}

void listToArray(list<double> li,double* & arr){
     int count=0;
     for (list<double>::iterator i=li.begin();i!=li.end();i++){
         arr[count]=(*i);
         count++;
     }
}
void listToArray(list<string> & inner,list<string> & leaves,string* & arr){ 
     int count=0;     
     for (list<string>::iterator i=inner.begin();i!=inner.end();i++){
              arr[inner.size()-count-1]=(*i);
              count++;
     }     
     for (list<string>::iterator i=leaves.begin();i!=leaves.end();i++){
         arr[count]=(*i);
         count++;
     }
}
int getPosition(string* & arr,string s,int n,int m){
    int i = n;
    int count=0;
    int k=0;
    while (i<m && count<2){
        if (arr[i].compare(s)==0){
           count++;
           k=i; 
        }
        i++;
    }
    if (count==0) return -1;
    else if (count>1){
        cout<<"There are at least two leaves that have the same label"<<endl;
        exit(EXIT_FAILURE);
    }
    else return k;
}

bool comparecs(char* c,string s){
     for (int i=0;i<s.size();i++) {
         if (c[i]!=s.at(i)) return false;
     }
     return true;
}
bool comparec(char* c,char* s,int n){
     for (int i=0;i<n;i++) {
         if (c[i]!=s[i]) return false;
     }
     return true;
}

double minV(double* & B,int n){
       int i=1;
       while (B[i]==0) i++;
       double M=B[i];
       while (i<n){
           if (B[i]!=0 && B[i]<M) M=B[i];
           i++;
       }
       return M;           
}

double Min(double a,double b){
       if (a<b) return a;
       else return b;
}

double Min(double a,double b,double c){
       if (a<b && a<c) return a;
       else if (b<a && b<c) return b;
       else return c;
}

double Max(double a,double b){
       if (a<b) return b;
       else return a;
}

double Maxi(double* & B,int n){
	   double max=0;
	   for (int i=0;i<n;i++){
	   	   if (max<B[i]) max=B[i];
	   }
	   return max;
}

double myabs(double a){
       if (a>=0) return a;
       else return -a;
}

bool isAncestor(int* & P,int i,int j){
     int x=j;
     while (x!=-1){
           if (x==i) return true;
           else x=P[x];
     }
     return false;
}

int lca(list<int> a, list<int> b){
    for (list<int>::iterator ia=a.begin();ia!=a.end();ia++){
        for (list<int>::iterator ib=b.begin();(ib!=b.end() && *ib>=*ia);ib++){
            if (*ia==*ib){return *ia;}
        }
    }
    return 0;
}

int mrca(int* & P,int i,int j){
    int c = i;
    while (!isAncestor(P,c,j)){
          c=P[c];
    }    
    return c;
}

bool decrease(int a,int b){
     if (a>=b) return true;
     else return false;
}

int index(list<int> & L, int e){
    int ind=0;
    for (list<int>::iterator i=L.begin();i!=L.end();i++){
        if (*i==e) return ind;
        ind++;
    }
    return -1;
}

int index(string s,string* & L,int n){
    for (int i=0;i<n;i++){
        if (s.compare(L[i])==0){ 
           return i;
        }
    }
    return -1;
}


bool contain(string s,list<string> l){
     for (list<string>::iterator iter = l.begin();iter!=l.end();iter++){
         if (s.compare(*iter)==0) return true;
     }
     return false;
}

string readLabel(FILE *f,FILE *w){       
       string s="";
       char c=readChar(f,"input tree");
       fprintf(w,"%c",c);
       while (c!=':'){
             s=s+c;
             fscanf(f,"%c",&c);
           fprintf(w,"%c",c);
       }
       return s.c_str();
}

string readLabel(char ch,FILE *f){       
       string s="";
       s=s+ch;
       char c=readChar(f,"input tree");
       while (c!=':'){
             s=s+c;
             fscanf(f,"%c",&c);
       }
       return s.c_str();
}


char readBracket(FILE *f,string fn){
	 char c=readChar(f,fn);
	 while (c!='('){
	 	   c=readChar(f,fn);
     }
     return c;
}
char read2P(FILE *f,string fn){
	 char c= readChar(f,fn);
	 if (c==';') return c;
	 while (c!=':' && c!=';') c=readChar(f,fn);
    return c;
}
/*
char readBracketC(FILE *f,string fn){
	 char c=readChar(f,fn);
	 while (c!=')'){
	 	   c=readChar(f,fn);
     }
     return c;
}




string readTaxon(FILE *f){
	   string s="";
	   char c=readChar(f);
	   while (c!=':') {
	   		 s+=c;
	   		 c=readChar(f);
	   }
	   return s.c_str();		 
}*/

void counting(string fn,int &n,int &m,int nb_tree){
     FILE * tree = fopen(fn.c_str(),"rt");
     if (tree==NULL) cout<<"Impossible to open tree file"<<endl;
     else{ 
          int n1,m1; 
          for (int i=0;i<nb_tree;i++){ 
              char c = readBracket(tree,fn); 
              n=0;//internal nodes
              int l=0;//leaves
              int a=0;
              while (c!=';'){
                 if (c=='(') {n++;a++;}
                 if (c==':' && a>0) l++;
                 if (c==')') {
                    a--;
                    if (a>0) read2P(tree,fn);
                 }       
                 c=readChar(tree,fn); 
              }
              if (a!=0) {
                  cout<<"Errors in the input trees"<<endl;
                  exit(EXIT_FAILURE);
              }
              m=n+l-1;//number of branches
              if (i>0 && n1!=n && m1!=m){ 
                 cout<<"The trees do not have the same size"<<endl;
                 exit (EXIT_FAILURE);
              }
              if (i==0){
                 n1=n;
                 m1=m;
              }
          }
          fclose(tree);
     }     
}


bool finish(bool* & flag,int size){
     for (int i=0;i<size;i++){
         if (flag[i]==false) {return false;}
     }
     return true;
}
int choose(bool* & flag,int* & Suc1,int* & Suc2,int size){
    for (int i=0;i<size;i++){
        if (!flag[i] && flag[Suc1[i]] && flag[Suc2[i]]) return i;
    }
    return -1;
}
int choose(bool* & flag,list<int>* SucI,int size){
    for (int i=0;i<size;i++){
        bool check = true;
        if (!flag[i]){
           for (list<int>::iterator iter = SucI[i].begin();iter!=SucI[i].end();iter++){
               check = check && flag[*iter];
           }
           if (check) return i;
        }
    }
    return -1;
}
void computeSuc(int* & Pre,int* & Suc1,int* & Suc2,int size,int n){
     for (int i=0;i<n;i++){
         Suc1[i]=-1;
     }
     for (int i=0;i<size;i++){
         if (Pre[i]!=-1){
            if (Suc1[Pre[i]]==-1) Suc1[Pre[i]]=i;
            else if (Suc1[Pre[i]]>i){
                 Suc2[Pre[i]]=Suc1[Pre[i]]; 
                 Suc1[Pre[i]]=i;
            }
            else Suc2[Pre[i]]=i;
         }
     }
}

int computeSuc_unrooted(int* & Pre,int* & Suc1,int* & Suc2,int size,int n){
     for (int i=0;i<n;i++){
         Suc1[i]=-1;
     }
     int s=-1;
     for (int i=0;i<size;i++){
         if (Pre[i]!=-1 && Pre[Pre[i]]==-1 && s==-1){
             s=i;
         }
         else if (Pre[i]!=-1){
              if (Suc1[Pre[i]]==-1) Suc1[Pre[i]]=i;
              else Suc2[Pre[i]]=i;
         }
     }
    return s;
}

void computeSuc(int* & P,list<int>* SucL,list<int>* SucI,int n,int m){
     for (int i=1;i<n;i++) SucI[P[i]].push_back(i);
     for (int i=n;i<=m;i++) SucL[P[i]].push_back(i);
}

void reduce(list<int>* A,list<int>* B,list<int> & Label,int size){
     for (int i=0;i<size;i++){
         B[i].push_back(A[i].front());
         Label.push_back(A[i].front());
         for (int j=0;j<size;j++){
             if (i!=j){
                int d = lca(A[i],A[j]);       
                B[i].push_back(d);
                Label.push_back(d);
             }
         }
     }
     Label.sort();
     Label.unique();
     for (int i=0;i<size;i++){
         B[i].sort();         
         B[i].unique();
     }
}

list<int> pos(int i,int* & P,int* & Suc1,int* & Suc2,int n){
          list<int> l;
          if (i>=n) return l;
          else{
               list<int> l1 = pos(Suc1[i],P,Suc1,Suc2,n);
               list<int> l2 = pos(Suc2[i],P,Suc1,Suc2,n);          
               for (list<int>::iterator iter=l1.begin();iter!=l1.end();iter++){             
                   l.push_back(*iter);
               }
               for (list<int>::iterator iter=l2.begin();iter!=l2.end();iter++){             
                   l.push_back(*iter);
               }          
               l.push_back(i);
               return l;
          }
}

list<int> pos_NB(int i,int* & P,list<int>* SucI,int n){
          list<int> l;
          if (i>=n) return l;
          else{
               for (list<int>::iterator iter = SucI[i].begin();iter!=SucI[i].end();iter++){
                   list<int> l1 = pos_NB(*iter,P,SucI,n);
                   for (list<int>::iterator iter1=l1.begin();iter1!=l1.end();iter1++){             
                       l.push_back(*iter1);
                   }
               }                         
               l.push_back(i);
               return l;
          }
}


list<int> postorder(int* & P,int* & Suc1,int* & Suc2,int n){
      int root=0;
      for (root=0;root<n;root++){
          if (P[root]==-1) break;
      }
      return pos(root,P,Suc1,Suc2,n);
}

list<int> postorder_NB(int* & P,list<int>* SucI,int n){
      int root=0;
      for (root=0;root<n;root++){
          if (P[root]==-1) break;
      }
      return pos_NB(root,P,SucI,n);
}

list<int> postorder_unrooted(int r,int* & P,int* & Suc1,int* & Suc2,int n){
          list<int> l1=pos(r,P,Suc1,Suc2,n);
          list<int> l2=pos(P[r],P,Suc1,Suc2,n);
          list<int> l;
          for (list<int>::iterator iter=l1.begin();iter!=l1.end();iter++){             
              l.push_back(*iter);
          }
          for (list<int>::iterator iter=l2.begin();iter!=l2.end();iter++){             
              l.push_back(*iter);
          }        
          return l;          
}

list<int> pre(int i,int* & P,int* & Suc1,int* & Suc2,int n){
          list<int> l;
          if (i>=n) return l;
          else{
               l.push_back(i);
               list<int> l1 = pre(Suc1[i],P,Suc1,Suc2,n);
               list<int> l2 = pre(Suc2[i],P,Suc1,Suc2,n);          
               for (list<int>::iterator iter=l1.begin();iter!=l1.end();iter++){             
                   l.push_back(*iter);
               }
               for (list<int>::iterator iter=l2.begin();iter!=l2.end();iter++){             
                   l.push_back(*iter);
               }          
               return l;
          }
}

list<int> preorder(int* & P,int* & Suc1,int* & Suc2,int n){
      int root=0;
      for (root=0;root<n;root++){
          if (P[root]==-1) break;
      }
      return pre(root,P,Suc1,Suc2,n);
}

list<int> preorder_unrooted(int r,int* & P,int* & Suc1,int* & Suc2,int n){
          list<int> l1=pre(P[r],P,Suc1,Suc2,n);
          list<int> l2=pre(r,P,Suc1,Suc2,n);
          list<int> l;
          for (list<int>::iterator iter=l1.begin();iter!=l1.end();iter++){             
              l.push_back(*iter);
          }
          for (list<int>::iterator iter=l2.begin();iter!=l2.end();iter++){             
              l.push_back(*iter);
          }        
          return l;          
}

void myExit( string msg, ... )
{
	va_list ptr;
	fprintf( stderr, "Error: " );
	va_start( ptr, msg );
	vfprintf( stderr, msg.c_str(), ptr );
	va_end( ptr );
	//fflush( NULL );
	exit( EXIT_FAILURE );
}


void myErrorMsg( string msg, ... )
{
	va_list ptr;
	fprintf( stderr, "Error: " );
	va_start( ptr, msg );
	vfprintf( stderr, msg.c_str(), ptr );
	va_end( ptr );

}

bool isReal( const char* str )
{
	while( *str!='\0' )
	{
		if( !( ('0'<=*str && *str<='9') || *str=='e' || *str=='E' || *str=='.' || *str=='-' ) )
			return false;
		str++;
	}
	return true;
}

bool isInteger( const char* str )
{
	if( *str=='-' )
		str++;
	while( *str!='\0' )
	{
		if( *str<'0' || '9'<*str )
			return false;
		str++;
	}
	return true;
}

void printSquareMatrixDouble( double** m, int size )
{
	for( int i=0 ; i<size ; i++ )
	{
		printf("%f",m[i][0]);
		for( int j=1 ; j<size ; j++ )
			printf("\t%f",m[i][j]);
		printf("\n");
	}
}

void newicktree(int n,int* & P,int* & Suc1,int* & Suc2,string* & Labels,double* & B,string* & Support,FILE *w){
    list<int> post_order=postorder(P,Suc1,Suc2,n);
    int r=0;
    for (list<int>::iterator iter = post_order.begin();iter!=post_order.end();iter++){
        int i = *iter;
        int s1 = Suc1[i];
        int s2 = Suc2[i];
        ostringstream st1,st2;
        st1<< B[s1];
        st2<< B[s2];
        Labels[i]="("+Labels[s1]+Support[s1]+":"+st1.str()+","+Labels[s2]+Support[s2]+":"+st2.str()+")";
        r=i;
    }
    fprintf(w,"%s;\n",Labels[r].c_str());
}

void nexustree(int n,int* P,int* Suc1,int* Suc2,string* Labels,double* B,string* Support,double* T,FILE *w){
    list<int> post_order=postorder(P,Suc1,Suc2,n);
    int r=0;
    for (list<int>::iterator iter = post_order.begin();iter!=post_order.end();iter++){
        int i = *iter;
        int s1 = Suc1[i];
        int s2 = Suc2[i];
        ostringstream st1,st2,date1,date2;
        st1<<B[s1];
        st2<<B[s2];
        date1<<T[s1];
        date2<<T[s2];
        if (s1<n && s2>=n) {
            Labels[i]="("+Labels[s1]+Support[s1]+"[&date="+date1.str()+"]:"+st1.str()+","+Labels[s2]+"[&date="+date2.str()+"]:"+st2.str()+")";
        }
        else if (s1>=n && s2<n){
            Labels[i]="("+Labels[s1]+"[&date="+date1.str()+"]:"+st1.str()+","+Labels[s2]+Support[s2]+"[&date="+date2.str()+"]:"+st2.str()+")";
        }
        else if (s1<n && s2<n){
            Labels[i]="("+Labels[s1]+Support[s1]+"[&date="+date1.str()+"]:"+st1.str()+","+Labels[s2]+Support[s2]+"[&date="+date2.str()+"]:"+st2.str()+")";
        }
        else Labels[i]="("+Labels[s1]+"[&date="+date1.str()+"]:"+st1.str()+","+Labels[s2]+"[&date="+date2.str()+"]:"+st2.str()+")";
        r=i;
    }
    ostringstream root;
    root<<T[r];
    fprintf(w,"%s;\n",(Labels[r]+"[&date="+root.str()+"]").c_str());
}

void nexustreeIC(int n,int* P,int* Suc1,int* Suc2,string* Labels,double* B,string* Support,double* T,double* T_min,double* T_max,FILE *w){
    list<int> post_order=postorder(P,Suc1,Suc2,n);
    int r=0;
    for (list<int>::iterator iter = post_order.begin();iter!=post_order.end();iter++){
        int i = *iter;
        int s1 = Suc1[i];
        int s2 = Suc2[i];
        ostringstream st1,st2,date1,date2,tmin1,tmax1,tmin2,tmax2;
        st1<<B[s1];
        st2<<B[s2];
        date1<<T[s1];
        date2<<T[s2];
        tmin1<<T_min[s1];tmax1<<T_max[s1];
        tmin2<<T_min[s2];tmax2<<T_max[s2];
        if (s1<n && s2>=n) {
            Labels[i]="("+Labels[s1]+Support[s1]+"[&date="+date1.str()+",IC=\""+date1.str()+"("+tmin1.str()+","+tmax1.str()+")\"]:"+st1.str()+","+Labels[s2]+"[&date="+Labels[s2]+"("+date2.str()+")]:"+st2.str()+")";
        }
        else if (s1>=n && s2<n){
            Labels[i]="("+Labels[s1]+"[&date="+Labels[s1]+"("+date1.str()+")]:"+st1.str()+","+Labels[s2]+Support[s2]+"[&date="+date2.str()+",IC=\""+date2.str()+"("+tmin2.str()+","+tmax2.str()+")\"]:"+st2.str()+")";
        }
        else if (s1<n && s2<n){
            Labels[i]="("+Labels[s1]+Support[s1]+"[&date="+date1.str()+",IC=\""+date1.str()+"("+tmin1.str()+","+tmax1.str()+")\"]:"+st1.str()+","+Labels[s2]+Support[s2]+"[&date="+date2.str()+",IC=\""+date2.str()+"("+tmin2.str()+","+tmax2.str()+")\"]:"+st2.str()+")";
        }
        else Labels[i]="("+Labels[s1]+"[&date="+Labels[s1]+"("+date1.str()+")]:"+st1.str()+","+Labels[s2]+"[&date="+Labels[s2]+"("+date2.str()+")]:"+st2.str()+")";
        r=i;
    }
    ostringstream root,rmin,rmax;
    root<<T[r];
    rmin<<T_min[r];
    rmax<<T_max[r];
    fprintf(w,"%s;\n",(Labels[r]+"[&date="+root.str()+",IC=\""+root.str()+"("+rmin.str()+","+rmax.str()+")\"]").c_str());
}

bool contain(int i,list<int> l){
     for (list<int>::iterator iter=l.begin();iter!=l.end();iter++){
         if (i==*iter) return true;
     }
     return false;
}

list<int> sub(int r,int n,int m,int* & P,int* & Suc1,int* & Suc2){
          list<int> result;
          result.push_back(r);
          if (r<n){
               list<int> l1 = sub(Suc1[r],n,m,P,Suc1,Suc2);
               list<int> l2 = sub(Suc2[r],n,m,P,Suc1,Suc2);
               for (list<int>::iterator iter = l1.begin();iter!=l1.end();iter++){
                   result.push_back(*iter);
               }
               for (list<int>::iterator iter = l2.begin();iter!=l2.end();iter++){
                   result.push_back(*iter);
               }
          }
          return result;
}
void sort(int* & tab,int size){
          for (int i=0;i<size;i++){
              for (int j=i;j<size;j++){
                  if (tab[i]>tab[j]){
                     int temp = tab[i];
                     tab[i]=tab[j];
                     tab[j]=temp;
                  }
              }
          }
}

int index(int* & tab,int value,int size){
    for (int i=0;i<size;i++){
        if (value==tab[i]) return i;
    }
    return -1;
}

void subTree(int n,int m,int r,int* & P,int* & Suc1,int* & Suc2,double* & B,string* & Support,string* & Labels,list<int> & P1,list<double> & B1,list<string> & Support1,list<string> &Labels1){
     list<int> nodes = sub(r,n,m,P,Suc1,Suc2);
     int *tab=new int[nodes.size()];
     listToArray(nodes,tab);
     sort(tab,(int)nodes.size());
    for (int i=0;i<nodes.size();i++){
         int j = index(tab,P[tab[i]],(int)nodes.size());
         P1.push_back(j);
         B1.push_back(B[tab[i]]);
         Support1.push_back(Support[tab[i]]);
         Labels1.push_back(Labels[tab[i]]);
     }
    delete[] tab;
}

int mrca(int* & P,list<int> taxa){
    int t=taxa.front();
    taxa.pop_front();
    bool flag = false;
    while (!flag && P[t]!=-1){
          t = P[t];
          flag=true;
          for (list<int>::iterator ia=taxa.begin();ia!=taxa.end();ia++){
              int j=*ia;
              if (!isAncestor(P,t,j)) {flag=false;break;}
          }
    }       
    return t;
}

void leaves(int n,int* & Suc1,int* & Suc2,int a,int & taxa1,int & taxa2){
     taxa1=Suc1[a];
     taxa2=Suc2[a];
     while (taxa1<n) taxa1 = Suc1[taxa1];
     while (taxa2<n) taxa2 = Suc2[taxa2];        
}


void computeSuc(int i,int n,int* & Suc1,int* & Suc2,bool* & flag,bool* & leaf,list<int> & SucL,list<int> & SucI){
     if (i<n){
        int s1=Suc1[i];
        int s2=Suc2[i];
        if (flag[s1]){
           if (leaf[s1]) SucL.push_back(s1);
           else SucI.push_back(s1);
        }
        else if (s1<n){
             list<int> SucL1,SucI1;
             computeSuc(s1,n,Suc1,Suc2,flag,leaf,SucL1,SucI1);
             concat(SucL,SucL1);
             concat(SucI,SucI1);
        }
        if (flag[s2]){
           if (leaf[s2]) SucL.push_back(s2);
           else SucI.push_back(s2);        
        }
        else if (s2<n){
             list<int> SucL2,SucI2;
             computeSuc(s2,n,Suc1,Suc2,flag,leaf,SucL2,SucI2);
             concat(SucL,SucL2);
             concat(SucI,SucI2);          
        }
     }
}

list<int> suc(int i,int n,int* & P,int* & Suc1,int* & Suc2,bool* & leaf,bool* & flag,double* & T){
//flag[i]=false
          list<int> result;
          result.push_back(i);
          leaf[i]=true;
          T[i]=T[P[i]];
          if (i<n){
                int s1=Suc1[i];
                int s2=Suc2[i];
                if (!flag[s1]){   
                   list<int> l1=suc(s1,n,P,Suc1,Suc2,leaf,flag,T);
                   concat(result,l1);
                }
                if (!flag[s2]){ 
                   list<int> l2=suc(s2,n,P,Suc1,Suc2,leaf,flag,T);
                   concat(result,l2);
                }             
          }
          return result;
}

list<int> suc(int i,int j,int n,int*  P,int*  Suc1,int*  Suc2,bool*  flag,bool*  leaf,int* & Pre,list<int> &sucL,list<int> &sucI){
          list<int> result;     
          int s1=Suc1[j];
          int s2=Suc2[j];
          if (flag[s1]){
             Pre[s1]=i;
             if (leaf[s1]) 
                sucL.push_back(s1);
             else
                 sucI.push_back(s1);             
          }
          else{
               list<int> l1=suc(i,s1,n,P,Suc1,Suc2,flag,leaf,Pre,sucL,sucI);
               concatPos(l1,result);
          }
          if (flag[s2]){
             Pre[s2]=i;
             if (leaf[s2])
                sucL.push_back(s2);
             else
                 sucI.push_back(s2);
          }
          else{
               list<int> l2=suc(i,s2,n,P,Suc1,Suc2,flag,leaf,Pre,sucL,sucI);
               concatPos(l2,result);
          }
          if (j!=i) result.push_back(j);
          return result;
}

list<int> suc(int i,int j,int n,int*  P,double* &T,int*  Suc1,int*  Suc2,bool*  flag,bool* & leaf,int* & Pre,list<int> &sucL,list<int> &sucI){
          list<int> result;     
          int s1=Suc1[j];
          int s2=Suc2[j];
          if (i==0 && j!=0){ leaf[j]=true; T[j]=T[0];}
          if (flag[s1]){
             Pre[s1]=i;
             if (leaf[s1]) 
                sucL.push_back(s1);
             else
                 sucI.push_back(s1);             
          }
          else{               
               list<int> l1=suc(i,s1,n,P,T,Suc1,Suc2,flag,leaf,Pre,sucL,sucI);
               concatPos(l1,result);
          }
          if (flag[s2]){
             Pre[s2]=i;
             if (leaf[s2])
                sucL.push_back(s2);
             else
                 sucI.push_back(s2);
          }
          else{
               list<int> l2=suc(i,s2,n,P,T,Suc1,Suc2,flag,leaf,Pre,sucL,sucI);
               concatPos(l2,result);
          }
          if (j!=i) result.push_back(j); 
          return result;
}

void computeFeuilles(int n,int m,int* & P,int* & Suc1,int* & Suc2,double* & T,bool* & flag,bool* & leaf,stack<int>* &feuilles,list<int> active_set){
     int count=0;
     for (list<int>::iterator iter=active_set.begin();iter!=active_set.end();iter++){
         int i = *iter;
         if (i>=n && i!=m+1 && i!=m+2){    
            int j=i;
            list<int> ai;
            while (!flag[j]){
                  feuilles[count].push(j);                                   
                  int k=j;
                  j=P[j];
                  leaf[j]=true;
                  T[j]=T[i];
                  if (Suc1[j]==k && !flag[Suc2[j]]){ 
                     ai.push_back(Suc2[j]);
                  }
                  if (Suc2[j]==k && !flag[Suc1[j]]){  
                     ai.push_back(Suc1[j]);
                  }
            }               
            for (list<int>::iterator it=ai.begin();it!=ai.end();it++){                
                concat(feuilles[count],suc(*it,n,P,Suc1,Suc2,leaf,flag,T));   
            }  
            count++;
         }
     }  
}

void reduceTree(int n,int m,int* & P,int* & Pre,int* & Suc1,int* & Suc2,list<int>* &SucL,list<int>* &SucI,bool*  flag,bool*  leaf,list<int>* &internal,list<int> active_set){
     int count=0;
     Pre[0]=-1;
     for (int i=0;i<n;i++){        
         int s1=Suc1[i];
         int s2=Suc2[i];
         if (flag[i]){            
            if ((!leaf[i]) && (!flag[s1] || !flag[s2])){
               internal[count]=suc(i,i,n,P,Suc1,Suc2,flag,leaf,Pre,SucL[i],SucI[i]);
               count++;
            }
            else {
                 Pre[s1]=i;
                 Pre[s2]=i;
                 if (leaf[s1]) 
                    SucL[i].push_back(s1);
                 else 
                      SucI[i].push_back(s1);
                 if (leaf[s2])
                    SucL[i].push_back(s2);
                 else
                     SucI[i].push_back(s2);
            }              
         }
         else if (leaf[i]){
              Pre[s1]=i;
              Pre[s2]=i;
              if (leaf[s1]) 
                 SucL[i].push_back(s1);
              else 
                   SucI[i].push_back(s1);
              if (leaf[s2])
                 SucL[i].push_back(s2);
              else
                  SucI[i].push_back(s2);
         }
     }
}

void reduceTree(int n,int m,int* & P,double* &T,int* & Pre,int* & Suc1,int* & Suc2,list<int>* &SucL,list<int>* &SucI,bool*  flag,bool*  &leaf,list<int>* &internal,list<int> active_set){
     int count=0;
     Pre[0]=-1;
     for (int i=0;i<n;i++){        
         int s1=Suc1[i];
         int s2=Suc2[i];
         if (flag[i]){
            if ((!leaf[i] || P[i]==-1) && (!flag[s1] || !flag[s2])){  
               internal[count]=suc(i,i,n,P,T,Suc1,Suc2,flag,leaf,Pre,SucL[i],SucI[i]);
               count++;
            }
            else {
                 Pre[s1]=i;
                 Pre[s2]=i;
                 if (leaf[s1]) 
                    SucL[i].push_back(s1);
                 else 
                      SucI[i].push_back(s1);
                 if (leaf[s2])
                    SucL[i].push_back(s2);
                 else
                     SucI[i].push_back(s2);
            }
         }
         else if (leaf[i]){
              Pre[s1]=i;
              Pre[s2]=i;
              if (leaf[s1]) 
                 SucL[i].push_back(s1);
              else 
                   SucI[i].push_back(s1);
              if (leaf[s2])
                 SucL[i].push_back(s2);
              else
                  SucI[i].push_back(s2);
         }

     } 

}


double root_mesure(int n,int m,int* & P1,double* & B1,string* & Labels1,int* & P2,double* & B2,string* & Labels2){
       double* R2T1=new double[m+1];
       double* R2T2=new double[m+1];
       for (int i=n;i<=m;i++){
           R2T1[i]=B1[i];
           int j1=i;
           while (P1[j1]!=-1){
                 j1=P1[j1];
                 if (P1[j1]!=-1) R2T1[i]+=B1[j1];
           }
           R2T2[i]=B2[i];
           int j2=i;
           while (P2[j2]!=-1){
                 j2=P2[j2];
                 if (P2[j2]!=-1) R2T2[i]+=B2[j2];
           }
       }
       double result=0;
       for (int i=n;i<=m;i++){
           int j = index(Labels1[i],Labels2,m+1);
           result += myabs(R2T1[i]-R2T2[j])/R2T1[i];           
       }
       return result/(m-n+1);
}

void figtree(FILE *w){
    fprintf(w,"begin figtree;\n");
	fprintf(w,"set appearance.backgroundColorAttribute=\"User Selection\";\n");
	fprintf(w,"set appearance.backgroundColour=#-1;\n");
	fprintf(w,"set appearance.branchColorAttribute=\"User Selection\";\n");
	fprintf(w,"set appearance.branchLineWidth=1.0;\n");
	fprintf(w,"set appearance.foregroundColour=#-16777216;\n");
	fprintf(w,"set appearance.selectionColour=#-2144520576;\n");
	fprintf(w,"set branchLabels.colorAttribute=\"User Selection\";\n");
	fprintf(w,"set branchLabels.displayAttribute=\"cv\";\n");
	fprintf(w,"set branchLabels.displayAttribute=\"ld\";\n");
	fprintf(w,"set branchLabels.fontName=\"Agency FB\";\n");
	fprintf(w,"set branchLabels.fontSize=10;\n");
	fprintf(w,"set branchLabels.fontStyle=0;\n");
	fprintf(w,"set branchLabels.isShown=true;\n");
	fprintf(w,"set branchLabels.significantDigits=4;\n");
	fprintf(w,"set layout.expansion=0;\n");
	fprintf(w,"set layout.layoutType=\"RECTILINEAR\";\n");
	fprintf(w,"set layout.zoom=0;\n");
	fprintf(w,"set nodeBars.barWidth=4.0;\n");
	fprintf(w,"set nodeLabels.colorAttribute=\"User Selection\";\n");
	fprintf(w,"set nodeLabels.displayAttribute=\"Node ages\";\n");
	fprintf(w,"set nodeLabels.fontName=\"Agency FB\";\n");
	fprintf(w,"set nodeLabels.displayAttribute=\"node\";\n");    	
	fprintf(w,"set nodeLabels.fontSize=8;\n");
	fprintf(w,"set nodeLabels.fontStyle=0;\n");
	fprintf(w,"set nodeLabels.isShown=false;\n");
	fprintf(w,"set nodeLabels.significantDigits=4;\n");
	fprintf(w,"set polarLayout.alignTipLabels=false;\n");
	fprintf(w,"set polarLayout.angularRange=0;\n");
	fprintf(w,"set polarLayout.rootAngle=0;\n");
	fprintf(w,"set polarLayout.rootLength=100;\n");
	fprintf(w,"set polarLayout.showRoot=true;\n");
	fprintf(w,"set radialLayout.spread=0.0;\n");
	fprintf(w,"set rectilinearLayout.alignTipLabels=false;\n");
	fprintf(w,"set rectilinearLayout.curvature=0;\n");
	fprintf(w,"set rectilinearLayout.rootLength=100;\n");
	fprintf(w,"set scale.offsetAge=0.0;\n");
	fprintf(w,"set scale.rootAge=1.0;\n");
	fprintf(w,"set scale.scaleFactor=1.0;\n");
	fprintf(w,"set scale.scaleRoot=false;\n");
	fprintf(w,"set scaleAxis.automaticScale=true;\n");
	fprintf(w,"set scaleAxis.fontSize=8.0;\n");
	fprintf(w,"set scaleAxis.isShown=false;\n");
	fprintf(w,"set scaleAxis.lineWidth=1.0;\n");
	fprintf(w,"set scaleAxis.majorTicks=1.0;\n");
	fprintf(w,"set scaleAxis.origin=0.0;\n");
	fprintf(w,"set scaleAxis.reverseAxis=false;\n");
	fprintf(w,"set scaleAxis.showGrid=true;\n");
	fprintf(w,"set scaleAxis.significantDigits=4;\n");
	fprintf(w,"set scaleBar.automaticScale=true;\n");
	fprintf(w,"set scaleBar.fontSize=10.0;\n");
	fprintf(w,"set scaleBar.isShown=false;\n");
	fprintf(w,"set scaleBar.lineWidth=1.0;\n");
	fprintf(w,"set scaleBar.scaleRange=0.06;\n");
	fprintf(w,"set scaleBar.significantDigits=4;\n");
	fprintf(w,"set tipLabels.colorAttribute=\"User Selection\";\n");
	fprintf(w,"set tipLabels.displayAttribute=\"Names\";\n");
	fprintf(w,"set tipLabels.fontName=\"Agency FB\";\n");
	fprintf(w,"set tipLabels.fontSize=10;\n");
	fprintf(w,"set tipLabels.fontStyle=0;\n");
	fprintf(w,"set tipLabels.isShown=true;\n");
	fprintf(w,"set tipLabels.significantDigits=4;\n");
	fprintf(w,"set trees.order=false;\n");
	fprintf(w,"set trees.orderType=\"increasing\";\n");
	fprintf(w,"set trees.rooting=false;\n");
	fprintf(w,"set trees.rootingType=\"User Selection\";\n");
	fprintf(w,"set trees.transform=false;\n");
	fprintf(w,"set trees.transformType=\"cladogram\";\n");
fprintf(w,"end;\n");
}


double reroot_unrootedtree(int m,int r,int* & P,double* & B,double* & T,int* & P_new,double* & B_new,double* & T_new){
//r neq 0       
//m: number of branches in the unrooted tree
       double br=B[r];
       //if (br<0) br=0;
       bool* bl = new bool[m+1];
       for (int i=0;i<=m;i++) bl[i]=true;
       P_new[0]=-1;
       P_new[r+1]=0;
       bl[r]=false;
       P_new[P[r]+1]=0;
       bl[P[r]]=false;
       B_new[r+1]=-1;
       B_new[P[r]+1]=-1;
       int i=P[r];
       int j=P[i];
       while (j!=-1){
             P_new[j+1]=i+1;
             B_new[j+1]=B[i];
             //if (B[i]<0) B_new[j+1]=0;
             i=j;
             j=P[i];
             bl[i]=false;
       }
       for (int i=0;i<=m;i++){
           if (bl[i]){
              P_new[i+1]=P[i]+1;
              B_new[i+1]=B[i];
              //if (B[i]<0) B_new[i+1]=0;
           }
           T_new[i+1]=T[i];
       }
       delete[] bl;
       return br;
}



void reroot_unrootedtree(int m,int r,double lambda,int* & P,double* & B,double* & T,int* & P_new,double* & B_new,double* & T_new){
//r neq 0       
//m: number of branches in the unrooted tree
       bool* bl = new bool[m+1];
       for (int i=0;i<=m;i++) bl[i]=true;
       P_new[0]=-1;
       P_new[r+1]=0;
       bl[r]=false;
       P_new[P[r]+1]=0;
       bl[P[r]]=false;
       B_new[r+1]=B[r]*lambda;
       B_new[P[r]+1]=(1-lambda)*B[r];
       int i=P[r];
       int j=P[i];
       while (j!=-1){
             P_new[j+1]=i+1;
             B_new[j+1]=B[i];
             //if (B[i]<0) B_new[j+1]=0;
             i=j;
             j=P[i];
             bl[i]=false;
       }
       for (int i=0;i<=m;i++){
           if (bl[i]){
              P_new[i+1]=P[i]+1;
              B_new[i+1]=B[i];
              //if (B[i]<0) B_new[i+1]=0;
           }
           T_new[i+1]=T[i];
       }
       delete[] bl;
}

void reroot_unrootedtree(int m,int r,double lambda,int* & P,double* & B,string* & Support,double* & T,int* & P_new,double* & B_new,string* & Support_new,double* & T_new){
    //r neq 0
    //m: number of branches in the unrooted tree
    bool* bl = new bool[m+1];
    for (int i=0;i<=m;i++) bl[i]=true;
    P_new[0]=-1;
    P_new[r+1]=0;
    bl[r]=false;
    P_new[P[r]+1]=0;
    bl[P[r]]=false;
    B_new[r+1]=B[r]*lambda;
    B_new[P[r]+1]=(1-lambda)*B[r];
    Support_new[r+1]="";
    Support_new[P[r]+1]="";
    int i=P[r];
    int j=P[i];
    while (j!=-1){
        P_new[j+1]=i+1;
        B_new[j+1]=B[i];
        Support_new[j+1]=Support[i];
        //if (B[i]<0) B_new[j+1]=0;
        i=j;
        j=P[i];
        bl[i]=false;
    }
    for (int i=0;i<=m;i++){
        if (bl[i]){
            P_new[i+1]=P[i]+1;
            B_new[i+1]=B[i];
            Support_new[i+1]=Support[i];
            //if (B[i]<0) B_new[i+1]=0;
        }
        T_new[i+1]=T[i];
    }
    delete[] bl;
}

double reroot_rootedtree(int m,int r,int* & P,int s10,int s20,double* & B,int* & P_new,double* & B_new){
       double br;
       if (r==s10 || r==s20){
          for (int i=0;i<=m;i++){
              P_new[i]=P[i];
              B_new[i]=B[i];
          } 
          B_new[s10]=-1;
          B_new[s20]=-1;
          br = B[s10]+B[s20];
       }
       else {
            bool* bl = new bool[m+1];
            for (int i=0;i<=m;i++) bl[i]=true;
            P_new[0]=-1;
            bl[0]=false;
            P_new[r]=0;
            bl[r]=false;
            P_new[P[r]]=0;
            bl[P[r]]=false;
            B_new[r]=-1;
            B_new[P[r]]=-1;
            int i=P[r];
            int j=P[i];     
            while (j!=0){
                  P_new[j]=i;
                  B_new[j]=B[i];
                  i=j;
                  j=P[i];
                  bl[i]=false;
            }
            int k=s10;
            if (k==i) k=s20;
            bl[k]=false;
            P_new[k]=i;
            br=B[r];
            B_new[k]=B[i]+B[k];
            for (int i=0;i<=m;i++){
                if (bl[i]){
                   P_new[i]=P[i];
                   B_new[i]=B[i];
                }
            }
            delete[] bl;
       }
       return br;
}

double reroot_rootedtree(int m,int r,int* & P,int s10,int s20,double* & B,int* & P_new,double* & B_new,int* & P_ref,int* & tab){
    double br;
    for (int i=0; i<=m; i++) {
        tab[i]=i;
    }
    if (r==s10 || r==s20){
        for (int i=0;i<=m;i++){
            P_new[i]=P[i];
            P_ref[i]=P[i];
            B_new[i]=B[i];
        }
        B_new[s10]=-1;
        B_new[s20]=-1;
        br = B[s10]+B[s20];
    }
    else {
        bool* bl = new bool[m+1];
        for (int i=0;i<=m;i++) bl[i]=true;
        P_new[0]=-1;
        P_ref[0]=-1;
        bl[0]=false;
        P_new[r]=0;
        P_ref[r]=0;
        bl[r]=false;
        P_new[P[r]]=0;
        P_ref[P[r]]=0;
        bl[P[r]]=false;
        B_new[r]=-1;
        B_new[P[r]]=-1;
        int i=P[r];
        int j=P[i];
        tab[i]=r;
        while (j!=0){
            tab[j]=i;
            P_new[j]=i;
            P_ref[j]=i;
            B_new[j]=B[i];
            i=j;
            j=P[i];
            bl[i]=false;
        }
        int k=s10;
        if (k==i) k=s20;
        bl[k]=false;
        P_new[k]=i;
        P_ref[k]=i;
        br=B[r];
        B_new[k]=B[i]+B[k];
        for (int i=0;i<=m;i++){
            if (bl[i]){
                P_new[i]=P[i];
                P_ref[i]=P[i];
                B_new[i]=B[i];
            }
        }
        delete[] bl;
    }
    return br;
}

void reroot_rootedtree(int m,int r,double lambda,int* & P,int s10,int s20,double* & B,int* & P_new,double* & B_new){
       double br=B[r];
       if (r==s10 || r==s20){
          for (int i=0;i<=m;i++){
              P_new[i]=P[i];
              B_new[i]=B[i];
          } 
          br=B[s10]+B[s20];
          B_new[s10]=(1-lambda)*br;                          
          B_new[s20]=lambda*br;
       }
       else {
            bool* bl = new bool[m+1];
            for (int i=0;i<=m;i++) bl[i]=true;
            P_new[0]=-1;
            bl[0]=false;
            P_new[r]=0;
            bl[r]=false;
            P_new[P[r]]=0;
            bl[P[r]]=false;
            int i=P[r];
            int j=P[i];     
            while (j!=0){
                  P_new[j]=i;
                  B_new[j]=B[i];
                  i=j;
                  j=P[i];
                  bl[i]=false;
            }
            int k=s10;
            if (k==i) k=s20;
            bl[k]=false;
            P_new[k]=i;
            B_new[k]=B[i]+B[k];
            for (int i=0;i<=m;i++){
                if (bl[i]){
                   P_new[i]=P[i];
                   B_new[i]=B[i];
                }
            }
            B_new[r]=lambda*br;
            B_new[P[r]]=(1-lambda)*br;
           if (r<P[r]) {
               double temp=B_new[r];
               B_new[r]=B_new[P[r]];
               B_new[P[r]]=temp;
           }
            delete[] bl;
       }    
}

void reroot_rootedtree(int m,int r,double lambda,int* & P,int s10,int s20,double* & B,string* & Support,int* & P_new,double* & B_new,string* & Support_new){
    double br=B[r];
    if (r==s10 || r==s20){
        for (int i=0;i<=m;i++){
            P_new[i]=P[i];
            B_new[i]=B[i];
            Support_new[i]=Support[i];
        }
        br=B[s10]+B[s20];
        B_new[s10]=(1-lambda)*br;
        B_new[s20]=lambda*br;
        Support_new[s10]="";
        Support_new[s20]="";
    }
    else {
        bool* bl = new bool[m+1];
        for (int i=0;i<=m;i++) bl[i]=true;
        P_new[0]=-1;
        bl[0]=false;
        P_new[r]=0;
        bl[r]=false;
        P_new[P[r]]=0;
        bl[P[r]]=false;
        int i=P[r];
        int j=P[i];
        while (j!=0){
            P_new[j]=i;
            B_new[j]=B[i];
            Support_new[j]=Support[i];
            i=j;
            j=P[i];
            bl[i]=false;
        }
        int k=s10;
        if (k==i) k=s20;
        bl[k]=false;
        P_new[k]=i;
        B_new[k]=B[i]+B[k];
        Support_new[k]="";
        for (int i=0;i<=m;i++){
            if (bl[i]){
                P_new[i]=P[i];
                B_new[i]=B[i];
                Support_new[i]=Support[i];
            }
        }
        B_new[r]=lambda*br;
        B_new[P[r]]=(1-lambda)*br;
        if (r<P[r]) {
            double temp=B_new[r];
            B_new[r]=B_new[P[r]];
            B_new[P[r]]=temp;
        }
        Support_new[s10]="";
        Support_new[s20]="";
        delete[] bl;
    }
}

void rooted2unrooted(int n,int m,int* & P,int* & Suc1,int* & Suc2,double* & B,double* & T,int* & Pur,double* & Bur,double* & Tur){
     int s1=Suc1[0];
     int s2=Suc2[0];
     if (s2<n){             
        Pur[s1-1]=s2-1;
        Pur[s2-1]=-1;
        Bur[s1-1]=B[s1]+B[s2];        
     }
     else{
          Pur[s2-1]=s1-1;
          Pur[s1-1]=-1;
          Bur[s2-1]=B[s1]+B[s2];
     }
     for (int i=0;i<m;i++){
         if (i!=(s1-1) && i!=(s2-1)){
            Pur[i]=P[i+1]-1;          
            Bur[i]=B[i+1];
         }        
         Tur[i]=T[i+1];
     }
}

void rooted2unrooted(int n,int m,int* & P,int* & Suc1,int* & Suc2,double* & B,string* & Labels){     
     int s1=Suc1[0];
     int s2=Suc2[0];
     int* Pur = new int[m];
     string* Labelsur = new string[m];
     double* Bur=new double[m];
     if (s2<n){             
        Pur[s1-1]=s2-1;
        Pur[s2-1]=-1;
        Bur[s1-1]=B[s1]+B[s2];        
     }
     else{
          Pur[s2-1]=s1-1;
          Pur[s1-1]=-1;
          Bur[s2-1]=B[s1]+B[s2];
     }
     for (int i=0;i<m;i++){
         if (i!=(s1-1) && i!=(s2-1)){
            Pur[i]=P[i+1]-1;          
            Bur[i]=B[i+1];
         }    
         Labelsur[i]=Labels[i+1];    
     }
     for (int i=0;i<m;i++){
         P[i]=Pur[i];
         B[i]=Bur[i];
         Labels[i]=Labelsur[i];
     }
     delete[] Labelsur;
     delete[] Bur;
     delete[] Pur;
}

void reroot(int n,int m,int r,int* & P,int &s,int* & Suc1,int* & Suc2,double* & B,string* & Support,int* & P_new,int* & Suc1_new,int* & Suc2_new,double* & B_new,string* & Support_new){
//reroot an unrooted tree to another unrooted tree
     for (int i=0;i<=m;i++){
         P_new[i]=P[i];
         B_new[i]=B[i];
         Support_new[i]=Support[i];
     }
     for (int i=0;i<n;i++){
         Suc1_new[i]=Suc1[i];
         Suc2_new[i]=Suc2[i];
     }
    P_new[r]=-1;
    int i=P[r];
     B_new[i]=B[r];
     Support_new[i]=Support[r];
     P_new[i]=r;
     int j=P[i];
    if (j==-1){
          if (Suc1[i]==r) Suc1_new[i]=s;
          else if (Suc2[i]==r) Suc2_new[i]=s;
          else s=i;
     }
     else{
            while (j!=-1){   
                  P_new[j]=i;
                  int s1=Suc1[i];
                  if (i==P[r]){
                     if (s1==r) Suc1_new[i]=j;
                     else Suc2_new[i]=j;
                  }
                  else{
                       if (s1==P_new[i]) Suc1_new[i]=j;
                       else Suc2_new[i]=j;
                  }
                  B_new[j]=B[i];
                  Support_new[j]=Support[i];
                  i=j;    
                  j=P[j];       
            }                
            if (Suc1[i]==P_new[i]) Suc1_new[i]=s;
            else if (Suc2[i]==P_new[i]) Suc2_new[i]=s;                        
       }
}

void unrooted2rooted(int n,int m,int r,double lambda,int* & P1,int* & Suc11,int* & Suc21,double* & B1,double* & T,int* & P_new,int* & Suc1_new,int* & Suc2_new,double* & B_new,double* & T_new){
//input: an unrooted tree, a branch r, and lambda
//output: the tree rooted at the position lambda of branch r         
     double b=B1[r];
     B1[r]=b*lambda; B1[P1[r]]=b*(1-lambda);
     //if (B1[r]<0) B1[r]=0;
     //if (B1[P1[r]]<0) B1[P1[r]]=0;
     for (int i=0;i<=m;i++){
         if (i==r || i==P1[r]) P_new[i+1]=0;
         else P_new[i+1]=P1[i]+1;
         B_new[i+1]=B1[i];
         T_new[i+1]=T[i];
     }
     for (int i=0;i<n;i++){
         Suc1_new[i+1]=Suc11[i]+1;
         Suc2_new[i+1]=Suc21[i]+1;
     }
     P_new[0]=-1;
     Suc1_new[0]=r+1;
     Suc2_new[0]=P1[r]+1;     
}

double unrooted2rooted(int n,int m,int r,int* & P1,int* & Suc11,int* & Suc21,double* & B1,double* & T,int* & P_new,int* & Suc1_new,int* & Suc2_new,double* & B_new,double* & T_new){
     double br = B1[r];
     B1[r]=-1;
     B1[P1[r]]=-1;
     for (int i=0;i<=m;i++){
         if (i==r || i==P1[r]) P_new[i+1]=0;
         else P_new[i+1]=P1[i]+1;
         B_new[i+1]=B1[i];
         T_new[i+1]=T[i];
     }
     for (int i=0;i<n;i++){
         Suc1_new[i+1]=Suc11[i]+1;
         Suc2_new[i+1]=Suc21[i]+1;
     }
     P_new[0]=-1;
     Suc1_new[0]=r+1;
     Suc2_new[0]=P1[r]+1;     
     return br;
}

double phi(int n,int m,int* & P,double* & B,double* & V,double* & T,double rho){
       double p = 0;
       for (int i=1;i<=m;i++){ 
           //p+=(B[i]-rho*T[i]+rho*T[P[i]])*(B[i]-rho*T[i]+rho*T[P[i]])/(2*V[i]) +log(2*M_PI*V[i])/2;
          p+=(B[i]-rho*T[i]+rho*T[P[i]])*(B[i]-rho*T[i]+rho*T[P[i]])/(V[i]); 
       }
       return p;//p = -(log likelihood)
}

double phi(int n,int m,double &rho,int* & P,double* & B,double* & V,double* & T){
       double a=0;
       double b=0;
       double c=0;       
       for (int i=1;i<=m;i++){ 
          a+=(T[P[i]]-T[i])*(T[P[i]]-T[i])/V[i];
          b+=2*(T[P[i]]-T[i])*B[i]/V[i];
          c+=B[i]*B[i]/V[i];
       }
       rho=-b/a/2;       
       return a*rho*rho+b*rho+c;
}
void output(bool constr,bool v,FILE* f,FILE* tree1,FILE* tree2,FILE* tree3,int m,int y,double rho,int* P,int* Suc1,int* Suc2,double* B,double* T,string* Labels,string* Support,double phi){
     /*if (v){
        fprintf(f,"rate %.6f, tMRCA %.3f, log_likelihood %.6f\n",rho,T[0],-phi);
        printf("rate %.6f, tMRCA %.3f, log_likelihood %.6f\n",rho,T[0],-phi);
     }
     else{
        fprintf(f,"rate %.6f, tMRCA %.3f\n\n",rho,T[0]);
        printf("rate %.6f, tMRCA %.3f\n\n",rho,T[0]);
     }*/
     fprintf(f,"rate %.6f, tMRCA %.3f, objective function %.6f\n",rho,T[0],phi);
     //printf("rate %.6f, tMRCA %.3f, objective function %.6f\n",rho,T[0],phi);
     fprintf(tree1,"tree %d = ",y);
     for (int i=1;i<=m;i++) B[i]=rho*(T[i]-T[P[i]]);
     nexustree(m/2,P,Suc1,Suc2,Labels,B,Support,T,tree1);
     newicktree(m/2,P,Suc1,Suc2,Labels,B,Support,tree2);
     for (int i=1;i<=m;i++) B[i]=B[i]/rho;
     newicktree(m/2,P,Suc1,Suc2,Labels,B,Support,tree3);   
     if (!constr){
        int count=0;
        for (int i=1;i<=m;i++){
            if (T[i]<T[P[i]]) count++;
        }
        fprintf(f,"Number of violated temporal constraints (nodes having date smaller than the one of its parent): %d (%.2f%%)\n",count,(count*100)/(double)m);
        printf("Number of violated temporal constraint (node having date smaller than the one of its parent): %d (%.2f%%)\n",count,(count*100)/(double)m);
     }
     fprintf(f,"\n");
     printf("\n");
}

double* variance(bool w,int m,double* B,int c,int s){
        double* V = new double[m+1];
        if (w){        
           for (int i=1;i<=m;i++){ 
               //V[i]=(B[i]+c/(c+s))/(c+s);
               V[i]=(B[i]+(double)c/s)/s;
           }
        }
        else{
             for (int i=1;i<=m;i++){ 
                 V[i]=1./(double)(s);
             }
        }
        return V;
}

list<string> getOutgroup(string fn){
             list<string> result;
             FILE * o = fopen(fn.c_str(),"rt");
             if (o==NULL) cout<<"Impossible to open outgroup file"<<endl;
             int ino=readInt(o,"Error in the outgroup file, the file should begin with an integer (the number of outgroups)");
             for (int i=0;i<ino;i++) result.push_back(readWord(o,fn));
             return result;
}
