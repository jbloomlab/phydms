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
#include "readData.h"

bool tree2data(FILE * tree,string date,int n,int m,int* &P,double* &B,string* &Support,double* &T,string* &leaves){
      int inode = 0;//number of internal nodes;
      stack<int> pileNode;    
      char c = readBracket(tree,"input tree");
      int a=1;
      int countleaf=n;
      P[0]=-1;
      B[0]=-1;
      bool all=false;
      do{
             c = readChar(tree,"input tree");
                   if (c==')'){    
                        a--;inode++;
                        int s1=pileNode.top();pileNode.pop();
                        int s2=pileNode.top();pileNode.pop();                   
                        P[s1]=n-inode;
                        P[s2]=n-inode;
                        if (a>0){
                           Support[n-inode]=readSupport(tree,"input tree");  
                           B[n-inode]=readdouble(tree,"input tree");
                        }
                        else
                            while (!pileNode.empty()){
                             int s=pileNode.top();pileNode.pop();
                             P[s]=n-inode;                                
                        }  
                        pileNode.push(n-inode);                  
                   }
                   else if (c!='(' && c!=',' && c!=-1 && c!=';' && c!='\n'){
                        string lb=readLabel(c,tree);  
                        pileNode.push(countleaf);                        
                        leaves[countleaf]=lb;
                        B[countleaf]=readdouble(tree,"input tree");
                        Support[countleaf]="";
                        countleaf++;
                   }                    
                   else if (c=='(') {a++;}
           }  while (a>0);
           FILE * dateFile = fopen(date.c_str(),"rt");
           if (dateFile == NULL)  cout << "Can not open the file " <<date<< endl;
           else{
               int ino=readInt(dateFile,"Error in the date file, the file should begin with an integer (the number of dated tips)");
               for (int i=0;i<=m;i++) T[i]=-1;
              int count=0;
              double first=-1;
              for (int i=0;i<ino;i++){
                     string s=readWord(dateFile,date);
                     double lt=readdouble(dateFile,date);                    
                     int k = getPosition(leaves,s,n,m+1);
                     if (k!=-1) {
                         if (count==0) first=lt;
                         else if (lt!=first) all=true;                                
                         T[k]=lt;
                         count++;
                     }
              }       
              if (count<countleaf-n){
                 cout<<"The dates of some taxa are missing. If the input tree contains some outgroups, use option -g to specify the file name of outgroups. "<<endl;
                 exit(EXIT_FAILURE);
              }  
              fclose(dateFile);              
           }
      return all;
}

void tree2dataS(FILE * tree,int n,int m,int* &P,double* &B,string* &Support,string* &leaves){
      int inode = 0;//number of internal nodes;
      stack<int> pileNode;    
      char c = readBracket(tree,"input tree");
      int a=1;
      int countleaf=n;
      P[0]=-1;
      B[0]=-1;
      do{
             c = readChar(tree,"input tree");
                   if (c==')'){    
                        a--;inode++;
                        int s1=pileNode.top();pileNode.pop();
                        int s2=pileNode.top();pileNode.pop();                   
                        P[s1]=n-inode;
                        P[s2]=n-inode;
                        if (a>0){                      
                           Support[n-inode]=readSupport(tree,"input tree");  
                           B[n-inode]=readdouble(tree,"input tree");                                       
                        }
                        else while (!pileNode.empty()){
                             int s=pileNode.top();pileNode.pop();
                             P[s]=n-inode;                                
                        }  
                        pileNode.push(n-inode);                  
                   }
                   else if (c!='(' && c!=',' && c!=-1 && c!=';' && c!='\n'){
                       string lb=readLabel(c,tree);
                        pileNode.push(countleaf);                        
                        leaves[countleaf]=lb;
                        B[countleaf]=readdouble(tree,"input tree");
                        Support[countleaf]="";
                        countleaf++;
                   }                    
                   else if (c=='(') {a++;}
           }  while (a>0);
}

void extrait_outgroup(string inFile,string outFile,list<string> &outgroups,int nb){
     int n;//the number of internal nodes
	 int m;//the number of branches
     int n1,m1;
    counting(inFile,n,m,nb);
    int *P = new int[m+1];//tableau de predecesseurs
    double *B = new double[m+1];//branch lengths
    string *Support = new string[m+1];
    string* Labels = new string[m+1];
	FILE * tree = fopen(inFile.c_str(),"rt");
	if (tree==NULL) cout<<"Can not open the tree file"<<endl;
	else{
         FILE * w = fopen(outFile.c_str(),"wt");
         bool rooted=false;
         for (int y=1;y<=nb;y++){
             if (rooted){
                n++;
                m++;
             }
             tree2dataS(tree,n,m,P,B,Support,Labels);
             int * Suc1= new int[n];
             int * Suc2= new int[n];
             computeSuc(P,Suc1,Suc2,m+1,n);
             list<int> P1;
             list<double> B1;
             list<string> Support1;
             list<string> Labels1;
             list<double> T1;
             int s;
             if (m==2*n){
                rooted=true;
                rooted2unrooted(n,m,P,Suc1,Suc2,B,Labels);
                m--;
                n--;
                computeSuc(P,Suc1,Suc2,m+1,n);
             }
             s= computeSuc_unrooted(P,Suc1,Suc2,m+1,n);
             list<int> out;
             list<int> in;
             for (int i=n;i<=m;i++){
                 bool flag = true;
                 for (list<string>::iterator iter=outgroups.begin();iter!=outgroups.end();iter++){
		             if (Labels[i].compare(*iter)==0){
                        out.push_back(i);                   
                        flag=false;
                        break;
                     }
                 }
                 if (flag) in.push_back(i);
             }
             if (out.size()==0){
                cout<<"Tree "<<y<<": the tree does not contain any outgroup"<<endl;  
                exit( EXIT_SUCCESS );
             }
             else{
                 int *P_new= new int[m+1];
                 int *Suc1_new= new int[n];
                 int *Suc2_new= new int[n];
                 double *B_new= new double[m+1];
                 string* Support_new = new string[m+1];
                 int t=out.front();
                 int t1=in.front();
                  if (out.size()==1){
                     reroot(n,m,t,P,s,Suc1,Suc2,B,Support,P_new,Suc1_new,Suc2_new,B_new,Support_new);
                     subTree(n,m,P[t],P_new,Suc1_new,Suc2_new,B_new,Support_new,Labels,P1,B1,Support1,Labels1);
                     int *Po=new int[m];
                     listToArray(P1,Po);
                     double *Bo=new double[m];
                     listToArray(B1,Bo);
                     string* Supporto=new string[m];
                     listToArray(Support1,Supporto);
                     string* Labelso = new string[m];
                     listToArray(Labels1,Labelso);
                     int *Suc1o=new int[n];
                     int *Suc2o=new int[n];
                     computeSuc(Po,Suc1o,Suc2o,m,n);
                     newicktree(n,Po,Suc1o,Suc2o,Labelso,Bo,Supporto,w);
                      delete[]  Po;
                      delete[]  Bo;
                      delete[] Labelso;
                      delete[] Suc1o;
                      delete[] Suc2o;
                      delete[] Supporto;
                  }
                  else{
	                   bool flag = false;
	                   while (!flag && P[t]!=-1){
                             t = P[t];
                              flag=true;
                              for (list<int>::iterator ia=out.begin();ia!=out.end();ia++){
                                  int j=*ia; 
                                  if (!isAncestor(P,t,j)) {
                                     flag=false;break;
                                  }
                              }
                       }//t is lca of outgroups
                       if (P[t]!=-1){//lca of outgroups is not the root
                          reroot(n,m,t,P,s,Suc1,Suc2,B,Support,P_new,Suc1_new,Suc2_new,B_new,Support_new);
                          subTree(n,m,P[t],P_new,Suc1_new,Suc2_new,B_new,Support_new,Labels,P1,B1,Support1,Labels1);
                       }
                       else{//lca of outgroups is the root
	                        bool flag = false;
	                        while (!flag && P[t1]!=-1){
                                  t1 = P[t1];
                                  flag=true;
                                  for (list<int>::iterator ia=in.begin();ia!=in.end();ia++){
                                      int j=*ia;
                                      if (!isAncestor(P,t1,j)) {flag=false;break;}
                                  }
                            }
                           //t1 is lca of ingroups
                            if (P[t1]==-1){//lca of ingroups is the root
                               cout<<"Tree "<<y<<": The outgroups are not separated from the ingroups"<<endl; 
                               exit( EXIT_FAILURE );
                            }
                            else{//lca of ingroups is not the root
                                reroot(n,m,t1,P,s,Suc1,Suc2,B,Support,P_new,Suc1_new,Suc2_new,B_new,Support_new);
                                subTree(n,m,t1,P_new,Suc1_new,Suc2_new,B_new,Support_new,Labels,P1,B1,Support1,Labels1);
                          }
                      }
                      m1=(int)P1.size()-1;
                      n1=m1/2;
                      if ((n1+out.size())==n+1){
                              int *Po = new int[m1+1];
                              listToArray(P1,Po);
                              double *Bo = new double[m1+1];
                              listToArray(B1,Bo);
                              string* Supporto = new string[m1+1];
                              listToArray(Support1,Supporto);
                              string* Labelso = new string[m1+1];
                              listToArray(Labels1,Labelso);
                              int *Suc1o = new int[n1];
                              int *Suc2o = new int[n1];
                              computeSuc(Po,Suc1o,Suc2o,m1+1,n1);
                              newicktree(n1,Po,Suc1o,Suc2o,Labelso,Bo,Supporto,w);
                              delete[]  Po;
                              delete[]  Bo;
                              delete[]  Supporto;
                              delete[] Labelso;
                              delete[] Suc1o;
                              delete[] Suc2o;
                      }
                      else {
                              cout<<"Tree "<<y<<": The outgroups are not separated from the ingroups"<<endl;
                              exit( EXIT_FAILURE );
                      }
                      delete[] P_new;
                      delete[] Suc1_new;
                      delete[] Suc2_new;
                      delete[] B_new;
                      delete[] Support_new;
                  }
             delete[] Suc1;
             delete[] Suc2;
             }
             
         }
         fclose(tree);
        fclose(w);
    }
    delete[]  P;
    delete[]  B;
    delete[]  Support;
    delete[] Labels;
}
