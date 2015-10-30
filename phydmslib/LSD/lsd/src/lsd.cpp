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
    
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <stack>
#include <list>
#include <math.h>
#include <ctime>
#include "options.h"
#include "readData.h"
#include "dating.h"
#include "dating_relative.h"
#include "estimate_root_relative_rooted.h"
#include "estimate_root_relative_unrooted.h"
#include "estimate_root_rooted.h"
#include "estimate_root_unrooted.h"
#include "utils.h"

using namespace std;

int main( int argc, char** argv ){
    bool rooted=true;
	options* opt = getOptions( argc, argv, rooted);
    FILE * result = fopen(opt->outFile,"wt");
    if (result==NULL){
        cout<<"Error: can not create the output file "<<opt->outFile<<endl;
        exit(EXIT_FAILURE);
    }
    printInterface( result, opt, rooted);
    fprintf(result,"\n");
	clock_t start = clock();
    double elapsed_time;
    if (opt->fnOutgroup!=NULL){
        printf("Removing outgroups ...\n\n");
        list<string> outgroup = getOutgroup(opt->fnOutgroup);
        opt->inFileOp=getDefaultInFileOpName((string)(opt->inFile));
        extrait_outgroup(opt->inFile,opt->inFileOp,outgroup,opt->nbData);
        opt->inFile=opt->inFileOp;
    }
    int n;//the number of internal nodes
	int m;//the number of branches
    counting(opt->inFile,n,m,opt->nbData);
    double rho;//the estimating rate
    double phi;//the value of the objective function
	FILE * tree = fopen(opt->inFile,"rt");
    if (tree==NULL){
        cout<<"Error: can not open the input file"<<endl;
        exit(EXIT_FAILURE);
    }
    opt->treeFile1=getDefaultOutputNexusTreeFileName((string)(opt->outFile));
    opt->treeFile2=getDefaultOutputNewick1TreeFileName((string)(opt->outFile));
    opt->treeFile3=getDefaultOutputNewick2TreeFileName((string)(opt->outFile));
    FILE * tree1 = fopen(opt->treeFile1,"wt");
    FILE * tree2 = fopen(opt->treeFile2,"wt");
    FILE * tree3 = fopen(opt->treeFile3,"wt");
    FILE * gr;//given rate
    if (tree1==NULL){
        cout<<"Error: can not create the output tree file "<<opt->treeFile1<<endl;
        exit(EXIT_FAILURE);
    }
    if (tree2==NULL){
        cout<<"Error: can not create the output tree file "<<opt->treeFile2<<endl;
        exit(EXIT_FAILURE);
    }
    if (tree2==NULL){
        cout<<"Error: can not create the output tree file "<<opt->treeFile3<<endl;
        exit(EXIT_FAILURE);
    }
    fprintf(tree1,"#NEXUS\nBegin trees;\n");
    if (opt->rate!=NULL)
        gr = fopen(opt->rate,"rt");
    int* P = new int[m+1];
    int* Suc1 = new int[n];
    int* Suc2 = new int[n];
    double* B = new double[m+1];
    double* V = new double[m+1];
    double* T = new double[m+1];
    string* Support= new string[m+1];
    string* Labels = new string[m+1];
    for (int y=1;y<=opt->nbData;y++){
        fprintf(result,"Tree %d ",y);
        cout<<"Reading the tree "<<y<< "..."<<endl;
        if (!opt->relative) opt->relative=!tree2data(tree,opt->inDateFile,n,m,P,B,Support,T,Labels);
        else tree2dataS(tree,n,m,P,B,Support,Labels);
        if (opt->relative){
           fprintf(result,"The results correspond to the estimation of relative dates when T[mrca]=%0.3f and T[tips]=%0.3f\n",opt->mrca,opt->leaves);
           printf("Estimating relative dates with T[mrca]=%0.3f and T[tips]=%0.3f\n",opt->mrca,opt->leaves);
           for (int i=n;i<=m;i++) T[i]=opt->leaves;
        }
        //////////////////////////////////rooted tree///////////////////////////////////////////////////////
        if (m==2*n){//rooted tree
            computeSuc(P,Suc1,Suc2,m+1,n);            
            if (!opt->constraint){//LD without constraints
                if (opt->estimate_root==NULL){//keep the given root
                   V=variance(opt->variance,m,B,opt->c,opt->seqLength);                   
                    if (opt->rate==NULL){//the rate is not given
                        //cout<<"Starting LD ..."<<endl;
                        if (!opt->relative) phi = without_constraint(n,m,P,Suc1,Suc2,B,V,T,rho,opt->delta);
                        else phi = without_constraint_relative(n,m,P,Suc1,Suc2,B,V,T,rho,opt->delta,opt->mrca);
                    }
                    else{//the rate is given
                        rho=readdouble(gr,(string)(opt->rate));
                        //cout<<"Starting LD ..."<<endl;
                        if (!opt->relative) phi = without_constraint_rate(rho,n,m,P,Suc1,Suc2,B,V,T);
                        else phi = without_constraint_rate_relative(rho,n,m,P,Suc1,Suc2,B,V,T,opt->mrca);
                    }                
                    output(opt->constraint,opt->variance,result,tree1,tree2,tree3,m,y,rho,P,Suc1,Suc2,B,T,Labels,Support,phi);
                }
                else{//estimate the root
                    int* P_new = new int[m+1];
                    double* B_new = new double[m+1];
                    string* Support_new = new string[m+1];                
                    int* Suc1_new = new int[n];
                    int* Suc2_new = new int[n];
                    V=variance(false,m,B_new,opt->c,opt->seqLength);
                    int r;
                    double lambda;
                    if (opt->rate==NULL){//the rate is not given                        
                        if (strcmp(opt->estimate_root,"l")==0){//improve locally the root around the given root                           
                            if (!opt->relative) phi=estimate_root_without_constraint_local_rooted(n,m,Suc1,Suc2,P,B,V,T,r,lambda,rho,opt->delta);
                            else phi=estimate_root_without_constraint_local_rooted_relative(n,m,Suc1,Suc2,P,B,V,T,r,lambda,rho,opt->delta,opt->mrca);
                        }
                        else{ //forget the given root and re-estimate the position of the root over all branhces
                            if (!opt->relative) phi=estimate_root_without_constraint_rooted(n,m,Suc1[0],Suc2[0],P,B,V,T,r,lambda,rho,opt->delta);
                            else phi=estimate_root_without_constraint_rooted_relative(n,m,Suc1[0],Suc2[0],P,B,V,T,r,lambda,rho,opt->delta,opt->mrca);
                        }
                        reroot_rootedtree(m,r,lambda,P,Suc1[0],Suc2[0],B,Support,P_new,B_new,Support_new);                        
                        V=variance(opt->variance,m,B_new,opt->c,opt->seqLength);                        
                        cout<<"Starting LD on the new rooted tree ..."<<endl;                        
                        computeSuc(P_new,Suc1_new,Suc2_new,m+1,n);
                        if (!opt->relative) phi = without_constraint(n,m,P_new,Suc1_new,Suc2_new,B_new,V,T,rho,opt->delta);    
                        else phi = without_constraint_relative(n,m,P_new,Suc1_new,Suc2_new,B_new,V,T,rho,opt->delta,opt->mrca);
                        output(opt->constraint,opt->variance,result,tree1,tree2,tree3,m,y,rho,P_new,Suc1_new,Suc2_new,B_new,T,Labels,Support_new,phi);
                    }
                    else{ //the rate is given
                        rho=readdouble(gr,(string)opt->rate);
                        if (strcmp(opt->estimate_root,"l")==0){//improve locally the root around the given root                          
                            if (!opt->relative) phi=estimate_root_without_constraint_local_rate_rooted(rho,n,m,Suc1,Suc2,P,V,B,T,r,lambda);
                            else phi=estimate_root_without_constraint_local_rate_rooted_relative(rho,n,m,Suc1,Suc2,P,B,V,T,r,lambda,opt->mrca);
                        }
                        else{ //forget the given root and re-estimate the position of the root over all branches                            
                            if (!opt->relative) phi=estimate_root_without_constraint_rate_rooted(rho,n,m,Suc1[0],Suc2[0],P,B,V,T,r,lambda);
                            else phi=estimate_root_without_constraint_rate_rooted_relative(rho,n,m,Suc1[0],Suc2[0],P,B,V,T,r,lambda,opt->mrca);
                        }                        
                        reroot_rootedtree(m,r,lambda,P,Suc1[0],Suc2[0],B,Support,P_new,B_new,Support_new);                        
                        cout<<"Starting LD on the new rooted tree ..."<<endl;
                        V=variance(opt->variance,m,B_new,opt->c,opt->seqLength);                                                
                        computeSuc(P_new,Suc1_new,Suc2_new,m+1,n);
                        if (!opt->relative) phi = without_constraint_rate(rho,n,m,P_new,Suc1_new,Suc2_new,B_new,V,T); 
                        else phi = without_constraint_rate_relative(rho,n,m,P_new,Suc1_new,Suc2_new,B_new,V,T,opt->mrca);      
                        output(opt->constraint,opt->variance,result,tree1,tree2,tree3,m,y,rho,P_new,Suc1_new,Suc2_new,B_new,T,Labels,Support_new,phi);
                    }
                    delete[] P_new;
                    delete[] Suc1_new;
                    delete[] Suc2_new;
                    delete[] B_new;
                    delete[] Support_new; 
                }                               
            }
            else {//QPD with temporal constrains
                if (opt->estimate_root==NULL){//keep the root
                    V=variance(opt->variance,m,B,opt->c,opt->seqLength);
                    if (opt->rate==NULL){//the rate is not given                        
                        if (!opt->relative) phi = with_constraint_active_set(n,m,P,Suc1,Suc2,B,V,T,rho,opt->delta);
                        else phi = with_constraint_active_set_relative(n,m,P,Suc1,Suc2,B,V,T,rho,opt->delta,opt->mrca);
                    }
                    else{//the rate is given
                        rho=readdouble(gr,(string)opt->rate);
                        if (!opt->relative) phi = with_constraint_active_set_rate(rho,n,m,P,Suc1,Suc2,B,V,T);
                        else phi = with_constraint_active_set_rate_relative(rho,n,m,P,Suc1,Suc2,B,V,T,opt->mrca);
                    }
                    output(opt->constraint,opt->variance,result,tree1,tree2,tree3,m,y,rho,P,Suc1,Suc2,B,T,Labels,Support,phi);
                }
                else{//estimate the root
                    int* P_new = new int[m+1];
                    double* B_new = new double[m+1];
                    string* Support_new = new string[m+1];            
                    int* Suc1_new = new int[n];
                    int* Suc2_new = new int[n];        
                    V=variance(false,m,B,opt->c,opt->seqLength);
                    int r;
                    double lambda;
                    if (opt->rate==NULL){//the rate is not given
                        if (strcmp(opt->estimate_root,"l")==0){//improve locally the root around the given root
                            if (!opt->relative) phi=estimate_root_with_constraint_local_rooted(n,m,Suc1,Suc2,P,B,V,T,r,lambda,rho,opt->delta);
                            else phi=estimate_root_with_constraint_local_rooted_relative(n,m,Suc1,Suc2,P,B,V,T,r,lambda,rho,opt->delta,opt->mrca);                            
                            reroot_rootedtree(m,r,lambda,P,Suc1[0],Suc2[0],B,Support,P_new,B_new,Support_new);                                                                                 
                            cout<<"Starting QPD on the new rooted tree ..."<<endl;
                            V=variance(opt->variance,m,B_new,opt->c,opt->seqLength);
                            computeSuc(P_new,Suc1_new,Suc2_new,m+1,n);                            
                            if (!opt->relative) phi = with_constraint_active_set(n,m,P_new,Suc1_new,Suc2_new,B_new,V,T,rho,opt->delta);      
                            else phi = with_constraint_active_set_relative(n,m,P_new,Suc1_new,Suc2_new,B_new,V,T,rho,opt->delta,opt->mrca); 
                            output(opt->constraint,opt->variance,result,tree1,tree2,tree3,m,y,rho,P_new,Suc1_new,Suc2_new,B_new,T,Labels,Support_new,phi);
                        }
                        else  if (strcmp(opt->estimate_root,"a")==0){//forget the given root and re-estimate the position of the root over all branhces using fast method
                            if (!opt->relative) phi=estimate_root_with_constraint_fast_rooted(n,m,Suc1[0],Suc2[0],P,B,V,T,r,lambda,rho,opt->delta);
                            else phi=estimate_root_with_constraint_fast_rooted_relative(n,m,Suc1[0],Suc2[0],P,B,V,T,r,lambda,rho,opt->delta,opt->mrca);                            
                            reroot_rootedtree(m,r,lambda,P,Suc1[0],Suc2[0],B,Support,P_new,B_new,Support_new);
                            cout<<"Starting QPD on the new rooted tree ..."<<endl;
                            V=variance(opt->variance,m,B_new,opt->c,opt->seqLength);
                            computeSuc(P_new,Suc1_new,Suc2_new,m+1,n);
                            if (!opt->relative) phi = with_constraint_active_set(n,m,P_new,Suc1_new,Suc2_new,B_new,V,T,rho,opt->delta);  
                            else phi = with_constraint_active_set_relative(n,m,P_new,Suc1_new,Suc2_new,B_new,V,T,rho,opt->delta,opt->mrca); 
                            output(opt->constraint,opt->variance,result,tree1,tree2,tree3,m,y,rho,P_new,Suc1_new,Suc2_new,B_new,T,Labels,Support_new,phi);
                        }
                        else{//forget the given root and re-estimate the position of the root over all branhces using complete method                                                    
                            if (!opt->relative) phi=estimate_root_with_constraint_rooted(n,m,Suc1[0],Suc2[0],P,B,V,T,r,lambda,rho,opt->delta);
                            else phi=estimate_root_with_constraint_rooted_relative(n,m,Suc1[0],Suc2[0],P,B,V,T,r,lambda,rho,opt->delta,opt->mrca);                            
                            reroot_rootedtree(m,r,lambda,P,Suc1[0],Suc2[0],B,Support,P_new,B_new,Support_new);
                            cout<<"Starting QPD on the new rooted tree ..."<<endl;
                            V=variance(opt->variance,m,B_new,opt->c,opt->seqLength);
                            computeSuc(P_new,Suc1_new,Suc2_new,m+1,n);
                            if (!opt->relative) phi = with_constraint_active_set(n,m,P_new,Suc1_new,Suc2_new,B_new,V,T,rho,opt->delta);  
                            else phi = with_constraint_active_set_relative(n,m,P_new,Suc1_new,Suc2_new,B_new,V,T,rho,opt->delta,opt->mrca); 
                            output(opt->constraint,opt->variance,result,tree1,tree2,tree3,m,y,rho,P_new,Suc1_new,Suc2_new,B_new,T,Labels,Support_new,phi);     
                        }
                    }
                    else{ //the rate is given
                        rho=readdouble(gr,(string)opt->rate);
                        if (strcmp(opt->estimate_root,"l")==0){//improve locally the root around the given root
                            if (!opt->relative) phi=estimate_root_with_constraint_local_rate_rooted(rho,n,m,Suc1,Suc2,P,B,V,T,r,lambda);
                            else  phi=estimate_root_with_constraint_local_rate_rooted_relative(rho,n,m,Suc1,Suc2,P,B,V,T,r,lambda,opt->mrca);
                            reroot_rootedtree(m,r,lambda,P,Suc1[0],Suc2[0],B,Support,P_new,B_new,Support_new);                                                        
                            cout<<"Starting QPD on the new rooted tree ..."<<endl;
                            V=variance(opt->variance,m,B_new,opt->c,opt->seqLength);
                            computeSuc(P_new,Suc1_new,Suc2_new,m+1,n);                            
                            if (!opt->relative) phi = with_constraint_active_set_rate(rho,n,m,P_new,Suc1_new,Suc2_new,B_new,V,T);  
                            else phi = with_constraint_active_set_rate_relative(rho,n,m,P_new,Suc1_new,Suc2_new,B_new,V,T,opt->mrca);
                            output(opt->constraint,opt->variance,result,tree1,tree2,tree3,m,y,rho,P_new,Suc1_new,Suc2_new,B_new,T,Labels,Support_new,phi);
                        }
                        else if (strcmp(opt->estimate_root,"a")==0){ //forget the given root and re-estimate the position of the root over all branches using fast method
                            if (!opt->relative) phi=estimate_root_with_constraint_fast_rate_rooted(rho,n,m,Suc1[0],Suc2[0],P,B,V,T,r,lambda);
                            else phi=estimate_root_with_constraint_fast_rate_rooted_relative(rho,n,m,Suc1[0],Suc2[0],P,B,V,T,r,lambda,opt->mrca);                            
                            reroot_rootedtree(m,r,lambda,P,Suc1[0],Suc2[0],B,Support,P_new,B_new,Support_new);
                            cout<<"Starting QPD on the new rooted tree ..."<<endl;
                            V=variance(opt->variance,m,B_new,opt->c,opt->seqLength);        
                            computeSuc(P_new,Suc1_new,Suc2_new,m+1,n);                                                
                            if (!opt->relative) phi = with_constraint_active_set_rate(rho,n,m,P_new,Suc1_new,Suc2_new,B_new,V,T);   
                            else  phi = with_constraint_active_set_rate_relative(rho,n,m,P_new,Suc1_new,Suc2_new,B_new,V,T,opt->mrca); 
                            output(opt->constraint,opt->variance,result,tree1,tree2,tree3,m,y,rho,P_new,Suc1_new,Suc2_new,B_new,T,Labels,Support_new,phi);
                        } 
                        else { //forget the root and re-estimate the position of the root over all branhces using complete method
                            if (!opt->relative) phi=estimate_root_with_constraint_rate_rooted(rho,n,m,Suc1[0],Suc2[0],P,B,V,T,r,lambda);
                            else phi=estimate_root_with_constraint_rate_rooted_relative(rho,n,m,Suc1[0],Suc2[0],P,B,V,T,r,lambda,opt->mrca);                            
                            reroot_rootedtree(m,r,lambda,P,Suc1[0],Suc2[0],B,Support,P_new,B_new,Support_new);
                            cout<<"Starting QPD on the new rooted tree ..."<<endl;
                            V=variance(opt->variance,m,B_new,opt->c,opt->seqLength);        
                            computeSuc(P_new,Suc1_new,Suc2_new,m+1,n);                                                
                            if (!opt->relative) phi = with_constraint_active_set_rate(rho,n,m,P_new,Suc1_new,Suc2_new,B_new,V,T);   
                            else  phi = with_constraint_active_set_rate_relative(rho,n,m,P_new,Suc1_new,Suc2_new,B_new,V,T,opt->mrca); 
                            output(opt->constraint,opt->variance,result,tree1,tree2,tree3,m,y,rho,P_new,Suc1_new,Suc2_new,B_new,T,Labels,Support_new,phi);
                        }                            
                    }                    
                    delete[] P_new;
                    delete[] Suc1_new;
                    delete[] Suc2_new;
                    delete[] B_new;
                    delete[] Support_new;
                }
            }            
        }
        ////////////////////////////unrooted tree////////////////////////////////////////
        else{ //unrooted tree
            int* P_new = new int[m+2];
            double* B_new = new double[m+2];
            double* T_new = new double[m+2];
            double* V_new = new double[m+2];
            string* Support_new = new string[m+2];
            string* Labels_new = new string[m+2];
            int* Suc1_new = new int[n+1];
            int* Suc2_new = new int[n+1];
            V_new=variance(false,m+1,B,opt->c,opt->seqLength);
            int r;
            double lambda;
            if (!opt->constraint){//LD without constraints
                if (opt->rate==NULL){//the rate is not given              
                    if (!opt->relative) phi=estimate_root_without_constraint_unrooted(n,m,P,B,V_new,T,r,lambda,rho,opt->delta);
                    else phi=estimate_root_without_constraint_unrooted_relative(n,m,P,B,V_new,T,r,lambda,rho,opt->delta,opt->mrca);  
                    reroot_unrootedtree(m,r,lambda,P,B,Support,T,P_new,B_new,Support_new,T_new);                     
                    cout<<"Starting LD on the new rooted tree ..."<<endl;
                    V_new=variance(opt->variance,m+1,B_new,opt->c,opt->seqLength);
                    computeSuc(P_new,Suc1_new,Suc2_new,m+2,n+1);                                        
                    if (!opt->relative) phi = without_constraint(n+1,m+1,P_new,Suc1_new,Suc2_new,B_new,V_new,T_new,rho,opt->delta); 
                    else phi = without_constraint_relative(n+1,m+1,P_new,Suc1_new,Suc2_new,B_new,V_new,T_new,rho,opt->delta,opt->mrca);
                    for (int i=n;i<=m;i++) Labels_new[i+1]=Labels[i];
                    output(opt->constraint,opt->variance,result,tree1,tree2,tree3,m+1,y,rho,P_new,Suc1_new,Suc2_new,B_new,T_new,Labels_new,Support_new,phi);
                }
                else {//the rate is given
                    rho=readdouble(gr,(string)opt->rate);
                    if (!opt->relative) phi=estimate_root_without_constraint_rate_unrooted(rho,n,m,P,B,V_new,T,r,lambda);
                    else phi=estimate_root_without_constraint_rate_unrooted_relative(rho,n,m,P,B,V_new,T,r,lambda,opt->mrca);
                    reroot_unrootedtree(m,r,lambda,P,B,Support,T,P_new,B_new,Support_new,T_new);
                    cout<<"Starting LD on the new rooted tree ..."<<endl;
                    V_new=variance(opt->variance,m+1,B_new,opt->c,opt->seqLength);
                    computeSuc(P_new,Suc1_new,Suc2_new,m+2,n+1);                    
                    if (!opt->relative) phi = without_constraint_rate(rho,n+1,m+1,P_new,Suc1_new,Suc2_new,B_new,V_new,T_new);     
                    else phi = without_constraint_rate_relative(rho,n+1,m+1,P_new,Suc1_new,Suc2_new,B_new,V_new,T_new,opt->mrca);                      
                    for (int i=n;i<=m;i++) Labels_new[i+1]=Labels[i];
                    output(opt->constraint,opt->variance,result,tree1,tree2,tree3,m+1,y,rho,P_new,Suc1_new,Suc2_new,B_new,T_new,Labels_new,Support_new,phi);
                }
            }
            else {//QPD
                if (opt->rate==NULL){//the rate is not given
                   if (strcmp(opt->estimate_root,"a")==0){//estimate the root using fast method
                      cout<<"Estimating the root using fast method."<<endl;;
                      if (!opt->relative) phi=estimate_root_without_constraint_unrooted(n,m,P,B,V_new,T,r,lambda,rho,opt->delta);
                      else phi=estimate_root_without_constraint_unrooted_relative(n,m,P,B,V_new,T,r,lambda,rho,opt->delta,opt->mrca);
                      reroot_unrootedtree(m,r,lambda,P,B,T,P_new,B_new,T_new);
                      computeSuc(P_new,Suc1_new,Suc2_new,m+2,n+1);
                      cout<<"The branches are re-enumarated."<<endl;
                      if (!opt->relative) phi=estimate_root_with_constraint_local_rooted(n+1,m+1,Suc1_new,Suc2_new,P_new,B_new,V_new,T_new,r,lambda,rho,opt->delta);
                      else phi=estimate_root_with_constraint_local_rooted_relative(n+1,m+1,Suc1_new,Suc2_new,P_new,B_new,V_new,T_new,r,lambda,rho,opt->delta,opt->mrca);                      
                      int* P_n = new int[m+2];
                      double* B_n = new double[m+2];
                      string* S_n = new string[m+2];
                      reroot_rootedtree(m+1,r,lambda,P_new,Suc1_new[0],Suc2_new[0],B_new,Support_new,P_n,B_n,S_n);
                      for (int i=0;i<m+2;i++){
                        P_new[i]=P_n[i];
                        B_new[i]=B_n[i];
                        Support_new[i]=S_n[i];
                      }
                      delete[] P_n;
                      delete[] B_n;
                      delete[] S_n;
                      computeSuc(P_new,Suc1_new,Suc2_new,m+2,n+1);
                      V_new=variance(opt->variance,m+1,B_new,opt->c,opt->seqLength);
                      cout<<"Starting QPD on the new rooted tree ..."<<endl;
                      if (!opt->relative) phi = with_constraint_active_set(n+1,m+1,P_new,Suc1_new,Suc2_new,B_new,V_new,T_new,rho,opt->delta); 
                      else phi = with_constraint_active_set_relative(n+1,m+1,P_new,Suc1_new,Suc2_new,B_new,V_new,T_new,rho,opt->delta,opt->mrca);                     
                      for (int i=n;i<=m;i++) Labels_new[i+1]=Labels[i];                    
                      output(opt->constraint,opt->variance,result,tree1,tree2,tree3,m+1,y,rho,P_new,Suc1_new,Suc2_new,B_new,T_new,Labels_new,Support_new,phi);
                   }
                   else {//estimate the root using complete method
                      if (!opt->relative) phi=estimate_root_with_constraint_unrooted(n,m,P,B,V_new,T,r,lambda,rho,opt->delta);
                      else phi=estimate_root_with_constraint_unrooted_relative(n,m,P,B,V_new,T,r,lambda,rho,opt->delta,opt->mrca);
                      reroot_unrootedtree(m,r,lambda,P,B,T,P_new,B_new,T_new);                                                               
                      V_new=variance(opt->variance,m+1,B_new,opt->c,opt->seqLength);
                      cout<<"Starting QPD on the new rooted tree ..."<<endl;
                      computeSuc(P_new,Suc1_new,Suc2_new,m+2,n+1);   
                      if (!opt->relative) phi = with_constraint_active_set(n+1,m+1,P_new,Suc1_new,Suc2_new,B_new,V_new,T_new,rho,opt->delta); 
                      else phi = with_constraint_active_set_relative(n+1,m+1,P_new,Suc1_new,Suc2_new,B_new,V_new,T_new,rho,opt->delta,opt->mrca);                     
                      for (int i=n;i<=m;i++) Labels_new[i+1]=Labels[i];                    
                      output(opt->constraint,opt->variance,result,tree1,tree2,tree3,m+1,y,rho,P_new,Suc1_new,Suc2_new,B_new,T_new,Labels_new,Support_new,phi);    
                   }
                }
                else {//the rate is known
                     if (strcmp(opt->estimate_root,"a")==0){//estimate the root using fast methode
                        rho=readdouble(gr,(string)(opt->rate));
                        cout<<"Estimating the root using fast method."<<endl;
                        if (!opt->relative) phi=estimate_root_without_constraint_rate_unrooted(rho,n,m,P,B,V_new,T,r,lambda);
                        else phi=estimate_root_without_constraint_rate_unrooted_relative(rho,n,m,P,B,V_new,T,r,lambda,opt->mrca);
                        reroot_unrootedtree(m,r,lambda,P,B,T,P_new,B_new,T_new);
                        computeSuc(P_new,Suc1_new,Suc2_new,m+2,n+1);
                        cout<<"The branches are re-enumarated."<<endl;
                        if (!opt->relative) phi=estimate_root_with_constraint_local_rate_rooted(rho,n+1,m+1,Suc1_new,Suc2_new,P_new,B_new,V_new,T_new,r,lambda);
                        else phi=estimate_root_with_constraint_local_rate_rooted_relative(rho,n+1,m+1,Suc1_new,Suc2_new,P_new,B_new,V_new,T_new,r,lambda,opt->mrca);                        
                        int* P_n = new int[m+2];
                        double* B_n = new double[m+2];
                        string* S_n = new string[m+2];
                        reroot_rootedtree(m+1,r,lambda,P_new,Suc1_new[0],Suc2_new[0],B_new,Support_new,P_n,B_n,S_n);
                        for (int i=0;i<m+2;i++){
                            P_new[i]=P_n[i];
                            B_new[i]=B_n[i];
                            Support_new[i]=S_n[i];
                        }
                        delete[] P_n;
                        delete[] B_n;
                        delete[] S_n;
                        computeSuc(P_new,Suc1_new,Suc2_new,m+2,n+1);
                        V_new=variance(opt->variance,m,B_new,opt->c,opt->seqLength);
                        cout<<"Starting QPD on the new rooted tree ..."<<endl;
                        if (!opt->relative) phi = with_constraint_active_set_rate(rho,n+1,m+1,P_new,Suc1_new,Suc2_new,B_new,V_new,T_new); 
                        else phi = with_constraint_active_set_rate_relative(rho,n+1,m+1,P_new,Suc1_new,Suc2_new,B_new,V_new,T_new,opt->mrca);                     
                        for (int i=n;i<=m;i++) Labels_new[i+1]=Labels[i];
                        output(opt->constraint,opt->variance,result,tree1,tree2,tree3,m+1,y,rho,P_new,Suc1_new,Suc2_new,B_new,T_new,Labels_new,Support_new,phi);
                     }
                     else{//estimate the root using complete methode
                        rho=readdouble(gr,(string)(opt->rate));
                        if (!opt->relative) phi=estimate_root_with_constraint_rate_unrooted(rho,n,m,P,B,V_new,T,r,lambda);
                        else phi=estimate_root_with_constraint_rate_unrooted_relative(rho,n,m,P,B,V_new,T,r,lambda,opt->mrca);
                        reroot_unrootedtree(m,r,lambda,P,B,T,P_new,B_new,T_new);                        
                        computeSuc(P_new,Suc1_new,Suc2_new,m+2,n+1);
                        V_new=variance(opt->variance,m,B_new,opt->c,opt->seqLength);
                        cout<<"Starting QPD on the new rooted tree ..."<<endl;
                        if (!opt->relative) phi = with_constraint_active_set_rate(rho,n+1,m+1,P_new,Suc1_new,Suc2_new,B_new,V_new,T_new); 
                        else phi = with_constraint_active_set_rate_relative(rho,n+1,m+1,P_new,Suc1_new,Suc2_new,B_new,V_new,T_new,opt->mrca);                     
                        for (int i=n;i<=m;i++) Labels_new[i+1]=Labels[i];
                        output(opt->constraint,opt->variance,result,tree1,tree2,tree3,m+1,y,rho,P_new,Suc1_new,Suc2_new,B_new,T_new,Labels_new,Support_new,phi);
                     }
                }
            }
            delete[] P_new;
            delete[] Suc1_new;
            delete[] Suc2_new;
            delete[] B_new;
            delete[] T_new;
            delete[] V_new;
            delete[] Support_new;
            delete[] Labels_new;
        }
    }
	elapsed_time = (double)(clock()-start)/CLOCKS_PER_SEC;
	fprintf(result, "\n*********************************************************\n" );
	fprintf(result, "Elapsed time: %.2f seconds\n", elapsed_time);
    fprintf(tree1,"End;");
	cout<<"Elapsed time: "<<elapsed_time<<" seconds"<<endl;;
	fclose(tree);
    fclose(result);
    fclose(tree1);
    fclose(tree2);
    fclose(tree3);
    if (opt->rate!=NULL) fclose(gr);
    delete[] P;
    delete[] B;
    delete[] V;
    delete[] T;
    delete[] Suc1;
    delete[] Suc2;
    delete[] Support;
    delete[] Labels;
	return EXIT_SUCCESS;
}

