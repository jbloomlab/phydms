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
#include "estimate_root_relative_unrooted.h"


double estimate_root_without_constraint_unrooted_relative(int n,int m,int* & Pur,double* & Bur,double* &V,double* & Tur,int &r,double &lambda,double &rho,double rho_min,double mrca){
    //Pur: unrooted tree
    //estimate root with LD algorithm for unrooted tree///////////////////////////
    cout<<"Estimating the root without constraints ... ";      
    double phi1=-1;
    int *P_new= new int[m+2];
    int *Suc1_new= new int[n+1];
    int *Suc2_new= new int[n+1];
    double *B_new= new double[m+2];
    double *T_new= new double[m+2];
    int y=1;
    double br=reroot_unrootedtree(m,y,Pur,Bur,Tur,P_new,B_new,T_new);
    computeSuc(P_new,Suc1_new,Suc2_new,m+2,n+1);
    double ld,rhor;
    //cout<<"Optimizing the root position on the branch "<<y<<" ... ";
    phi1=without_constraint_lambda_relative(n+1,m+1,br,P_new,Suc1_new,Suc2_new,B_new,V,T_new,rhor,rho_min,ld,mrca);
    //printf("%.10f\n",phi1);
    lambda=ld;
    rho=rhor;
    r=y;
    y++;
    double phi;
    while (y<=m){
        br=reroot_unrootedtree(m,y,Pur,Bur,Tur,P_new,B_new,T_new);
        computeSuc(P_new,Suc1_new,Suc2_new,m+2,n+1);
        //cout<<"Optimizing the root position on the branch "<<y<<" ... ";
        phi=without_constraint_lambda_relative(n+1,m+1,br,P_new,Suc1_new,Suc2_new,B_new,V,T_new,rhor,rho_min,ld,mrca);
        //printf("%.10f\n",phi);
        if (phi1>phi){
            phi1=phi;
            rho = rhor;
            r=y;
            lambda=ld;
        }
        y++;
    }
    cout<<"The tree is rooted on the branch "<<r<<endl;    
    delete[] P_new;
    delete[] Suc1_new;
    delete[] Suc2_new;
    delete[] B_new;
    delete[] T_new;
    return phi1;
}

double estimate_root_without_constraint_rate_unrooted_relative(double rho,int n,int m,int* & Pur,double* & Bur,double* &V,double* & Tur,int &r,double &lambda,double mrca){
    //Pur: unrooted tree
    cout<<"Estimating the root without constraints ... ";    
    double phi1=-1;
    int *P_new= new int[m+2];
    int *Suc1_new= new int[n+1];
    int *Suc2_new= new int[n+1];
    double *B_new= new double[m+2];
    double *T_new= new double[m+2];
    int y=1;
    double br=reroot_unrootedtree(m,y,Pur,Bur,Tur,P_new,B_new,T_new);
    computeSuc(P_new,Suc1_new,Suc2_new,m+2,n+1);
    double ld;
    //cout<<"Optimizing the root position on the branch "<<y<<" ... ";
    phi1=without_constraint_lambda_rate_relative(rho,n+1,m+1,br,P_new,Suc1_new,Suc2_new,B_new,V,T_new,ld,mrca);
    //printf("%.10f\n",phi1); 
    lambda=ld;
    r=y;
    y++;
    double phi;
    while (y<=m){
        br=reroot_unrootedtree(m,y,Pur,Bur,Tur,P_new,B_new,T_new);
        computeSuc(P_new,Suc1_new,Suc2_new,m+2,n+1);
        //cout<<"Optimizing the root position on the branch "<<y<<" ... ";
        phi=without_constraint_lambda_rate_relative(rho,n+1,m+1,br,P_new,Suc1_new,Suc2_new,B_new,V,T_new,ld,mrca);
        //printf("%.10f\n",phi);   
        if (phi1>phi){
            phi1=phi;
            r=y;
            lambda=ld;
        }
        y++;
    }
    cout<<"The tree is rooted on the branch "<<r<<endl;    
    delete[] P_new;
    delete[] Suc1_new;
    delete[] Suc2_new;
    delete[] B_new;
    delete[] T_new;
    return phi1;
}


double estimate_root_with_constraint_unrooted_relative(int n,int m,int* & Pur,double* & Bur,double* &V,double* & Tur,int &r,double &lambda,double &rho,double rho_min,double mrca){
    //Pur: unrooted tree
    //estimate root with QPD algorithm for unrooted tree////////////////////////////
    cout<<"Estimating the root using constrained mode on all branches ... ";    
    double phi1=-1;
    int *P_new= new int[m+2];
    int *Suc1_new= new int[n+1];
    int *Suc2_new= new int[n+1];
    double *B_new= new double[m+2];
    double *T_new= new double[m+2];
    int y=1;
    double br=reroot_unrootedtree(m,y,Pur,Bur,Tur,P_new,B_new,T_new);
    computeSuc(P_new,Suc1_new,Suc2_new,m+2,n+1);
    double ld,rhor;
    //cout<<"Optimizing the root position on the branch "<<y<<" ... ";
    phi1=with_constraint_active_set_lambda_relative(n+1,m+1,br,P_new,Suc1_new,Suc2_new,B_new,V,T_new,rhor,rho_min,ld,mrca);
    //printf("%.10f\n",phi1);
    lambda=ld;
    rho=rhor;
    r=y;
    y++;
    double phi;
    while (y<=m){
        br=reroot_unrootedtree(m,y,Pur,Bur,Tur,P_new,B_new,T_new);
        computeSuc(P_new,Suc1_new,Suc2_new,m+2,n+1);
        //cout<<"Optimizing the root position on the branch "<<y<<" ... ";
        phi=with_constraint_active_set_lambda_relative(n+1,m+1,br,P_new,Suc1_new,Suc2_new,B_new,V,T_new,rhor,rho_min,ld,mrca);
        //printf("%.10f\n",phi);
        if (phi1>phi){
            phi1=phi;
            r=y;
            lambda=ld;
        }
        y++;
    }
    cout<<"The tree is rooted on the branch "<<r<<endl;    
    delete[] P_new;
    delete[] Suc1_new;
    delete[] Suc2_new;
    delete[] B_new;
    delete[] T_new;
    return phi1;
}


double estimate_root_with_constraint_rate_unrooted_relative(double rho,int n,int m,int* & Pur,double* & Bur,double* &V,double* & Tur,int &r,double &lambda,double mrca){
    //Pur: unrooted tree
    cout<<"Estimating the root using contrained mode on all branches ... ";    
    double phi1=-1;
    int *P_new= new int[m+2];
    int *Suc1_new= new int[n+1];
    int *Suc2_new= new int[n+1];
    double *B_new= new double[m+2];
    double *T_new= new double[m+2];
    int y=1;
    double br=reroot_unrootedtree(m,y,Pur,Bur,Tur,P_new,B_new,T_new);
    computeSuc(P_new,Suc1_new,Suc2_new,m+2,n+1);
    double ld;
    //cout<<"Optimizing the root position on the branch "<<y<<" ... ";
    phi1=with_constraint_active_set_lambda_rate_relative(rho,n+1,m+1,br,P_new,Suc1_new,Suc2_new,B_new,V,T_new,ld,mrca);
    //printf("%.10f\n",phi1);
    lambda=ld;
    r=y;
    y++;
    double phi;
    while (y<=m){
        br=reroot_unrootedtree(m,y,Pur,Bur,Tur,P_new,B_new,T_new);
        computeSuc(P_new,Suc1_new,Suc2_new,m+2,n+1);
        //cout<<"Optimizing the root position on the branch "<<y<<" ... ";
        phi=with_constraint_active_set_lambda_rate_relative(rho,n+1,m+1,br,P_new,Suc1_new,Suc2_new,B_new,V,T_new,ld,mrca);
        //printf("%.10f\n",phi);
        if (phi1>phi){
            phi1=phi;
            r=y;
            lambda=ld;
        }
        y++;
    }
    cout<<"The tree is rooted on the branch "<<r<<endl;    
    delete[] P_new;
    delete[] Suc1_new;
    delete[] Suc2_new;
    delete[] B_new;
    delete[] T_new;
    return phi1;
}

