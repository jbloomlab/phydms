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
#include "dating.h"


double without_constraint(int n,int m,int* & P,int* & Suc1,int* & Suc2,double* & B,double* & V,double* & T,double &rho,double rho_min){
    //This function implements LD algorithm.
    //The estimated substitution rate and dates of all internal nodes are returned in variable paramaters.
    //It returns the value of the objective function
    list<int> pos = postorder(P,Suc1,Suc2,n);
    double *W= new double[n];
    double *C = new double[n];    
    double *D= new double[n];//T[i]=D[i]+E[i]/rho
    double * E = new double[n];
    double *X =new double[n];//T[i]=W[i].T[a(i)]+C[i]+X[i]/rho    
    for (list<int>::iterator iter=pos.begin();iter!=pos.end();iter++){
        int i =  *iter;
        int s1 = Suc1[i];
        int s2 = Suc2[i];
        if (i==0){//(1/V[s1]+1/V[s2])t[i]=t[s1]/V[s1]+t[s2]/V[s2]-(b[s1]/V[s1]+b[s2]/V[s2])/rho
            if (s1>=n && s2>=n){
                D[0]=(T[s1]*V[s2]+T[s2]*V[s1])/(V[s1]+V[s2]);
                E[0]=-(B[s1]*V[s2]+B[s2]*V[s1])/(V[s1]+V[s2]);
            }
            else if (s1>=n){
                D[0]=(V[s2]*T[s1]+V[s1]*C[s2])/(V[s1]+V[s2]-W[s2]*V[s1]);
                E[0]=(V[s1]*X[s2]-B[s1]*V[s2]-B[s2]*V[s1])/(V[s1]+V[s2]-W[s2]*V[s1]);
            }
            else if (s2>=n){
                D[0]=(V[s1]*T[s2]+V[s2]*C[s1])/(V[s1]+V[s2]-W[s1]*V[s2]);
                E[0]=(V[s2]*X[s1]-B[s1]*V[s2]-B[s2]*V[s1])/(V[s1]+V[s2]-W[s1]*V[s2]);
            }
            else{
                D[0]=(C[s1]*V[s2]+C[s2]*V[s1])/(V[s1]+V[s2]-W[s1]*V[s2]-W[s2]*V[s1]);
                E[0]=(X[s1]*V[s2]+X[s2]*V[s1]-B[s1]*V[s2]-B[s2]*V[s1])/(V[s1]+V[s2]-W[s1]*V[s2]-W[s2]*V[s1]);
            }
            
        }
        else{//(1/V[i]+1/V[s1]+1/V[s2])t[i]=t[p]/V[i]+t[s1]/V[s1]+t[s2]/V[s2]+(B[i]/V[i]-B[s1]/V[s1]-B[s2]/V[s2])/rho
            if (s1>=n && s2>=n){
                W[i]=V[s1]*V[s2]/(V[s1]*V[s2]+V[s1]*V[i]+V[s2]*V[i]);
                C[i]=(V[s2]*T[s1]+V[s1]*T[s2])*V[i]/(V[s1]*V[s2]+V[s1]*V[i]+V[s2]*V[i]);
                X[i]=(B[i]*V[s1]*V[s2]-B[s1]*V[i]*V[s2]-B[s2]*V[i]*V[s1])/(V[s1]*V[s2]+V[s1]*V[i]+V[s2]*V[i]);
            }
            else if (s1>=n){
                W[i]=V[s1]*V[s2]/(V[s1]*V[s2]+V[i]*V[s2]+(1-W[s2])*V[i]*V[s1]);
                C[i]=(T[s1]*V[s2]+C[s2]*V[s1])*V[i]/(V[s1]*V[s2]+V[i]*V[s2]+(1-W[s2])*V[i]*V[s1]);
                X[i]=(X[s2]*V[i]*V[s1]+B[i]*V[s1]*V[s2]-B[s1]*V[i]*V[s2]-B[s2]*V[i]*V[s1])/(V[s1]*V[s2]+V[i]*V[s2]+(1-W[s2])*V[i]*V[s1]);
            }
            else if (s2>=n){
                W[i]=V[s1]*V[s2]/(V[s1]*V[s2]+V[i]*V[s1]+(1-W[s1])*V[i]*V[s2]);
                C[i]=(T[s2]*V[s1]+C[s1]*V[s2])*V[i]/(V[s1]*V[s2]+V[i]*V[s1]+(1-W[s1])*V[i]*V[s2]);
                X[i]=(X[s1]*V[i]*V[s2]+B[i]*V[s1]*V[s2]-B[s1]*V[i]*V[s2]-B[s2]*V[i]*V[s1])/(V[s1]*V[s2]+V[i]*V[s1]+(1-W[s1])*V[i]*V[s2]);
            }
            else{
                W[i]=V[s1]*V[s2]/(V[s1]*V[s2]+V[i]*V[s2]*(1-W[s1])+V[i]*V[s1]*(1-W[s2]));
                C[i]=(C[s1]*V[s2]+C[s2]*V[s1])*V[i]/(V[s1]*V[s2]+V[i]*V[s2]*(1-W[s1])+V[i]*V[s1]*(1-W[s2]));
                X[i]=(X[s1]*V[i]*V[s2]+X[s2]*V[i]*V[s1]+B[i]*V[s1]*V[s2]-B[s1]*V[i]*V[s2]-B[s2]*V[i]*V[s1])/(V[s1]*V[s2]+V[i]*V[s2]*(1-W[s1])+V[i]*V[s1]*(1-W[s2]));
	        }
        }
    }
    list<int> pr = preorder(P,Suc1,Suc2,n);
    for (list<int>::iterator iter = pr.begin();iter!=pr.end();iter++){
        int i = *iter;
        if (i!=0){
            D[i]=W[i]*D[P[i]]+C[i];
            E[i]=W[i]*E[P[i]]+X[i];
        }
    }
    double *F = new double[m+1];
    double *G = new double[m+1];
    // T[P[i]]-T[i]=F[i]+G[i]/rho
    F[0]=0;
    G[0]=0;
    for (int i = 1;i<n;i++){
        F[i] = D[P[i]]-D[i];
        G[i] = E[P[i]]-E[i];
    }
    
    for (int i = n;i<=m;i++){
        F[i] = D[P[i]]-T[i];
        G[i] = E[P[i]];
    }
	
    double a = 0;
    double b = 0;
    double c = 0;//phi = a*rho^2+b*rho+c
    for (int i=1;i<=m;i++){
        a += F[i]*F[i]/V[i];
        b += 2*(B[i]+G[i])*F[i]/V[i];
        c += (B[i]+G[i])*(B[i]+G[i])/V[i];
    }
    rho = -b/(2*a);
    if (rho<rho_min) rho=rho_min;    
    for (int i =0;i<n;i++) T[i] = D[i]+E[i]/rho;
    delete[] W;
    delete[] C;
    delete[] D;
    delete[] E;
    delete[] F;
    delete[] G;
    delete[] X;
    return phi(n,m,P,B,V,T,rho);//a*rho*rho+b*rho+c;
}

double without_constraint_rate(double rho,int n,int m,int* & P,int* & Suc1,int* & Suc2,double* & B,double* & V,double* & T){
    //This function implements LD algorithm in the case that the rate is given
    list<int> pos = postorder(P,Suc1,Suc2,n);
    double *W= new double[n];
    double *C = new double[n];
    double *X =new double[n];//T[i]=W[i].T[a(i)]+C[i]+X[i]/rho
    double *D= new double[n];//T[i]=D[i]+E[i]/rho
    double * E = new double[n];
    for (list<int>::iterator iter=pos.begin();iter!=pos.end();iter++){
        int i =  *iter;
        int s1 = Suc1[i];
        int s2 = Suc2[i];
        if (i==0){//(1/V[s1]+1/V[s2])t[i]=t[s1]/V[s1]+t[s2]/V[s2]-(b[s1]/V[s1]+b[s2]/V[s2])/rho
            if (s1>=n && s2>=n){
                D[0]=(T[s1]*V[s2]+T[s2]*V[s1])/(V[s1]+V[s2]);
                E[0]=-(B[s1]*V[s2]+B[s2]*V[s1])/(V[s1]+V[s2]);
            }
            else if (s1>=n){
                D[0]=(V[s2]*T[s1]+V[s1]*C[s2])/(V[s1]+V[s2]-W[s2]*V[s1]);
                E[0]=(V[s1]*X[s2]-B[s1]*V[s2]-B[s2]*V[s1])/(V[s1]+V[s2]-W[s2]*V[s1]);
            }
            else if (s2>=n){
                D[0]=(V[s1]*T[s2]+V[s2]*C[s1])/(V[s1]+V[s2]-W[s1]*V[s2]);
                E[0]=(V[s2]*X[s1]-B[s1]*V[s2]-B[s2]*V[s1])/(V[s1]+V[s2]-W[s1]*V[s2]);
            }
            else{
                D[0]=(C[s1]*V[s2]+C[s2]*V[s1])/(V[s1]+V[s2]-W[s1]*V[s2]-W[s2]*V[s1]);
                E[0]=(X[s1]*V[s2]+X[s2]*V[s1]-B[s1]*V[s2]-B[s2]*V[s1])/(V[s1]+V[s2]-W[s1]*V[s2]-W[s2]*V[s1]);
            }
        }
        else{//(1/V[i]+1/V[s1]+1/V[s2])t[i]=t[p]/V[i]+t[s1]/V[s1]+t[s2]/V[s2]+(B[i]/V[i]-B[s1]/V[s1]-B[s2]/V[s2])/rho
            if (s1>=n && s2>=n){
                W[i]=V[s1]*V[s2]/(V[s1]*V[s2]+V[s1]*V[i]+V[s2]*V[i]);
                C[i]=(V[s2]*T[s1]+V[s1]*T[s2])*V[i]/(V[s1]*V[s2]+V[s1]*V[i]+V[s2]*V[i]);
                X[i]=(B[i]*V[s1]*V[s2]-B[s1]*V[i]*V[s2]-B[s2]*V[i]*V[s1])/(V[s1]*V[s2]+V[s1]*V[i]+V[s2]*V[i]);
            }
            else if (s1>=n){
                W[i]=V[s1]*V[s2]/(V[s1]*V[s2]+V[i]*V[s2]+(1-W[s2])*V[i]*V[s1]);
                C[i]=(T[s1]*V[s2]+C[s2]*V[s1])*V[i]/(V[s1]*V[s2]+V[i]*V[s2]+(1-W[s2])*V[i]*V[s1]);
                X[i]=(X[s2]*V[i]*V[s1]+B[i]*V[s1]*V[s2]-B[s1]*V[i]*V[s2]-B[s2]*V[i]*V[s1])/(V[s1]*V[s2]+V[i]*V[s2]+(1-W[s2])*V[i]*V[s1]);
            }
            else if (s2>=n){
                W[i]=V[s1]*V[s2]/(V[s1]*V[s2]+V[i]*V[s1]+(1-W[s1])*V[i]*V[s2]);
                C[i]=(T[s2]*V[s1]+C[s1]*V[s2])*V[i]/(V[s1]*V[s2]+V[i]*V[s1]+(1-W[s1])*V[i]*V[s2]);
                X[i]=(X[s1]*V[i]*V[s2]+B[i]*V[s1]*V[s2]-B[s1]*V[i]*V[s2]-B[s2]*V[i]*V[s1])/(V[s1]*V[s2]+V[i]*V[s1]+(1-W[s1])*V[i]*V[s2]);
            }
            else{
                W[i]=V[s1]*V[s2]/(V[s1]*V[s2]+V[i]*V[s2]*(1-W[s1])+V[i]*V[s1]*(1-W[s2]));
                C[i]=(C[s1]*V[s2]+C[s2]*V[s1])*V[i]/(V[s1]*V[s2]+V[i]*V[s2]*(1-W[s1])+V[i]*V[s1]*(1-W[s2]));
                X[i]=(X[s1]*V[i]*V[s2]+X[s2]*V[i]*V[s1]+B[i]*V[s1]*V[s2]-B[s1]*V[i]*V[s2]-B[s2]*V[i]*V[s1])/(V[s1]*V[s2]+V[i]*V[s2]*(1-W[s1])+V[i]*V[s1]*(1-W[s2]));
	        }
        }
    }
    list<int> pr = preorder(P,Suc1,Suc2,n);
    for (list<int>::iterator iter = pr.begin();iter!=pr.end();iter++){
        int i = *iter;
        if (i!=0){
            D[i]=W[i]*D[P[i]]+C[i];
            E[i]=W[i]*E[P[i]]+X[i];
        }
    }
    for (int i =0;i<n;i++) T[i] = D[i]+E[i]/rho;
    //double phi=0;
    //for (int i = 1;i<=m;i++) phi += (B[i]-rho*(T[i]-T[P[i]]))*(B[i]-rho*(T[i]-T[P[i]]))/V[i];
    delete[] W;
    delete[] C;
    delete[] X;
    delete[] D;
    delete[] E;
    return phi(n,m,P,B,V,T,rho);
}

/////////////////////////////////////////////////////////////////////////////////////////////////
//////////////This part implements QPD algorithm (active set methods)////////////////////////////

int remove_ne_lambda(list<double> & lambda,list<int> & active_set){
    //remove the active set corresponding to the most negative lagrange multiplier
    double ne_lambda=0;
    double min = 0;
    int as=-1;
    list<int>::iterator iteras=active_set.begin();
    for (list<double>::iterator iter = lambda.begin();iter!=lambda.end();iter++){
        ne_lambda = *iter;
        if (ne_lambda<min){
            min = ne_lambda;
            as=*iteras;
        }
        iteras++;
    }
    return as;
    
}

bool conditions(int m,list<double>& ldLagrange,int* & P,double* & T){
    //the stop condition of active set method
    for (list<double>::iterator iter = ldLagrange.begin();iter!=ldLagrange.end();iter++){
        if (*iter < 0) {
            return false;
        }
    }
    for (int i=1;i<=m;i++)
        if (T[i]<T[P[i]]-(10e-15)) {
            return false;
        }
    return true;
}

double starting_point(int n,int m,int* & P,int* & Suc1,int* & Suc2,double* & B,double* & V,double* & T,double &rho,double rho_min,list<int> & active_set){
    //compute a starting feasible point, which is the solution from without constraint and then collapse all the branches that violate the constraints.
    without_constraint(n,m,P,Suc1,Suc2,B,V,T,rho,rho_min);
    list<int> pos = postorder(P,Suc1,Suc2,n);
    for (list<int>::iterator iter = pos.begin();iter!=pos.end();iter++){
        int i = *iter;
        int s1 = Suc1[i];
        int s2 = Suc2[i];
        if (T[i]>T[s1] || T[i]>T[s2]){
            if (T[s1]<T[s2]){
                T[i]=T[s1];
                active_set.push_back(s1);
            }
            else{
                T[i]=T[s2];
                active_set.push_back(s2);
            }
        }
    }
    return phi(n,m,P,B,V,T,rho);
}

list<double> without_constraint_active_set(int n,int m,bool* & flag,int* & P,int* & Suc1,int* & Suc2,double* & B,double* & V,double* & T,double &rho,double rho_min,list<int>& active_set,double &phi1){
    //this method computes the optimized solution of the Lagrange function.
    //it returns the list of Lagrange multipliers
    list<double> ld;
    list<int>* SucL = new list<int>[n];
    list<int>* SucI = new list<int>[n];
    int* Pre=new int[m+1];
    for (int i=0;i<=m;i++) Pre[i]=-1;
    list<int> ls;
    for (int i=0;i<=m;i++) flag[i]=true;
    bool* leaf=new bool[m+1];
    for (int i=0;i<n;i++) leaf[i]=false;
    for (int i=n;i<=m;i++) leaf[i]=true;
    //flag
    
    for (list<int>::iterator iter=active_set.begin();iter!=active_set.end();iter++){
        int i = *iter;
        flag[i]=false;
        if (i>=n && i!=m+1 && i!=m+2) ls.push_back(i);
    }
    
    stack<int>* feuilles = new stack<int>[ls.size()];
    
    list<int> top;
    //leaf
    computeFeuilles(n,m,P,Suc1,Suc2,T,flag,leaf,feuilles,active_set);
    for (int i=0;i<n;i++){
        int s1= Suc1[i];
        int s2= Suc2[i];
        if (flag[i] && !leaf[i] && (!flag[s1] || !flag[s2]))
            top.push_back(i);
    }
    list<int>* internal = new list<int>[top.size()];
    
    reduceTree(n,m,P,Pre,Suc1,Suc2,SucL,SucI,flag,leaf,internal,active_set);
    list<int> pos = postorder(P,Suc1,Suc2,n);
    list<int> pre = preorder(P,Suc1,Suc2,n);
    double *W= new double[n];
    double *C = new double[n];
    double *X =new double[n];
    double *D= new double[n];
    double * E = new double[n];
    for (int i=0;i<n;i++){
        W[i]=0;C[i]=0;X[i]=0;D[i]=0;E[i]=0;
    }
    for (list<int>::iterator it = pos.begin();it!=pos.end();it++){
        int i = *it;
        if (flag[i] && !leaf[i]){
            int p = Pre[i];
            if (p==-1){
                double coefs=0;
                double xtemp=0;
                double ctemp=0;
                for (list<int>::iterator iter = SucL[i].begin();iter!=SucL[i].end();iter++){
                    coefs+=1/V[*iter];
                    ctemp+=T[*iter]/V[*iter];
                    xtemp-=B[*iter]/V[*iter];
                }
                for (list<int>::iterator iter = SucI[i].begin();iter!=SucI[i].end();iter++){
                    coefs+=(1-W[*iter])/V[*iter];
                    ctemp+=C[*iter]/V[*iter];
                    xtemp+=(X[*iter]-B[*iter])/V[*iter];
                }
                D[i]=ctemp/coefs;
                E[i]=xtemp/coefs;
            }
            else{
                if (!leaf[p]){
                    double coefs=1/V[i];
                    double wtemp=1/V[i];
                    double ctemp=0;
                    double xtemp=B[i]/V[i];
                    for (list<int>::iterator iter = SucL[i].begin();iter!=SucL[i].end();iter++){
                        coefs+=1/V[*iter];
                        ctemp+=T[*iter]/V[*iter];
                        xtemp-=B[*iter]/V[*iter];
                    }
                    for (list<int>::iterator iter = SucI[i].begin();iter!=SucI[i].end();iter++){
                        coefs+=(1-W[*iter])/V[*iter];
                        ctemp+=C[*iter]/V[*iter];
                        xtemp+=(X[*iter]-B[*iter])/V[*iter];
                    }
                    W[i]=wtemp/coefs;
                    C[i]=ctemp/coefs;
                    X[i]=xtemp/coefs;
                }
                else{
                    double coefs=1/V[i];
                    double dtemp=T[p]/V[i];
                    double etemp=B[i]/V[i];
                    for (list<int>::iterator iter = SucL[i].begin();iter!=SucL[i].end();iter++){
                        coefs+=1/V[*iter];
                        dtemp+=T[*iter]/V[*iter];
                        etemp-=B[*iter]/V[*iter];
                    }
                    for (list<int>::iterator iter = SucI[i].begin();iter!=SucI[i].end();iter++){
                        coefs+=(1-W[*iter])/V[*iter];
                        dtemp+=C[*iter]/V[*iter];
                        etemp+=(X[*iter]-B[*iter])/V[*iter];
                    }
                    D[i]=dtemp/coefs;
                    E[i]=etemp/coefs;
                }
            }
        }
    }
    delete[] SucL;
    delete[] SucI;
    for (list<int>::iterator iter = pre.begin();iter!=pre.end();iter++){
        int i = *iter;
        if (Pre[i]!=-1 && flag[i] && !leaf[i] && !leaf[Pre[i]]){
            D[i]=W[i]*D[Pre[i]]+C[i];
            E[i]=W[i]*E[Pre[i]]+X[i];
        }
        int s1 = Suc1[i];
        int s2 = Suc2[i];
        if (!flag[s1] && !leaf[s1]) {
            D[s1]=D[i];
            E[s1]=E[i];
        }
        if (!flag[s2] && !leaf[s2]) {
            D[s2]=D[i];
            E[s2]=E[i];
        }
    }
    double *F=new double[m+1];
    double *G=new double[m+1];//T[P[i]]-T[i]=F[i]+G[i]/rho;
    for (int i=0;i<=m;i++){
        F[i]=0;G[i]=0;
    }
  /*  for (int i=1;i<=m;i++){
        if (leaf[P[i]]) F[i]+=T[P[i]];
        else {
            F[i]+=D[P[i]];
            G[i]+=E[P[i]];
        }
        if (leaf[i]) F[i]-=T[i];
        else {
            F[i]-=D[i];
            G[i]-=E[i];
        }
    }
    double a = 0;
    double b = 0;
    double c = 0; //a*rho*rho+b*rho+c=0
    for (int i=1;i<=m;i++){
        a+=F[i]*F[i]/V[i];
        b+=(2*F[i]*G[i]+B[i]*F[i])/V[i];
        c+=(G[i]*G[i]+B[i]*G[i])/V[i];
    }
    double delta = b*b-4*a*c;
    if (delta>=0){
        rho = (-b+sqrt(delta))/(2*a);
        if (rho<rho_min) rho=rho_min;
    }
    else {
        rho = rho_min;
    }*/
    
    
    for (int i=1;i<=m;i++){
        if (leaf[P[i]] && leaf[i]){
           F[i]=T[P[i]]-T[i];
           G[i]=0;
        }
        else if (leaf[P[i]]){
             F[i]=T[P[i]]-D[i];
             G[i]=-E[i];
        }
        else if (leaf[i]){
             F[i]=D[P[i]]-T[i];
             G[i]=E[P[i]];
        }        
        else {
            F[i]=D[P[i]]-D[i];
            G[i]=E[P[i]]-E[i];            
        }
    }
    double a = 0;
    double b = 0;
    double c = 0;//phi = a*rho^2+b*rho+c
    for (int i=1;i<=m;i++){
        a += F[i]*F[i]/V[i];
        b += 2*(B[i]+G[i])*F[i]/V[i];
        c += (B[i]+G[i])*(B[i]+G[i])/V[i];
    }
    rho = -b/(2*a);
    if (rho<rho_min) rho=rho_min;
    for (list<int>::iterator iter = pre.begin();iter!=pre.end();iter++){
        int i = *iter;
        if (flag[i] && !leaf[i]){
            T[i] = D[i]+E[i]/rho;
        }
        int s1 = Suc1[i];
        int s2 = Suc2[i];
        if (!flag[s1])
            T[s1]=T[i];
        if (!flag[s2])
            T[s2]=T[i];
    }
    phi1 = phi(n,m,P,B,V,T,rho);
    delete[] W;
    delete[] C;
    delete[] X;
    delete[] D;
    delete[] E;
    delete[] F;
    delete[] G;
    
    delete[] Pre;
    int *as=new int[m+1];
    for (int i=0;i<=m;i++) as[i]=-1;
    int count = 0;
    for (list<int>::iterator iter = active_set.begin();iter!=active_set.end();iter++){
        int a = *iter;
        if (a<=m){
            as[a]=count;
            count++;
        }
    }
    bool *bl=new bool[m+1];
    for (int i=0;i<=m;i++) bl[i]=false;
    double* lambda = new double[count];
    for (int y=0;y<ls.size();y++){
        while (!feuilles[y].empty()){
            int i=feuilles[y].top();
            if (!flag[i]){
                feuilles[y].pop();
                int p=P[i];
                int s1=Suc1[i];
                int s2=Suc2[i];
                if (i<n &&(flag[s1]||bl[s1]) && (flag[s2]||bl[s2])){
                        lambda[as[i]]=2*rho*(B[s1]-rho*T[s1]+rho*T[i])/V[s1]+2*rho*(B[s2]-rho*T[s2]+rho*T[i])/V[s2]-2*rho*(B[i]-rho*T[i]+rho*T[p])/V[i];
                    if (flag[s1] && flag[s2]) bl[i]=true;
                    if (!flag[s1] && bl[s1] && flag[s2]) {lambda[as[i]]+=lambda[as[s1]];bl[i]=true;}
                    if (!flag[s2] && bl[s2] && flag[s1]) {lambda[as[i]]+=lambda[as[s2]];bl[i]=true;}
                    if (!flag[s1] && bl[s1] && !flag[s2] && bl[s2]) {lambda[as[i]]+=lambda[as[s1]]+lambda[as[s2]];bl[i]=true;}
                }
                else{
                    int j=Suc1[p];
                    if (j==i)
                        j=Suc2[p];
                    lambda[as[i]]=-2*rho*(B[i]-rho*T[i]+rho*T[p])/V[i]-2*rho*(B[j]-rho*T[j]+rho*T[p])/V[j];
                    if (P[p]!=-1)
                        lambda[as[i]]+=2*rho*(B[p]-rho*T[p]+rho*T[P[p]])/V[p];
                    if (flag[j] && flag[p]) bl[i]=true;
                    if (flag[p] && !flag[j] && bl[j]) {lambda[as[i]]-=lambda[as[j]];bl[i]=true;}
                    if (!flag[p] && bl[p] && flag[j]) {lambda[as[i]]+=lambda[as[p]];bl[i]=true;}
                    if (!flag[j] && bl[j] && !flag[p] && bl[p]){lambda[as[i]]+=lambda[as[p]]-lambda[as[j]];bl[i]=true;}
                }
            }
        }
    }
    for (int y=0;y<top.size();y++){
        for (list<int>::iterator iter = internal[y].begin();iter!=internal[y].end();iter++){
            int i= *iter;
            int p = P[i];
            int s1=Suc1[i];
            int s2=Suc2[i];
            if (i<n &&(flag[s1]||bl[s1]) && (flag[s2]||bl[s2])){
                    lambda[as[i]]=2*rho*(B[s1]-rho*T[s1]+rho*T[i])/V[s1]+2*rho*(B[s2]-rho*T[s2]+rho*T[i])/V[s2]-2*rho*(B[i]-rho*T[i]+rho*T[p])/V[i];
                if (flag[s1] && flag[s2]) bl[i]=true;
                if (!flag[s1] && bl[s1] && flag[s2]) {lambda[as[i]]+=lambda[as[s1]];bl[i]=true;}
                if (!flag[s2] && bl[s2] && flag[s1]) {lambda[as[i]]+=lambda[as[s2]];bl[i]=true;}
                if (!flag[s1] && bl[s1] && !flag[s2] && bl[s2]) {lambda[as[i]]+=lambda[as[s1]]+lambda[as[s2]];bl[i]=true;}
            }
        }
    }
    for (int i=0;i<count;i++){
        if (myabs(lambda[i])<(10e-15))
            lambda[i]=0;
        ld.push_back(lambda[i]);
    }
    delete[] as;
    delete[] feuilles;
    delete[] internal;
    delete[] leaf;
    delete[] bl;
    delete[] lambda;
    return ld;
}

double with_constraint_active_set(int n,int m,int* & P,int* & Suc1,int* & Suc2,double* & B,double* & V,double* & T,double &rho,double rho_min){
    //this methods implements the QPD algorithm (active set method)
    list<int> active_set;
  //  cout<<"Computing starting feasible point ..."<<endl;
    double phi1 = starting_point(n,m,P,Suc1,Suc2,B,V,T,rho,rho_min,active_set);
    double *T_new=new double[m+1];
    double *dir=new double[n];
    double alpha;
    bool *flag=new bool[m+1];
    for (int i=n;i<=m;i++) T_new[i]=T[i];
    list<double> lambda = without_constraint_active_set(n,m,flag,P,Suc1,Suc2,B,V,T_new,rho,rho_min,active_set,phi1);
    int nb_iter=1;
    while (!conditions(m,lambda,P,T_new) && nb_iter<=maxIter){
     //   cout<<nb_iter<<": ";
        for (int i=0;i<n;i++)  dir[i]=T_new[i]-T[i];
        alpha=1;
        int as=0;
        double a;
        for (int i=n;i<=m;i++){
            if (flag[i]){
                int p = P[i];
                if (dir[p]>0){
                    a = (T[i]-T[p])/dir[p];
                    if (a<alpha){
                        alpha = a;
                        as = i;
                    }
                }
            }
        }
        for (int i=0;i<n;i++){
            if (flag[i]){
                int p = P[i];
                if (p!=-1){
                    if (dir[p]>dir[i]){
                        a = (T[i]-T[p])/(dir[p]-dir[i]);
                        if (a<alpha){
                            alpha = a;
                            as = i;
                        }
                    }
                }
            }
        }
        for (int i=0;i<n;i++) T[i]=T[i]+alpha*dir[i];
        int asrm= remove_ne_lambda(lambda,active_set);
        if (asrm!=-1) {
            active_set.remove(asrm);
        }
        if (as!=0) {
            active_set.push_back(as);
        }
        lambda = without_constraint_active_set(n,m,flag,P,Suc1,Suc2,B,V,T_new,rho,rho_min,active_set,phi1);        
        nb_iter++;
    }
    if (nb_iter<=maxIter) for (int i=0;i<n;i++) T[i]=T_new[i];
    else phi1=phi(n,m,rho,P,B,V,T);
    delete[] T_new;
    delete[] dir;
    delete[] flag;
    return phi1;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////this part implements QPD algorithm in the case that the rate is given/////////////////////////

double starting_point_rate(double rho,int n,int m,int* &P,int* & Suc1,int* & Suc2,double* & B,double* & V,double* & T,list<int> & active_set){
    without_constraint_rate(rho,n,m,P,Suc1,Suc2,B,V,T);
    list<int> pos = postorder(P,Suc1,Suc2,n);
    for (list<int>::iterator iter = pos.begin();iter!=pos.end();iter++){
        int i = *iter;
        int s1 = Suc1[i];
        int s2 = Suc2[i];
        if (T[i]>T[s1] || T[i]>T[s2]){
            if (T[s1]<T[s2]){
                T[i]=T[s1];
                active_set.push_back(s1);
            }
            else{
                T[i]=T[s2];
                active_set.push_back(s2);
            }
        }
    }
    return phi(n,m,P,B,V,T,rho);
}


list<double> without_constraint_active_set_rate(double rho,int n,int m,bool* & flag,int* &P,int* & Suc1,int* & Suc2,double* & B,double* & V,double* & T,list<int>& active_set,double &phi1){
    list<double> ld;
    
    // FIX FOR CLANG
    list<int>* SucL = new list<int>[n];
    list<int>* SucI = new list<int>[n];
    
    int *Pre=new int[m+1];
    for (int i=0;i<=m;i++) Pre[i]=-1;
    list<int> ls;
    for (int i=0;i<=m;i++) flag[i]=true;
    bool *leaf=new bool[m+1];
    for (int i=0;i<n;i++) leaf[i]=false;
    for (int i=n;i<=m;i++) leaf[i]=true;
    //flag
    for (list<int>::iterator iter=active_set.begin();iter!=active_set.end();iter++){
        int i = *iter;
        flag[i]=false;
        if (i>=n && i!=m+1 && i!=m+2) ls.push_back(i);
    }
    
    // stack<int> feuilles[ls.size()];
    stack<int>* feuilles = new stack<int>[ls.size()];
    
    list<int> top;
    //leaf
    computeFeuilles(n,m,P,Suc1,Suc2,T,flag,leaf,feuilles,active_set);
    for (int i=0;i<n;i++){
        int s1= Suc1[i];
        int s2= Suc2[i];
        if (flag[i] && !leaf[i] && (!flag[s1] || !flag[s2]))
            top.push_back(i);
    }
    
    // list<int> internal[top.size()];
    list<int>* internal = new list<int>[top.size()];
    
    reduceTree(n,m,P,Pre,Suc1,Suc2,SucL,SucI,flag,leaf,internal,active_set);
    
    list<int> pos = postorder(P,Suc1,Suc2,n);
    list<int> pre = preorder(P,Suc1,Suc2,n);
    double *W= new double[n];
    double *C = new double[n];
    double *X =new double[n];
    double *D= new double[n];
    double * E = new double[n];
    for (int i=0;i<n;i++){
        W[i]=0;C[i]=0;X[i]=0;D[i]=0;E[i]=0;
    }
    
    for (list<int>::iterator it = pos.begin();it!=pos.end();it++){
        int i = *it;
        if (flag[i] && !leaf[i]){
            int p = Pre[i];
            if (p==-1){
                double coefs=0;
                double xtemp=0;
                double ctemp=0;
                for (list<int>::iterator iter = SucL[i].begin();iter!=SucL[i].end();iter++){
                    coefs+=1/V[*iter];
                    ctemp+=T[*iter]/V[*iter];
                    xtemp-=B[*iter]/V[*iter];
                }
                for (list<int>::iterator iter = SucI[i].begin();iter!=SucI[i].end();iter++){
                    coefs+=(1-W[*iter])/V[*iter];
                    ctemp+=C[*iter]/V[*iter];
                    xtemp+=(X[*iter]-B[*iter])/V[*iter];
                }
                D[i]=ctemp/coefs;
                E[i]=xtemp/coefs;
            }
            else{
                if (!leaf[p]){
                    double coefs=1/V[i];
                    double wtemp=1/V[i];
                    double ctemp=0;
                    double xtemp=B[i]/V[i];
                    for (list<int>::iterator iter = SucL[i].begin();iter!=SucL[i].end();iter++){
                        coefs+=1/V[*iter];
                        ctemp+=T[*iter]/V[*iter];
                        xtemp-=B[*iter]/V[*iter];
                    }
                    for (list<int>::iterator iter = SucI[i].begin();iter!=SucI[i].end();iter++){
                        coefs+=(1-W[*iter])/V[*iter];
                        ctemp+=C[*iter]/V[*iter];
                        xtemp+=(X[*iter]-B[*iter])/V[*iter];
                    }
                    W[i]=wtemp/coefs;
                    C[i]=ctemp/coefs;
                    X[i]=xtemp/coefs;
                }
                else{
                    double coefs=1/V[i];
                    double dtemp=T[p]/V[i];
                    double etemp=B[i]/V[i];
                    for (list<int>::iterator iter = SucL[i].begin();iter!=SucL[i].end();iter++){
                        coefs+=1/V[*iter];
                        dtemp+=T[*iter]/V[*iter];
                        etemp-=B[*iter]/V[*iter];
                    }
                    for (list<int>::iterator iter = SucI[i].begin();iter!=SucI[i].end();iter++){
                        coefs+=(1-W[*iter])/V[*iter];
                        dtemp+=C[*iter]/V[*iter];
                        etemp+=(X[*iter]-B[*iter])/V[*iter];
                    }
                    D[i]=dtemp/coefs;
                    E[i]=etemp/coefs;
                }
            }
        }
    }
    for (list<int>::iterator iter = pre.begin();iter!=pre.end();iter++){
        int i = *iter;
        if (Pre[i]!=-1 && flag[i] && !leaf[i] && !leaf[Pre[i]]){
            D[i]=W[i]*D[Pre[i]]+C[i];
            E[i]=W[i]*E[Pre[i]]+X[i];
        }
        int s1 = Suc1[i];
        int s2 = Suc2[i];
        if (!flag[s1] && !leaf[s1]) {
            D[s1]=D[i];
            E[s1]=E[i];
        }
        if (!flag[s2] && !leaf[s2]) {
            D[s2]=D[i];
            E[s2]=E[i];
        }
    }
    
    for (list<int>::iterator iter = pre.begin();iter!=pre.end();iter++){
        int i = *iter;
        if (!leaf[i]){
            if (flag[i]){
                T[i] = D[i]+E[i]/rho;
            }
            int s1 = Suc1[i];
            int s2 = Suc2[i];
            if (!flag[s1]) T[s1]=T[i];
            if (!flag[s2]) T[s2]=T[i];
        }
    }
    delete[] W;
    delete[] C;
    delete[] X;
    delete[] D;
    delete[] E;
    phi1 = phi(n,m,P,B,V,T,rho);
    
    int *as=new int[m+1];
    for (int i=0;i<=m;i++) as[i]=-1;
    int count = 0;
    for (list<int>::iterator iter = active_set.begin();iter!=active_set.end();iter++){
        int a =*iter;
        if (a<=m){
            as[*iter]=count;
            count++;
        }
    }
    bool *bl=new bool[m+1];
    for (int i=0;i<=m;i++) bl[i]=false;
    
    double *lambda=new double[count];
    for (int y=0;y<ls.size();y++){
        while (!feuilles[y].empty()){
            int i=feuilles[y].top();
            if (!flag[i]){
                feuilles[y].pop();
                int p=P[i];
                int s1=Suc1[i];
                int s2=Suc2[i];
                if (i<n &&(flag[s1]||bl[s1]) && (flag[s2]||bl[s2])){
                        lambda[as[i]]=2*rho*(B[s1]-rho*T[s1]+rho*T[i])/V[s1]+2*rho*(B[s2]-rho*T[s2]+rho*T[i])/V[s2]-2*rho*(B[i]-rho*T[i]+rho*T[p])/V[i];
                    if (flag[s1] && flag[s2]) bl[i]=true;
                    if (!flag[s1] && bl[s1] && flag[s2]) {lambda[as[i]]+=lambda[as[s1]];bl[i]=true;}
                    if (!flag[s2] && bl[s2] && flag[s1]) {lambda[as[i]]+=lambda[as[s2]];bl[i]=true;}
                    if (!flag[s1] && bl[s1] && !flag[s2] && bl[s2]) {lambda[as[i]]+=lambda[as[s1]]+lambda[as[s2]];bl[i]=true;}
                }
                else{
                    int j=Suc1[p];
                    if (j==i)
                        j=Suc2[p];
                    lambda[as[i]]=-2*rho*(B[i]-rho*T[i]+rho*T[p])/V[i]-2*rho*(B[j]-rho*T[j]+rho*T[p])/V[j];
                    if (P[p]!=-1)
                        lambda[as[i]]+=2*rho*(B[p]-rho*T[p]+rho*T[P[p]])/V[p];
                    if (flag[j] && flag[p]) bl[i]=true;
                    if (flag[p] && !flag[j] && bl[j]) {lambda[as[i]]-=lambda[as[j]];bl[i]=true;}
                    if (!flag[p] && bl[p] && flag[j]) {lambda[as[i]]+=lambda[as[p]];bl[i]=true;}
                    if (!flag[j] && bl[j] && !flag[p] && bl[p]){lambda[as[i]]+=lambda[as[p]]-lambda[as[j]];bl[i]=true;}
                }
            }
        }
    }
    for (int y=0;y<top.size();y++){
        for (list<int>::iterator iter = internal[y].begin();iter!=internal[y].end();iter++){
            int i= *iter;
            int p = P[i];
            int s1=Suc1[i];
            int s2=Suc2[i];
            if (i<n &&(flag[s1]||bl[s1]) && (flag[s2]||bl[s2])){
                    lambda[as[i]]=2*rho*(B[s1]-rho*T[s1]+rho*T[i])/V[s1]+2*rho*(B[s2]-rho*T[s2]+rho*T[i])/V[s2]-2*rho*(B[i]-rho*T[i]+rho*T[p])/V[i];
                if (flag[s1] && flag[s2]) bl[i]=true;
                if (!flag[s1] && bl[s1] && flag[s2]) {lambda[as[i]]+=lambda[as[s1]];bl[i]=true;}
                if (!flag[s2] && bl[s2] && flag[s1]) {lambda[as[i]]+=lambda[as[s2]];bl[i]=true;}
                if (!flag[s1] && bl[s1] && !flag[s2] && bl[s2]) {lambda[as[i]]+=lambda[as[s1]]+lambda[as[s2]];bl[i]=true;}
            }
        }
    }
    
    for (int i=0;i<count;i++){
        if (myabs(lambda[i])<(10e-15))
            lambda[i]=0;
        ld.push_back(lambda[i]);
    }
    delete[] Pre;
    delete[] leaf;
    delete[] SucL;
    delete[] SucI;
    delete[] internal;
    delete[] feuilles;
    delete[] as;
    delete[] bl;
    delete[] lambda;
    return ld;
}

double with_constraint_active_set_rate(double rho,int n,int m,int* &P,int* & Suc1,int* & Suc2,double* & B,double* & V,double* & T){
    list<int> active_set;
    double phi1;
   // cout<<"Computing the starting feasible point"<<endl;
    phi1 = starting_point_rate(rho,n,m,P,Suc1,Suc2,B,V,T,active_set);
    double *T_new=new double[m+1];
    double *dir=new double[n];
    double alpha;
    for (int i=n;i<=m;i++) T_new[i]=T[i];
    bool *flag=new bool[m+1];
    list<double> lambda = without_constraint_active_set_rate(rho,n,m,flag,P,Suc1,Suc2,B,V,T_new,active_set,phi1);
    int nb_iter=1;
    while (!conditions(m,lambda,P,T_new) && nb_iter<=maxIter){
       // cout<<nb_iter<<": ";
        for (int i=0;i<n;i++)  dir[i]=T_new[i]-T[i];
        alpha=1;
        int as=0;
        double a;
        for (int i=n;i<=m;i++){
            if (flag[i]){
                int p = P[i];
                if (dir[p]>0){
                    a = (T[i]-T[p])/dir[p];
                    if (a<alpha){
                        alpha = a;
                        as = i;
                    }
                }
            }
        }
        for (int i=0;i<n;i++){
            if (flag[i]){
                int p = P[i];
                if (p!=-1){
                    if (dir[p]>dir[i]){
                        a = (T[i]-T[p])/(dir[p]-dir[i]);
                        if (a<alpha){
                            alpha = a;
                            as = i;
                        }
                    }
                }
            }
        }
        for (int i=0;i<n;i++) T[i]=T[i]+alpha*dir[i];
        int asrm= remove_ne_lambda(lambda,active_set);
        if (asrm!=-1) {
            active_set.remove(asrm);
        }
        if (as!=0) {
            active_set.push_back(as);
        }
        lambda = without_constraint_active_set_rate(rho,n,m,flag,P,Suc1,Suc2,B,V,T_new,active_set,phi1);
        nb_iter++;
    }
    if (nb_iter<=maxIter) for (int i=0;i<n;i++) T[i]=T_new[i];
    else phi1=phi(n,m,P,B,V,T,rho);
    delete[] T_new;
    delete[] flag;
    delete[] dir;
    return phi1;
}
