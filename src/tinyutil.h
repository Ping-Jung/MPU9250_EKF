#ifndef _TINY_UTIL_H
#define _TINY_UTIL_H
#endif

#define FTYPE float
#define PI 3.1415926

// variable of state: 4 q0 qx qy qz
// variable of observation 6 bax bay baz emx emy emz 

#include <stdint.h>
#include <stdio.h>
#include "Arduino.h"
#include <math.h>

FTYPE InvSqrt(FTYPE x);
/*matrix operation*/

void mat_multiply(FTYPE *A, FTYPE *B, FTYPE *C, int m, int n, int k); // m:row of A, n: column of A, k-column of B  (v)
void mat_inverse(FTYPE *A, FTYPE *B, int n);                        //symmetric matrix inverse (v by 4x4)


void mat_add(FTYPE *A, FTYPE *B, FTYPE *C, int k);      
void mat_sub(FTYPE *A, FTYPE *B, FTYPE *C, int k);       
void mat_zero(FTYPE *A, int k);                         
void mat_trans(FTYPE *A, FTYPE *C,int m, int n);     // m:row of A, n: column of A  
void mat_prod(FTYPE *A, FTYPE *B ,FTYPE *C, int s ,int t);   // single production, ex C=ax*q0, t size of b, s represent size of A 

void print_mat(FTYPE *A, int m, int n);   //

FTYPE InvSqrt(FTYPE x)
{
   uint32_t i = 0x5F1F1412 - (*(uint32_t*)&x >> 1);
   float tmp = *(float*)&i;
   return tmp * (1.69000231f - 0.714158168f * x * tmp * tmp);
}

// C=AB
void mat_multiply(FTYPE *A, FTYPE *B, FTYPE *C, int m, int n, int k){
   int i,j,r;

   for(i=0;i<m;i++){
    for(r=0;r<k;r++){
      FTYPE temp=0;
      for(j=0;j<n;j++){
       temp+=A[i*n+j]*B[j*k+r];
      }
      C[i*k+r]=temp;
    }  
   }   
}

void mat_inverse(FTYPE *A, FTYPE *B, int dim){         // use gaussian elimination
    int i,j,k;
    float temp;
    for(i=0;i<(dim*dim);i++){
        if(i%(dim+1)==0){
            B[i]=1;
        }
        else{
          B[i]=0;
        }
    }

    for(i=1;i<dim;i++){
      for(j=i;j<dim;j++){
        if(A[j*dim+(i-1)]!=0){
           temp=-1*A[j*dim+(i-1)]/A[(i-1)*dim+(i-1)];
           for(k=0;k<dim;k++){
            A[j*dim+k]+=temp*A[(i-1)*dim+k];
           }
           for(k=i;k>0;k--){
            B[j*dim+(k-1)]+=temp*B[(i-1)*dim+(k-1)];
           }
        }
      }
    }

    for(i=(dim-1);i>0;i--){
      for(j=i;j>0;j--){
        if(A[(j-1)*dim+i]!=0){
          temp=(-1*A[(j-1)*dim+i]/A[i*dim+i]);
          for(k=0;k<dim;k++){
            B[(j-1)*dim+k]+=temp*B[i*dim+k];
          }
        }
      }
    }
   
    for(i=0;i<dim;i++){
      temp=A[i*dim+i];
      for(j=0;j<dim;j++){
       B[i*dim+j]/=temp;
      }
    }
}

// C=A+B
void mat_add(FTYPE *A, FTYPE *B, FTYPE *C, int k){
  for(int i=0;i<k;i++){
    C[i]=A[i]+B[i];
  }
}

//C=A-B
void mat_zero(FTYPE *A, int k){
  for(int i=0;i<k;i++){
    A[i]=0;
  }
}

void mat_sub(FTYPE *A, FTYPE *B, FTYPE *C, int k){
  for(int i=0;i<k;i++){
    C[i]=A[i]-B[i];
  }
}

// C=A^T
void mat_trans(FTYPE *A, FTYPE *C,int m, int n){
   int i,j;
   for(i=0;i<m;i++){
    for(j=0;j<n;j++){
      C[j*m+i]=A[i*n+j];
    }
   }
}

void print_mat(FTYPE *A, int m, int n){
  int i,j;
  for (i=0;i<m;i++){
    for(j=0;j<n;j++){
      Serial.print(A[i*n+j]);
      Serial.print(" ");
    }
    Serial.print("\n");
  }

}

void mat_prod(FTYPE *A, FTYPE *B ,FTYPE *C, int s ,int t){

  for(int i=0;i<t;i++){
    for(int j=0;j<s; j++){
      C[i*s+j]=A[j]*B[i];
   }   
  } 
}

void mat_copy(FTYPE *A, FTYPE *B ,int size){
  for(int i=0;i<size;i++){
    B[i]=A[i];
  }
}

void norm(FTYPE *A){
  FTYPE magni=sqrt(A[0]*A[0]+A[1]*A[1]+A[2]*A[2]);
  A[0]=A[0]/magni;
  A[1]=A[1]/magni;
  A[2]=A[2]/magni;
}

