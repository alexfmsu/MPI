#ifndef MATRIX_H
#define MATRIX_H 

#include <typeinfo>
#include <cstring>
#include <limits>
#include <cmath>
#include <mpi.h>
#include <iomanip>

#include "../src/local_matrix.cpp"

using namespace std;

template <typename T> 
class Matrix{
public:
   //BEGIN-------------------------------------------------CONSTRUCTOR-------------------------------------------------
   Matrix<T>(void){
      m = -1;
      n = -1;
      
      type = -1;
      
      init_mpi();
   }   
   
   Matrix<T>(long _m, long _n){
      init_mpi();
      
      m = _m;
      n = _n;
      
      set_type();      
      
      alloc_memory();
      empty();
   }
   
   Matrix<T>(const Matrix<T> &tmp){
      (*this) = Matrix<T>(tmp.m, tmp.n);
      
      for(long i = 0; i < m; i++){
         for(long j = 0; j < n; j++){
            A[i][j] = tmp.A[i][j];            
         }      
      }
      
      loc = tmp.loc;
   }
   //END---------------------------------------------------CONSTRUCTOR-------------------------------------------------
   
   void init_mpi();
   
   void alloc_memory();
   
   void empty();
   
   void distrib_dimensions(long& m, long& n, int& type, bool mpiroot);
   
   void free_T();
   
   Matrix<T> operator*(Matrix<T>& tmp);  
   
   //BEGIN-------------------------------------------------SET---------------------------------------------------------
   void set_type();
   //END---------------------------------------------------SET---------------------------------------------------------
   
   //BEGIN-------------------------------------------------GET---------------------------------------------------------
   long get_m() const{ return m; }
   long get_n() const{ return n; }
   
   int get_type() const{ return type; }
   //END---------------------------------------------------GET---------------------------------------------------------
    
   //BEGIN-------------------------------------------------FILE READ/WRITE---------------------------------------------
   void read_from_file(const char* filename);
   
   void write_to_file(const char* filename);
   //END---------------------------------------------------FILE READ/WRITE---------------------------------------------
   
   //BEGIN-------------------------------------------------SCATTER/GATHER----------------------------------------------
   void scatter();
   void gather();
   //END---------------------------------------------------SCATTER/GATHER----------------------------------------------
   
   //BEGIN-------------------------------------------------PRINT-------------------------------------------------------
   void print() const;
   void print_nodes();
   //END---------------------------------------------------PRINT-------------------------------------------------------
   
   //BEGIN-------------------------------------------------DESTRUCTOR--------------------------------------------------
   ~Matrix<T>(){   
      free_T();
      
      loc.free_data();
   }
   //END---------------------------------------------------DESTRUCTOR--------------------------------------------------
   
   LocalMatrix<T> loc;	
   
private:
   T** A;
   //-----------
   long m, n;
   
   int type;
   MPI_Datatype mpitype;
   
   bool mem_alloced;
   //-----------
   int rank;
   int size;
   
   bool mpiroot;
   //-----------
};
#endif
