#ifndef MATRIX_H
#define MATRIX_H 

#include <typeinfo>
#include <cstring>
#include <limits>
#include <cmath>
#include <vector>
#include <mpi.h>

enum strip_type{ROWS,COLS,NONE};

using namespace std;

template<typename T> class Matrix{
public:
   //BEGIN------------------------CONSTRUCTOR----------------------------------
   Matrix<T>(long _m, long _n, long _dm = 0, long _dn = 0){
      init_mpi();
      				
      m = _m;
      n = _n;

      dm = _dm;
      dn = _dn;
      
      set_type();      
      
      allocate();      

      set_zero();
   }
   
   Matrix<T>(const Matrix<T> &tmp){
      init_mpi();
      		
      m = tmp.m;
      n = tmp.n;
      
      //dm = tmp.dm;
      //dn = tmp.dn;
           
      set_type();      
      	
      long _m, _n;

      if(dm && dn){
         _m = dm;
         _n = dn;
      }else{
         _m = m;
         _n = n;
      }

      allocate();
      
      std::copy(tmp.A.begin(),tmp.A.end(),A.begin());
   }
   //END--------------------------CONSTRUCTOR----------------------------------
   
   //BEGIN------------------------DESTRUCTOR-----------------------------------
   ~Matrix(){    
      A.clear();      
   }
   //END--------------------------DESTRUCTOR-----------------------------------
   
   //BEGIN------------------------INIT-MPI-------------------------------------
   void init_mpi(){
      MPI_Comm_rank(MPI_COMM_WORLD, &rank);
      MPI_Comm_size(MPI_COMM_WORLD, &size);
      
      int sq = (int)sqrt(size);
		nrows = (int)sqrt(size);
		while(size % nrows)
			nrows--;
		
		ncols = size/nrows;
   }   
   //END--------------------------INIT-MPI-------------------------------------
  
   //BEGIN------------------------ALLOCATE-------------------------------------
   void allocate(){
      long _m, _n;
   
      if(dm && dn){
         _m = dm;
         _n = dn;
      }else{
         _m = m;
         _n = n;
      }

      A.resize(_m*_n);
   }		
   //END--------------------------ALLOCATE-------------------------------------
   
   //BEGIN------------------------SET-ZERO-------------------------------------
   void set_zero(){
      long _m, _n;

      if(dm && dn){
         _m = dm;
         _n = dn;
      }else{
         _m = m;
         _n = n;
      }

      std::fill(A.begin(),A.end(),0);  
   } 
   //END--------------------------SET-ZERO-------------------------------------
   
   //BEGIN------------------------FILE-HANDLERS--------------------------------
   void read_file_header(const char* filename);
   
   Matrix(const char* filename, strip_type);
   //END--------------------------FILE-HANDLERS--------------------------------
   
   //BEGIN------------------------GET------------------------------------------
   long get_m() const{ return m; }
   long get_n() const{ return n; }
   
   long get_dm() const{ return dm; }
   long get_dn() const{ return dn; }
     	
   int get_type() const{ return type; }
   //END--------------------------GET------------------------------------------
    
   //BEGIN------------------------SET------------------------------------------
   void set_type();
   //END--------------------------SET------------------------------------------
   	
   Matrix<T> operator*(Matrix<T>& tmp);  
   
   //BEGIN------------------------PRINT----------------------------------------
   void print(strip_type);
   //END--------------------------PRINT----------------------------------------
   	
   //BEGIN------------------------WRITE-TO-FILE--------------------------------
   void writeToFile(const char* filename);
   //END--------------------------WRITE-TO-FILE--------------------------------
   
   //BEGIN------------------------TRANSPOSE------------------------------------
   void local_transpose(T *a, int n);
   void transpose();
   //END--------------------------TRANSPOSE------------------------------------
private:
   std::vector<T> A;
   //----------
   long m, n;
   long dm, dn;
   			 
   int type;
   //----------
   int rank;
   int size;

	int nrows;
	int ncols;

   MPI_Datatype mpi_type;
};
#endif
