#ifndef MATRIX_H
#define MATRIX_H 

#include <typeinfo>
#include <cstring>
#include <limits>
#include <cmath>
#include <mpi.h>
#include <iomanip>

#include "block.h"

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

      block.set_size((long)(m/sqrt(size)), (long)(m/sqrt(size)));
   }
   
   Matrix<T>(const Matrix<T> &tmp){
      (*this) = Matrix<T>(tmp.m, tmp.n);
      
      long k = 0;
   
      for(long i = 0; i < m; i++){
         for(long j = 0; j < n; j++){
            A[k] = tmp.A[k];
   
            k++;            
         }      
      }
   }
   //END---------------------------------------------------CONSTRUCTOR-------------------------------------------------
   
   void init_mpi();
   
   void alloc_memory();
   
   void empty();
   
   void distrib_dimensions(long& m, long& n, int& type, bool mpiroot);
   
   void free_T();
   
   //BEGIN-------------------------------------------------SET---------------------------------------------------------
   void set_type();
   //END---------------------------------------------------SET---------------------------------------------------------
   
   //BEGIN-------------------------------------------------GET---------------------------------------------------------
   long get_m() const{ return m; }
   long get_n() const{ return n; }
   
   int get_type() const{ return type; }
   MPI_Datatype get_mpitype() const{ return mpitype; }
   
   int get_nprocs() const{ return size; }
   int get_rank() const{ return rank; }
   
   bool get_mpiroot() const{ return mpiroot; }
   //END---------------------------------------------------GET---------------------------------------------------------
    
   //BEGIN-------------------------------------------------FILE READ/WRITE---------------------------------------------
   void read_from_file(const char* filename);
   void write_to_file(const char* filename);
   //END---------------------------------------------------FILE READ/WRITE---------------------------------------------
   
   //BEGIN-------------------------------------------------GATTER/GATHER----------------------------------------------
   void scatter(int*, MPI_Comm, MPI_Comm);
   void gather(int*, MPI_Comm, MPI_Comm);
   //END---------------------------------------------------SCATTER/GATHER----------------------------------------------
   
   //BEGIN-------------------------------------------------PRINT-------------------------------------------------------
   void print() const;
   //END---------------------------------------------------PRINT-------------------------------------------------------
   
   //BEGIN-------------------------------------------------DESTRUCTOR--------------------------------------------------
   ~Matrix<T>(){  
      free_T();
   }
   //END---------------------------------------------------DESTRUCTOR--------------------------------------------------
   
   T* A;
   Block<T> block;
   
private:
   //-------------------
   long m, n;
   long block_size;
   
   int type;
   MPI_Datatype mpitype;
   
   bool mem_alloced;
   //-------------------
   int rank;
   int size;
   
   bool mpiroot;
   //-------------------
};
#endif

int proc_pos(int row, int col, int sq){
   return ((row + sq) % sq)*sq + (col + sq) % sq;
}

template<typename T>
void cannon_multiply(Matrix<T>* A, Matrix<T>* B, Matrix<T>* C){
   //------------------------------------------------------------------------------------------------------------------
   MPI_Barrier(MPI_COMM_WORLD);
   
   double t0 = MPI_Wtime();
	//------------------------------------------------------------------------------------------------------------------
   T* a = (A->block).data;
   T* b = (B->block).data;
   T* c = (C->block).data;
   
   int dm = (A->block).get_m();
   
   MPI_Datatype mpitype = A->get_mpitype();
   //------------------------------------------------------------------------------------------------------------------
   int rank = A->get_rank();
   int size = A->get_nprocs();
   
   int sq = (int)sqrt(size);
   
   int myrow = rank / sq;
   int mycol = rank % sq;
   //------------------------------------------------------------------------------------------------------------------   
   MPI_Status status;
   
   T tmp[dm * dm];
   //------------------------------------------------------------------------------------------------------------------
   MPI_Sendrecv(a, dm * dm, mpitype, proc_pos(myrow, mycol - myrow, sq), 1, tmp, dm * dm, mpitype, 
   proc_pos(myrow, mycol + myrow, sq), 1, MPI_COMM_WORLD, &status);
   memcpy(a, tmp, dm * dm * sizeof(T));
   
   MPI_Sendrecv(b, dm * dm, mpitype, proc_pos(myrow - mycol, mycol, sq), 1, tmp, dm * dm, mpitype, 
   proc_pos(myrow + mycol, mycol, sq), 1, MPI_COMM_WORLD, &status); 
   memcpy(b, tmp, dm * dm * sizeof(T));
   //------------------------------------------------------------------------------------------------------------------
   for(int s = 0; s < sq; s++){
      for(int i = 0; i < dm; i++){
         for(int k = 0; k < dm; k++){
            for(int j = 0; j < dm; j++){
               c[i*dm + j] += a[i*dm + k] * b[k*dm + j];
            }
         }
      }
      
      // left shift
      MPI_Sendrecv(a, dm * dm, mpitype, proc_pos(myrow, mycol - 1, sq), 1, tmp, dm * dm, mpitype, 
      proc_pos(myrow, mycol + 1, sq), 1, MPI_COMM_WORLD, &status);
      memcpy(a, tmp, dm * dm * sizeof(T));
      
      // upper shift  
      MPI_Sendrecv(b, dm * dm, mpitype, proc_pos(myrow - 1, mycol, sq), 1, tmp, dm * dm, mpitype, 
      proc_pos(myrow + 1, mycol, sq), 1, MPI_COMM_WORLD, &status);
      memcpy(b, tmp, dm * dm * sizeof(T));
   }
   //------------------------------------------------------------------------------------------------------------------
   MPI_Sendrecv(a, dm * dm, mpitype, proc_pos(myrow, mycol + myrow, sq), 1, tmp, dm * dm, mpitype, 
   proc_pos(myrow, mycol - myrow, sq), 1, MPI_COMM_WORLD, &status); 
   memcpy(a, tmp, dm * dm * sizeof(T));
   
   MPI_Sendrecv(b, dm * dm, mpitype, proc_pos(myrow + mycol, mycol, sq), 1, tmp, dm * dm, mpitype, 
   proc_pos(myrow - mycol, mycol, sq), 1, MPI_COMM_WORLD, &status);
   memcpy(b, tmp, dm * dm * sizeof(T));
   //------------------------------------------------------------------------------------------------------------------
   MPI_Barrier(MPI_COMM_WORLD);
   
   double t1 = MPI_Wtime();
   //------------------------------------------------------------------------------------------------------------------	
   if(A->get_mpiroot())cout << "Multiply time: " << (t1 - t0) << endl << endl;
}


