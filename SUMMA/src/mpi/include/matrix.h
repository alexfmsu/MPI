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

template<typename T>
void outer_product(int m, int s, T* a, T* b, T* c, int rank){
   for(int i = 0; i < m; i++){
      for(int j = 0; j < m; j++){
         c[i*m+j] += a[s+m*i] * b[j+s*m];  
      }
   }
}

template<typename T>
void SUMMA(Matrix<T>* A, Matrix<T>* B, Matrix<T>* C, int* grid_coords, MPI_Comm row_comm, MPI_Comm col_comm){
   MPI_Barrier(MPI_COMM_WORLD);

   double t0 = MPI_Wtime();   
   //------------------------------------------------------------------------------------------------------------------
   T* a = (A->block).data;
   T* b = (B->block).data;
   T* c = (C->block).data;
   
   int dm = (A->block).get_m();
   MPI_Datatype mpitype = A->get_mpitype();
   int rank = A->get_rank();
   int size = A->get_nprocs();

   int sq = (int)sqrt(size);

   for(int s = 0; s < sq * dm; s++){   
      MPI_Datatype row, col;

      MPI_Type_contiguous(dm, mpitype, &row);
      MPI_Type_commit(&row);
      
      MPI_Type_vector(dm, 1, dm, mpitype, &col);
      MPI_Type_commit(&col);

      MPI_Status status_a[size], status_b[size];
      MPI_Request req_s_a[size], req_s_b[size], req_r_a[size], req_r_b[size];
      //------------------------------------------------------------------------------------------------------------------
      int t = s % dm;

      if(rank==0 || rank == size-1)outer_product(dm, t, a, b, c, rank);
      
      MPI_Isend(&b[(t) * dm], 1, row, ( ((rank + sq)/sq >= sq) ? (rank + sq)%sq : (rank + sq) ), 1, MPI_COMM_WORLD, &req_s_b[rank]);
      MPI_Irecv(&b[(t) * dm], 1, row, ((rank - sq < 0) ? (size + rank - sq) : (rank - sq)), 1, MPI_COMM_WORLD, &req_r_b[rank]);
      
      MPI_Wait(&req_r_b[rank], &status_b[rank]);
      
      if(rank!=0 && rank != size-1)outer_product(dm, t, a, b, c, rank);
      
      MPI_Isend(&a[t], 1, col, rank/sq*sq+(rank+1)%sq, 0, MPI_COMM_WORLD, &req_s_a[rank]);
      MPI_Irecv(&a[t], 1, col, (rank%sq - 1 < 0) ? (rank+sq-1) : (rank - 1), 0, MPI_COMM_WORLD, &req_r_a[rank]);   
      
      MPI_Wait(&req_r_a[rank], &status_a[rank]);
      //------------------------------------------------------------------------------------------------------------------
   }
   //------------------------------
   MPI_Barrier(MPI_COMM_WORLD);

   double t1 = MPI_Wtime();

   if(A->get_rank() == 0)cout << "Multiply time: " << (t1 - t0) << endl << endl;
   //-----------------------------
}


