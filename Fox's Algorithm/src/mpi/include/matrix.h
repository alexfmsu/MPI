#ifndef MATRIX_H
#define MATRIX_H 

#include <typeinfo>
#include <cstring>
#include <limits>
#include <cmath>
#include <mpi.h>
#include <iomanip>

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
   
   void distrib_dimensions(long& m, long& n, int& type, long& block_size, bool mpiroot);
   
   void free_T();
   
   //BEGIN-------------------------------------------------SET---------------------------------------------------------
   void set_block_size(long _block_size){ block_size = _block_size; }
   
   void set_type();
   //END---------------------------------------------------SET---------------------------------------------------------
   
   //BEGIN-------------------------------------------------GET---------------------------------------------------------
   long get_m() const{ return m; }
   long get_n() const{ return n; }
   long get_block_size() const{ return block_size; }
   
   int get_type() const{ return type; }
   MPI_Datatype get_mpitype() const{ return mpitype; }
   
   int get_nprocs() const{ return size; }

   bool get_mpiroot() const{ return mpiroot; }
   //END---------------------------------------------------GET---------------------------------------------------------
    
   //BEGIN-------------------------------------------------FILE READ/WRITE---------------------------------------------
   void read_from_file(const char* filename);
   void write_to_file(const char* filename);
   //END---------------------------------------------------FILE READ/WRITE---------------------------------------------
   
   //BEGIN-------------------------------------------------SCATTER/GATHER----------------------------------------------
   void scatter(const int* GridCoords,  T* pMatrixBlock, MPI_Comm row_comm, MPI_Comm col_comm);
   void gather(const int* GridCoords, MPI_Comm row_comm, MPI_Comm col_comm);
   //END---------------------------------------------------SCATTER/GATHER----------------------------------------------
   
   //BEGIN-------------------------------------------------PRINT-------------------------------------------------------
   void print() const;
   void print_nodes();
   //END---------------------------------------------------PRINT-------------------------------------------------------
   
   //BEGIN-------------------------------------------------DESTRUCTOR--------------------------------------------------
   ~Matrix<T>(){   
      free_T();
   }
   //END---------------------------------------------------DESTRUCTOR--------------------------------------------------
   
   T* A_block;
   
private:
   T* A;
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
void mult(Matrix<T>& A, Matrix<T>& B, Matrix<T>& C, int* GridCoords, T* pMatrixAblock, MPI_Comm row_comm, MPI_Comm col_comm){
   //------------------------------------------------------------------------------------------------------------------
   int np = A.get_nprocs();   
   int grid_size = (int)sqrt(np);
   
   long block_size = A.get_block_size();
   
   MPI_Datatype mpitype = A.get_mpitype();   
   //------------------------------------------------------------------------------------------------------------------
   C.set_block_size(block_size);
   
   C.A_block = new T[block_size * block_size];
   for(long i = 0; i < block_size * block_size; i++){
      C.A_block[i] = 0;
   }
   //------------------------------------------------------------------------------------------------------------------
   MPI_Barrier(MPI_COMM_WORLD);
 
   double t0 = MPI_Wtime();
   //------------------------------------------------------------------------------------------------------------------
   for(int i = 0; i < grid_size; i++) {
      int Pivot = (GridCoords[0] + i) % grid_size;
            
      if(GridCoords[1] == Pivot){ // копирование передаваемого блока в отдельный буфер
         for(long j = 0; j < block_size * block_size; j++){
            A.A_block[j] = pMatrixAblock[j];
         }
      }
      //---------------------------------------------------------------------------------------------------------------
      long I = 0;
      
      for(long i1 = 0; i1 < block_size; i1++){
          for(long k1 = 0; k1 < block_size; k1++){
            for(long j1 = 0; j1 < block_size; j1++){
               C.A_block[i1 * block_size + j1] += A.A_block[i1 * block_size + k1] * B.A_block[k1 * block_size + j1];
            }
         }
      }
      //---------------------------------------------------------------------------------------------------------------
      MPI_Status status;
      
      int prev_proc = GridCoords[0] - 1;
      int next_proc = GridCoords[0] + 1;
      
      if(GridCoords[0] == (grid_size - 1)){
          next_proc = 0;
      }
      
      if(GridCoords[0] == 0){
         prev_proc = grid_size - 1;
      }
       
      MPI_Sendrecv_replace(B.A_block, block_size * block_size, mpitype, next_proc, 0, prev_proc, 0, col_comm, &status); 
   }
   //------------------------------------------------------------------------------------------------------------------
   MPI_Barrier(MPI_COMM_WORLD);

   double t1 = MPI_Wtime();
   //------------------------------------------------------------------------------------------------------------------
   if(A.get_mpiroot()){ 
      cout << "Multiply time: " << (t1 - t0) << endl << endl; 
   }
   //------------------------------------------------------------------------------------------------------------------
   C.gather(GridCoords, row_comm, col_comm);
}


