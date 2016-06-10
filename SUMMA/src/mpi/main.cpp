#include <iostream>
#include <cstdlib>
#include <cmath>
#include "mpi.h"

#include "src/matrix.cpp"

int rank, size;
bool mpiroot;

int grid_size;      	// размер решетки процессов
int grid_coords[2];  // координаты процесса в решетке

MPI_Comm grid_comm; //коммуникатор в виде квадратной решетки
MPI_Comm col_comm;  //коммуникатор – столбец решетки
MPI_Comm row_comm;  //коммуникатор – строка решетки



void create_grid_comms() {
   int dim_size[2]; // количество процессов в каждом измерении решетки
   int periodic[2];	    
   int sub_dims[2];  	
  
   dim_size[0] = grid_size; 
   dim_size[1] = grid_size;
  
   periodic[0] = 0;
   periodic[1] = 0;
   
   MPI_Cart_create(MPI_COMM_WORLD, 2, dim_size, periodic, 1, &grid_comm); //коммуникатор в виде квадратной решетки 
   
   MPI_Cart_coords(grid_comm, rank, 2, grid_coords); //координаты процесса в решетке 
   
   //коммуникаторы для строк решетки
   sub_dims[0] = 0;
   sub_dims[1] = 1;
   
   MPI_Cart_sub(grid_comm, sub_dims, &row_comm);
   //------------------------------------------------------------------------------------------------------------------
   //коммуникаторы для столбцов решетки
   sub_dims[0] = 1;
   sub_dims[1] = 0;
   
   MPI_Cart_sub(grid_comm, sub_dims, &col_comm);
   //------------------------------------------------------------------------------------------------------------------
}

int readType(const char* filename){
    int type = -1;   
   //------------------------------------------------------------------------------------------------------------------
   FILE* f;
   
   if((f = fopen(filename, "rb")) == NULL){
      if(mpiroot){
         cout << endl;
         cout << "Error: Cannot open file '" << filename << "'" << endl;
         cout << endl;	
      }
      			
      MPI_Finalize();
      
      exit(1);
   }
   //------------------------------------------------------------------------------------------------------------------
   fseek(f, sizeof(long) * 2, SEEK_SET);
   
   if(fread(&type, sizeof(int), 1, f) != 1){
      MPI_Abort(MPI_COMM_WORLD, 1);	
      
      exit(1);	
   }
   
   if(type < 1){
      if(mpiroot){
         cout << endl;
         cout << "Incorrect value: 'type'" << filename << "'" << endl;
         cout << endl;	
      }
      
      MPI_Finalize();
      
      exit(1);
   }
   
   fclose(f);
   //------------------------------------------------------------------------------------------------------------------
   return type;
}

void multiply(const char* file1, const char* file2, const char* file3){
   int type1 = readType(file1);
   int type2 = readType(file2);
   //------------------------------------------------------------------------------------------------------------------   
   if(type1 != type2){
      if(mpiroot){
         cout << endl;
         cout << "Type mismatch" << endl;
         cout << endl;
      }
      
      MPI_Finalize();
      
      exit(1);
   }
   
   if(mpiroot){
      cout << "np: " << size << endl;
      cout << "-------------------------" << endl;
   }	
   //------------------------------------------------------------------------------------------------------------------
   switch(type1){
      case 1:{
         //------------------------------------------------------------------------------------------------------------         
         if(mpiroot)cout << "Matrix A:" << endl << endl;
         Matrix<int> A;
         A.read_from_file(file1);
         A.scatter(grid_coords, row_comm, col_comm);
         if(mpiroot)cout << "-------------------------" << endl;
         //------------------------------------------------------------------------------------------------------------         
         if(mpiroot)cout << "Matrix B:" << endl << endl;
         Matrix<int> B;
         B.read_from_file(file2);
         B.scatter(grid_coords, row_comm, col_comm);
         if(mpiroot)cout << "-------------------------" << endl;
         //------------------------------------------------------------------------------------------------------------
         if(mpiroot)cout << "Matrix C:" << endl << endl;
         Matrix<int> C(A.get_m(), B.get_n());
         
         int sq = (int)sqrt(size);
         
         //cout << "send: from " << rank << " to " << ( ((rank + sq)/sq >= sq) ? (rank + sq)%sq : (rank + sq) ) << endl;
         //cout << "recv:  " << rank << " from " << ((rank - sq < 0) ? (size + rank - sq) : (rank - sq)) << endl;

         SUMMA(&A, &B, &C,grid_coords, row_comm, col_comm);
      
         C.gather(grid_coords, row_comm, col_comm);
         
         //C.print();
         C.write_to_file(file3);
         //------------------------------------------------------------------------------------------------------------         
         break;
      }
      case 2:{
         //------------------------------------------------------------------------------------------------------------         
         if(mpiroot)cout << "Matrix A:" << endl << endl;
         Matrix<long> A;
         A.read_from_file(file1);
         A.scatter(grid_coords, row_comm, col_comm);
         if(mpiroot)cout << "-------------------------" << endl;
         //------------------------------------------------------------------------------------------------------------         
         if(mpiroot)cout << "Matrix B:" << endl << endl;
         Matrix<long> B;
         B.read_from_file(file2);
         B.scatter(grid_coords, row_comm, col_comm);
         if(mpiroot)cout << "-------------------------" << endl;
         //------------------------------------------------------------------------------------------------------------
         if(mpiroot)cout << "Matrix C:" << endl << endl;
         Matrix<long> C(A.get_m(), B.get_n());
         
         int sq = (int)sqrt(size);
         
         //cout << "send: from " << rank << " to " << ( ((rank + sq)/sq >= sq) ? (rank + sq)%sq : (rank + sq) ) << endl;
         //cout << "recv:  " << rank << " from " << ((rank - sq < 0) ? (size + rank - sq) : (rank - sq)) << endl;

         SUMMA(&A, &B, &C,grid_coords, row_comm, col_comm);
      
         C.gather(grid_coords, row_comm, col_comm);
         
         //C.print();
         C.write_to_file(file3);
         //------------------------------------------------------------------------------------------------------------         
         break;
      }
      case 3:{
         //------------------------------------------------------------------------------------------------------------         
         if(mpiroot)cout << "Matrix A:" << endl << endl;
         Matrix<float> A;
         A.read_from_file(file1);
         A.scatter(grid_coords, row_comm, col_comm);
         if(mpiroot)cout << "-------------------------" << endl;
         //------------------------------------------------------------------------------------------------------------         
         if(mpiroot)cout << "Matrix B:" << endl << endl;
         Matrix<float> B;
         B.read_from_file(file2);
         B.scatter(grid_coords, row_comm, col_comm);
         if(mpiroot)cout << "-------------------------" << endl;
         //------------------------------------------------------------------------------------------------------------
         if(mpiroot)cout << "Matrix C:" << endl << endl;
         Matrix<float> C(A.get_m(), B.get_n());
         
         int sq = (int)sqrt(size);
         
         //cout << "send: from " << rank << " to " << ( ((rank + sq)/sq >= sq) ? (rank + sq)%sq : (rank + sq) ) << endl;
         //cout << "recv:  " << rank << " from " << ((rank - sq < 0) ? (size + rank - sq) : (rank - sq)) << endl;

         SUMMA(&A, &B, &C,grid_coords, row_comm, col_comm);
      
         C.gather(grid_coords, row_comm, col_comm);
         
         //C.print();
         C.write_to_file(file3);
         //------------------------------------------------------------------------------------------------------------         
         break;
      }
      case 4:{
         //------------------------------------------------------------------------------------------------------------         
         if(mpiroot)cout << "Matrix A:" << endl << endl;
         Matrix<double> A;
         A.read_from_file(file1);
         A.scatter(grid_coords, row_comm, col_comm);
         if(mpiroot)cout << "-------------------------" << endl;
         //------------------------------------------------------------------------------------------------------------         
         if(mpiroot)cout << "Matrix B:" << endl << endl;
         Matrix<double> B;
         B.read_from_file(file2);
         B.scatter(grid_coords, row_comm, col_comm);
         if(mpiroot)cout << "-------------------------" << endl;
         //------------------------------------------------------------------------------------------------------------
         if(mpiroot)cout << "Matrix C:" << endl << endl;
         Matrix<double> C(A.get_m(), B.get_n());
         
         int sq = (int)sqrt(size);
         
         //cout << "send: from " << rank << " to " << ( ((rank + sq)/sq >= sq) ? (rank + sq)%sq : (rank + sq) ) << endl;
         //cout << "recv:  " << rank << " from " << ((rank - sq < 0) ? (size + rank - sq) : (rank - sq)) << endl;
         
         SUMMA(&A, &B, &C,grid_coords, row_comm, col_comm);
         
         C.gather(grid_coords, row_comm, col_comm);
         
         //C.print();
         C.write_to_file(file3);
         //------------------------------------------------------------------------------------------------------------         
         break;
      }
      default:{
         break;
      }
   }
}

using namespace std;

int main(int argc, char** argv){
   //------------------------------------------------------------------------------------------------------------------
   MPI_Init(&argc, &argv);
   
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   MPI_Comm_size(MPI_COMM_WORLD, &size);
   
   mpiroot = (rank == 0);
   
   grid_size = (int)sqrt(size);
   //------------------------------------------------------------------------------------------------------------------
   if(size != (grid_size * grid_size)) {
      if(mpiroot == 0) {
         printf ("Error: 'np' must be a perfect square\n");
      }
      
      MPI_Finalize();		
      
      return 1;
   }else{
      create_grid_comms();
   }
   //------------------------------------------------------------------------------------------------------------------
   MPI_Barrier(MPI_COMM_WORLD);	
   
   double t1 = MPI_Wtime(); 	
   //------------------------------------------------------------------------------------------------------------------
   multiply(argv[1], argv[2], argv[3]);
   //------------------------------------------------------------------------------------------------------------------
   MPI_Barrier(MPI_COMM_WORLD);	
   
   double t2 = MPI_Wtime(); 	
   //------------------------------------------------------------------------------------------------------------------
   if(mpiroot){
      cout << "-------------------------" << endl;
      cout << "Total time: " << (t2 - t1) << endl;
   }
   //------------------------------------------------------------------------------------------------------------------
   
   MPI_Finalize();
   
   //------------------------------------------------------------------------------------------------------------------
   
   return 0;
}

