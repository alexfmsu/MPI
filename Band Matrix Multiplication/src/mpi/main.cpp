#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <mpi.h>

#include "src/matrix.cpp"

using namespace std;

int rank, size;
bool mpiroot;

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
         if(mpiroot)cout << "Matrix A:" << endl << endl;
         Matrix<int> A;
         A.read_from_file(file1);
         if(mpiroot)cout << "-------------------------" << endl;
         
         if(mpiroot)cout << "Matrix B:" << endl << endl;
         Matrix<int> B;
         B.read_from_file(file2); 
         if(mpiroot)cout << "-------------------------" << endl;
         
         if(mpiroot)cout << "Matrix C:" << endl << endl;
         Matrix<int> C = A*B;
         C.write_to_file(file3);
          
         break;
      }
      case 2:{
         if(mpiroot)cout << "Matrix A:" << endl << endl;
         Matrix<long> A;
         A.read_from_file(file1);
         if(mpiroot)cout << "-------------------------" << endl;
         
         if(mpiroot)cout << "Matrix B:" << endl << endl;
         Matrix<long> B;
         B.read_from_file(file2);
         if(mpiroot)cout << "-------------------------" << endl;
         
         if(mpiroot)cout << "Matrix C:" << endl << endl;
         Matrix<long> C = A*B;   
         C.write_to_file(file3);
         
         break;
      }
      case 3:{
         if(mpiroot)cout << "Matrix A:" << endl << endl;
         Matrix<float> A;
         A.read_from_file(file1); 
         if(mpiroot)cout << "-------------------------" << endl;
         
         if(mpiroot)cout << "Matrix B:" << endl << endl;
         Matrix<float> B;
         B.read_from_file(file2); 
         if(mpiroot)cout << "-------------------------" << endl;
         
         if(mpiroot)cout << "Matrix C:" << endl << endl;
         Matrix<float> C = A*B;  
         C.write_to_file(file3);
         
         break;
      }
      case 4:{
         if(mpiroot)cout << "Matrix A:" << endl << endl;
         Matrix<double> A;
         A.read_from_file(file1); 
         if(mpiroot)cout << "-------------------------" << endl;
         
         if(mpiroot)cout << "Matrix B:" << endl << endl;
         Matrix<double> B;
         B.read_from_file(file2); 
         if(mpiroot)cout << "-------------------------" << endl;
         
         if(mpiroot)cout << "Matrix C:" << endl << endl;
         Matrix<double> C = A*B;  
         C.write_to_file(file3);
         
         break;
      }
      default:{
         break;
      }
   } 
}

int main(int argc, char** argv){
   MPI_Init(&argc, &argv);
   
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   MPI_Comm_size(MPI_COMM_WORLD, &size);
   
   mpiroot = (rank == 0);
   //------------------------------------------------------------------------------------------------------------------	
   MPI_Barrier(MPI_COMM_WORLD);	
   
   double t1 = MPI_Wtime(); 	
   //------------------------------------------------------------------------------------------------------------------	
   if(argc != 4){
      if(mpiroot){
         cout << endl;
         cout << "Usage: ./a.out <file1_in> <file2_in> <file3_out>" << endl;
         cout << endl;
      }
      
      MPI_Finalize();		
      
      return 1;
   }
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
