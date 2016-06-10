#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <mpi.h>

#include "src/matrix.cpp"

int readType(const char* filename){
   //---------------------------------------------------------------------------
   int rank, size;
   
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   MPI_Comm_size(MPI_COMM_WORLD, &size);
   //---------------------------------------------------------------------------
   FILE* f;
   
   if((f = fopen(filename, "rb")) == NULL){
      if(rank == 0){
         printf("\nError: Cannot open file '%s'\n\n", filename);
      }
       
      MPI_Finalize();
      exit(1);
   }
   
   fseek(f, sizeof(long)*2, SEEK_SET);
   
   int type = -1;   
   	
   if(fread(&type, sizeof(int), 1, f) != 1 || type<1 || type>4){
      if(rank == 0){
         printf("\nIncorrect value: 'type'\n\n");      
      }
      
      MPI_Finalize();
      exit(1);
   }
   
   fclose(f);
   
   return type;
}

void transpose(const char* file_in, const char* file_out){
   //---------------------------------------------------------------------------
   int rank, size;
   
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   MPI_Comm_size(MPI_COMM_WORLD, &size);
   //---------------------------------------------------------------------------
   int type = readType(file_in);
   
   switch(type){
      case 1:{
         if(rank == 0){
            cout << "file_in: " << "'" << file_in << "'" << endl;
         }
         
         Matrix<int> A(file_in, COLS); //A.print(COLS);
         
         if(rank == 0){
            cout << "m = " << A.get_m() << ", n = " << A.get_n() << "\n";
         }
         
         A.transpose(); //AT.print(COLS);
         A.writeToFile(file_out);
     		break;
      }
      case 2:{
         if(rank == 0){
            cout << "file_in: " << "'" << file_in << "'" << endl;
         }
         
         Matrix<long> A(file_in, COLS); //A.print(COLS);
         
         if(rank == 0){
            cout << "m = " << A.get_m() << ", n = " << A.get_n() << "\n";
         }
         
         A.transpose(); //AT.print(COLS);
         A.writeToFile(file_out);
     		
         break;
      }
      case 3:{
         if(rank == 0){
            cout << "file_in" << "'" << file_in << "'" << endl;
         }
         
         Matrix<float> A(file_in, COLS); //A.print(COLS);
         
         if(rank == 0){
            cout << "m = " << A.get_m() << ", n = " << A.get_n() << "\n";
         }
         
         A.transpose(); //AT.print(COLS);
         A.writeToFile(file_out);
     		
         break;
      }
      case 4:{
         if(rank == 0){
            cout << "file_in" << "'" << file_in << "'" << endl;
         }
         
         Matrix<double> A(file_in, COLS); //A.print(COLS);
         
         if(rank == 0){
            cout << "m = " << A.get_m() << ", n = " << A.get_n() << "\n";
         }
         
         A.transpose(); //AT.print(COLS);
         A.writeToFile(file_out);
     		
         break;
      }
      default:{
         break;
      }
   }
}

int main(int argc, char** argv){ 	
   //---------------------------------------------------------------------------
   MPI_Init(&argc, &argv);
   
   int rank,size;
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   MPI_Comm_size(MPI_COMM_WORLD, &size);
   //---------------------------------------------------------------------------
   if(argc != 3){
      if(rank == 0){
         printf("\nUsage: ./a.out <file_in> <file_out>\n\n");
      }
    
      MPI_Finalize();
      return 1;
   }
   //---------------------------------------------------------------------------
   if(rank == 0){
      cout << endl << "nproc: " << size << endl << endl;
   }
    
   transpose(argv[1], argv[2]);
   
   if(rank == 0){
      cout << endl;
   }
   //---------------------------------------------------------------------------
   MPI_Finalize();	
   
   return 0;
}
