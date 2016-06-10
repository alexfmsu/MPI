#include "../include/matrix.h"

//BEGIN----------------------------------------------------INIT-MPI----------------------------------------------------
template <typename T>
void Matrix<T>::init_mpi(){
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   MPI_Comm_size(MPI_COMM_WORLD, &size);
   
   mpiroot = (rank == 0);
}   
//END------------------------------------------------------INIT-MPI----------------------------------------------------

//BEGIN----------------------------------------------------ALLOC-MEMORY------------------------------------------------
template <typename T>
void Matrix<T>::alloc_memory(){
   if(mpiroot){
      A = new T[m*n];
      
      empty();
      
      mem_alloced = true;
   }
}		
//END------------------------------------------------------ALLOC-MEMORY------------------------------------------------
   
//BEGIN----------------------------------------------------EMPTY-------------------------------------------------------
template <typename T>
void Matrix<T>::empty(){
   if(mpiroot){
      long k = 0;

      for(long i = 0; i < m; i++){
         for(long j = 0; j < n; j++){
            A[k++] = 0;         
         }
      }
   }
} 
//END------------------------------------------------------EMPTY-------------------------------------------------------

//BEGIN----------------------------------------------------DISTRIB-DIMENSIONS------------------------------------------
template <typename T>
void Matrix<T>::distrib_dimensions(long& m, long& n, int& type, bool mpiroot){
   int dimensions[3];
   
   if (mpiroot){
      dimensions[0] = (int)m;
      dimensions[1] = (int)n;
      dimensions[2] = type;
   }
   
   MPI_Bcast(dimensions, 3, MPI_INT, 0, MPI_COMM_WORLD);
   
   m = (long)dimensions[0];
   n = (long)dimensions[1];
   type = dimensions[2];
   
   MPI_Barrier(MPI_COMM_WORLD);
}
//END------------------------------------------------------DISTRIB-DIMENSIONS------------------------------------------

//BEGIN----------------------------------------------------FREE-T------------------------------------------------------
template <typename T>
void Matrix<T>::free_T(){
   if(mpiroot){
      if(mem_alloced){
         delete[] A;
         
         A = NULL;
         
         mem_alloced = false;
      }
   }
} 
//END------------------------------------------------------FREE-T------------------------------------------------------

//BEGIN----------------------------------------------------SET---------------------------------------------------------
template <typename T>
void Matrix<T>::set_type(){
   char s[3];
   
   strcpy(s, typeid(T).name());
   
   if(strcmp(s, "i") == 0){
      type = 1;
      mpitype = MPI_INT;
   }else if(strcmp(s, "l") == 0){
      type = 2;
      mpitype = MPI_LONG;
   }else if(strcmp(s, "f") == 0){
      type = 3;
      mpitype = MPI_FLOAT;
   }else if(strcmp(s, "d") == 0){
      type = 4;
      mpitype = MPI_DOUBLE;
   }else type = -1;
}
//END------------------------------------------------------SET---------------------------------------------------------

//BEGIN----------------------------------------------------SCATTER/GATHER----------------------------------------------
template <typename T>
void Matrix<T>::scatter(int* GridCoords, MPI_Comm row_comm, MPI_Comm col_comm){
   int block_size = block.get_m();

   T* pMatrixRow = new T[block_size * m];
   //------------------------------------------------------------------------------------------------------------------	
   MPI_Barrier(col_comm);
   
   if(GridCoords[1] == 0){
      MPI_Scatter(A, block_size * m, mpitype, pMatrixRow, block_size * m, mpitype, 0, col_comm);
   }
   
   MPI_Barrier(col_comm);
   //------------------------------------------------------------------------------------------------------------------
   MPI_Barrier(row_comm);   

   for(long i = 0; i < block_size; i++){
      MPI_Scatter(&pMatrixRow[i * m], block_size, mpitype, &(block.data[i * block_size]), block_size, mpitype, 0, row_comm);
   }

   MPI_Barrier(row_comm);   

   delete[] pMatrixRow;
}

template <typename T>
void Matrix<T>::gather(int* GridCoords, MPI_Comm row_comm, MPI_Comm col_comm){
   int block_size = block.get_m();
   
   T* pResultRow = new T[m * block_size];
   //----------------------------------------------------------------------------------------------------------------- 
   MPI_Barrier(row_comm);   
   
   for(long i = 0; i < block_size; i++){
      MPI_Gather(&block.data[i * block_size], block_size, mpitype, &pResultRow[i * m], block_size, mpitype, 0, row_comm);
   }
   //-----------------------------------------------------------------------------------------------------------------
   MPI_Barrier(col_comm);
   
   if(GridCoords[1] == 0){
      MPI_Gather(pResultRow, block_size * m, mpitype, A, block_size * m, mpitype, 0, col_comm);
   }
   
   MPI_Barrier(col_comm);   
   //-----------------------------------------------------------------------------------------------------------------
   delete[] pResultRow;
}
//END------------------------------------------------------SCATTER/GATHER----------------------------------------------

//BEGIN----------------------------------------------------READ-FROM-FILE----------------------------------------------
template <typename T>
void Matrix<T>::read_from_file(const char* filename){
   MPI_Barrier(MPI_COMM_WORLD);
   
   double time1 = MPI_Wtime();      
   //------------------------------------------------------------------------------------------------------------------
   FILE* f;   
   
   if(mpiroot){
      if((f = fopen(filename, "rb")) == NULL){
         cout << endl;
         cout << "Error: Cannot open file " << filename << endl;
         cout << endl;	
         
         MPI_Abort(MPI_COMM_WORLD, 1);
         
         exit(1);	
      }
      //---------------------------------------------------------------------------------------------------------------
      if(fread(&m, sizeof(long), 1, f) != 1){
         MPI_Abort(MPI_COMM_WORLD, 1);	
         
         exit(1);
      }
      
      if(m < 1){
         cout << endl;
         cout << "Error: parameter 'm' not found" << endl;
         cout << endl;
         
         MPI_Abort(MPI_COMM_WORLD, 1);	
         
         exit(1);	
      }
      
      if(m % size){
         cout << endl;
         cout << "Error: Matrix's size must be divisible by 'np'" << endl;
         cout << endl;	
         
         MPI_Abort(MPI_COMM_WORLD, 1);	
         exit(1);	
      }
      //---------------------------------------------------------------------------------------------------------------
      if(fread(&n, sizeof(long), 1, f) != 1){
         MPI_Abort(MPI_COMM_WORLD, 1);	
         
         exit(1);
      }
      
      if(n < 1){
         cout << endl;
         cout << "Error: parameter 'n' not found" << endl;
         cout << endl;
         
         MPI_Abort(MPI_COMM_WORLD, 1);	
         
         exit(1);	
      }
      
      block_size = (long)(m/sqrt(size));
   }
  
   MPI_Barrier(MPI_COMM_WORLD);  
 
   distrib_dimensions(m, n, type, mpiroot);
   
   set_type();	
   
   alloc_memory();
   
   if(mpiroot){
      fseek(f, sizeof(long) * 2 + sizeof(int), SEEK_SET);
      
      long k = 0;      

      for(long i = 0; i < m; i++){
         for(long j = 0; j < n; j++){
            if(fread(&A[k], sizeof(T), 1, f) != 1){
               MPI_Abort(MPI_COMM_WORLD, 1);	
               
               exit(1);	
            }      
            
            if(sizeof(A[k]) != sizeof(T)){
               if(rank == 0){
                  cout << endl;
                  cout << "Error: Cannot read matrix from file " << filename << endl;
                  cout << endl;
               }
               
               MPI_Abort(MPI_COMM_WORLD, 1);
               
               exit(1);	
            }         

            k++;
         }      
      }
      
      fclose(f);			
   }
   //------------------------------------------------------------------------------------------------------------------
   MPI_Barrier(MPI_COMM_WORLD);
   
   double time2 = MPI_Wtime();      
   //------------------------------------------------------------------------------------------------------------------
   if(mpiroot){
      cout << "Read time: " << (time2 - time1) << endl;
      cout << "type: 'int'" << endl;
      cout << "size: " << get_m() << endl;
      cout << endl;
   }
   //------------------------------------------------------------------------------------------------------------------	
   block.set_size((long)(m/sqrt(size)), (long)(n/sqrt(size)));
}
//END------------------------------------------------------READ-FROM-FILE----------------------------------------------

//BEGIN----------------------------------------------------PRINT-------------------------------------------------------
template <typename T>
void Matrix<T>::print() const{
   using namespace std;

   cout.precision(10);   

   if(rank == 0) {
      for(int i = 0; i < m; i++){
         for(int j = 0; j < n; j++) {
            cout << (T)A[i*m+j] << " ";
         }

         cout << endl;                                                                                                                       
      }

		cout << endl;
   }
}
//END------------------------------------------------------PRINT-------------------------------------------------------

//BEGIN----------------------------------------------------WRITE-TO-FILE-----------------------------------------------
template <typename T>
void Matrix<T>::write_to_file(const char* filename){
   if(mpiroot){
      double time1 = MPI_Wtime();      
      //---------------------------------------------------------------------------------------------------------------	
      FILE* f;
         
      if((f = fopen(filename, "wb")) == NULL){
         if(mpiroot){
            cout << endl;
            cout << "Error: Cannot open file " << filename << endl;
            cout << endl;	
         }
         
         MPI_Finalize();
         
         exit(1);
      }
      
      fwrite(&m, sizeof(long), 1, f);
      fwrite(&n, sizeof(long), 1, f);
      fwrite(&type, sizeof(int), 1, f);
      
      long k = 0;
      
      for(long i = 0; i < m; i++){
         for(long j = 0; j < n; j++){
            fwrite(&A[k++], sizeof(T), 1, f);
         }
      }
      
      fclose(f);      
      //---------------------------------------------------------------------------------------------------------------
      double time2 = MPI_Wtime();      
      
      cout << "Write time: " << (time2 - time1) << endl;
   }
}
//END------------------------------------------------------WRITE-TO-FILE-----------------------------------------------

