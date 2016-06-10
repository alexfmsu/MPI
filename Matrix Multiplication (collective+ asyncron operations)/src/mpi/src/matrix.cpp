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
      A = new T*[m];
      
      for(long i = 0; i < m; i++){
         A[i] = new T[n];      
      }
      
      empty();
      
      mem_alloced = true;
   }
}		
//END------------------------------------------------------ALLOC-MEMORY------------------------------------------------
   
//BEGIN----------------------------------------------------EMPTY-------------------------------------------------------
template <typename T>
void Matrix<T>::empty(){
   if(mpiroot){
      for(long i = 0; i < m; i++){
         for(long j = 0; j < n; j++){
            A[i][j] = 0;         
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
         for(long i = 0; i < m; i++){
            delete[] A[i];      
         }
         
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
void Matrix<T>::scatter(){
   long dm = m/size; 
   
   loc.set_size(dm, n);
   //------------------------------------------------------------------------------------------------------------------    
   MPI_Barrier(MPI_COMM_WORLD);
   
   double time1 = MPI_Wtime();   	
   //------------------------------------------------------------------------------------------------------------------
   MPI_Request req_s[dm];
   MPI_Request req_r[dm];
   
   MPI_Status st_r[dm];
   //---------------------------------------------------------------------------------------------------------------
   for(long j = 0; j < dm; j++){
      if(mpiroot){
         for(int i = 0; i < size; i++){
            MPI_Isend(A[i*dm+j], n, mpitype, i, (i*dm+j), MPI_COMM_WORLD, &req_s[j]);            
         }      
      }
      
      MPI_Irecv(loc.data[j], n, mpitype, 0, (rank*dm+j), MPI_COMM_WORLD, &req_r[j]);
   }
   
   MPI_Waitall(dm, req_r, st_r);     
   //---------------------------------------------------------------------------------------------------------------
   MPI_Barrier(MPI_COMM_WORLD);
   
   double time2 = MPI_Wtime();   	
   //------------------------------------------------------------------------------------------------------------------
   if(mpiroot){
      cout << "Scatter time: " << (time2 - time1) << endl; 
   }
   //------------------------------------------------------------------------------------------------------------------	
   free_T();
}

template <typename T>
void Matrix<T>::gather(){
   long dm = loc.get_m();
   //------------------------------------------------------------------------------------------------------------------	
   MPI_Barrier(MPI_COMM_WORLD);
   
   double time1 = MPI_Wtime();   	
   //------------------------------------------------------------------------------------------------------------------
   MPI_Request req_s[dm];
   MPI_Request req_r[dm];
   
   MPI_Status st_r[dm];
   //---------------------------------------------------------------------------------------------------------------
   for(long j = 0; j < dm; j++){
      if(mpiroot){
         for(int i = 0; i < size; i++){
            MPI_Irecv(A[i*dm+j], n, mpitype, i, (i*dm+j), MPI_COMM_WORLD, &req_s[j]);            
         }      
      }
      
      MPI_Isend(loc.data[j], n, mpitype, 0, (rank*dm+j), MPI_COMM_WORLD, &req_r[j]);
   }

   MPI_Waitall(dm, req_r, st_r);
   //---------------------------------------------------------------------------------------------------------------	
   MPI_Barrier(MPI_COMM_WORLD);
   
   double time2 = MPI_Wtime();   	
   //------------------------------------------------------------------------------------------------------------------	
   if(mpiroot){
      cout << "Gather time: " << (time2 - time1) << endl;
      cout << endl;   
   }
   //------------------------------------------------------------------------------------------------------------------
   loc.free_data();
}
//END------------------------------------------------------SCATTER/GATHER----------------------------------------------

//BEGIN----------------------------------------------------OPERATOR * -------------------------------------------------
template <typename T>
Matrix<T> Matrix<T>::operator*(Matrix<T>& tmp){   
   //------------------------------------------------------------------------------------------------------------------
   long m1 = loc.get_m();
   long n1 = loc.get_n();
   
   long m2 = tmp.loc.get_m();
   long n2 = tmp.loc.get_n();
   //------------------------------------------------------------------------------------------------------------------
   if(n != tmp.get_m()){
      if(rank == 0){
         cout << endl;
         cout << "Error: size mismatch" << endl;
         cout << endl;
      }
      
      MPI_Finalize();
      
      exit(1);
   }
   //------------------------------------------------------------------------------------------------------------------
   Matrix<T> C(m, n2);
   
   C.loc.set_size(m1, n2);
   
   long dj = m2*rank;
   //------------------------------------------------------------------------------------------------------------------
   MPI_Barrier(MPI_COMM_WORLD);  	
   
   double time1 = MPI_Wtime();      
   //------------------------------------------------------------------------------------------------------------------
   for(long i2 = 0; i2 < m1; i2++){
      for(long j1 = 0; j1 < m2; j1++){
         for(long i1 = 0; i1 < n2; i1++){
            C.loc.data[i2][i1] += loc.data[i2][j1+dj] * tmp.loc.data[j1][i1];
         }
      }
   }
   //------------------------------------------------------------------------------------------------------------------
   MPI_Barrier(MPI_COMM_WORLD);
   //------------------------------------------------------------------------------------------------------------------
   for(int i = 0; i < size-1; i++){	
      MPI_Request req_r[m2];
      MPI_Status status[m2];
      
      for(long j = 0; j < m2; j++){
         MPI_Isend(tmp.loc.data[j], n, mpitype, (rank) ? (rank-1) : (size-1), j, MPI_COMM_WORLD, &req_r[j]);
         MPI_Irecv(tmp.loc.data[j], n, mpitype, (rank == size-1) ? (0) : (rank+1), j, MPI_COMM_WORLD, &req_r[j]);
      }

      MPI_Waitall(m2, req_r, status);
      
      dj = (dj + m2) % n1;
      //---------------------------------------------------------------------------------------------------------------
      for(long i2 = 0; i2 < m1; i2++){
         for(long j1 = 0; j1 < m2; j1++){
            for(long i1 = 0; i1 < n2; i1++){
               C.loc.data[i2][i1] += loc.data[i2][j1+dj] * tmp.loc.data[j1][i1];
           }
         }
      }
      //---------------------------------------------------------------------------------------------------------------
      MPI_Barrier(MPI_COMM_WORLD);
   }
   //------------------------------------------------------------------------------------------------------------------
   MPI_Barrier(MPI_COMM_WORLD);  	

   double time2 = MPI_Wtime();   	
   //------------------------------------------------------------------------------------------------------------------
   if(mpiroot){
      cout << "Multiply time: " << (time2 - time1) << endl;
      cout << endl;   
   }
   //------------------------------------------------------------------------------------------------------------------
   C.gather();
   //------------------------------------------------------------------------------------------------------------------ 
   return C;
}
//END------------------------------------------------------OPERATOR * -------------------------------------------------

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
   }
  
   MPI_Barrier(MPI_COMM_WORLD);  
 
   distrib_dimensions(m, n, type, mpiroot);
   
   set_type();	
   
   alloc_memory();
   
   if(mpiroot){
      fseek(f, sizeof(long) * 2 + sizeof(int), SEEK_SET);
      
      for(long i = 0; i < m; i++){
         for(long j = 0; j < n; j++){
            if(fread(&A[i][j], sizeof(T), 1, f) != 1){
               MPI_Abort(MPI_COMM_WORLD, 1);	
               
               exit(1);	
            }      
            
            if(sizeof(A[i][j]) != sizeof(T)){
               if(rank == 0){
                  cout << endl;
                  cout << "Error: Cannot read matrix from file " << filename << endl;
                  cout << endl;
               }
               
               MPI_Abort(MPI_COMM_WORLD, 1);
               
               exit(1);	
            }         
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
   scatter();
}
//END------------------------------------------------------READ-FROM-FILE----------------------------------------------

//BEGIN----------------------------------------------------PRINT-------------------------------------------------------
template <typename T>
void Matrix<T>::print() const{
   using namespace std;
   
   if(mpiroot){
      if(!mem_alloced){
         cout << "Matrix is empty" << endl;
         cout << endl;
         
         return;
      }
      
      cout.precision(10);		
      
      for(long i = 0; i < m; i++){
         for(long j = 0; j < n; j++){
            cout << setw(10) << A[i][j] << " ";
         }
         
         cout << endl;
      }
      
      cout << endl;
   }
}

template <typename T>
void Matrix<T>::print_nodes(){
   for(int r = 0; r < size; r++){
      MPI_Barrier(MPI_COMM_WORLD);
      
      if(rank == r){
         cout << "node " << rank << ":" << endl;
         
         loc.print();
      }
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
      
      for(long i = 0; i < m; i++){
         for(long j = 0; j < n; j++){
            fwrite(&A[i][j], sizeof(T), 1, f);
         }
      }
      
      fclose(f);      
      //---------------------------------------------------------------------------------------------------------------
      double time2 = MPI_Wtime();      
      
      cout << "Write time: " << (time2 - time1) << endl;
   }
}
//END------------------------------------------------------WRITE-TO-FILE-----------------------------------------------

