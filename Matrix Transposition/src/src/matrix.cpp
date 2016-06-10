#include "../include/matrix.h"

//BEGIN---------------------------READ-FILE'S-HEADER----------------------------
template <class T>
void Matrix<T>::read_file_header(const char* filename){
   char* file = const_cast<char *>(filename);
   	
   MPI_File f;	
   //----------------------------------------------------------------------
   m = -1;
   n = -1;
   type = -1;
   //----------------------------------------------------------------------
   MPI_Status status;
   
   MPI_File_open(MPI_COMM_WORLD, file, MPI_MODE_RDONLY, MPI_INFO_NULL, &f);
   	
   MPI_File_read(f, &m, 1, MPI_LONG, &status);
   MPI_File_read(f, &n, 1, MPI_LONG, &status);
   MPI_File_read(f, &type, 1, MPI_INT, &status);
   //----------------------------------------------------------------------	
   if(m < 1){
      printf("Error: incorrect parameter 'm'\n");
      MPI_Finalize();
   }
   	
   if(n < 1){
      printf("Error: incorrect parameter 'n'\n");
      MPI_Finalize();
   }
   	
   if(type < 1 || type > 4){
      printf("Error: incorrect parameter 'type'\n");
      MPI_Finalize();
   }
   	
   set_type();
   //----------------------------------------------------------------------
   MPI_Barrier(MPI_COMM_WORLD);
   	
   MPI_File_close(&f);
}
//END-----------------------------READ-FILE'S-HEADER----------------------------

//BEGIN---------------------------READ-FROM-FILE--------------------------------
template <class T>
Matrix<T>::Matrix(const char* filename, strip_type s){
   init_mpi();
   //----------------------------------------------------------------------
   char* file = const_cast<char *>(filename);
   read_file_header(filename);
   //----------------------------------------------------------------------
   MPI_File f;
   MPI_Status status;
   		
   int header_shift = sizeof(long)*2+sizeof(int);			
   
   double t0;
   
   switch(s){
      case ROWS:{
         int blocksize_i = m/nrows;
         			
         dm = ((rank/ncols)+1 < nrows) ? (blocksize_i) : (m-(nrows-1)*blocksize_i);
         
         long bufsize = blocksize_i*n*sizeof(T);
	 	   
         A.resize(dm*n);
	 		
         t0 = MPI_Wtime();		
         MPI_File_open(MPI_COMM_WORLD, file, MPI_MODE_RDONLY, MPI_INFO_NULL,&f);
         
         int i = rank/ncols;
         
         MPI_File_seek(f, header_shift+i*bufsize, MPI_SEEK_SET);
         
         for(long i = 0; i < dm; i++){
	   MPI_File_read_all(f, &A[i*m], n, mpi_type, &status);
         }
         
         break;
      }
      
      case COLS:{
         int blocksize_j = n/size;
	 		
         dm = m;
         dn = blocksize_j;
	 long bufsize = blocksize_j*dm*sizeof(T);
         	
         A.resize(m*dn);
         
         t0 = MPI_Wtime();
         
         MPI_File_open(MPI_COMM_WORLD, file, MPI_MODE_RDONLY, MPI_INFO_NULL,&f);
	 		
         MPI_File_seek(f, header_shift+rank*blocksize_j*sizeof(T), MPI_SEEK_SET);
	 	 
         for(long i = 0; i < m; i++){
            MPI_File_read(f, &A[i*dn], dn, mpi_type, &status);	 		
            MPI_File_seek(f, (n-dn)*sizeof(T), MPI_SEEK_CUR);			
         }
         
         break;
      }
      default:{
         break;
      }
   }	
   //---------------------------------------------------------------------------
   MPI_Barrier(MPI_COMM_WORLD);
   double t1 = MPI_Wtime();
   
   MPI_File_close(&f);
}

//END-----------------------------READ-FROM-FILE--------------------------------

//BEGIN---------------------------SET-------------------------------------------
template <class T>
void Matrix<T>::set_type(){
   char s[3];
   
   strcpy(s, typeid(T).name());
   	
   if(strcmp(s, "i") == 0){
      type = 1;
      mpi_type = MPI_INT;
   }else if(strcmp(s, "l") == 0){
      type = 2;
      mpi_type = MPI_LONG;
   }else if(strcmp(s, "f") == 0){
      type = 3;
      mpi_type = MPI_FLOAT;
   }else if(strcmp(s, "d") == 0){
      type = 4;
      mpi_type = MPI_DOUBLE;
   }else type = 0;
}
//END-----------------------------SET-------------------------------------------

//BEGIN---------------------------PRINT-----------------------------------------
template <class T>
void Matrix<T>::print(strip_type s = NONE){
   long M, N;
   
   switch(s){
      case ROWS:
         M = dm;
         N = n;
         break;
      case COLS:
         M = m;
         N = dn;
         break;
      default:
         M = dm;
         N = dn;
         break;	
   }
   
   for(int k = 0; k < size; k++){
      if(k == rank){
         cout<<"rank: "<<rank<<endl;
         for(long i = 0; i < M; i++){
            for(long j = 0; j < N; j++){
               cout << A[i*N+j] << " ";         
            }
            cout << endl;
         }
         cout << endl << flush;
      }   
	
      MPI_Barrier(MPI_COMM_WORLD);
   }
}
//END-----------------------------PRINT-----------------------------------------

//BEGIN---------------------------OPERATOR * -----------------------------------
template <class T>
Matrix<T> Matrix<T>::operator*(Matrix<T>& tmp){   
   long m1 = dm;
   long n1 = n;
   
   long m2 = tmp.get_m();
   long n2 = tmp.get_dn();
   
   if(n1 != m2){
      if(rank == 0){
         printf("\nSize mismatch\n\n");
      }

      MPI_Finalize();
      exit(1);
   }
   
   long M = m;
   long N = tmp.get_n();
   Matrix<T> C(M,N,m1,n2);
   
   double t0 = MPI_Wtime();   
   for(long i = 0; i < m1; i++){
      for(long j = 0; j < n2; j++){
         for(long k = 0; k < n1; k++){
            C.A[i][j] += A[i][k]*tmp.A[k][j];
         }
      }
   }
   MPI_Barrier(MPI_COMM_WORLD);
   double t1 =  MPI_Wtime();   
   	
   if(rank == 0){
      cout << "Multiply time: " << (t1-t0) << endl;   
   }
   MPI_Barrier(MPI_COMM_WORLD);
   
   return C;
}
//END-----------------------------OPERATOR * -----------------------------------

//BEGIN---------------------------TRANSPOSE-------------------------------------
template <class T>
void Matrix<T>::local_transpose(T *a, int n){
   int ij = 0, ji, l = -1;
   T tmp;
      
   for(int i = 0; i < n; i++){
      l += n + 1;
      ji = l;
      ij += i + 1;
      
      for(int j = i+1; j < n; j++){
         tmp = a[ij];
      	 a[ij] = a[ji];
      	 a[ji] = tmp;
      	      
         ij++;
      	 ji += n;
      }
   }
}
//-----------------------------------------------------------------------------
template <class T>
void Matrix<T>::transpose(){    
   T a[m][dn];
   
   MPI_Barrier(MPI_COMM_WORLD);
   //-------------------------------------
   double t0 = MPI_Wtime();           
   
   MPI_Alltoall(A.data(), dn*dn, mpi_type,
   	             &a[0][0], dn*dn, mpi_type, 
                MPI_COMM_WORLD);      
      
   for(int i = 0; i < size; i++){
      local_transpose(&a[i*dn][0], dn);   
   }
   
   A.assign(&a[0][0],&a[m][dn]);
   //-------------------------------------
   MPI_Barrier(MPI_COMM_WORLD);
   double t1 = MPI_Wtime();
   
   if(rank == 0){
      cout<<"\nTranspose time: "<<(t1-t0)<<endl;
   }
}	
//END-----------------------------TRANSPOSE------------------------------------

//BEGIN---------------------------WRITE-TO-FILE--------------------------------
template <class T>
void Matrix<T>::writeToFile(const char* filename){
   char* file = const_cast<char *>(filename);
   char native[10] = "native";
   			
   long M = m;
   long N = n;
   
   double t0 = MPI_Wtime();
   //---------------------------------------------------------------------------
   if(rank == 0){
      FILE* f;
      		
      if((f = fopen(filename,"wb")) == NULL) {
         printf("Cannot open file: %s\n",file);
         MPI_Finalize();
         exit(1);
      }
      	
      fwrite(&M,sizeof(long),1,f);
      fwrite(&N,sizeof(long),1,f);
      fwrite(&type,sizeof(int),1,f);
		
      fclose(f);    			
   }
   	
   MPI_Barrier(MPI_COMM_WORLD);
   //---------------------------------------------------------------------------
   MPI_File f;
   MPI_File_open(MPI_COMM_WORLD, file, MPI_MODE_CREATE | MPI_MODE_WRONLY | MPI_MODE_APPEND, MPI_INFO_NULL, &f);
   int block_i = m;
   int block_j = n/size;
   long displ = sizeof(int)+sizeof(long)*2;
   		
   MPI_File_set_view(f, displ, mpi_type, mpi_type, native, MPI_INFO_NULL);
   MPI_Status status;		
   	
   long i = rank;
   
   displ = 0;
   displ += i*block_j;
      			
   long displ_j = displ;			
   
   for(long j = 0; j < dm; j++){
      MPI_File_write_at_all(f, displ_j, &A[j*dn], dn, mpi_type, &status);		   
      displ_j += n;
   }

   MPI_File_close(&f);

   MPI_Barrier(MPI_COMM_WORLD);
   double t1 = MPI_Wtime();
}
//END-----------------------------WRITE-TO-FILE---------------------------------
