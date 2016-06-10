#include "../include/matrix.h"

//BEGIN----------------------------------------------------ALLOC-MEMORY------------------------------------------------
template <typename T>
void Matrix<T>::alloc_memory(){
   A = new T*[m];
   			   	
   for(long i = 0; i < m; i++){
      A[i] = new T[n];      
   }
   				
   empty();
}		
//END------------------------------------------------------ALLOC-MEMORY------------------------------------------------
   
//BEGIN----------------------------------------------------EMPTY-------------------------------------------------------
template <typename T>
void Matrix<T>::empty(){
   for(long i = 0; i < m; i++){
      for(long j = 0; j < n; j++){
         A[i][j] = 0;         
      }
   }
} 
//END------------------------------------------------------EMPTY-------------------------------------------------------

//BEGIN----------------------------------------------------SET---------------------------------------------------------
template <typename T>
void Matrix<T>::set_type(){
   char s[3];
   	
   strcpy(s, typeid(T).name());
   		
   if(strcmp(s, "i") == 0){
      type = 1;
   }else if(strcmp(s, "l") == 0){
      type = 2;
   }else if(strcmp(s, "f") == 0){
      type = 3;
   }else if(strcmp(s, "d") == 0){
      type = 4;
   }else type = -1;
}
//END------------------------------------------------------SET---------------------------------------------------------

//BEGIN----------------------------------------------------OPERATOR * -------------------------------------------------
template <typename T>
Matrix<T> Matrix<T>::operator*(Matrix<T>& tmp){   
   //------------------------------------------------------------------------------------------------------------------
   long m1 = m;
   long n1 = n;
   	
   long m2 = tmp.get_m();
   long n2 = tmp.get_n();
   //------------------------------------------------------------------------------------------------------------------
   if(n1 != m2){
      printf("Error: size mismatch\n");
      exit(1);
   }
   //------------------------------------------------------------------------------------------------------------------
   Matrix<T> C(m1, n2);
   //------------------------------------------------------------------------------------------------------------------
   clock_t time1 = clock();      
   
   for(long i = 0;i < m1; i++){
      for(long j = 0; j < n2; j++){
         for(long k = 0; k < n1; k++){
            C.A[i][j] += A[i][k] * tmp.A[k][j];
         }
      }
   }
   clock_t time2 = clock();   
   //------------------------------------------------------------------------------------------------------------------
   cout << "Multiply time: " << (time2 - time1) / (double)(CLOCKS_PER_SEC) << endl;
   cout << endl;   
   //------------------------------------------------------------------------------------------------------------------ 
   return C;
}
//END------------------------------------------------------OPERATOR * -------------------------------------------------

//BEGIN----------------------------------------------------READ-FROM-FILE----------------------------------------------
template <typename T>
void Matrix<T>::read_from_file(const char* filename){
   clock_t time1 = clock();      
   //------------------------------------------------------------------------------------------------------------------
   FILE* f;   
   	
   if((f = fopen(filename, "rb")) == NULL){
      cout << endl;
      cout << "Error: Cannot open file " << filename << endl;
      cout << endl;	
      		
      exit(1);
   }
   //------------------------------------------------------------------------------------------------------------------
   if(fread(&m, sizeof(long), 1, f) != 1){
      exit(1);
   }
   	   
   if(m < 1){
      cout << endl;
      cout << "Error: parameter 'm' not found" << endl;
      cout << endl;
      		
      exit(1);
   }
   //------------------------------------------------------------------------------------------------------------------
   if(fread(&n, sizeof(long), 1, f) != 1){
      exit(1);
   }
   	   
   if(n < 1){
      cout << endl;
      cout << "Error: parameter 'n' not found" << endl;
      cout << endl;
      	
      exit(1);
   }
   //------------------------------------------------------------------------------------------------------------------
   set_type();	
   //------------------------------------------------------------------------------------------------------------------
   alloc_memory();
   	
   fseek(f, sizeof(long) * 2 + sizeof(int), SEEK_SET);
   		
   for(long i = 0; i < m; i++){
      for(long j = 0; j < n; j++){
         if(fread(&A[i][j], sizeof(T), 1, f) != 1){
            exit(1);
         }     
          		
         if(sizeof(A[i][j]) != sizeof(T)){
            cout << endl;
            cout << "Error: Cannot read matrix from file " << filename << endl;
            cout << endl;
            
            exit(1);
         }               
      }	
   }
   	
   fclose(f);	
   //------------------------------------------------------------------------------------------------------------------
   clock_t time2 = clock();      
   	
   cout << "size: " << get_m() << endl << endl;
   cout << "Read time: " << (time2 - time1) / (double)(CLOCKS_PER_SEC) << endl;
   //------------------------------------------------------------------------------------------------------------------	
}
//END------------------------------------------------------READ-FROM-FILE----------------------------------------------

//BEGIN----------------------------------------------------PRINT-------------------------------------------------------
template <typename T>
void Matrix<T>::print() const{
   using namespace std;
   	
   cout.precision(10);		
   		
   for(long i = 0; i < m; i++){
      for(long j = 0; j < n; j++){
         cout << setw(10) << A[i][j] << " ";
      }
      
      cout << endl;
   }
		
   cout << endl;
}
//END------------------------------------------------------PRINT-------------------------------------------------------

//BEGIN----------------------------------------------------WRITE-TO-FILE-----------------------------------------------
template <typename T>
void Matrix<T>::write_to_file(const char* filename){
   clock_t time1 = clock();      
   		
   FILE* f;
   	        
   if((f = fopen(filename, "wb")) == NULL){
      cout << endl;
      cout << "Error: Cannot open file " << filename << endl;
      cout << endl;	
      	
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
   clock_t time2 = clock();      
   cout << "Write time: " << (time2 - time1) / (double)(CLOCKS_PER_SEC) << endl;
}
//END------------------------------------------------------WRITE-TO-FILE-----------------------------------------------

