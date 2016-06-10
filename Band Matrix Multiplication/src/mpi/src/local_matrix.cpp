#include "../include/local_matrix.h"

//BEGIN----------------------------------------------------ALLOC-MEMORY------------------------------------------------
template <typename T> 
void LocalMatrix<T>::alloc_memory(){
   data = new T*[m];
   
   for(long i = 0; i < m; i++){
      data[i] = new T[n];
      
      for(long j = 0; j < n; j++){
         data[i][j] = 0;
      }      
   }
   
   mem_alloced = true;
}
//END------------------------------------------------------ALLOC-MEMORY------------------------------------------------

//BEGIN----------------------------------------------------FREE-DATA---------------------------------------------------
template <typename T> 
void LocalMatrix<T>::free_data(){
   if(mem_alloced){
      for(long i = 0; i < m; i++){
         delete[] data[i];      
      }
      
      delete[] data;
      
      data = NULL;
      
      mem_alloced = false;
   }
}
//END------------------------------------------------------FREE-DATA---------------------------------------------------

//BEGIN----------------------------------------------------SET-SIZE----------------------------------------------------
template <typename T> 
void LocalMatrix<T>::set_size(long _m, long _n){
   m = _m;
   n = _n;
   
   alloc_memory();
}
//END------------------------------------------------------SET-SIZE----------------------------------------------------

