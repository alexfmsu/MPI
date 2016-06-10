#ifndef LOCALMATRIX_H
#define LOCALMATRIX_H

#include <iomanip>

template <typename T> 
class LocalMatrix{

public:
   //BEGIN-------------------------------------------------CONSTRUCTOR-------------------------------------------------	
   LocalMatrix<T>(){ 
      m = -1;
      n = -1;
      
      mem_alloced = false;
   };
   //END---------------------------------------------------CONSTRUCTOR-------------------------------------------------	
   
   void set_size(long _m, long _n);
   
   void alloc_memory();
   
   void free_data();
   
   //BEGIN-------------------------------------------------GET---------------------------------------------------------
   long get_m() const{ return m; }
   long get_n() const{ return n; }
   //END---------------------------------------------------GET---------------------------------------------------------
   
   //BEGIN-------------------------------------------------PRINT-------------------------------------------------------
   void print(){
      using namespace std;
      
      if(!mem_alloced){
         cout << "Matrix is empty" << endl;
         cout << endl;

         return;
      }
      
      cout.precision(10);		
      
      for(long i = 0; i < m; i++){
         for(long j = 0; j < n; j++){
            cout << setw(10) << data[i][j] << " ";
         }
         
         cout << endl;
      }
      
      cout << endl;
   }
   //END---------------------------------------------------PRINT-------------------------------------------------------
   
   //BEGIN-------------------------------------------------DESTRUCTOR--------------------------------------------------
   ~LocalMatrix<T>(){
      free_data();
   }
   //END---------------------------------------------------DESTRUCTOR--------------------------------------------------
   
   T** data;	
   
private:
   long m, n;
   int type;
   
   bool mem_alloced;
};
#endif
