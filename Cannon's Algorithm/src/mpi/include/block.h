#ifndef _BLOCK_H_
#define _BLOCK_H_

template <typename T> 
class Block{

public:
   Block<T>(){}
   
   void set_size(int _m, int _n){
      m = _m;
      n = _n;

      alloc_memory();
      empty();
   }
   
   void alloc_memory(){
      data = new T[m * n];
   }
   
   void empty(){
      for(long i = 0; i < m*n; i++){
         data[i] = 0;
      }
   }
   
   int get_m() const{ return m; }

   ~Block<T>(){
      delete[] data;

      data = NULL;
   }   

   T* data;

private:
   int m;
   int n;
};

#endif

