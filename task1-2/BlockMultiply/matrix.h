#include <typeinfo>
#include <cstring>
#include <limits>

template <typename T> class Matrix{
public:
   Matrix<T>(long _m,long _n){
      m = _m;
      n = _n;

      setType();      
            
      A = new T*[m];
      for(long i=0;i<m;i++){
         A[i] = new T[n];      
         for(long j=0;j<n;j++){
            A[i][j] = 0;         
         }
      }
   }
   
   Matrix<T>(const Matrix<T> &tmp){
      m = tmp.m;
      n = tmp.n;
      
      A = new T*[m];
      for(long i=0;i<m;i++){
         A[i] = new T[n];    
         for(long j=0;j<n;j++){
            A[i][j] = tmp.A[i][j];            
         }      
      }
   }
   
   ~Matrix(){      
      if(A){
         for(long i=0;i<m;i++){
            delete[] A[i];      
         }

         delete[] A;
         A = NULL;
      }      
   }

   void generate(){
      if(!A){
         printf("Error\n");
         exit(1);
      }

      srand((unsigned)time(0));
      
      T max = std::numeric_limits<T>::max();
      T min = -max+1; 
      
      for(long i=0;i<m;i++){
         for(long j=0;j<n;j++){
            double r = (double)rand()/(double)RAND_MAX;
            A[i][j] = min+(T)(max*r)-(T)(min*r);
         }
      }      
   }  
   
   long getM(){
      return m;
   }
   
   long getN(){
      return n;
   }
   
   void setType(){
      char s[10];
      
      strcpy(s,typeid(T).name());
      
      if(strcmp(s,"i")==0){
         type = 1;
      }else if(strcmp(s,"l")==0){
         type = 2;
      }else if(strcmp(s,"f")==0){
         type = 3;
      }else if(strcmp(s,"d")==0){
         type = 4;
      }else type = 0;
   }

   int getType(){
      return type;
   }

   Matrix(const char* filename){
      FILE* f;
      
      if((f = fopen(filename,"rb")) == NULL) {
         printf("Error: Cannot open file %s\n",filename);
         exit(1);
      }
      
      m = -1;
      fread(&m,sizeof(long),1,f);
      if(m<1){
         printf("Error: parameter 'm' not found\n");
         exit(1);
      }
      
      n = -1;
      fread(&n,sizeof(long),1,f);
      if(n<1){
         printf("Error: parameter 'n' not found\n");
         exit(1);
      }
      
    
      A = new T*[m];
      for(long i=0;i<m;i++){
         A[i] = new T[n];      
      }

      fseek(f,sizeof(long)*2+sizeof(int),SEEK_SET);

      for(long i=0;i<m;i++){
         for(long j=0;j<n;j++){
            fread(&A[i][j],sizeof(T),1,f);      
            if(sizeof(A[i][j])!=sizeof(T)){
               printf("Error: Cannot read matrix from file '%s'",filename);               
               fclose(f);
               exit(1);
            }         
         }      
      }

      fclose(f);   
   }

   void writeToFile(const char* filename){
      FILE* f;
         
      if((f = fopen(filename,"wb")) == NULL) {
         printf("Cannot open file: %s\n",filename);
         exit(1);
      }
      
      fwrite(&m,sizeof(long),1,f);
      fwrite(&n,sizeof(long),1,f);
      fwrite(&type,sizeof(int),1,f);
      
      for(long i=0;i<m;i++){
         for(long j=0;j<n;j++){
            fwrite(&A[i][j],sizeof(long),1,f);
         }
      }
      
      fclose(f);      
   }
  
   void print(){
      using namespace std;

      for(long i=0;i<m;i++){
         for(long j=0;j<n;j++){
            cout<<A[i][j]<<" ";         
         }
         cout<<endl;
      }   
   }
   
   void println(){
      using namespace std;
      print();
      cout<<endl;      
   }
   
   Matrix<T> operator*(Matrix<T>& tmp);  
   
protected:
   T** A;

private:
   long m,n;
   int type;
};
