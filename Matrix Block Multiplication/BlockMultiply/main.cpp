#include <iostream>
#include <cstdlib>
#include <stdio.h>
#include "matrix.cpp"
#include <papi.h>

int readType(const char* filename){
   int type = -1;   
   
   FILE* f;

   if((f = fopen(filename,"rb")) == NULL) {
      printf("Error: Cannot open file %s\n",filename);
      exit(1);
   }

   fseek(f,sizeof(long)*2,SEEK_SET);
  
   fread(&type,sizeof(int),1,f);
   if(type<1){
      printf("Incorrect value: 'type'\n");      
      exit(1);
   }

   fclose(f);
   
   return type;
}

void multiply(const char* file1,const char* file2,const char* file3){
   int type1 = readType(file1);
   int type2 = readType(file2);
      
   if(type1!=type2){
      printf("Type mismatch\n");      
      exit(1);
   }

   switch(type1){
      case 1:{
            Matrix<int> A(file1);
            Matrix<int> B(file2);
            Matrix<int> C = A*B; 
            C.writeToFile(file3);
            break;
      }
      case 2:{
            Matrix<long> A(file1);
            Matrix<long> B(file2);
            Matrix<long> C = A*B;               
            C.writeToFile(file3);
            break;
      }
      case 3:{
            Matrix<float> A(file1);
            Matrix<float> B(file2);
            Matrix<float> C = A*B;               
            C.writeToFile(file3);
            break;      
      }
      case 4:{
            Matrix<double> A(file1);
            Matrix<double> B(file2);          
            Matrix<double> C = A*B;               
            C.writeToFile(file3);
            break;
      }
      default:{
            break;
      }
   }   
}

int main(int argc,char** argv){
  if(argc!=4){
      printf("\nUsage: ./a.out <file1_in> <file2_in> <file3_out>\n\n");
      return 1;
   }
 
	if(PAPI_num_counters()<2) {
   	printf("No hardware counters here, or PAPI not supported.\n");
   	exit(1);
	}

   multiply(argv[1],argv[2],argv[3]);
   
   return 0;

}
