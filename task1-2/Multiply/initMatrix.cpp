#include <iostream>
#include <cstdlib>
#include <cstdio>
#include "matrix.cpp"

using namespace std;

int main(int argc,char** argv){
   long m,n;
   int type;

   cout<<endl;

   if(argc!=5){
      cout<<"Usage: ./init <m> <n> <type> <filename>"<<endl<<endl;
      cout<<"       <type>: 1 - int, 2 - long, 3 - float, 4 - double"<<endl<<endl;
      exit(1);
   }

   m = -1;
   m = (long)atoi(argv[1]);
   if(m<1){
      cout<<"Error: Incorrect value \"m\""<<endl<<endl;
      exit(1);
   }
   
   n = -1;
   n = (long)atoi(argv[2]);
   if(n<1){
      cout<<"Error: Incorrect value \"n\""<<endl<<endl;
      exit(1);
   }

   type = -1;
   type = atoi(argv[3]);
   if(type<1){
      cout<<"Error: Incorrect value \"type\""<<endl<<endl;
      exit(1);
   }
   
   switch(type){
      case 1:{
         Matrix<int> A(m,n);
         A.generate();
         A.writeToFile(argv[4]);
         if(m<256 && n<256)A.println();
         break;
      }
      case 2:{
         Matrix<long> A(m,n);
         A.generate();
         A.writeToFile(argv[4]);
         if(m<256 && n<256)A.println();
         break;
      }
      case 3:{
         Matrix<float> A(m,n);
         A.generate();
         A.writeToFile(argv[4]);
         if(m<256 && n<256)A.println();
         break;
      }
      case 4:{
         Matrix<double> A(m,n);
         A.generate();
         A.writeToFile(argv[4]);
         if(m<256 && n<256)A.println();
         break;
      }
      default:{
         break;
      }
   }
   
   return 0;
}
