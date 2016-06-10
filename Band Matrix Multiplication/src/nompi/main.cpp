#include <iostream>
#include <cstdlib>
#include <cstdio>

#include "src/matrix.cpp"

using namespace std;

int readType(const char* filename){
	int type = -1;   
	//------------------------------------------------------------------------------------------------------------------
	FILE* f;
			
	if((f = fopen(filename, "rb")) == NULL){
		cout << endl;
		cout << "Error: Cannot open file '" << filename << "'" << endl;
		cout << endl;	
			
		exit(1);
	}
	//------------------------------------------------------------------------------------------------------------------
	fseek(f, sizeof(long) * 2, SEEK_SET);
	
	if(fread(&type, sizeof(int), 1, f) != 1){
		exit(1);
	}
	
	if(type < 1){
		cout << endl;
		cout << "Incorrect value: 'type'" << filename << "'" << endl;
		cout << endl;	
			
		exit(1);
	}
	
   fclose(f);
   //------------------------------------------------------------------------------------------------------------------
   return type;
}

void multiply(const char* file1, const char* file2, const char* file3){
   int type1 = readType(file1);
   int type2 = readType(file2);
   //------------------------------------------------------------------------------------------------------------------   
   if(type1 != type2){
      cout << endl;
		cout << "Error: type mismatch" << endl;
		cout << endl;
		
      exit(1);
   }
	//------------------------------------------------------------------------------------------------------------------
   switch(type1){
      case 1:{
			cout << "type: 'int'" << endl;  
			cout << "-------------------------" << endl;
			cout << "Matrix A:" << endl << endl;
			Matrix<int> A;
			A.read_from_file(file1);
			cout << "-------------------------" << endl;
			
			cout << "Matrix B:" << endl << endl;
			Matrix<int> B;
			B.read_from_file(file2);
			cout << "-------------------------" << endl;
			
			cout << "Matrix C:" << endl << endl;
			Matrix<int> C = A*B; 
			C.write_to_file(file3);
						
			break;
		}
		case 2:{
			cout << "type: 'long'" << endl;  
			cout << "-------------------------" << endl;
			cout << "Matrix A:" << endl << endl;
			Matrix<long> A;
			A.read_from_file(file1);
			cout << "-------------------------" << endl;
			
			cout << "Matrix B:" << endl << endl;
			Matrix<long> B;
			B.read_from_file(file2);
			cout << "-------------------------" << endl;
			
			cout << "Matrix C:" << endl << endl;
			Matrix<long> C = A*B;
			C.write_to_file(file3);
			
			break;
		}
		case 3:{
			cout << "type: 'float'" << endl;  
			cout << "-------------------------" << endl;
			cout << "Matrix A:" << endl << endl;
			Matrix<float> A;
			A.read_from_file(file1);
			cout << "-------------------------" << endl;
			
			cout << "Matrix B:" << endl << endl;
			Matrix<float> B;
			B.read_from_file(file2);
			cout << "-------------------------" << endl;
			
			cout << "Matrix C:" << endl << endl;
			Matrix<float> C = A*B; 
			C.write_to_file(file3);
			
			break;
		}
		case 4:{
			cout << "type: 'double'" << endl;  
			cout << "-------------------------" << endl;
			cout << "Matrix A:" << endl << endl;
			Matrix<double> A;
			A.read_from_file(file1);
			cout << "-------------------------" << endl;
																			
			cout << "Matrix B:" << endl << endl;
			Matrix<double> B;
			B.read_from_file(file2);
			cout << "-------------------------" << endl;
			
			cout << "Matrix C:" << endl << endl;
			Matrix<double> C = A*B;
			C.write_to_file(file3);
			
			break;
		}
		default:{
			break;
		}
	} 
}

int main(int argc, char** argv){
	clock_t t1 = clock(); 	
	//------------------------------------------------------------------------------------------------------------------	
	if(argc != 4){
		cout << endl;
		cout << "Usage: ./a.out <file1_in> <file2_in> <file3_out>" << endl;
		cout << endl;
		
		return 1;
	}
	//------------------------------------------------------------------------------------------------------------------	
	multiply(argv[1], argv[2], argv[3]);
	//------------------------------------------------------------------------------------------------------------------
	clock_t t2 = clock(); 	
	
	cout << "-------------------------" << endl;
	cout << "Total time: " << (t2 - t1) / double(CLOCKS_PER_SEC) << endl;
	//------------------------------------------------------------------------------------------------------------------
	return 0;
}
