#ifndef MATRIX_H
#define MATRIX_H 

#include <typeinfo>
#include <cstring>
#include <limits>
#include <cmath>
#include <iomanip>

using namespace std;

template <typename T> 
class Matrix{
public:
	//BEGIN-------------------------------------------------CONSTRUCTOR-------------------------------------------------
	Matrix<T>(){
		m = -1;
		n = -1;
		
		type = -1;
	}   
	
	Matrix<T>(long _m, long _n){
		m = _m;
		n = _n;
		
		set_type();      
		
		alloc_memory();
		empty();
	}
	
	Matrix<T>(const Matrix<T> &tmp){
		(*this) = Matrix<T>(tmp.m, tmp.n);
		
		for(long i = 0; i < m; i++){
			for(long j = 0; j < n; j++){
				A[i][j] = tmp.A[i][j];            
			}      
		}
	}
	//END---------------------------------------------------CONSTRUCTOR-------------------------------------------------
	
	void alloc_memory();
	
	void empty();
	
	Matrix<T> operator*(Matrix<T>& tmp);  
	
	//BEGIN-------------------------------------------------SET---------------------------------------------------------
	void set_type();
	//END---------------------------------------------------SET---------------------------------------------------------
		
	//BEGIN-------------------------------------------------GET---------------------------------------------------------
	long get_m() const{ return m; }
	long get_n() const{ return n; }
	   
	int get_type() const{ return type; }
	//END---------------------------------------------------GET---------------------------------------------------------
	   
	//BEGIN-------------------------------------------------FILE READ/WRITE---------------------------------------------
	void read_from_file(const char* filename);
	
	void write_to_file(const char* filename);
	//END---------------------------------------------------FILE READ/WRITE---------------------------------------------
		
	//BEGIN-------------------------------------------------PRINT-------------------------------------------------------
	void print() const;
	//END---------------------------------------------------PRINT-------------------------------------------------------
	
	//BEGIN-------------------------------------------------DESTRUCTOR--------------------------------------------------
	~Matrix<T>(){   
		if((m > 0) && (n > 0)){
			for(int i = 0; i < m; i++){
				delete[] A[i];      
			}
				
			delete[] A;
			
			A = NULL;
		}
	}
	//END---------------------------------------------------DESTRUCTOR--------------------------------------------------
	
	T** A;
	
private:
	//-----------
	long m, n;
				 
	int type;
	//-----------
};
#endif
