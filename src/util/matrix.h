/*
	Dense Matrix and Vectors

	Author:   Steffen Rendle, http://www.libfm.org/
	modified: 2012-06-08

	Copyright 2010-2012 Steffen Rendle, see license.txt for more information
*/

#ifndef MATRIX_H_
#define MATRIX_H_

#include <vector>
#include <assert.h>
#include <iostream>
#include <fstream>
#include "../util/memory.h"
#include "../util/random.h"

using namespace std;
const unsigned DVECTOR_EXPECTED_FILE_ID = 1;
const unsigned DMATRIX_EXPECTED_FILE_ID = 1001;

struct dmatrix_file_header {
	unsigned id;
	unsigned type_size;
	unsigned num_rows;
	unsigned num_cols;
}; 

template <typename T> class DMatrix {
	public:
        vector<T> i_value;
		std::vector<std::string> col_names;
		unsigned dim1, dim2;
		
		T get(unsigned x, unsigned y) {
   			//assert((x < dim1) && (y < dim2));
			return i_value[x*dim2 + y];	
		}

		DMatrix(unsigned p_dim1, unsigned p_dim2) {
			dim1 = 0;
			dim2 = 0;	
			setSize(p_dim1, p_dim2);
		}
		
		DMatrix() {
			dim1 = 0;
			dim2 = 0;	
		}
		
		~DMatrix() {
		}
		
		void init(const T& v) {
            i_value.assign(dim1*dim2, v);
   		}
		void setSize(unsigned p_dim1, unsigned p_dim2) {
			if ((p_dim1 == dim1) && (p_dim2 == dim2)) {
				return;
			}
			dim1 = p_dim1;
			dim2 = p_dim2;
            i_value.reserve(dim1*dim2);
			col_names.assign(dim2, "");					
		}

		T& operator() (unsigned x, unsigned y) {
   		//	assert((x < dim1) && (y < dim2));
			return i_value[x*dim2 + y];	
		}
   		T operator() (unsigned x, unsigned y) const {
   		//	assert((x < dim1) && (y < dim2));
   			return i_value[x*dim2 + y];	
   		}
   		
   		void save(std::string filename, bool has_header = false) {
		   	std::ofstream out_file (filename.c_str());
			if (out_file.is_open())	{
				if (has_header) {
					for (unsigned i_2 = 0; i_2 < dim2; i_2++) {
						if (i_2 > 0) {
							out_file << "\t";
						}
						out_file << col_names[i_2];					
					}	
					out_file << std::endl;
				}
				for (unsigned i_1 = 0; i_1 < dim1; i_1++) {
					for (unsigned i_2 = 0; i_2 < dim2; i_2++) {
						if (i_2 > 0) {
							out_file << "\t";
						}
						out_file << value[i_1][i_2];
					}
					out_file << std::endl;
				}
				out_file.close();
			} else {
				std::cout << "Unable to open file " << filename;
			}   			
   		}

		void saveToBinaryFile(std::string filename) {
			std::cout << "writing to " << filename << std::endl; std::cout.flush();
			std::ofstream out(filename.c_str(), std::ios_base::out | std::ios_base::binary);
			if (out.is_open()) {
				dmatrix_file_header fh;
				fh.id = DMATRIX_EXPECTED_FILE_ID;
				fh.num_rows = dim1;
				fh.num_cols = dim2;
				fh.type_size = sizeof(T);
				out.write(reinterpret_cast<char*>(&fh), sizeof(fh));
				for (unsigned i = 0; i < dim1; i++) {
					out.write(reinterpret_cast<char*>(value[i]), sizeof(T)*dim2);
				}
				out.close();
			} else {
				throw "could not open " + filename;
			}
		}

		void loadFromBinaryFile(std::string filename) {
			std::cout << "reading " << filename << std::endl; std::cout.flush();
			std::ifstream in(filename.c_str(), std::ios_base::in | std::ios_base::binary);
			if (in.is_open()) {
				dmatrix_file_header fh;
				in.read(reinterpret_cast<char*>(&fh), sizeof(fh));
				assert(fh.id == DMATRIX_EXPECTED_FILE_ID);
				assert(fh.type_size == sizeof(T));
				setSize(fh.num_rows, fh.num_cols);
				for (unsigned i = 0; i < dim1; i++) {
					in.read(reinterpret_cast<char*>(value[i]), sizeof(T)*dim2);
				}
				in.close();
			} else {
				throw "could not open " + filename;
			}
		}

		void load(std::string filename) {
			std::ifstream in_file (filename.c_str());
			if (! in_file.is_open()) {
				throw "Unable to open file " + filename;
			}	
			for (unsigned i_1 = 0; i_1 < dim1; i_1++) {
				for (unsigned i_2 = 0; i_2 < dim2; i_2++) {
					T v;
					in_file >> v;
					i_value[i_1*dim2 + i_2] = v;
				}
			}
			in_file.close();
		} 		
   		
};

template <typename T> class DVector {
	public:
		unsigned dim;
        vector<T> i_value;
		DVector() {
			dim = 0;
		}
		DVector(unsigned p_dim) {
			dim = 0;
			setSize(p_dim);
		}
		~DVector() {
				vector<T> tmp_vec_2_swap;
                i_value.swap(tmp_vec_2_swap);
		}
		T get(unsigned x) {
			return i_value[x];
		}
		void setSize(unsigned p_dim) {
            dim = p_dim;
            i_value.reserve(p_dim);		
		}
		T& operator[] (unsigned x) {
            return i_value[x];
		}

        T operator() (unsigned x) const {
            return i_value[x];	
        }

   		void init(const T &v) {
            i_value.assign(dim, v);
   		}

   		void save(std::string filename) {
		   	std::ofstream out_file (filename.c_str());
			if (out_file.is_open())	{
				for (auto itr: i_value) {
					out_file << itr << std::endl;
				}
				out_file.close();
			} else {
				std::cout << "Unable to open file " << filename;
			}   			
   		}

		void saveToBinaryFile(std::string filename) {
		   	std::ofstream out (filename.c_str(), std::ios_base::out | std::ios_base::binary);
			if (out.is_open())	{
				unsigned file_version = DVECTOR_EXPECTED_FILE_ID;
				unsigned data_size = sizeof(T);
				unsigned num_rows = dim;
				out.write(reinterpret_cast<char*>(&file_version), sizeof(file_version));
				out.write(reinterpret_cast<char*>(&data_size), sizeof(data_size));
				out.write(reinterpret_cast<char*>(&num_rows), sizeof(num_rows));
				out.write(reinterpret_cast<char*>(value), sizeof(T)*dim);
				out.close();
			} else {
				std::cout << "Unable to open file " << filename;
			}   			
   		}


		void load(std::string filename) {
			std::ifstream in_file (filename.c_str());
			if (! in_file.is_open()) {
				throw "Unable to open file " + filename;
			}	
			for (unsigned i = 0; i < dim; i++) {
				T v;
				in_file >> v;
				i_value[i] = v;
			}
			in_file.close();
		}


		void loadFromBinaryFile(std::string filename) {
			std::ifstream in (filename.c_str(), std::ios_base::in | std::ios_base::binary);
			if (in.is_open())	{
				unsigned file_version;
				unsigned data_size;
				unsigned num_rows;
				in.read(reinterpret_cast<char*>(&file_version), sizeof(file_version));
				in.read(reinterpret_cast<char*>(&data_size), sizeof(data_size));
				in.read(reinterpret_cast<char*>(&num_rows), sizeof(num_rows));
				assert(file_version == DVECTOR_EXPECTED_FILE_ID);
				assert(data_size == sizeof(T));
				T* value = new T(num_rows);
				in.read(reinterpret_cast<char*>(value), sizeof(T)*dim);
				in.close();
                i_value.assign(value, value + num_rows);
                delete[] value;
			} else {
				std::cout << "Unable to open file " << filename;
			}   			
		}
};


class DVectorDouble : public DVector<double> {
	public:
		void init_normal(double mean, double stdev) {	
			for (unsigned i_2 = 0; i_2 < dim; i_2++) {
				i_value[i_2] = ran_gaussian(mean, stdev);
			}
		}
};

class DMatrixDouble : public DMatrix<double> {
	public:
		void init(double mean, double stdev) {	
			for (unsigned i_1 = 0; i_1 < dim1; i_1++) {
				for (unsigned i_2 = 0; i_2 < dim2; i_2++) {
					i_value[i_1*dim2 + i_2] = ran_gaussian(mean, stdev);
				}
			}
		}
		void init_column(double mean, double stdev, int column) {	
			for (unsigned i_1 = 0; i_1 < dim1; i_1++) {
				i_value[i_1*dim2 + column] = ran_gaussian(mean, stdev);
			}
		}
};


#endif /*MATRIX_H_*/
