/**
 * @brief Shortcuts to handle vectors
 * 
 * @file vectorops.cpp
 * @author Cristobal Rodero
 * @date 2018
 */

#include "../headers/vectorops.h"

/**
 * @brief Sum between two vectors.
 * 
 * @param v1 (<CODE>std::vector<double></CODE>) First vector to sum.
 * @param v2 (<CODE>std::vector<double></CODE>) Second vector to sum.
 * @return <CODE>std::vector<double></CODE> Result vector.
 */
std::vector<double> operator+(const std::vector<double> &v1,const std::vector<double> &v2){
	
	std::vector<double> res(v1.size());
	
	for(short int i = 0; i<v1.size(); i++){
		res[i] = v1[i]+v2[i];
	}
		
	return res;
}

/**
 * @brief Substraction between two vectors.
 * 
 * @param v1 (<CODE>std::vector<double></CODE>) Vector to substract from.
 * @param v2 (<CODE>std::vector<double></CODE>) Vector to substract.
 * @return <CODE>std::vector<double></CODE> Result vector.
 */
std::vector<double> operator-(const std::vector<double> &v1,const std::vector<double> &v2){

	std::vector<double> res(v1.size());
	
	for(short int i = 0; i<v1.size(); i++){
		res[i] = v1[i]-v2[i];
	}
		
	return res;
}
/**
 * @brief Product of a vector with an scalar.
 * 
 * @param alpha (<CODE>double</CODE>) Scalar.
 * @param v1 (<CODE>std::vector<double></CODE>) Vector.
 * @return <CODE>std::vector<double></CODE> Result vector.
 */
std::vector<double> operator*(const double &alpha, const std::vector<double> &v1){
	
	std::vector<double> res(v1.size());
	
	#pragma omp parallel
	for(short int i = 0; i<v1.size(); i++)
		res[i] = alpha*v1[i];
		
	return res;
}
/**
 * @brief Product of a vector with an scalar.
 * 
 * @param v1 (<CODE>std::vector<double></CODE>) Vector.
 * @param alpha (<CODE>double</CODE>) Scalar.
 * @return <CODE>std::vector<double></CODE> Result vector.
 */
std::vector<double> operator*(const std::vector<double> &v1, const double &alpha){
	return alpha*v1;
}
/**
 * @brief Division of a vector by an scalar.
 * 
 * @param v1 (<CODE>std::vector<double></CODE>) Vector.
 * @param alpha (<CODE>double</CODE>) Scalar.
 * @return <CODE>std::vector<double></CODE> Result vector.
 */
std::vector<double> operator/(const std::vector<double> &v1, const double &alpha){
	return v1*(1/alpha);
}
/**
 * @brief Euclidean norm of a vector.
 * 
 * @param v (<CODE>std::vector<double></CODE>) Vector.
 * @return (<CODE>double</CODE>) Norm of the vector.
 */
double norm(const std::vector<double> &v){

	double res = 0;
	
	for(short int i = 0; i<v.size(); i++){
		res += pow(v[i],2);
	}
	
	res = sqrt(res);
	
	return res;
}
/**
 * @brief Extract a vector from a <CODE>points</CODE> structure.
 * 
 * @param pts (<CODE>points</CODE>) Points structure.
 * @param pos (<CODE>int</CODE>) Component of the structure to extract.
 * @return (<CODE>std::vector<double></CODE>) Extracted vector. 
 */
std::vector<double> pts2vec(const points &pts, const int &pos){
	
	std::vector<double> res(3);
	
	res[0] = pts.x[pos];
	res[1] = pts.y[pos];
	res[2] = pts.z[pos];
	
	return res;
}
/**
 * @brief Absolute value component-wise of a vector.
 * 
 * @param v (<CODE>std::vector<double></CODE>) Vector.
 * @return (<CODE>std::vector<double></CODE>) Vector with positive components.
 */
std::vector<double> abs(const std::vector<double> &v){

	std::vector<double> res(v.size());
	
	for(short int i = 0; i < v.size(); i++){
		res[i] = abs(v[i]);
	}
	
	return res;
}
/**
 * @brief Computes where the minimum value of a vector is.
 * 
 * @param v (<CODE>std::vector<double></CODE>) Vector.
 * @return (<CODE>int</CODE>) Index of the minimum value of the vector.
 */
int min_pos(const std::vector<double> &v){
    int minPos = 0;
    
    for(int i = 0; i < v.size(); i++){
        if (v[i] < v[minPos]) 
            minPos = i;
    }
    
    return minPos;    
}
/**
 * @brief Computes the factorial of an integer.
 * 
 * @param x (<CODE>int</CODE>) Positive integer.
 * @return (<CODE>int</CODE>) Factorial of <VAR>x</CODE>.
 */
int factorial(const int &x){
	if(x == 1)
		return 1;
	else
		return x*factorial(x-1);
}
/**
 * @brief Finds the maximum value of a vector.
 * 
 * @param v (<CODE>std::vector<double></CODE>) Vector.
 * @return (<CODE>double</CODE>) Maximum value of the vector.
 */
double max_value_vec(const std::vector<double> &v){
	auto max_iterator = std::max_element(v.begin(),v.end());

	return *max_iterator;
}
/**
 * @brief Finds the minimum <B>positive</B> value of a vector. If all the components are negative, returns -1.
 * 
 * @param v (<CODE>std::vector<double></CODE>) Vector.
 * @return (<CODE>double</CODE>) Minimum positive value. -1 if all the components are negative.
 */
double min_positive(const std::vector<double> &v){
	double res;

	for(int i=0; i < v.size(); i++){
		if(v[i] > 0){
			res = v[i];
			break;
		}
	}

	for(int i = 0; i < v.size(); i++){
		if(v[i] > 0 && v[i] < res)
			res = v[i];
    }

	return res;
}
/**
 * @brief Computes a vector formed by the minimum value of two vectors component-wise.
 * 
 * @param v1 (<CODE>std::vector<double></CODE>) First vector.
 * @param v2 (<CODE>std::vector<double></CODE>) Second vector.
 * @return (<CODE>std::vector<double></CODE>) Vector whose <VAR>i</VAR>-th component is the minimum between <CODE>v1[i]</CODE>
 * and <CODE>v2[i]</CODE>.
 */
std::vector<double> min_vec(const std::vector<double> &v1, const std::vector<double> &v2){
	std::vector<double> res(v1.size());
	
	#pragma omp parallel
	for(int i=0; i < v1.size(); i++)
		res[i] = std::min(v1.at(i),v2.at(i));
	
	return res;
}
/**
 * @brief Computes a vector with each component being the minimum of the positive components between two vectors.
 * 
 * @bug If a component is negative sometimes takes the component of the first vector and sometimes from the second one.
 * 
 * @param v1 (<CODE>std::vector<double></CODE>) First vector.
 * @param v2 (<CODE>std::vector<double></CODE>) Second vector.
 * @return (<CODE>std::vector<double></CODE>) Result vector. 
 */
std::vector<double> min_positive(const std::vector<double> &v1, const std::vector<double> &v2){
	
	std::vector<double> res = v1;

	for(int i = 0; i < v1.size(); i++){
        if(v1[i] < 0 || (v2[i] < v1[i] && v2[i] >= 0))
			res[i] = v2[i];
    }

	return res;
}
/**
 * @brief Computes the cross product between two 3D vectors.
 * 
 * @param u (<CODE>std::vector<double></CODE>) First vector.
 * @param v (<CODE>std::vector<double></CODE>) Second vector.
 * @return std::vector<double> (<CODE>std::vector<double></CODE>) Result vector.
 */
std::vector<double> cross(const std::vector<double> &u,const std::vector<double> &v){
	
	std::vector<double> res(3);

	res.at(0) = u.at(1)*v.at(2) - u.at(2)*v.at(1);
	res.at(1) = u.at(2)*v.at(0) - u.at(0)*v.at(2);
	res.at(2) = u.at(0)*v.at(1) - u.at(1)*v.at(0);

	return res;

}
/**
 * @brief Computes the sum of the components of a vector.
 * 
 * @param v (<CODE>std::vector<double></CODE>) Vector.
 * @return (<CODE>std::vector<double></CODE>) Sum of the components of the vector.
 */
double sum(const std::vector<double> &v){

	double res = 0;
	#pragma parallel for
	for(int i = 0; i < v.size(); i++)
		res += v.at(i);

	return res;
}
/**
 * @brief Computes the dot product between two 3D vectors.
 * 
 * @param u (<CODE>std::vector<double></CODE>) First vector.
 * @param v (<CODE>std::vector<double></CODE>) Second vector.
 * @return std::vector<double> (<CODE>std::vector<double></CODE>) Result vector.
 */
double dot(const std::vector<double> &u, const std::vector<double> &v){

	return u[0]*v[0]+u[1]*v[1]+u[2]*v[2];
}

double average(const std::vector<double> &v){
	double res;

	res = sum(v);
	res /= v.size();

	return res;
}