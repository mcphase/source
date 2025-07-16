/******************************************************************

Introduction:
mathparser is a simple c++ program to parse math expressions.

The program is a modified version of math expression parser 
presented in the book : "C++ The Complete Reference" by H.Schildt.

It supports operators: + - * / ^ ( )

It supports math functions : SIN, COS, TAN, ASIN, ACOS, ATAN, SINH, 
COSH, TANH, ASINH, ACOSH, ATANH, LN, LOG, EXP, SQRT, SQR, ROUND, INT.

It supports variables A to Z.

Sample:
25 * 3 + 1.5*(-2 ^ 4 * log(30) / 3)
x = 3
y = 4
r = sqrt(x ^ 2 + y ^ 2)
t = atan(y / x)

mathparser version 1.0 by Hamid Soltani. (gmail: hsoltanim)
Last modified: Aug. 2016.

*******************************************************************/
#ifndef MATHPARSE
#define MATHPARSE
//#include "stdafx.h"
#include "mathparser.hpp"
#include <iostream>
#include <cstdlib>
#include <cctype>
#include <cstring>
#include <math.h> 
#include <map>
#ifndef PI
#define PI 3.14159265358979323846 
#endif

using namespace std;

enum types { DELIMITER = 1, VARIABLE, NUMBER, FUNCTION };
//const int NUMVARS = 58; 
class parser {
	char *exp_ptr; // points to the expression
	char token[256]; // holds current token
	char tok_type; // holds token's type
        std::map<std::string, double> vars;
//	double vars[NUMVARS]; // holds variable's values
	void eval_exp1(double &result);
	void eval_exp2(double &result);
	void eval_exp3(double &result);
	void eval_exp4(double &result);
	void eval_exp5(double &result);
	void eval_exp6(double &result);
	void get_token();
public:
	parser();
	double eval_exp(char *exp);
	char errormsg[64];
};
#endif