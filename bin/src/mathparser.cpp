/******************************************************************

Introduction:
mathparser is a simple c++ program to parse math expressions.

The program is a modified version of math expression parser 
presented in the book : "C++ The Complete Reference" by H.Schildt.

It supports operators: + - * / ^ ( )

It supports math functions : sin, cos, tan, asin, acos, atan, sinh, 
cosh, tanh, asinh, acosh, atanh, ln, log, exp, sqrt, sqr, round, int.

It supports variables A to Z, a-z.

Sample:
25 * 3 + 1.5*(-2 ^ 4 * log(30) / 3)
x = 3
y = 4
r = sqrt(x ^ 2 + y ^ 2)
t = atan(y / x)

mathparser version 1.0 by Hamid Soltani. (gmail: hsoltanim)
Last modified: Aug. 2016.

*******************************************************************/

//#include "stdafx.h"
#include "mathparser.hpp"

char *chngChar (char *str, char oldChar, char newChar) {
    char *strPtr = str;
    while ((strPtr = strchr (strPtr, oldChar)) != NULL)
        *strPtr++ = newChar;
    return str;
}
// Parser constructor.
parser::parser()
{
	int i;
	exp_ptr = NULL;
//	for (i = 0; i < NUMVARS; i++)vars[i] = 0.0;
	errormsg[0] = '\0';
}
// Parser entry point.
double parser::eval_exp(char *exp)
{//cout << "Expression:" << exp<< "\n";
	errormsg[0] = '\0';
	double result;
	exp_ptr = exp;
while (isspace(*exp_ptr))  // skip over white space at beginning
		++exp_ptr; 
	get_token();
         if (*errormsg){
			cout << "Parse Variable Error: expression " << exp << "\n" << errormsg << "\n\n";
                        exit(EXIT_FAILURE);
                       }
	if (!*token) 
	{
		strcpy(errormsg, "No Expression Present"); // no expression present
		return (double)0;
	}
	eval_exp1(result);
	if (*token) // last token must be null or space or newline or =
		if(!isspace(*token)&&!strchr("=)\n\t\r",*token))sprintf(errormsg, " Last char %c not space or newline",*token);
         if (*errormsg){
			cout << "Parse Variable Error: expression " << exp << "\n" << errormsg << "\n\n";
                        exit(EXIT_FAILURE);
                       }
       //  else {cout  << "result: " << result << "\n"; } // for debug
	return result;
}
// Process an assignment.
void parser::eval_exp1(double &result)
{
	int slot;
	char temp_token[80];
	if (tok_type == VARIABLE) 
	{
		// save old token
		char *t_ptr = exp_ptr;
		strcpy(temp_token, token);
		// compute the index of the variable
//		slot = *token - 'A';
                std::string tstr(token);
		get_token();
		if (*token != '=') 
		{
			exp_ptr = t_ptr; // return current token
			strcpy(token, temp_token); // restore old token
			tok_type = VARIABLE;
		}
		else {
			get_token(); // get next part of exp
			eval_exp2(result);
                        vars[tstr]=result;
			return;
		}
	}
	eval_exp2(result);
}
// Add or subtract two terms.
void parser::eval_exp2(double &result)
{
	register char op;
	double temp;
	eval_exp3(result);
	while ((op = *token) == '+' || op == '-')
	{
		get_token();
		eval_exp3(temp);
		switch (op) 
		{
		case '-':
			result = result - temp;
			break;
		case '+':
			result = result + temp;
			break;
		}
	}
}
// Multiply or divide two factors.
void parser::eval_exp3(double &result)
{
	register char op;
	double temp;
	eval_exp4(result);
	while ((op = *token) == '*' || op == '/') 
	{
		get_token();
		eval_exp4(temp);
		switch (op) 
		{
		case '*':
			result = result * temp;
			break;
		case '/':
			result = result / temp;
			break;
		}
	}
}
// Process an exponent.
void parser::eval_exp4(double &result)
{
	double temp;
	eval_exp5(result);
	while (*token == '^')
	{
		get_token();
		eval_exp5(temp);
		result = pow(result, temp);
	}
}
// Evaluate a unary + or -.
void parser::eval_exp5(double &result)
{
	register char op;
	op = 0;
	if ((tok_type == DELIMITER) && *token == '+' || *token == '-')
	{
		op = *token;
		get_token();
	}
	eval_exp6(result);
	if (op == '-')
		result = -result;
}
// Process a function, a parenthesized expression, a value or a variable
void parser::eval_exp6(double &result)
{
	bool isfunc = (tok_type == FUNCTION);
	char temp_token[80];
	if (isfunc)
	{
		strcpy(temp_token, token);
		get_token();
	} 
	if ((*token == '(')) 
	{
		get_token();
		eval_exp2(result);
		if (*token != ')')
			strcpy(errormsg, "Unbalanced Parentheses");
		if (isfunc)
		{
			if (!strcmp(temp_token, "sin"))
				result = sin(PI / 180 * result);
			else if (!strcmp(temp_token, "cos"))
				result = cos(PI / 180 * result);
			else if (!strcmp(temp_token, "tan"))
				result = tan(PI / 180 * result);
			else if (!strcmp(temp_token, "asin"))
				result = 180 / PI*asin(result);
			else if (!strcmp(temp_token, "acos"))
				result = 180 / PI*acos(result);
			else if (!strcmp(temp_token, "atan"))
				result = 180 / PI*atan(result);
			else if (!strcmp(temp_token, "sinh"))
				result = sinh(result);
			else if (!strcmp(temp_token, "cosh"))
				result = cosh(result);
			else if (!strcmp(temp_token, "tanh"))
				result = tanh(result);
			else if (!strcmp(temp_token, "asinh"))
				result = asinh(result);
			else if (!strcmp(temp_token, "acosh"))
				result = acosh(result);
			else if (!strcmp(temp_token, "atanh"))
				result = atanh(result);
			else if (!strcmp(temp_token, "ln"))
				result = log(result);
			else if (!strcmp(temp_token, "log"))
				result = log10(result);
			else if (!strcmp(temp_token, "exp"))
				result = exp(result);
			else if (!strcmp(temp_token, "sqrt"))
				result = sqrt(result);
			else if (!strcmp(temp_token, "sqr"))
				result = result*result;
			else if (!strcmp(temp_token, "round"))
				result = round(result);
			else if (!strcmp(temp_token, "int"))
				result = floor(result);
			else
				strcpy(errormsg, "Unknown Function");
		}
		get_token();
	}
	else    {std::string tstr(token);
                 switch (tok_type)
		{
		case VARIABLE:
			result = vars[tstr];
			get_token();
			return;
		case NUMBER:
			chngChar (token, 'D', 'E'); // for exponential number format with D and d ...
                       chngChar (token, 'd', 'E');
	               result = atof(token);
			get_token();
			return;
		default: ;
			//strcpy(errormsg, "Syntax Error - not Variable or number");
		}
               }
}
// Obtain the next token.
void parser::get_token()
{
	register char *temp;
	tok_type = 0;
	temp = token;
	*temp = '\0';
	if (!*exp_ptr||isspace(*exp_ptr)) return;// at end of expression
         
        if (strchr("+-*/%^()=", *exp_ptr)) 
	{
		tok_type = DELIMITER;
          if(*exp_ptr=='=') {
          *temp++ = *exp_ptr++;  // advance to next char
           while (isspace(*exp_ptr))  // skip over white space
          *temp++ = *exp_ptr++;  // advance to next char
             }                  
            else 
            {
                 	*temp++ = *exp_ptr++;  // advance to next char
            }
	}
	else if (isalpha(*exp_ptr)) 
	{

		while (!strchr(" +-/*%^=()\t\r\n", *exp_ptr) && (*exp_ptr)) // advance until end of alphabetic characters
//			*temp++ = toupper(*exp_ptr++);
			*temp++ = *exp_ptr++;
		tok_type = (*exp_ptr == '(') ? FUNCTION : VARIABLE;
	}
	else if (isdigit(*exp_ptr) || *exp_ptr == '.')
	{while (!strchr(" +-/*%^=()\t\rEeDd", *exp_ptr) && (*exp_ptr))
//			*temp++ = toupper(*exp_ptr++);
			*temp++ = *exp_ptr++;
             if(strchr("EeDd", *exp_ptr)&& (*exp_ptr))  // for exponential number format detection
             {
//                      *temp++ = toupper(*exp_ptr++);
			*temp++ = *exp_ptr++;
              if(strchr("+-", *exp_ptr) && (*exp_ptr))
//                      *temp++ = toupper(*exp_ptr++);
			*temp++ = *exp_ptr++;
               while (!strchr(" +-/*%^=()\t\r", *exp_ptr) && (*exp_ptr))
//			*temp++ = toupper(*exp_ptr++);
			*temp++ = *exp_ptr++;
             }
		
		tok_type = NUMBER;
	}
	*temp = '\0';
	//if ((tok_type == VARIABLE) && (token[1]))
	//	sprintf(errormsg, "Variable %c too long - only variables with one letter possible",*token);
}
/*
int main()
{
	char expstr[256];
	parser ob;
	cout << "Math expression parser. Enter a blank line to stop.\n\n";
	do
	{
		cout << "Enter expression: ";
		cin.getline(expstr, 255);
		double ans = ob.eval_exp(expstr);
		if (*ob.errormsg)
			cout << "Error: " << ob.errormsg << "\n\n";
		else
			cout << "Answer: " << ans << "\n\n";
	} while (*expstr);
	return 0;
}
*/
