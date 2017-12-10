

/***********************************************************************
 *
 *	maxmin .h
 *
 *	return maximum value, minimum value.
 *	
 *
 *	2015 July 14 : NOT using "if" sentence in max function.
 *  2012 Oct. 08 : coded by Sho Nakamura (Tohoku Univ).
 *
 **********************************************************************/


double max(double variable1, double variable2)
{
	if(variable1>variable2) return variable1;
	else return variable2;
}

double min(double variable1, double variable2)
{
	if(variable1>variable2) return variable2;
	else return variable1;
}





