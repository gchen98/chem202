
/*
* Author: Teresa Breyer, Gary Chen
* Created: Thursday, November 20, 2003 8:18:56 PM
* Modified: Thursday, November 20, 2003 8:18:56 PM
*/
package edu.ucla.chem202.algorithm.hmm;


/**
* a utility class of useful math subroutines
*/

public class MathUtil
{

	// to make our debugging statements more readable.  Tune precision as you see fit.

	public static double round(double value) {
		return (double)Math.round(value*100)/100;
	}

	public static double addLog(double logProbA, double probA, double logProbB, 
	double probB, double logProbC, double probC) {
		if ( logProbA+Math.log(probA) >= logProbB + Math.log(probB) && logProbA+Math.log(probA) >= logProbC + Math.log(probC)) {
			return logProbA + Math.log(probA) + Math.log(1 + Math.exp(logProbB + Math.log(probB) - logProbA - Math.log(probA)) +
				Math.exp(logProbC + Math.log(probC) - logProbA - Math.log(probA)));
		} else if ( logProbC+Math.log(probC) >= logProbB + Math.log(probB) && logProbA+Math.log(probA) <= logProbC + Math.log(probC)) {
			return logProbB + Math.log(probB) + Math.log(1 + Math.exp(logProbA + Math.log(probA) - logProbB - Math.log(probB)) +
				Math.exp(logProbC + Math.log(probC) - logProbB - Math.log(probB)));
		} else {
		return logProbC + Math.log(probC) + Math.log(1 + Math.exp(logProbB + Math.log(probB) - logProbC - Math.log(probC)) +
				Math.exp(logProbA + Math.log(probA) - logProbC - Math.log(probC)));
		}
	}

	public static double addLog(double logProbA, double probA, double logProbB, double probB) {
		if (  logProbA+Math.log(probA) >= logProbB + Math.log(probB) ) {
			return logProbA + Math.log(probA) + Math.log(1 + Math.exp(logProbB + Math.log(probB) - logProbA - Math.log(probA)));
		} else {
			return logProbB + Math.log(probB) + Math.log(1 + Math.exp(logProbA + Math.log(probA) - logProbB - Math.log(probB)));
		}
	}

	/**
	* addLogLog: logProbA = logProbA + log(probA) from the addLog function
	*/
	public static double addLogLog(double logProbA, double logProbB) {
		if (  logProbA  >= logProbB )  {
			return logProbA + Math.log(1 + Math.exp(logProbB - logProbA));
		} else {
			return logProbB + Math.log(1 + Math.exp(logProbA - logProbB));
		}
	}

	public static double addLogLog(double logProbA, double logProbB, double 
logProbC) {
		if (  logProbA  >= logProbB && logProbA  >= logProbC ) {
			return logProbA + Math.log(1 + Math.exp(logProbB - logProbA)) +  
Math.log(1 + Math.exp(logProbC - logProbA));
		} else if (  logProbB  >= logProbA  &&  logProbB  >= logProbC ) {
			return logProbB + Math.log(1 + Math.exp(logProbA - logProbB)) +  
Math.log(1 + Math.exp(logProbC - logProbB));
		} else {
			return logProbC + Math.log(1 + Math.exp(logProbA - logProbC)) +  
Math.log(1 + Math.exp(logProbB - logProbC));
		}
	}

	public static double addLog(double logProbA, double probA) {
		return logProbA + Math.log(probA);
	}

	public static double getMaxPairSummed(double logProbA, double probA, double 
logProbB, double probB) {
		if (  logProbA+Math.log(probA) >= logProbB + Math.log(probB) ) {
			return logProbA + Math.log(probA);
		} else {
			return logProbB + Math.log(probB);
		}
	}

	public static double getMaxPairSummed(double logProbA, double probA, double 
logProbB, double probB, double logProbC, double probC) {
		double max = logProbA+Math.log(probA);
		if ( max < logProbB + Math.log(probB) ) max = logProbB + Math.log(probB);
		if ( max < logProbC + Math.log(probC) ) max = logProbC + Math.log(probC);
		return max;
	}
}
