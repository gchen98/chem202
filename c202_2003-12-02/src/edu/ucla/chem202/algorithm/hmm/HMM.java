package edu.ucla.chem202.algorithm.hmm;

import edu.ucla.chem202.DataFormatException;

public interface HMM {
	
	public static final byte UNDEF = 0;
	public static final byte MATCH = 1;
	public static final byte DELETION = 2;
	public static final byte INSERTION = 3;
	public static final byte BEGIN = 4;
	public static final byte END = 5;
	
	/**
	* returns the length of all sequences in the alignment profile
	*/

	public int getLength();
	
	/**
	* emission probabilities
	*
	* @param state is any of the state constants
	* @param residue is the residue interested
	* @param column is column # of the profile indexed from 0 to length-1
	* @returns the emission probability for that particular state
	*/
	
	public double getEmissionProb(byte state, char residue, int column)
	throws DataFormatException;
	
	/**
	* transition probabilities
	*
	* @param fromState is any of the state constants
	* @param toState is any of the state constants	
	* @param toCol is the destination column of the transition edge indexed from 0 to length-1
	* @returns the transition probability for that particular state pair
	*/	
	
	public double getTransitionProb(byte fromState, byte toState, int toCol)
	throws DataFormatException;
		
	// add/modify any methods you need in this interface.  email me
	// when you are done and I will adjust the implementation.
	
}
