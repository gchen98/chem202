package edu.ucla.chem202.algorithm.hmm;

import edu.ucla.chem202.io.Logger;
import edu.ucla.chem202.algorithm.*;
import edu.ucla.chem202.DataFormatException;
import java.io.*;
import java.util.List;
import java.util.ArrayList;

/**
* Encapsulates finished alignment. Should be able to marshall out into various formats
*/

public class HMMAlignment
extends Alignment
implements HMM{
	
	//private byte[] colStates;
	private List residueDefs = new ArrayList();
	private double[][] matchProbs;
	
	//private int totalTransitionMoves;
	//private final short TOTAL_MOVES = 3;
	
	/**
	*  The input stream should be a fileinput stream for instance
	* that will read in a residue (DNA or protein or anything else) definition file
	*/
	
	public HMMAlignment(List sequences, InputStream in)
	throws DataFormatException, IOException{
		super(sequences);
		BufferedReader bufRead=new BufferedReader(new InputStreamReader(in));		
		String newLine = null;
		while ((newLine = bufRead.readLine())!=null) {			
			residueDefs.add(new Character(newLine.charAt(0)));	// add permitted residues to our list
		}
		createProbModel();
	}
	
	/**
	* Compute probabilities at each column of alignment
	*/
	
	private void createProbModel(){
		matchProbs = new double[residueDefs.size()][seqLength];
		int numSequences = sequences.size();
		for (int colCount = 0; colCount<seqLength; colCount++) {
			int gapCount = 0;
			for (int rowCount = 0; rowCount<sequences.size(); rowCount++ ) { // count gaps
				char observedChar = getCharAt(rowCount,colCount);
				if (observedChar == GAP_CHAR) {
					gapCount++;
				}
			}
			for (int resCount = 0;resCount<residueDefs.size();resCount++) {  // count each res
				char curResidue = ((Character)residueDefs.get(resCount)).charValue();
				int curResidueCount = 1;				
				for (int rowCount = 0; rowCount<sequences.size(); rowCount++ ) {
					char observedChar = getCharAt(rowCount,colCount);
					if (observedChar == curResidue) {
						curResidueCount++;
					}
				}
				double residueProb = (double)curResidueCount/
				(sequences.size()+residueDefs.size()-gapCount);
				matchProbs[resCount][colCount] = residueProb;
				//Logger.debugln("column "+colCount+" residue "+curResidue+" has prob "+residueProb);
			}			
		}		
	}
	
	private double getMatchEmission(char residue, int column)
	throws DataFormatException {
		int charCol = column - 1;
		if (!residueDefs.contains(new Character(residue))) {
			throw new DataFormatException("The residue "+residue+" was not defined in the residue file.");
		}
		if (charCol == -1) {
			throw new DataFormatException("The first column of the profile must be 1");
		}
		int resIndex = residueDefs.indexOf(new Character(residue));
		return matchProbs[resIndex][charCol];
	}	
	
	public double getEmissionProb(byte state, char residue, int column)
	throws DataFormatException {
		switch (state) {
		case MATCH:
			return getMatchEmission(residue, column);			
		case INSERTION:
			return getMatchEmission(residue, column);			
		case DELETION:
			return 0;			
		default:
			throw new DataFormatException("Invalid state queried");
		}		
	}
	
	public int getLength() {
		return seqLength;
	}
	
	private int getResidueCounts(int col) {
		int residueCount = 1;
		for (int rowCount = 0; rowCount<sequences.size(); rowCount++ ) { // count gaps
			char observedChar = getCharAt(rowCount,col);
			if (observedChar != GAP_CHAR) {
				residueCount++;
			}
		}
		return residueCount;
	}
	
	private int getResidueCounts(int col, char ch) {
		int residueCount = 1;
		for (int rowCount = 0; rowCount<sequences.size(); rowCount++ ) { // count gaps
			char observedChar = getCharAt(rowCount,col);
			if (observedChar == ch) {
				residueCount++;
			}
		}
		return residueCount;
	}
	
	private boolean isCharGap(char ch) {
		return ch==GAP_CHAR?true:false;
	}
	
	private int getCol2ColCounts(boolean isFirstGap, boolean isSecondGap, int startColumn) {
		int c2cCount = 1;
		for (int rowCount = 0; rowCount<sequences.size(); rowCount++) {
			char ch1 = getCharAt(rowCount,startColumn);
			char ch2 = getCharAt(rowCount,startColumn+1);
			if ((isFirstGap==isCharGap(ch1)) && (isSecondGap==isCharGap(ch2))) {				
				c2cCount++;
			}
		}
		return c2cCount;
	}	
	
	private double getMatchInsertTransitionProb(byte toState, int startColumn) {
		int m2m = getCol2ColCounts(false,false,startColumn);
		int m2d = getCol2ColCounts(false,false,startColumn);
		int m2i = getResidueCounts(startColumn,GAP_CHAR);
		int total = m2m+m2d+m2i;
		switch (toState) {
		case MATCH:
			return (double)m2m/total;
		case DELETION:
			return (double)m2d/total;
		case INSERTION:
			return (double)m2i/total;		
		}
		return -1;		
	}
	
	private double getInsertionTransitionProb(byte toState, int startColumn) {
		int i2m = getCol2ColCounts(true,false,startColumn);
		int i2d = getCol2ColCounts(true,true,startColumn);
		int i2i = getResidueCounts(startColumn,GAP_CHAR);
		int total = i2m + i2d + i2i;
		switch (toState) {
		case MATCH:
			return (double)i2m/total;
		case DELETION:
			return (double)i2d/total;
		case INSERTION:
			return (double)i2i/total;		
		}
		return -1;
	}
	
	private double getDeletionTransitionProb(byte toState, int startColumn) {
		int d2m = getCol2ColCounts(true,false,startColumn);
		int d2d = getCol2ColCounts(true,true,startColumn);
		int d2i = getResidueCounts(startColumn,GAP_CHAR);
		int total = d2m + d2d + d2i;
		switch (toState) {
		case MATCH:
			return (double)d2m/total;
		case DELETION:
			return (double)d2d/total;
		case INSERTION:
			return (double)d2i/total;		
		}
		return -1;
	}
	
	public double getTransitionProb(byte fromState, byte toState, int startCol) {
		switch (fromState) {
		case BEGIN:
			return getMatchInsertTransitionProb(toState, startCol);
		case MATCH:
			return getMatchInsertTransitionProb(toState, startCol);
		case DELETION:
			return getDeletionTransitionProb(toState, startCol);
		case INSERTION:
			return getMatchInsertTransitionProb(toState, startCol);
		}
		return -1;
	}
	
	
	public double getMtoMProb(int startColumn){ // same as I to M		
		throw new java.lang.UnsupportedOperationException("Use getTranstionProb() instead");
	}	

	public double getMtoIProb(int startColumn){ // same as I to I
		throw new java.lang.UnsupportedOperationException("Use getTranstionProb() instead");
	}
	
	public double getDtoMProb(int startColumn){
		throw new java.lang.UnsupportedOperationException("Use getTranstionProb() instead");
	}
	
	// here 1 <= startColumn <= // seqLength, so return 0 when startColumn == 0	
	public double getDtoDProb(int startColumn){ 
		throw new java.lang.UnsupportedOperationException("Use getTranstionProb() instead");
	}
												
	public double getMtoDProb(int startColumn){  // same as I to D
		throw new java.lang.UnsupportedOperationException("Use getTranstionProb() instead");
	}

	public double getRandomStateProb(char residue, int column){	
		throw new java.lang.UnsupportedOperationException("How do we compute the random state?");
	}
	
	public double getDeleteProb(int column) {
		throw new java.lang.UnsupportedOperationException("Deletion emissions are silent?");
	}
	
	public double getInsertProb(int column) {
		throw new java.lang.UnsupportedOperationException("Use random distribution model?");
	}	
	
	public byte getState(int column) {
		throw new java.lang.UnsupportedOperationException("Do we need this?");
	}
	
	public void setState(int column, byte newState) {
		throw new java.lang.UnsupportedOperationException("Do we need this?");
	}	
	
	public double getMatchStateProb(char residue, int column)
	throws DataFormatException {
		throw new java.lang.UnsupportedOperationException("Do we need this?");
	}
}
