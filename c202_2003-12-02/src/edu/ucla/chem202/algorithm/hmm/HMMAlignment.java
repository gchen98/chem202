package edu.ucla.chem202.algorithm.hmm;

import edu.ucla.chem202.io.Logger;
import edu.ucla.chem202.algorithm.*;
import edu.ucla.chem202.DataFormatException;
import java.io.*;
import java.util.List;
import java.util.ArrayList;

/**
* represents an HMM profile, probabilitic model based on a MSA.
* Author: Gary Chen
* Created: Wednesday, November 12, 2003 11:06:55 PM
* Modified: Saturday, November 15, 2003 11:06:55 PM
*/

public class HMMAlignment
extends Alignment
implements HMM{
	
	private List residueDefs = new ArrayList();
	private double[][] matchProbs;
	
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
		initProbModel();
	}
	
	public HMMAlignment(List sequences)
	throws DataFormatException{
		super(sequences);
	}
	
	/**
	* Compute probabilities at each column of alignment
	*/
	
	private void initProbModel(){
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
		if (!residueDefs.contains(new Character(residue))) {
			throw new DataFormatException("The residue "+residue+" was not defined in the residue file.");
		}
		if (column <0) {
			throw new DataFormatException("Invalid column of "+column+" queried. Column indices begin at 0.");
		}
		int resIndex = residueDefs.indexOf(new Character(residue));
		return matchProbs[resIndex][column];
	}	
	
	public double getEmissionProb(byte state, char residue, int column)
	throws DataFormatException {
		switch (state) {
		case MATCH:
			return getMatchEmission(residue, column);			
		case INSERTION:	
			return getMatchEmission(residue, column);// we are using the same state for insert
		case DELETION:
			return 0;
		default:
			throw new DataFormatException("Invalid state "+state+" queried");
		}		
	}
	
	public int getLength() {
		return seqLength;
	}
	
	private int getResidueCounts(int col) {
		if ((col<0) || (col>=seqLength)) return 0;
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
		if ((col<0) || (col>=seqLength)) return 0;
		int residueCount = 1;
		for (int rowCount = 0; rowCount<sequences.size(); rowCount++ ) { // count gaps
			char observedChar = getCharAt(rowCount,col);
			if (observedChar == ch) {
				residueCount++;
			}
		}
		return residueCount;
	}

	private double getDeletionTransitionProb(byte toState, int toCol)
	throws DataFormatException{
		int residues = getResidueCounts(toCol);
		int gaps = getResidueCounts(toCol,GAP_CHAR);
		int total = residues + gaps;
		switch (toState) {
		case DELETION:
			return (double)gaps/total;
		case MATCH:
			return (double)residues/total;
		case END:
			return 1;
		default: // we do not support delete to insert transitions
			throw new DataFormatException("Invalid destination state: "+toState);
		}
	}
	
	private double getBeginTransitionProb(byte toState)
	throws DataFormatException{
		return getDeletionTransitionProb(toState,0);		
	}
	
	private double getMatchInsertTransitionProb(byte toState, int toCol)
	throws DataFormatException{
		double prevInsertionProb = (toCol>0)?getDeletionTransitionProb(DELETION, toCol-1):0;
		double d2d = getDeletionTransitionProb(DELETION,toCol);
		double d2m = getDeletionTransitionProb(MATCH,toCol);
		//int residues = getResidueCounts(toCol);		
		//int gaps = getResidueCounts(toCol,GAP_CHAR);
		//int total = residues + gaps;				
		switch (toState) {
		case INSERTION:
			return getDeletionTransitionProb(DELETION, toCol);
		case DELETION:
			return (double)(1-prevInsertionProb)*d2d;
		case MATCH:
			return (double)(1-prevInsertionProb)*d2m;
		case END:
			return 1-prevInsertionProb;
		default:
			throw new DataFormatException("Invalid destination state: "+toState);
		}
	}
	
	public double getTransitionProb(byte fromState, byte toState, int toCol)
	throws DataFormatException{		
		switch (fromState) {
		case MATCH:
			return getMatchInsertTransitionProb(toState, toCol);
		case INSERTION: // match and insert are the same state really
			return getMatchInsertTransitionProb(toState, toCol);
		case BEGIN:
			if (toCol!=0) throw new DataFormatException("begin state must transition to 0");
			return getDeletionTransitionProb(toState, toCol);
		case DELETION:
			return getDeletionTransitionProb(toState, toCol);
		default:		
			throw new DataFormatException("Invalid from state "+fromState+" queried");
		}
	}
	
	private Sequence plugGapChar(Sequence seq, int col) {
		//Logger.debugln("before: seqlength: "+seqLength);
		String rawSeq = seq.getRawSeq();
		rawSeq = rawSeq.substring(0,col+1) + GAP_CHAR + rawSeq.substring(col+1,rawSeq.length());
		//Logger.debugln(rawSeq);
		if (rawSeq.length()>seqLength) {
			seqLength = rawSeq.length();
			//Logger.debugln("seqlength: "+seqLength);
		}
		return new Sequence(seq.getName(),seq.getDesc(),rawSeq);
	}
	
	public void addDeleteAtCol(int col) {
		Sequence seq = (Sequence)sequences.get(sequences.size()-1);
		sequences.set(sequences.size()-1,plugGapChar(seq,col));
	}
	
	public void addInsertAtCol(int col) {
		for (int count=0;count<sequences.size()-1;count++) {
			sequences.set(count,plugGapChar((Sequence)sequences.get(count),col));
		}
	}
	
	public Object clone() {
	 	List newList = new ArrayList();
		for (int count=0;count<sequences.size();count++) {
			Sequence oldSeq = (Sequence)sequences.get(count);
			Sequence seq = new Sequence(oldSeq.getName(),oldSeq.getDesc(),oldSeq.getRawSeq());
			newList.add(seq);
		}
		HMMAlignment align = null;
		try {
			align = new HMMAlignment(newList);		 
		}
		catch (Exception ex) {
			Logger.debugln(ex.getMessage());
		}
		return align;
	}	


	


}
