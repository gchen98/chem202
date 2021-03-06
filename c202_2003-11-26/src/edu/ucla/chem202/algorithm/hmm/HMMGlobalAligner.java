package edu.ucla.chem202.algorithm.hmm;

/*
 * Author: Teresa Breyer, Gary Chen, Sul-min Kim
 * Created: Saturday, November 15, 2003 11:06:55 PM
 * Modified: Saturday, November 15, 2003 11:06:55 PM
 */

import edu.ucla.chem202.*;
import edu.ucla.chem202.io.*;
import edu.ucla.chem202.algorithm.*;
import edu.ucla.chem202.algorithm.hmm.*;
import java.util.List;
import java.lang.Math;
import java.util.*;

import java.io.FileInputStream;
import java.io.FileOutputStream;

/**
* Align two sequences using dynamic programming.
*/
	
public class HMMGlobalAligner
	implements HMMAligner
{
	private HMM myHMM;
	private Sequence curSeq;
	private double probability;
	private int numberResidues;
	private int numberColumns;
	private Cell[][] matrix;
	// the following variables are deprecated. remove asap.
	private	double[][] D; // move to DELETION
	private double[][] M; // move to MATCH
	
	
	/**
	* @author: Teresa Breyer
	*/
	
	
	public void init(HMM hiddenMarkovModel)
	{
		myHMM = hiddenMarkovModel;
	}
	
	/**
	* @author: Teresa Breyer
	*/

	public void computeBackwardProbs()
	throws DataFormatException	{
		StringBuffer sequence = new StringBuffer(curSeq.getRawSeq());

		numberResidues = curSeq.getLength();
		numberColumns = myHMM.getLength();

		matrix[numberColumns - 1][numberResidues].bMatch = Math.log(myHMM.getTransitionProb(myHMM.MATCH, myHMM.END, numberColumns));

		matrix[numberColumns - 1][numberResidues].bDelete = Math.log(myHMM.getTransitionProb(myHMM.DELETION, myHMM.END, numberColumns));

		// last row of M
		for (int j = numberResidues - 1; j >=0 ; j--)
			matrix[numberColumns - 1][j].bMatch = matrix[numberColumns - 1][j+1].bMatch +				
			Math.log(myHMM.getEmissionProb(myHMM.MATCH,sequence.charAt(j), 
			numberColumns - 1)) + Math.log(myHMM.getTransitionProb(myHMM.MATCH, myHMM.INSERTION,  
			numberColumns - 1));

		// last row of D
		for (int j = numberResidues - 2; j >=0 ; j--)
			matrix[numberColumns - 1][j].bDelete = 0; // -INF	// can't go to INSERTION directly, so not possible

		// last column of D
		for (int i = numberColumns - 2; i >= 0;  i--)
		{
			//Logger.debugln("i=" + i + " " + myHMM.getTransitionProb(myHMM.DELETION, myHMM.DELETION,i+1));
			matrix[i][numberResidues].bDelete = matrix[i+1][numberResidues].bDelete +
				Math.log(myHMM.getTransitionProb(myHMM.DELETION, myHMM.DELETION,i+1));
		}

		// last column of M
		for (int i = numberColumns - 2; i >= 0; i--)
		{	
			//Logger.debugln("i=" + i);
			matrix[i][numberResidues].bMatch = matrix[i+1][numberResidues].bDelete +
				Math.log(myHMM.getTransitionProb(myHMM.MATCH, myHMM.DELETION,i+1));
		}

		// need the second last row and columns because of the -INF entries in last row and column

		// second last row of M
		for (int j = numberResidues - 1; j >=0 ; j--)
			matrix[numberColumns - 2][j].bMatch = MathUtil.addLogLog(matrix[numberColumns - 1][j+1].bMatch +
					Math.log(myHMM.getEmissionProb(myHMM.MATCH,sequence.charAt(j), numberColumns - 1)) +
					Math.log(myHMM.getTransitionProb(myHMM.MATCH, myHMM.MATCH, numberColumns - 1)),
					matrix[numberColumns - 2][j+1].bMatch + Math.log(myHMM.getEmissionProb(myHMM.MATCH, sequence.charAt(j), 
					numberColumns - 2)) + Math.log(myHMM.getTransitionProb(myHMM.MATCH, myHMM.INSERTION, numberColumns - 2)));

		// second last row of D
		for (int j = numberResidues - 1; j >=0 ; j--) {
			//Logger.debugln("j=" + j + " " +  matrix[numberColumns - 1][j+1].bMatch + " " + myHMM.getEmissionProb(myHMM.MATCH,sequence.charAt(j), numberColumns - 1) + " " +myHMM.getTransitionProb(myHMM.DELETION, myHMM.MATCH, numberColumns - 1));
			matrix[numberColumns - 2][j].bDelete = matrix[numberColumns - 1][j+1].bMatch + Math.log(myHMM.getEmissionProb(myHMM.MATCH,sequence.charAt(j), numberColumns - 1)) +
				Math.log(myHMM.getTransitionProb(myHMM.DELETION, myHMM.MATCH, numberColumns - 1));
		}

		// general case
		for (int i = numberColumns - 3; i >= 0; i--) {
			for (int j = numberResidues - 1; j >= 0; j--) {

				/*matrix[i][j].bMatch = matrix[i+1][j+1].bMatch * 
				myHMM.getEmissionProb(myHMM.MATCH,sequence.charAt(j+1), i+1) * 
				myHMM.getTransissionProb(myHMM.MATCH, myHMM.MATCH, i+1) +
					matrix[i][j+1].bMatch * myHMM.getEmissionProb(myHMM.MATCH, 
				sequence.charAt(j+1), i) * myHMM.getTransitionProb(myHMM.MATCH, 
				myHMM.INSERTION, i) +
					matrix[i+1][j].bDelete * myHMM.getTransitionProb(myHMM.MATCH, 
				myHMM.DELETION, i+1); */
				/*Logger.debugln("Match " + matrix[i+1][j+1].bMatch + " " + myHMM.getEmissionProb(myHMM.MATCH,sequence.charAt(j), i+1) + " " +
						myHMM.getTransitionProb(myHMM.MATCH, myHMM.MATCH, i+1));
				Logger.debugln("Insertion " + matrix[i][j+1].bMatch + " " + myHMM.getEmissionProb(myHMM.MATCH, sequence.charAt(j), i) + " " +
					myHMM.getTransitionProb(myHMM.MATCH, myHMM.INSERTION, i));
				Logger.debugln("Deletion " + matrix[i+1][j].bDelete + " " + myHMM.getTransitionProb(myHMM.MATCH, myHMM.DELETION, i+1) + " " );*/
				matrix[i][j].bMatch = MathUtil.addLogLog(matrix[i+1][j+1].bMatch + 
					Math.log(myHMM.getEmissionProb(myHMM.MATCH,sequence.charAt(j), i+1)) +
					Math.log(myHMM.getTransitionProb(myHMM.MATCH, myHMM.MATCH, i+1)),
					matrix[i][j+1].bMatch + Math.log(myHMM.getEmissionProb(myHMM.MATCH, sequence.charAt(j), i)) +
					Math.log(myHMM.getTransitionProb(myHMM.MATCH, myHMM.INSERTION, i)),
					matrix[i+1][j].bDelete + Math.log(myHMM.getTransitionProb(myHMM.MATCH, myHMM.DELETION, i+1)));

				/*matrix[i][j].bDelete = matrix[i+1][j].bDelete *  
				myHMM.getTransitionProb(myHMM.DELETION, myHMM.DELETION, i+1) +
					matrix[i+1][j+1].bMatch * 
				myHMM.getEmissionProb(myHMM.MATCH,sequence.charAt(j+1), i+1) * 
				myHMM.getTransitionProb(myHMM.DELETION, myHMM.MATCH, i+1);*/
				/*Logger.debugln("Match " + matrix[i+1][j+1].bMatch + " " + myHMM.getEmissionProb(myHMM.MATCH,sequence.charAt(j), i+1) + " "
					+ myHMM.getTransitionProb(myHMM.DELETION, myHMM.MATCH, i+1));
				Logger.debugln("Delete " + matrix[i+1][j].bDelete + " " + myHMM.getTransitionProb(myHMM.DELETION, myHMM.DELETION, i+1));*/
				matrix[i][j].bDelete = MathUtil.addLogLog(matrix[i+1][j].bDelete + 
					Math.log(myHMM.getTransitionProb(myHMM.DELETION, myHMM.DELETION, i+1)),
					matrix[i+1][j+1].bMatch + Math.log(myHMM.getEmissionProb(myHMM.MATCH,sequence.charAt(j), i+1)) +
					Math.log(myHMM.getTransitionProb(myHMM.DELETION, myHMM.MATCH, i+1)));

				}
		}
		// first row - first colum - last row
		// first column
		for (int i = 0; i < numberColumns; i++)
			matrix[i][0].pMatch =  0;//-INF
		for (int i = 0; i < numberColumns-1; i++)
			matrix[i][0].pDelete =  matrix[i][0].bDelete + matrix[i][0].fDelete;
		
		// first row
		for (int  j = 1; j < numberResidues + 1 ; j++)
			matrix[0][j].pMatch = matrix[0][j].bMatch + matrix[0][j].fMatch;	
		for (int  j = 1; j < numberResidues + 1 ; j++)
			matrix[0][j].pDelete = 0;//-INF
		
		//last row
		for (int  j = 1; j < numberResidues + 1 ; j++)
			matrix[numberColumns-1][j].pMatch = matrix[numberColumns-1][j].bMatch + matrix[numberColumns-1][j].fMatch;	
		for (int  j = 0; j < numberResidues ; j++)
			matrix[numberColumns-1][j].pDelete = 0;//-INF
		matrix[numberColumns-1][numberResidues].pDelete = matrix[numberColumns-1][numberResidues].bDelete + matrix[numberColumns-1][numberResidues].fDelete;

		// general case
		for (int i = 1; i < numberColumns - 1; i++) {
			for (int j = 1; j < numberResidues + 1; j++) {
				matrix[i][j].pMatch = matrix[i][j].bMatch + matrix[i][j].fMatch;
				matrix[i][j].pDelete = matrix[i][j].bDelete + matrix[i][j].fDelete;
			}
		}
	}

	/**
	* @author: Teresa Breyer
	*/

	public void computeForwardProbs()
	throws DataFormatException
	{
		StringBuffer sequence = new StringBuffer(curSeq.getRawSeq());

		numberResidues = curSeq.getLength();
		numberColumns = myHMM.getLength();

		matrix[0][0].fDelete = Math.log(myHMM.getTransitionProb(myHMM.BEGIN, 
			myHMM.DELETION, 0));
		matrix[0][1].fMatch = Math.log(myHMM.getEmissionProb(myHMM.MATCH, 
			sequence.charAt(0), 0) ) +
			Math.log(myHMM.getTransitionProb(myHMM.BEGIN, myHMM.MATCH, 0));

		// first row of M
		for (int j = 2; j < numberResidues + 1; j++) {
			matrix[0][j].fMatch = Math.log(myHMM.getEmissionProb(myHMM.MATCH, 
				sequence.charAt(j-1), 0)) +
				matrix[0][j-1].fMatch + Math.log( myHMM.getTransitionProb(myHMM.MATCH, 
				myHMM.INSERTION, 0)) ;
		}

		// first column of D
		for (int i = 1; i < numberColumns; i++) {
			matrix[i][0].fDelete = matrix[i-1][0].fDelete + 
			Math.log(myHMM.getTransitionProb(myHMM.DELETION, myHMM.DELETION, i));
		}

		matrix[1][1].fDelete =  matrix[0][1].fMatch + Math.log(myHMM.getTransitionProb(myHMM.MATCH, myHMM.DELETION, 1));
		Logger.debugln("1 1");
		matrix[1][1].fMatch = Math.log(myHMM.getEmissionProb(myHMM.MATCH, sequence.charAt(0), 1)) +
		matrix[0][0].fDelete + Math.log(myHMM.getTransitionProb(myHMM.DELETION, myHMM.MATCH, 1)) ;

		//second column both
		for (int i = 2; i < numberColumns; i++) {
			matrix[i][1].fDelete = MathUtil.addLogLog(matrix[i-1][1].fMatch + Math.log(myHMM.getTransitionProb(myHMM.MATCH, myHMM.DELETION, i)),
				matrix[i-1][1].fDelete + Math.log(myHMM.getTransitionProb(myHMM.DELETION, myHMM.DELETION, i))); // second column of D
			matrix[i][1].fMatch = Math.log(myHMM.getEmissionProb(myHMM.MATCH, sequence.charAt(0), i)) +
				matrix[i-1][0].fDelete + Math.log(myHMM.getTransitionProb(myHMM.DELETION, myHMM.MATCH, i)) ;
		}

		// second row both
		for (int j = 2; j < numberResidues + 1; j++) {
			matrix[1][j].fMatch = Math.log(myHMM.getEmissionProb(myHMM.MATCH, sequence.charAt(j-1), 1)) +
				MathUtil.addLogLog(matrix[0][j-1].fMatch + Math.log(myHMM.getTransitionProb(myHMM.MATCH, myHMM.MATCH, 1)),
					matrix[1][j-1].fMatch + Math.log(myHMM.getTransitionProb(myHMM.MATCH, myHMM.INSERTION, 1)));
			matrix[1][j].fDelete = matrix[0][j].fMatch + Math.log(myHMM.getTransitionProb(myHMM.MATCH, myHMM.DELETION, 1));
		}

		// general case
		for (int i = 2; i < numberColumns; i++) {
			for (int j = 2; j < numberResidues + 1; j++) {

				matrix[i][j].fMatch = Math.log(myHMM.getEmissionProb(myHMM.MATCH, sequence.charAt(j-1), i)) +
					MathUtil.addLogLog(matrix[i-1][j-1].fMatch + Math.log(myHMM.getTransitionProb(myHMM.MATCH, myHMM.MATCH, i)),
					matrix[i][j-1].fMatch + Math.log(myHMM.getTransitionProb(myHMM.MATCH, myHMM.INSERTION, i)),
					matrix[i-1][j-1].fDelete + Math.log(myHMM.getTransitionProb(myHMM.DELETION, myHMM.MATCH, i))) ;
				matrix[i][j].fDelete = MathUtil.addLogLog(matrix[i-1][j].fMatch + Math.log(myHMM.getTransitionProb(myHMM.MATCH, myHMM.DELETION, i)),
					matrix[i-1][j].fDelete + Math.log(myHMM.getTransitionProb(myHMM.DELETION, myHMM.DELETION, i)));
			}
		}

		if (matrix[numberColumns-1][numberResidues].fMatch >= matrix[numberColumns-1][numberResidues].fDelete) {
			Logger.debugln("prob for match greater than delete");
			probability = (matrix[numberColumns-1][numberResidues].fMatch +	Math.log(1 + Math.exp(matrix[numberColumns-1][numberResidues].fDelete- 
				matrix[numberColumns-1][numberResidues].fMatch))) / Math.log(10);
		} else {
			Logger.debugln("prob for match less than delete");
			probability = (matrix[numberColumns-1][numberResidues].fDelete +
				Math.log(1 + Math.exp(matrix[numberColumns-1][numberResidues].fMatch- 
				matrix[numberColumns-1][numberResidues].fDelete))) / Math.log(10);
		}
	}

	/**
	* @author: Teresa Breyer
	*/

	public void printDataStructure()
	{
		Logger.debugln(" ");
		Logger.debugln("M= ");
		for (int i = 0; i < matrix.length; i++) {
			for (int j = 0; j < matrix[0].length; j++) {
				Logger.debug(MathUtil.round(matrix[i][j].fMatch) + " ");
			}
			Logger.debugln(" ");
		}
		Logger.debugln(" ");
		Logger.debugln("D= ");
		for (int i = 0; i < matrix.length; i++) {
			for (int j = 0; j < matrix[0].length; j++) {
				Logger.debug(MathUtil.round(matrix[i][j].fDelete) + " ");
			}
			Logger.debugln(" ");
		}
		Logger.debugln(" ");
		Logger.debugln("backward:\nM= ");
		for (int i = 0; i < matrix.length; i++) {
			for (int j = 0; j < matrix[0].length; j++) {
				Logger.debug(MathUtil.round(matrix[i][j].bMatch) + " ");
			}
			Logger.debugln(" ");
		}
		Logger.debugln(" ");
		Logger.debugln("D= ");
		for (int i = 0; i < matrix.length; i++) {
			for (int j = 0; j < matrix[0].length; j++) {
				Logger.debug(MathUtil.round(matrix[i][j].bDelete) + " ");
			}
			Logger.debugln(" ");
		}
		Logger.debugln(" ");
		Logger.debugln("viterbi:\nM= ");
		for (int i = 0; i < matrix.length; i++) {
			for (int j = 0; j < matrix[0].length; j++) {
				Logger.debug(MathUtil.round(matrix[i][j].vMatch) + " ");
			}
			Logger.debugln(" ");
		}
		Logger.debugln(" ");
		Logger.debugln("D= ");
		for (int i = 0; i < matrix.length; i++) {
			for (int j = 0; j < matrix[0].length; j++) {
				Logger.debug(MathUtil.round(matrix[i][j].vDelete) + " ");
			}
			Logger.debugln(" ");
		}
		Logger.debugln(" ");
		Logger.debugln("traceBackState:\nM= ");
		for (int i = 0; i < matrix.length; i++) {
			for (int j = 0; j < matrix[0].length; j++) {
				Logger.debug(MathUtil.round(matrix[i][j].traceBackState) + " ");
			}
			Logger.debugln(" ");
		}
	}
	
	/**
	* @author: Gary Chen
	* All this does is populate the grid	
	*/
	
	public void align(Sequence seq)
	throws DataFormatException{
		curSeq = seq;
		int hmmLengh = myHMM.getLength();
		int seqLength = seq.getLength()+1;
		matrix = new Cell[hmmLengh][seqLength];
		for (int y=0;y<hmmLengh;y++) {
			for (int x=0;x<seqLength;x++) {
				matrix[y][x] = new Cell();
				matrix[y][x].yLoc = y;
				matrix[y][x].xLoc = x;
			}
		}
		computeForwardProbs();
		//computeBackwardProbs();
		//computeViterbiProbs();
	}
	
	/**
	* @author: Gary Chen
	* Compute best path and generate traceback
	* Code is handling boundary conditions first before general conditions for performance reasons.
	*/
	
	public void computeViterbiProbs()
	throws DataFormatException {
		String sequence = curSeq.getRawSeq();
		// *** FIRST ROW AND COLUMN ***
		// initial condition for delete state from begin
		matrix[0][0].vDelete = Math.log(myHMM.getTransitionProb(HMM.BEGIN, HMM.DELETION, 0));		
		// first column of D.  DEGENERATE condition where only previous move available which is deletion.
		for (int i = 1; i < numberColumns; i++) {	
			matrix[i][0].vDelete = matrix[i-1][0].vDelete + Math.log(myHMM.getTransitionProb(HMM.DELETION, HMM.DELETION, i));			
			matrix[i][0].traceBack = matrix[i-1][0]; // only one traceback possible. point up for deletion
			matrix[i-1][0].traceBackState = HMM.DELETION;	// previous move was a deletion
		}
		// initial condition for match state from begin		
		matrix[0][1].vMatch = Math.log(myHMM.getEmissionProb(HMM.MATCH, sequence.charAt(0), 0) ) + 
			Math.log(myHMM.getTransitionProb(HMM.BEGIN, HMM.MATCH, 0));
		// first row of M. DEGENERATE condition where only previous move available which is insertion.
		for (int j = 2; j < numberResidues + 1; j++) {
			matrix[0][j].vMatch = Math.log(myHMM.getEmissionProb(HMM.MATCH, sequence.charAt(j-1), 0)) + 
				matrix[0][j-1].vMatch + Math.log( myHMM.getTransitionProb(HMM.MATCH, HMM.INSERTION, 0)) ;
			matrix[0][j].traceBack = matrix[0][j-1];  // only one traceback possible. point left for insertion
			matrix[0][j-1].traceBackState = HMM.INSERTION; // previous move was a insertion
		}
		// *** SECOND ROW AND COLUMN ***
		// DEGENERATE condition for delete state where the ONLY previous state can be a match
		matrix[1][1].vDelete =  matrix[0][1].vMatch + Math.log(myHMM.getTransitionProb(HMM.MATCH, HMM.DELETION, 1));
		// DEGENERATE condition for match state where the ONLY previous state can be a delete
		matrix[1][1].vMatch = Math.log(myHMM.getEmissionProb(HMM.MATCH, sequence.charAt(0), 1)) + 
				matrix[0][0].vDelete + Math.log(myHMM.getTransitionProb(HMM.DELETION, HMM.MATCH, 1)) ;
		// determine traceback for 1,1
		if (matrix[1][1].vMatch>=matrix[1][1].vDelete) { // optimistically choose the to-match move
			matrix[1][1].traceBack = matrix[0][0];
			matrix[0][0].traceBackState = HMM.DELETION; // represents delete to match
		}
		else { // choose the to-delete move
			matrix[1][1].traceBack = matrix[0][1];
			matrix[0][1].traceBackState = HMM.MATCH; // represents match to delete
		}		
		//second column of match and delete		
		for (int i = 2; i < numberColumns; i++) {
			// normal condition for delete state where the ONLY two previous states are a match and delete
			double m2d = MathUtil.addLog(matrix[i-1][1].vMatch, myHMM.getTransitionProb(HMM.MATCH, HMM.DELETION, i));
			double d2d = MathUtil.addLog(matrix[i-1][1].vDelete, myHMM.getTransitionProb(HMM.DELETION, HMM.DELETION, i));			
			matrix[i][1].vDelete = (d2d>=m2d)?d2d:m2d;  // discourage matches from going to deletes
			// DEGENERATE condition for match state where the ONLY previous state can be a delete
			matrix[i][1].vMatch = Math.log(myHMM.getEmissionProb(HMM.MATCH, sequence.charAt(0), i)) + 
				matrix[i-1][0].vDelete + Math.log(myHMM.getTransitionProb(HMM.DELETION, HMM.MATCH, i));
			// determine traceback for i,1
			if (matrix[i][1].vMatch>=matrix[i][1].vDelete) { // optimistically choose the to-match move
				matrix[i][1].traceBack = matrix[i-1][0];
				matrix[i-1][0].traceBackState = HMM.DELETION;  // represents delete to match move
			}
			else { // choose the to-delete move
				matrix[i][1].traceBack = matrix[i-1][1];
				matrix[i-1][1].traceBackState = (m2d>=d2d)?HMM.MATCH:HMM.DELETION;  // represents either match to delete or delete to delete move
			}			
		}		
		// second row of match and delete.  
		for (int j = 2; j < numberResidues + 1; j++) {
			// DEGENERATE condition for match state where the ONLY two previous state can be a match or insertion
			double m2m = MathUtil.addLog(matrix[0][j-1].vMatch, myHMM.getTransitionProb(HMM.MATCH, HMM.MATCH, 1));
			double m2i = MathUtil.addLog(matrix[1][j-1].vMatch, myHMM.getTransitionProb(HMM.MATCH, HMM.INSERTION, 1));
			double max = (m2m>=m2i)?m2m:m2i; // encourage entering match state
			matrix[1][j].vMatch = Math.log(myHMM.getEmissionProb(HMM.MATCH, sequence.charAt(j-1), 1)) + max;
			// DEGENERATE condition for delete state where the ONLY previous state can be a match
			matrix[1][j].vDelete = matrix[0][j].vMatch + Math.log(myHMM.getTransitionProb(HMM.MATCH, HMM.DELETION, 1));
			if (matrix[1][j].vMatch>=matrix[1][j].vDelete) { // optimistically choose the to-match move
				if (m2m>=m2i) {  // encourage to-match moves
					matrix[1][j].traceBack = matrix[0][j-1];
					matrix[0][j-1].traceBackState = HMM.MATCH; // match to match
				}
				else {
					matrix[1][j].traceBack = matrix[1][j-1];
					matrix[1][j-1].traceBackState = HMM.INSERTION; // match to insert
				}
			}
			else { // choose the to-delete move
				matrix[1][j] = matrix[0][j];
				matrix[0][j].traceBackState = HMM.MATCH; // match to delete
			}
		}		
		// *** GENERAL CASE ***
		for (int i = 2; i < numberColumns; i++) {
			for (int j = 2; j < numberResidues + 1; j++) {				
				// normal condition for match state where the ONLY three previous states are a match, insert, and delete
				double m2m = MathUtil.addLog(matrix[i-1][j-1].vMatch, myHMM.getTransitionProb(HMM.MATCH, HMM.MATCH, i));
				double m2i = MathUtil.addLog(matrix[i][j-1].vMatch, myHMM.getTransitionProb(HMM.MATCH, HMM.INSERTION, i));
				double d2m = MathUtil.addLog(matrix[i-1][j-1].vDelete, myHMM.getTransitionProb(HMM.DELETION, HMM.MATCH, i));
				double mMax = m2m; // encourage m2m moves
				if (m2i>mMax) mMax = m2i;
				if (d2m>mMax) mMax = d2m;	// compare three to-match moves and choose the best
				matrix[i][j].vMatch = Math.log(myHMM.getEmissionProb(HMM.MATCH, sequence.charAt(j-1), i)) + mMax;
				// normal condition for delete state where the ONLY two previous states are a match and delete
				double m2d = MathUtil.addLog(matrix[i-1][j].vMatch,  myHMM.getTransitionProb(HMM.MATCH, HMM.DELETION, i));
				double d2d = MathUtil.addLog(matrix[i-1][j].vDelete, myHMM.getTransitionProb(HMM.DELETION, HMM.DELETION, i));				
				matrix[i][j].vDelete = (d2d>=m2d)?d2d:m2d; // discourage match to delete moves
				// determine traceback
				if (matrix[i][j].vMatch>=matrix[i][j].vDelete) { // choose the to-match move
					if ((m2i>m2m) && (m2i>d2m)) {  // if match to insert is the biggest value. if ties, then go to else
						matrix[i][j].traceBack = matrix[i][j-1];
						matrix[i][j-1].traceBackState = HMM.INSERTION; // match to insert
					}
					else { // otherwise one of the other two is the biggest. find it.
						matrix[i][j].traceBack =  matrix[i-1][j-1];
						matrix[i-1][j-1].traceBackState = (m2m>=d2m) ? HMM.MATCH : HMM.DELETION; // either match to match or delete to match
					}					
				}
				else { // choose the to-delete move
					matrix[i][j].traceBack = matrix[i-1][j];
					matrix[i-1][j].traceBackState = (d2d>=m2d)?HMM.DELETION:HMM.MATCH; // discourage match to deletes
				}
			}
		}
	}

	/**
	* @author: Teresa Breyer
	*/
	
	
	public double getTotalProb()
	{
		return probability;
	}
	
	/**
	* @author: Sul-min Kim
	*/	
	
	public void printPosteriorProbs() {
	    // start from bottom right cell
	    // get its traceback cell
	    // get the traceback state of the traceback cell to determine whether to emit something or not
	    // print the posterior prob at the cell
	    // repeat until traceback is null	
	    
	    String sequence = curSeq.getRawSeq();
	    int i = matrix.length-1; //index of last row element;
	    int j = matrix[0].length-1; //index of last col element -> numberResidues-1;
	    Cell currCell = matrix[i][j]; //start traceback at last array element
	    byte currTBState;
	    List output = new ArrayList(); //hold output for stdout

	    Logger.debugln("*****SMK DEBUG*****");
	    Logger.debugln("(j,i) "+j+" "+i+" seq: "+sequence+ " TotalProb: "+getTotalProb());

	    while(currCell != null){
		Logger.debugln("while x,y = "+currCell.xLoc+","+currCell.yLoc+" pMatch:"+currCell.pMatch );
		//Logger.debugln("\tcurrCell tb = "+currCell.traceBack.xLoc+" "+currCell.traceBack.yLoc);
		if(currCell.traceBackState != HMM.DELETION){
		    //output.add("COLUMN:"+j+" LETTER:"+sequence.charAt(j-1)+" PROB="+ MathUtil.round( Math.exp( currCell.fMatch + currCell.bMatch - getTotalProb() ) ) ); //push into the array the 
		    output.add("COLUMN:"+j+" LETTER:"+sequence.charAt(j-1)+" PROB="+ MathUtil.round( Math.exp( currCell.pMatch - getTotalProb() ) * 100 ) + 
			       " x,y: " +currCell.xLoc+","+currCell.yLoc);
		}
		if(currCell.traceBackState == HMM.DELETION){
		    output.add("COLUMN:"+j+" LETTER:- PROB="+ MathUtil.round( Math.exp( currCell.pMatch - getTotalProb() ) * 100 ) +
			       " x,y: " +currCell.xLoc+","+currCell.yLoc); //push into the array the 
		}
		if(matrix[i][j].traceBack != null){
		    currCell = matrix[i][j].traceBack;
		    i = currCell.yLoc;
		    j = currCell.xLoc;
		}
		else{
		    break;
		}
	    }
	    if(output.size() > 0){
		for(int k=output.size()-1; k>=0; k--){
		    Logger.println( output.get(k).toString() );
		}
	    }
	    Logger.debugln("*****SMK DEBUG DONE*****");	    
	}
		
	class Cell {
		double fMatch;
		double fDelete;
		double bMatch;
		double bDelete;
		double vMatch;
		double vDelete;
		double pMatch;
		double pDelete;
		int xLoc;
		int yLoc;
		byte traceBackState;	// use this to query what kind of state this was when the traceback was generated.
		Cell traceBack;			// points to a previous cell as the traceback.
	}

}
