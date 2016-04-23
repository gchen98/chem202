/*
 * Author: gchen
 * Created: Sunday, October 12, 2003 3:06:59 PM
 * Modified: Sunday, October 12, 2003 3:06:59 PM
 */


package edu.ucla.chem202.gchen.algorithm.dp;
	

import edu.ucla.chem202.gchen.*;
import edu.ucla.chem202.gchen.io.*;
import edu.ucla.chem202.gchen.algorithm.*;

import java.io.InputStream;
import java.util.Arrays;
import java.util.List;
import java.util.ArrayList;


public class TwoSequenceLoopAligner
extends TwoSequenceDPAligner{
	
	private List subSequences = new ArrayList();
	
	class SubSequence {
		int startIndex;
		int endIndex;
		int currentScore;		
	}
	
	// need to override parent as this is really different
	
	
	protected void initEdges(){
		grid[0][0] = 0;
		matches[0][0] = false;
		for (int x=1;x<xDim;x++) {
			grid[0][x]=0;
			matches[0][x] = false;
			pointers[0][x] = WEST;
		}		
		for (int y=1;y<yDim;y++) {
			grid[y][0]=0;
			matches[y][0] = false;
			pointers[y][0] = NORTH;
		}
	}
	
	protected boolean useDefaultCellValue() {
		return true;
	}
	
	protected int getDefaultCellValue(){
		return 0;
	}
	
	protected Cell getTracebackHead(){
		Cell cell = new Cell();
		cell.yCoor = maxCell.yCoor;
		cell.xCoor = maxCell.xCoor;
		cell.score = maxCell.score;		
		return cell;
	}
	
	protected boolean isTerminate(Cell cell){
		if (cell.score==0) return true;
		return false;
	}
	
	public int getScore() {
		return maxCell.score;
	}
	
	public Alignment getAlignment() {
		Alignment alignment = new Alignment(seq1.getName(),seq2.getName());
		Cell cell = getTracebackHead();
		int yCoor = cell.yCoor;
		int xCoor = cell.xCoor;
		int score = cell.score;
		boolean terminate=false;		
		while (!terminate){				
			Logger.debugln("Traceback coordinates are: "+yCoor+" and "+xCoor+" with score "+score);
			String residue1 = null,residue2 = null;
			int prevY = yCoor;
			int prevX = xCoor;
			Logger.debugln("Pointer is pointing "+pointers[yCoor][xCoor]);
			switch (pointers[yCoor][xCoor]) {
			case NORTH:
				prevY--;
				break;
			case NW:
				prevY--;
				prevX--;				
				break;
			case WEST:
				prevX--;
			}
			if (xCoor==prevX) {
				residue1=seq1.seqString().substring(yCoor-1,yCoor);
				residue2=Alignment.GAP_CHAR;			
			}
			else if (yCoor==prevY) {
				residue1=Alignment.GAP_CHAR;				
				residue2=Sequence.bases.get(seq2.seqString().substring(xCoor-1,xCoor)).toString(); 
				
			}
			else {
				residue1=seq1.seqString().substring(yCoor-1,yCoor);
				//residue2=seq2.seqString().substring(xCoor-1,xCoor);
				residue2=Sequence.bases.get(residue1).toString();
			}
			alignment.addPosition(residue1,residue2);
			yCoor = cell.yCoor = prevY;
			xCoor = cell.xCoor = prevX;
			score = cell.score = grid[yCoor][xCoor];
			//currCell = currCell.getPrevious();
			terminate=isTerminate(cell);
		}
		return alignment;
	}
	
	public void align()
	throws DataFormatException{
		maxCell.score = 0;
		maxCell.xCoor = 0;
		maxCell.yCoor = 0;
		initEdges();
		int yDimTemp = yDim-1;		
		int xDimTemp = xDim-1;
		for (int y=1;y<yDimTemp;y++) {
			for (int x=1;x<xDimTemp;x++) {
				//Logger.debugln("y and x coordinates: "+y+" and "+x);
				int yPenalty = gapOpen, xPenalty = gapOpen;		
				
				int prevYScore = grid[y-1][x];
				if (matches[y-1][x]) {				
					yPenalty = gapExtend;					
				}
				int proposeYScore = prevYScore-yPenalty;
					
				int prevXScore = grid[y][x-1];
				if (matches[y][x-1]) {
					xPenalty = gapExtend;					
				}
				int proposeXScore = prevXScore-xPenalty;					
									
				char yResidue = seq1.seqString().charAt(y-1);
				char xResidue = seq2.seqString().charAt(x-1);
				int prevDiagScore = grid[y-1][x-1];
				int proposeMatchScore = prevDiagScore + matrix.getScore(yResidue,xResidue);					
									
				int prevYCoor = y-1;
				int prevXCoor = x-1;
				byte pointer = NW;
				int proposedMax = -1;				
				
				if (useDefaultCellValue()) { 
					proposedMax = getDefaultCellValue(); // in algorithms such as Smith Waterman
				} 
				else {
					proposedMax = proposeMatchScore; // in algorithms such as Needleman-Wunsch
				}
				
				if (proposeMatchScore>=proposedMax) {
					proposedMax = proposeMatchScore; // let's default ties to the previous diagonal					
				}
				if (proposeYScore>proposedMax) {
					proposedMax = proposeYScore;
					prevYCoor = prevYCoor;
					prevXCoor = prevXCoor;
					pointer = NORTH;
					//Logger.debugln("Arrow up at "+y+" and "+x);
				}
				if (proposeXScore>proposedMax) {
					proposedMax = proposeXScore;
					prevYCoor = prevYCoor;
					prevXCoor = prevXCoor;
					pointer = WEST;
					//Logger.debugln("Arrow left at "+y+" and "+x);
				}
				
				//Cell cell=new Cell(proposedMax,y,x);
				if (proposedMax == proposeMatchScore) {
					matches[y][x] = true;					
				}
				pointers[y][x] = pointer;
				//cell.setPrevious(prevCell);
				// for local alignment
				if (proposedMax>maxCell.score) {
					maxCell.score = proposedMax;					
					maxCell.yCoor = y;
					maxCell.xCoor = x;
				}
				grid[y][x]=proposedMax;
				
			}
			xDimTemp--;
			
		}
		yDimTemp--;
	}	
}
