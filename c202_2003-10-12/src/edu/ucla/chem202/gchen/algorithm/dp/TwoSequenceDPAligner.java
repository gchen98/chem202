/*
 * Author: gchen
 * Created: Sunday, October 12, 2003 3:06:59 PM
 * Modified: Sunday, October 12, 2003 3:06:59 PM
 */


/**
*
* incomplete: used for benchmarking memory usage of arrays vs objects.
* Suspected objects might be more expensive but this didn't seem to be the case.
*/

package edu.ucla.chem202.gchen.algorithm.dp;

import edu.ucla.chem202.gchen.*;
import edu.ucla.chem202.gchen.algorithm.*;
import edu.ucla.chem202.gchen.io.Logger;
	
import java.io.InputStream;
import java.util.Arrays;
import java.util.List;

public abstract class TwoSequenceDPAligner
implements Aligner{
	
	protected int gapOpen, gapExtend;
	protected int yDim,xDim;
	//protected int maxScore, maxX, maxY;
	protected Cell maxCell;
		
	protected int [][] grid;
	protected byte [][] pointers;
	protected boolean [][] matches;
	
	protected Sequence seq1, seq2;
	protected ScoringMatrix matrix;
	protected Cell cell;
	
	public static final byte UNDEF = 0;
	public static final byte NORTH = 1;
	public static final byte NW = 2;
	public static final byte WEST = 3;
	
	// Let's assume sequence 1 to rest on the y axis and sequence 2 on the 
	
	public void init(List sequences, int gapOpen, int gapExtend, ScoringMatrix matrix)	{		
		Sequence seq1 = (Sequence)sequences.get(0);
		Sequence seq2 = (Sequence)sequences.get(1);
		yDim = seq1.seqString().length()+1;
		xDim = seq2.seqString().length()+1;
		grid = new int[yDim][xDim];
		pointers = new byte[yDim][xDim];
		matches = new boolean[yDim][xDim];
		this.seq1 = seq1;
		this.seq2 = seq2;
		this.gapOpen = gapOpen;
		this.gapExtend = gapExtend;
		this.matrix = matrix;
		maxCell = new Cell();
	}
	
	// do traceback
	
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
				residue2=seq2.seqString().substring(xCoor-1,xCoor);
			}
			else {
				residue1=seq1.seqString().substring(yCoor-1,yCoor);
				residue2=seq2.seqString().substring(xCoor-1,xCoor);				
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
	
	protected abstract void initEdges();
	protected abstract boolean useDefaultCellValue();
	protected abstract int getDefaultCellValue();
	protected abstract Cell getTracebackHead();	
	protected abstract boolean isTerminate(Cell cell);
	
		
	// benchmarking purposes
	
	public void align()
	throws DataFormatException{
		maxCell.score = 0;
		maxCell.xCoor = 0;
		maxCell.yCoor = 0;
		initEdges();
		for (int y=1;y<yDim;y++) {
			for (int x=1;x<xDim;x++) {			
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
		}
	}
	
	
	
	private void pad(String str, String direction) {
		int cellsize = 4;
		int padSize = cellsize - str.length();
		
		for (short count=0;count<padSize;count++) {
			Logger.debug(" ");
		}
		Logger.debug(direction);
		Logger.debug(str);
	}	
	
	public void dump() {
		pad(""," ");
		pad(""," ");		
		String seq2Str = seq2.seqString();
		for (int count=0;count<seq2Str.length();count++) {
		  pad(seq2Str.substring(count,count+1)," ");
		}
		Logger.debugln("");
		
		for (int y=0;y<yDim;y++) {
			if (y>0) pad(seq1.seqString().substring(y-1,y)," ");
			else pad(""," ");
			for (int x=0;x<xDim;x++) {				
				//Logger.debugln("This cell match type is: "+cell.isMatch());
				//Logger.debugln("Grid printing "+y+" and "+x+" and xDim is "+xDim);
				String dir = null;
				switch (pointers[y][x]) {
				case NORTH:
					dir = "|";
					break;
				case NW:
					dir = "\\";				
					break;
				case WEST:
					dir = "-";
					break;
				case UNDEF:
					dir="*";				
				}
				pad(Integer.toString(grid[y][x]),dir);
			}
			Logger.debugln("");
		}
	}
	
	protected class Cell {
		protected int xCoor;
		protected int yCoor;
		protected int score;
	}
	
}
