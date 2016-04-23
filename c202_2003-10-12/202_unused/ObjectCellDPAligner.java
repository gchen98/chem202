/*
 * Author: gchen
 * Created: Sunday, October 12, 2003 3:06:59 PM
 * Modified: Sunday, October 12, 2003 3:06:59 PM
 */

/**
* This class does the heavy lifting of the dynamic programming algorithm
* Variants such as Needleman-Wunsch and Smith-Waterman extend this class.
*/

package edu.ucla.chem202.gchen.algorithm.dp;
import edu.ucla.chem202.gchen.*;
import edu.ucla.chem202.gchen.io.Logger;
	
import java.io.InputStream;
import java.util.Arrays;
import java.util.List;

public abstract class ObjectCellDPAligner
implements DPAligner{
	
	protected int gapOpen, gapExtend;
	protected int yDim,xDim;	
	protected Cell[][] grid;
	protected Cell maxCell = null;
	protected Sequence seq1, seq2;
	protected ScoringMatrix matrix;
	
	// algorithm specific subclasses implement these
	
	protected abstract void initEdges();	
	protected abstract boolean useDefaultCellValue();
	protected abstract int getDefaultCellValue();
	protected abstract Cell getTracebackHead();	
	protected abstract boolean isTerminate(Cell cell);
	
	// Let's assume sequence 1 to rest on the y axis and sequence 2 on the x axis
	
	public void init(Sequence seq1, Sequence seq2, int gapOpen, int gapExtend, ScoringMatrix matrix)	{
		yDim = seq1.seqString().length()+1;
		xDim = seq2.seqString().length()+1;
		grid = new Cell[yDim][xDim];
		this.seq1 = seq1;
		this.seq2 = seq2;
		this.gapOpen = gapOpen;
		this.gapExtend = gapExtend;
		this.matrix = matrix;
	}
	
	public int getScore() {
		return maxCell.getScore();
	}
	
	// do traceback
	
	public Alignment getAlignment() {
		Alignment alignment = new Alignment(seq1.getName(),seq2.getName());
		Cell currCell = getTracebackHead();
		boolean terminate=false;		
		while (!terminate){
			
			int yCoor = currCell.getYCoor();
			int xCoor = currCell.getXCoor();
			Logger.debugln("Traceback coordinates are: "+yCoor+" and "+xCoor);
			String residue1 = null,residue2 = null;
			if (xCoor==currCell.getPrevious().getXCoor()) {
				residue1=seq1.seqString().substring(yCoor-1,yCoor);
				residue2=Alignment.GAP_CHAR;				
			}
			else if (yCoor==currCell.getPrevious().getYCoor()) {
				residue1=Alignment.GAP_CHAR;				
				residue2=seq2.seqString().substring(xCoor-1,xCoor);
			}
			else {
				residue1=seq1.seqString().substring(yCoor-1,yCoor);
				residue2=seq2.seqString().substring(xCoor-1,xCoor);				
			}
			alignment.addPosition(residue1,residue2);
			currCell = currCell.getPrevious();
			terminate=isTerminate(currCell);			
		}
		return alignment;
	}
	
	// populate grid
	
	public void align()
	throws DataFormatException{
		maxCell = new Cell(0,0,0);
		initEdges();		
		for (int y=1;y<yDim;y++) {
			for (int x=1;x<xDim;x++) {			
				int yPenalty = gapOpen, xPenalty = gapOpen;		
				
				Cell prevY = grid[y-1][x];
				if (!prevY.isMatch()) {
					yPenalty = gapExtend;					
				}
				int proposeYScore = prevY.getScore()-yPenalty;
					
				Cell prevX = grid[y][x-1];
				if (!prevX.isMatch()) {
					xPenalty = gapExtend;					
				}
				int proposeXScore = prevX.getScore()-xPenalty;					
									
				String yResidue = seq1.seqString().substring(y-1,y);
				String xResidue = seq2.seqString().substring(x-1,x);
				Cell prevDiag = grid[y-1][x-1];
				int proposeMatchScore = prevDiag.getScore() + matrix.getScore(yResidue,xResidue);					
									
				Cell prevCell = prevDiag; 
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
					prevCell = prevY;
					//Logger.debugln("Arrow up at "+y+" and "+x);
				}
				if (proposeXScore>proposedMax) {
					proposedMax = proposeXScore;
					prevCell = prevX;
					//Logger.debugln("Arrow left at "+y+" and "+x);
				}
				
				Cell cell=new Cell(proposedMax,y,x);
				if (proposedMax == proposeMatchScore) {
					cell.setMatch();
				}
				cell.setPrevious(prevCell);
				// for local alignment
				if (cell.getScore()>maxCell.getScore()) {
					maxCell = cell;
				}
				grid[y][x]=cell;
			}
		}
	}
	
	private void pad(String str) {
		int cellsize = 5;
		int padSize = cellsize - str.length();
		Logger.debug(str);
		for (short count=0;count<padSize;count++) {
			Logger.debug(" ");
		}
	}
	
	public void dump() {
		pad("");
		pad("");		
		String seq2Str = seq2.seqString();
		for (int count=0;count<seq2Str.length();count++) {
		  pad(seq2Str.substring(count,count+1));
		}
		Logger.debugln("");
		
		for (int y=0;y<yDim;y++) {
			if (y>0) pad(seq1.seqString().substring(y-1,y));
			else pad("");
			for (int x=0;x<xDim;x++) {
				Cell cell = grid[y][x];
				//Logger.debugln("This cell match type is: "+cell.isMatch());
				pad(Integer.toString(cell.getScore()));
			}
			Logger.debugln("");
		}
	}
	
	protected class Cell{
	
		private boolean isMatch = false;	
		private int score;
		private int yCoor, xCoor;	
		private Cell previous;
		
		public Cell(int score, int yCoor, int xCoor) {
			this.score = score;
			this.yCoor = yCoor;
			this.xCoor = xCoor;
		}
		
		public void setPrevious(Cell previous) {
			this.previous = previous;
		}
		
		public Cell getPrevious() {
			return previous;
		}
		
		public int getScore() {
			return this.score;
		}		
		
		public void setMatch() {
			this.isMatch = true;
		}
		
		// if this is true, then a pointer to this cell from a gap 
		// will be considered a gap opening, otherwise a gap extension
		
		public boolean isMatch() {
			return isMatch;
		}
		
		public int getXCoor() {
			return xCoor;
		}
		
		public int getYCoor() {
			return yCoor;
		}
	}
	
}
