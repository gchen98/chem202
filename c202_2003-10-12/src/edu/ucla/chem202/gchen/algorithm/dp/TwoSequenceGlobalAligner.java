/*
 * Author: gchen
 * Created: Sunday, October 12, 2003 3:06:59 PM
 * Modified: Sunday, October 12, 2003 3:06:59 PM
 */

/**
* Incomplete: not gonna use primive based cells
*/


package edu.ucla.chem202.gchen.algorithm.dp;
	
import java.io.InputStream;
import java.util.Arrays;
import java.util.List;

public class TwoSequenceGlobalAligner
extends TwoSequenceDPAligner{
	
	protected void initEdges(){
		grid[0][0]=0;
		
		int extension=0-gapExtend;		
		grid[0][1] = 0-gapOpen;
		pointers[0][1] = WEST;		
		for (int x=2;x<xDim;x++) {
			grid[0][x]=extension;			
			pointers[0][x]= WEST;
			extension-=gapExtend;
		}
		
		extension=0-gapExtend;
		grid[1][0]=0-gapOpen;
		pointers[1][0]=NORTH;
		for (int y=2;y<yDim;y++) {
			grid[y][0]=extension;
			pointers[y][0]=NORTH;
			extension-=gapExtend;
		}
	}
	
	protected boolean useDefaultCellValue() {
		return false;
	}
	
	protected int getDefaultCellValue(){
		return 0;
	}
	
	protected Cell getTracebackHead(){
		Cell cell = new Cell();		
		cell.yCoor = yDim-1;
		cell.xCoor = xDim-1;
		cell.score = grid[yDim-1][xDim-1];
		return cell;
	}
	
	protected boolean isTerminate(Cell cell){
		if ((cell.xCoor==0)&&(cell.yCoor==0)) return true;
		return false;		
	}
	
	public int getScore() {
		return getTracebackHead().score;
	}
}
