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

public class TwoSequenceLocalAligner
extends TwoSequenceDPAligner{
	
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
}
