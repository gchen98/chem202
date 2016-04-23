/*
 * Author: gchen
 * Created: Sunday, October 12, 2003 3:06:59 PM
 * Modified: Sunday, October 12, 2003 3:06:59 PM
 */

/**
* Smith Waterman variant
*/

package edu.ucla.chem202.gchen.algorithm.dp;
	
import java.io.InputStream;
import java.util.Arrays;
import java.util.List;

public class ObjectCellLocalAligner
extends ObjectCellDPAligner{	
	
	// Let's assume sequence 1 to rest on the y axis and sequence 2 on the x axis	
	
	
	protected boolean isTerminate(Cell cell) {
		if (cell.getScore()==0) return true;
		return false;
	}
	
	protected boolean useDefaultCellValue() {
		return true;
	}
	
	protected int getDefaultCellValue() {
		return 0;
	}
	
	protected Cell getTracebackHead() {
		return maxCell;
	}
	
	protected void initEdges() {
		grid[0][0] = new Cell(0,0,0);
		for (int x=1;x<xDim;x++) {
			grid[0][x]=new Cell(0,0,x);
			grid[0][x].setPrevious(grid[0][x-1]);
		}		
		for (int y=1;y<yDim;y++) {
			grid[y][0]=new Cell(0,y,0);
			grid[y][0].setPrevious(grid[y-1][0]);
		}
	}	
	
}
