/*
 * Author: gchen
 * Created: Sunday, October 12, 2003 3:06:59 PM
 * Modified: Sunday, October 12, 2003 3:06:59 PM
 */

/**
* Needleman-Wunsch variant
*/

package edu.ucla.chem202.gchen.algorithm.dp;
	
import java.io.InputStream;
import java.util.Arrays;
import java.util.List;

public class ObjectCellGlobalAligner
extends ObjectCellDPAligner{	
	
	// Let's assume sequence 1 to rest on the y axis and sequence 2 on the x axis	
	
	
	protected boolean isTerminate(Cell cell) {
		if ((cell.getXCoor()==0)&&(cell.getYCoor()==0)) return true;
		return false;
	}
	
	protected boolean useDefaultCellValue() {
		return false;
	}
	
	protected int getDefaultCellValue() {
		return -1;
	}
	
	protected Cell getTracebackHead() {
		return grid[yDim-1][xDim-1];
	}
	
	protected void initEdges() {
		grid[0][0]=new Cell(0,0,0);
		
		int extension=0-gapExtend;
		grid[0][1]=new Cell(0-gapOpen,0,1);
		grid[0][1].setPrevious(grid[0][0]);
		for (int x=2;x<xDim;x++) {
			grid[0][x]=new Cell(extension,0,x);
			grid[0][x].setPrevious(grid[0][x-1]);
			extension-=gapExtend;
		}
		
		extension=0-gapExtend;
		grid[1][0]=new Cell(0-gapOpen,1,0);
		grid[1][0].setPrevious(grid[0][0]);
		for (int y=2;y<yDim;y++) {
			grid[y][0]=new Cell(extension,y,0);
			grid[y][0].setPrevious(grid[y-1][0]);
			extension-=gapExtend;
		}
	}	
	
}
