/*
 * Author: gchen
 * Created: Sunday, October 12, 2003 3:06:59 PM
 * Modified: Sunday, October 12, 2003 3:06:59 PM
 */

package edu.ucla.chem202.gchen.algorithm.poa;

import edu.ucla.chem202.gchen.io.*;
import edu.ucla.chem202.gchen.*;
import edu.ucla.chem202.gchen.algorithm.*;
import java.util.Stack;
import java.util.List;
import java.util.ArrayList;
import java.util.Iterator;
import java.io.*;

public class POALocalAligner
implements Aligner, Constants{
	
	private ScoringMatrix matrix;
	private int gapPenalty;
	private Cell[][] cells;
	private Cell maxCell;
	private List nodeList;
	private List sortedList;
	private List sequences;
	private int currentSeqId, headExtensionStart, tailExtensionStart;
	
	// encapsulates a cell grid
	
	class Cell {
		TraceBack traceBack;
		int col,row;
		Cell(int row, int col, TraceBack traceBack) {
			this.col = col;
			this.row = row;			
			this.traceBack = traceBack;
		}		
	}
	
	// encapsulates info about the score at a position, and which cell to traceback to
	
	class TraceBack{
		static final char MATCH='M';
		static final char INSERT='I';
		static final char DELETE='D';
		Cell cell;
		int propScore;
		int type;
		TraceBack(Cell cell, int propScore, char type) {
			this.cell = cell;
			this.propScore = propScore;
			this.type = type;
		}
	}
	
	class Node
	implements Constants{		
		char id; // represents sequence character
		int index = UNDEF; // the position in the graph
		List prevNodeIndices = new ArrayList();
		List successorNodeIndices = new ArrayList();
		List containedSeqs = new ArrayList();
				
		void addPrevNodeIndex(int newNodeIndex) {
			if (!prevNodeIndices.contains(new Integer(newNodeIndex))) {
				prevNodeIndices.add(new Integer(newNodeIndex));
			}
		}
		void addSuccessorNodeIndex(int newNodeIndex) {			
			if (!successorNodeIndices.contains(new Integer(newNodeIndex))) {
				successorNodeIndices.add(new Integer(newNodeIndex));
			}
		}
		void addContainedSeqId(int seqId) {
			if (!containedSeqs.contains(new Integer(seqId))) {
				containedSeqs.add(new Integer(seqId));
			}
		}
		/**
		* convenience methods for adding a node to the graph
		*/		
		void addPrevNode(int nodeIndex) {
			this.addPrevNodeIndex(nodeIndex);
			Node node = (Node)nodeList.get(nodeIndex);
			node.addSuccessorNodeIndex(index);						
		}		
		void addSuccessorNode(int nodeIndex) {
			this.addSuccessorNodeIndex(nodeIndex);
			Node node = (Node)nodeList.get(nodeIndex);
			node.addPrevNodeIndex(index);						
		}	
		Node(char id, int index){
			this.id = id;
			this.index = index;
		}
	}	
	
	public void init(List sequences, int gapOpen, int gapExtend, ScoringMatrix matrix) {
		nodeList = new ArrayList();		// init list of nodes
		sortedList = new ArrayList();	// init ordering list of nodes
		this.sequences = sequences;
		this.gapPenalty = gapOpen;
		this.matrix = matrix;
		currentSeqId = 0;				// variable to track what sequence u'r aligning
		Sequence seq = (Sequence)sequences.get(currentSeqId);
		String seqString = seq.seqString();
		// create our initial graph
		for (int charCount = 0;charCount<seqString.length();charCount++) {			
			char curChar = seqString.charAt(charCount);			
			Node node = new Node(curChar,charCount);
			nodeList.add(node);
			node.addContainedSeqId(currentSeqId);
			if (charCount>0) {
				node.addPrevNode(charCount-1);
			}
		}
	}	
	
	private TraceBack score(int seqIndex, int nodeIndex)
	throws DataFormatException{
		TraceBack maxTraceBack = null;		
		TraceBack deletion = delete(seqIndex,nodeIndex);
		TraceBack insertion = insert(seqIndex,nodeIndex);
		TraceBack match = match(seqIndex,nodeIndex);
		if (match.propScore<0) {
			match.propScore=0;// negatives don't count for local align so let it be 0
		}
		maxTraceBack = match;		// deletion and insertion have to beat match
		if (deletion.propScore>maxTraceBack.propScore) {
			maxTraceBack = deletion;
		}
		if (insertion.propScore>maxTraceBack.propScore) {
			maxTraceBack = insertion;
		}
		return maxTraceBack;
	}
	
	private TraceBack delete(int seqIndex, int nodeIndex)
	throws DataFormatException{
		Node node = (Node)nodeList.get(nodeIndex);
		Iterator it = node.prevNodeIndices.iterator();		
		TraceBack traceBack = null, maxTraceBack = null;
		int rowIndex = seqIndex+1;		
		if (!it.hasNext()) {
			int colIndex = 0;	// if no predecessors, your only candidate traceback is the 0 column
			traceBack=new TraceBack(cells[rowIndex][colIndex],cells[rowIndex][colIndex].traceBack.propScore - gapPenalty,TraceBack.DELETE);
			maxTraceBack = traceBack;			
		}
		else {
			int propScore = 0;
			while (it.hasNext()) {
				int prevNodeIndex = ((Integer)it.next()).intValue();			
				int colIndex = prevNodeIndex+1;				
				traceBack = new TraceBack(cells[rowIndex][colIndex],cells[rowIndex][colIndex].traceBack.propScore - gapPenalty,TraceBack.DELETE);				
				if (traceBack.propScore>=propScore){
					maxTraceBack = traceBack;
					propScore=traceBack.propScore;
				}
			}
			if (maxTraceBack==null) { // if nothing is 0 or greater then just get the last one seen
				maxTraceBack = traceBack;
			}
		}
		return maxTraceBack;
	}

	private TraceBack insert(int seqIndex, int nodeIndex)
	throws DataFormatException{
		Node node = (Node)nodeList.get(nodeIndex);
		int rowIndex = seqIndex+1;
		int colIndex = nodeIndex+1;
		TraceBack traceBack = new TraceBack(cells[rowIndex-1][colIndex],cells[rowIndex-1][colIndex].traceBack.propScore-gapPenalty,TraceBack.INSERT);
		return traceBack;
	}

	
	private TraceBack match(int seqIndex, int nodeIndex)
	throws DataFormatException{		
		Node node = (Node)nodeList.get(nodeIndex);
		String currentSeqString = ((Sequence)sequences.get(currentSeqId)).seqString();
		int alignScore = matrix.getScore(currentSeqString.charAt(seqIndex),node.id);
		TraceBack maxTraceBack = null,traceBack = null;
		Iterator it = node.prevNodeIndices.iterator();		
		int rowIndex = seqIndex+1;
		if (!it.hasNext()) {  
			int colIndex = 0;	// if no predecessors, your only candidate traceback is the 0 column
			traceBack = new TraceBack(cells[rowIndex-1][colIndex],cells[rowIndex-1][colIndex].traceBack.propScore + alignScore,TraceBack.MATCH);
			maxTraceBack = traceBack;			
		}
		else {
			int propScore = 0;
			while (it.hasNext()) {
				//Logger.debugln("In match iterator");
				int prevNodeIndex = ((Integer)it.next()).intValue();				
				int colIndex = prevNodeIndex+1;
				traceBack = new TraceBack(cells[rowIndex-1][colIndex],cells[rowIndex-1][colIndex].traceBack.propScore + alignScore,TraceBack.MATCH);
				if (traceBack.propScore>=propScore){
					maxTraceBack = traceBack;
					propScore=traceBack.propScore;
				}
			}
			if (maxTraceBack==null) { // if nothing is 0 or greater then just get the last one seen
				maxTraceBack = traceBack;
			}
		}
		return maxTraceBack;
	}
	
	private int getSortedIndex(int pos) {
		return ((Integer)sortedList.get(pos)).intValue();
	}
	
	/**
	* Uses zero in-degree sorting
	*/

	private void sortGraph() {		
		sortedList = new ArrayList();
		int[] incomingEdges = new int[nodeList.size()]; 
 		Stack stack = new Stack();
 		int counter = 0;
 		for (int count=0;count<nodeList.size();count++) {
			Node node = (Node)nodeList.get(count);
			incomingEdges[count] = node.prevNodeIndices.size();
			if (incomingEdges[count]==0){
				stack.push(new Integer(count));
			}
		}
		while (!stack.empty()) {
			int headIndex = ((Integer)stack.pop()).intValue();
			sortedList.add(new Integer(headIndex));
			Node node = (Node)nodeList.get(headIndex);
			List successors = node.successorNodeIndices;			
			for (int count=0;count<successors.size();count++) {
				int nodeIndex = ((Integer)successors.get(count)).intValue();
				if (--incomingEdges[nodeIndex]==0) {
					stack.push(new Integer(nodeIndex));
				}
			}
		}
		Logger.debug("Sorted List: ");
		for (int count=0;count<sortedList.size();count++) {
			Logger.debug(sortedList.get(count)+" ");			
		}
		Logger.debugln("");
		
	}

	public int getScore(){
		throw new UnsupportedOperationException("getScore() is not supported");
	}
	
	public Alignment getAlignment() {
		throw new UnsupportedOperationException("getAlignment() is not supported");
	}
	
	// really not necessary but assumes a initialization row 
	// and col like standard local alignment
	
	private void initGrid(int y, int x) {
		int initScore = 0;
		Logger.debugln("Building a grid of dim "+y+" by "+x);
		cells = new Cell[y][x];
		//cells[0][0] = new Cell(0,0,new TraceBack(cells[0][0],initScore,TraceBack.MATCH));
		cells[0][0] = new Cell(0,0,null);
		TraceBack tb = new TraceBack(cells[0][0],initScore,TraceBack.MATCH);
		cells[0][0].traceBack = tb;
		maxCell = cells[0][0];
		for (int count=1;count<x;count++) {
			cells[0][count] = new Cell(0,count,new TraceBack(cells[0][count-1],initScore,TraceBack.DELETE));			
		}
		for (int count=1;count<y;count++) {
			cells[count][0] = new Cell(count,0,new TraceBack(cells[count-1][0],initScore,TraceBack.INSERT));			
		}
		sortGraph();
	}
	
	// fills the dynamic programming grid
	
	private void populateGrid(int maxY,int maxX)
	throws DataFormatException{
		initGrid(maxY,maxX);			
		for (int rowCount=1;rowCount<maxY;rowCount++) {			
			for (int count=0;count<sortedList.size();count++) {				
				int colCount = getSortedIndex(count) + 1;				
				int seqIndex = rowCount-1;
				int colIndex = colCount-1;
				TraceBack traceBack = score(seqIndex,colIndex);
				cells[rowCount][colCount] = new Cell(rowCount,colCount,traceBack);
				if (cells[rowCount][colCount].traceBack.propScore>maxCell.traceBack.propScore) {
					maxCell = cells[rowCount][colCount];
				}
			}
		}		
	}
	
	// insert some of the characters into the graph from the newly aligned sequence
	// that are not in the graph
	
	private int insertSeqEnds(int high, int low, int lastNodeSeen) {
		String seqString = ((Sequence)sequences.get(currentSeqId)).seqString();
		for (int curRow=high;curRow>low;curRow--) {
			int charPos = curRow-1;				
			char curChar = seqString.charAt(charPos);
			int nodePosition = nodeList.size();
			Node node = new Node(curChar,nodePosition);
			node.addContainedSeqId(currentSeqId);
			nodeList.add(node);				
			Logger.debugln("Appending inserted nodes at "+nodePosition+" with char "+curChar);
			node.addContainedSeqId(currentSeqId);
			if (curRow<high) {
				Logger.debugln("Connecting successor node of position "+(lastNodeSeen));
				node.addSuccessorNode(lastNodeSeen);
			}
			else if (lastNodeSeen!=UNDEF) {
				Logger.debugln("Chaining front node set to position "+(lastNodeSeen));
				node.addSuccessorNode(lastNodeSeen);
			}
			lastNodeSeen = nodePosition;
		}
		return lastNodeSeen;
	}
	
	// collapse the nodes to the graph
	
	private void fuseNodes(int yDim, int xDim){
		int lastNodeSeen=UNDEF;
		String seqString = ((Sequence)sequences.get(currentSeqId)).seqString();
		if (maxCell.row<yDim-1) { // if we didn't reach the end, add some tail nodes
			Logger.debugln("Inserting tail of sequence...");
			lastNodeSeen = insertSeqEnds(yDim-1,maxCell.row,lastNodeSeen);
		}
		Cell curCell = maxCell;
		int lastTraceBackType = UNDEF;
		int insertSize = 0;		
		while (curCell.traceBack.propScore>0) {
			int nodePos = curCell.col-1;
			int charPos = curCell.row-1;
			char yChar = seqString.charAt(charPos);
			char xChar = ((Node)nodeList.get(nodePos)).id;
			Logger.debug("The current cell has score "+curCell.traceBack.propScore+
				" with coordinates "+curCell.row+","+curCell.col+
				" and characters "+yChar+","+xChar+" and traceback type ");
			switch(curCell.traceBack.type) {
			case TraceBack.MATCH:
				Logger.debugln("match");
				Node node = (Node)nodeList.get(nodePos);
				node.addContainedSeqId(currentSeqId);
				if (lastNodeSeen!=UNDEF) {
					node.addSuccessorNode(lastNodeSeen);
					Logger.debugln("Connecting successor node at position "+(lastNodeSeen));
				}				
				lastNodeSeen = nodePos;
				break;			
			case TraceBack.DELETE:
				Logger.debugln("delete");
				break;			
			case TraceBack.INSERT:
				Logger.debugln("insert");
				int lastPos = nodeList.size();
				node = new Node(yChar,lastPos);
				node.addContainedSeqId(currentSeqId);
				if (lastNodeSeen!=UNDEF) {
					node.addSuccessorNode(lastNodeSeen);
					Logger.debugln("Connecting successor node at position "+(lastNodeSeen));
				}				
				nodeList.add(node);
				lastNodeSeen = lastPos;
				break;
			}
			lastTraceBackType = curCell.traceBack.type;
			curCell = curCell.traceBack.cell;			
		}
		if (curCell.row>0) { // if we didn't reach the end, add some tail nodes
			Logger.debugln("Extending head of sequence...");
			lastNodeSeen = insertSeqEnds(curCell.row,0,lastNodeSeen);
		}		
	}
	
	// marshall the clustal output
	
	public void clustalOut(OutputStream out)
	throws IOException{
		BufferedWriter writer = new BufferedWriter(new OutputStreamWriter(out));		
		sortGraph();
		int sequenceLength = sortedList.size();		
		int newStart = 0;
		int newEnd = (MAX_SEQCLUSTAL_LENGTH>sequenceLength) ? sequenceLength: MAX_SEQCLUSTAL_LENGTH;		
		do {
			for (int seqNum=0;seqNum<sequences.size();seqNum++) {
				Sequence curSeq = (Sequence)sequences.get(seqNum);
				writer.write(curSeq.getName());
				for (int count=0;count<MAX_SEQNAME_LENGTH-curSeq.getName().length();count++) {
					writer.write(" ");	// pad with spaces
				}
				for (int count=newStart;count<newEnd;count++) {
					int nodeId = getSortedIndex(count);
					Node node = (Node)nodeList.get(nodeId);
					if (node.containedSeqs.contains(new Integer(seqNum))){
						writer.write(node.id);					
					}
					else {
						writer.write(GAP_CHAR);
					}
				}
				writer.newLine();
			}
			writer.newLine();
			writer.newLine();
			newStart+=MAX_SEQCLUSTAL_LENGTH;
			newEnd=(newStart+MAX_SEQCLUSTAL_LENGTH)>sequenceLength?sequenceLength:MAX_SEQCLUSTAL_LENGTH;			
	} while ((newEnd-newStart>0) );
		writer.flush();
	}
	
	// heart of algorithm
	
	public void align()
	throws DataFormatException {
		for (int count=1;count<sequences.size();count++) {
			currentSeqId = count;
			Sequence curSeq = (Sequence)sequences.get(count);
			int yDim = curSeq.seqString().length()+1;
			int xDim = nodeList.size()+1;
			populateGrid(yDim,xDim);
			Logger.debugln("The max cell is "+maxCell.row+","+maxCell.col+
				" with score: "+maxCell.traceBack.propScore); // print to std out
			Logger.println("SCORE:"+maxCell.traceBack.propScore+" ("+curSeq.getName()+")"); // print to std out
			fuseNodes(yDim,xDim);
		}		
	}
	
	private void pad(String str) {
		int cellsize = 9;
		int padSize = cellsize - str.length();		
		for (short count=0;count<padSize;count++) {
			Logger.debug(" ");
		}
		Logger.debug(str);
	}
	
	// dump the grid
	
	public void dump() {
		pad("");
		pad("");		
		int xDim = sortedList.size()+1;
		String curStr = ((Sequence)sequences.get(currentSeqId)).seqString();
		int yDim = curStr.length();
		for (int count=0;count<sortedList.size();count++) {			
			Node node = (Node)nodeList.get(getSortedIndex(count));
		  	pad(node.id+"");
		}
		Logger.debugln("");
		
		for (int y=0;y<yDim;y++) {
			if (y>0) pad(curStr.charAt(y-1)+"");
			else pad("");
			for (int x=0;x<xDim;x++) {				
				TraceBack curTraceBack = cells[y][x].traceBack;		
				pad(curTraceBack.propScore+";"+curTraceBack.cell.col+","+curTraceBack.cell.row);
			}
			Logger.debugln("");
		}
		
	}

}
