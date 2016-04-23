/*
 * Author: gchen
 * Created: Sunday, October 12, 2003 3:06:59 PM
 * Modified: Sunday, October 12, 2003 3:06:59 PM
 */

package edu.ucla.chem202.gchen.algorithm.poa;

import edu.ucla.chem202.gchen.*;
import edu.ucla.chem202.gchen.io.*;
import java.util.*;
import java.io.*;

public class POABundler
implements Constants{
	
	private int seqLength;
	private int numSequences;
	private double percIdentity;
	private List sequenceList;
	private List consensusList;	
	private Node[] nodeList;
	
	public POABundler(double percIdentity) {		
		this.percIdentity = percIdentity;
	}
	
	/**
	* lazy initialization of variables
	*/
	
	private void initVars(List sequenceList)
	throws DataFormatException{
		this.sequenceList = sequenceList;
		this.consensusList = new ArrayList();
		this.numSequences = sequenceList.size();
		Iterator iterator = sequenceList.iterator();
		if (iterator.hasNext()) {
			Sequence seq = (Sequence)iterator.next();	
			seqLength = seq.seqString().length();
		}
		else {
			throw new DataFormatException("No sequences were provided.");
		}		
		nodeList = new Node[seqLength];		
	}	
	
	/**
	* Adding a node to our graph
	*/
	
	private void addNode(char curChar,int prevPosition, int curPosition)
	throws DataFormatException{
		Node node = nodeList[curPosition];

		if (node==null) { // if we haven't seen it before, add it to our node list
			if (nodeList==null) {
				Logger.debugln("It's null 2!!");
			}
			// FIX THIS LATER!!!!!
			node = new Node(curChar,curPosition, nodeList);
			nodeList[curPosition] = node; 
		}
		else if (node.getId()!=curChar) { // if we break our assumption that only one non gap character will be seen for each column
			throw new DataFormatException("Multiple non-gap characters were seen for column: "+curPosition);
		}
		if (prevPosition!=UNDEF) {
			node.addPrevNode(prevPosition);		
		}		
	}
	
	/**
	* algorithm for building a consensus sequence
	*/
	
	private Sequence buildConsensus(String name) {
		// STEP 1: traverse the node list from left to right to find the highest scoring node
		Node maxNode = null;		
		Logger.debugln("Positions;Scores:");		
		for (int charCount = 0;charCount<seqLength;charCount++) {
			Node curNode = nodeList[charCount];			
			if (curNode!=null) {
				double prevScore = curNode.getPrevScore();				
				if (prevScore!=UNDEF) { // if this current node does indeed have a predecessor, proceed
					curNode.setScore(prevScore);
					if ((maxNode==null) || (curNode.getScore()>maxNode.getScore())) {
						maxNode = curNode;	// track the top scoring node we encounter
					}					
					Logger.debug(charCount+";"+curNode.getScore()+" ");
				}
			}
		}
		// STEP 2: begin traceback from highest scoring node
		Node curNode = maxNode;
		Logger.debugln("Traceback will begin at position "+curNode.getIndex());
		char[] consensus = new char[seqLength];
		for (int count=0;count<seqLength;count++) {
			consensus[count]=GAP_CHAR; // initialize the consensus with gaps for now.
		}
		while (curNode.getPrevScore()!=UNDEF){ // while current node has a predecessor
			consensus[curNode.getIndex()] = curNode.getId();
			curNode = nodeList[curNode.getTraceBackIndex()];  // follow the traceback path
		}
		consensus[curNode.getIndex()] = curNode.getId(); // fill in the appropriate column of the consensus array
		Logger.debugln("The consensus is "+new String(consensus));
		return new Sequence(name,name,new String(consensus));
	}
	
	/**
	* high level function that populates the graph
	*/
	
	private void buildGraph()
	throws DataFormatException{
		for (int seqcount=0;seqcount<numSequences;seqcount++) {
			Sequence seq = (Sequence)sequenceList.get(seqcount);
			Logger.debugln("Building sequence: "+seq.getName());
			String seqString = seq.seqString();
			int lastCharPosition=UNDEF; // we want to remember where our last non gap character was which will represent our previous node
			if (seqString.length()!=seqLength) {
				throw new DataFormatException("The sequences in the RC-MSA are not the same length");
			}
			for (int charCount = 0;charCount<seqLength;charCount++) {
				char curChar = seqString.charAt(charCount);				
				if (curChar!=GAP_CHAR) {										
					addNode(curChar,lastCharPosition,charCount);
					lastCharPosition = charCount;
				}				
			}
		}
	}
	
	/**
	* determine whether or not this sequence fits the consensus
	*/
		
	private boolean include(Sequence seq, Sequence cons) {		
		String seqString = seq.seqString();
		String consString = cons.seqString();
		int matches = 0; // tally the number of matches detected
		int consChars = 0; // number of non-gap characters in consensus
		for (int letterCount=0;letterCount<seqLength;letterCount++) {
			char queryChar = seqString.charAt(letterCount);
			char consChar = consString.charAt(letterCount);
			//if (consChar!=GAP_CHAR) {
			//	consChars++;
			//}
			//if ((queryChar!=GAP_CHAR) && (queryChar==consChar)) {
			//	matches++;
			if (queryChar==consChar) {
				matches++;
			}
			consChars++;
		}
		double percMatch = (double)matches/consChars*100;
		Logger.debugln("The sequence to consensus identity value is "+percMatch);		
		if ((percMatch) >= percIdentity) {
			seq.setIdentityToConsensus(percMatch); // track this for our desired output later
			Logger.debugln("Sequence: "+seq.getName()+" bundled to consensus "+cons.getDesc());
			return true;
		}
		else {
			seq.setIdentityToConsensus(0); // track this for our desired output later
		}
		return false;
	}
	
	/**
	* removing a sequence from the graph by downsizing/removing predecessor edges.
	* this is performed on sequences that have been bundled.
	*/
	
	private void removeFromGraph(Sequence seq) {
		int lastCharIndex = UNDEF; // we want to remember where we last saw a character since this represents the previous node
		String seqString = seq.seqString();			
		for (int charCount = 0; charCount<seqLength; charCount++) {
			char queryChar = seqString.charAt(charCount);			
			if ((charCount>0) && (queryChar!=GAP_CHAR)) { // don't bother for node 0. it will never have a predecessor to remove
				Node node = nodeList[charCount];
				node.removePrev(lastCharIndex);
			}
			if (queryChar!=GAP_CHAR) {
				lastCharIndex = charCount;
			}
		}
	}
	
	/**
	* this labels the appropriate sequences with a bundle ID (consensus name)
	* it returns the number of sequences "bundled"
	*/
	
	private int makeBundle(Sequence consensus) {
		int bundleSize=0;
		// step one: run through the sequences and label sequences that are to be bundled.
		Iterator it = sequenceList.iterator();
		while (sequenceList.iterator().hasNext()) {
			Sequence seq = (Sequence)sequenceList.iterator().next();			
			if ((seq.getDesc()==null) && (include(seq,consensus))){
				seq.setDesc(consensus.getDesc());
				bundleSize++;
			}
		}
		// step two: run through the newly found bundled sequences and remove these from the graph
		it = sequenceList.iterator();
		while (it.hasNext()) {
			Sequence seq = (Sequence)it.next();			
			if ((seq.getDesc()!=null) && (seq.getDesc().equals(consensus.getDesc()))){				
				removeFromGraph(seq);				
			}
		}
		return bundleSize;
	}
	
	/**
	* This returns the number of sequences in our input sequence list that have 
	* not been assigned a consensus
	*/
	
	private int getUnbundledSeqCount() {
		int count = 0;
		for (int seqNum=0;seqNum<numSequences;seqNum++) {
			Sequence curSeq = (Sequence)sequenceList.get(seqNum);
			if (curSeq.getDesc()==null) {
				count++;
			}
		}
		return count;
	}
	
	/**
	* The heart of algorithm
	*/
	
	public void build(List sequenceList) 
		throws DataFormatException{
		initVars(sequenceList);
		buildGraph();
		boolean bundled = false;
		int consensCount = 0;
		int unbundledSeqCount = UNDEF;
		do { // keep looping until at least one of two conditions is met:
			// 1) no sequences could be bundled in this iteration
			// 2) no sequences are left to be bundled.
			Sequence cons = buildConsensus(CONSENS_PREFIX+consensCount++);
			consensusList.add(cons); // store the new consensus in a separate list for output later
			if (makeBundle(cons)>0) { // make a bundle if possible
				bundled = true;
				Logger.debugln("Made a bundle this round");
			}
			else {
				bundled = false;
				Logger.debugln("Did not make a bundle this round");
			}
			unbundledSeqCount = getUnbundledSeqCount();
			Logger.debugln("Remaining sequences to join to consensus: "+unbundledSeqCount);			
		} while ((bundled) && (unbundledSeqCount>0));		
	}
	
	private void printClustalLists(BufferedWriter writer, List seqList)
	throws IOException{
		Iterator it = seqList.iterator();
		while (it.hasNext()) {
			Sequence seq = (Sequence)it.next();
			String seqName = seq.getName();
			writer.write(seqName);
			int padSize = MAX_SEQNAME_LENGTH-seqName.length();
			for (int count=0;count<padSize;count++) {
				writer.write(" ");
			}
			writer.write(seq.seqString());
			writer.newLine();
		}	
	}
	
	public void clustalOut(OutputStream out)
	throws IOException{
		BufferedWriter writer = new BufferedWriter(new OutputStreamWriter(out));
		printClustalLists(writer,sequenceList);
		printClustalLists(writer,consensusList);
		writer.flush();
	}
	
	public void outputBundleClassification() {
		Iterator it = sequenceList.iterator();
		while (it.hasNext()) {
			Sequence seq = (Sequence)it.next();
			Logger.println(seq.getName()+" "+seq.getDesc()+" "+seq.getIdentityToConsensus() );
		}
	}


	class Node
	implements Constants{
		
		
		
		private Node[] nodeList;
		private char id; // represents sequence character
		private double score = 0; // node score
		private int index = UNDEF; // the position in the sequence
		private int traceBackIndex; // the previous node that this node should trace back to
		private List prevWeights = new ArrayList();
		private List prevNodeIndices = new ArrayList();
		
		// the following setters and getters convert Objects to primitives.
		// using a List allows us to have a dynamic array at the expense of some memory
		private double getPrevWeight(int listIndex) {
			return ((Double)prevWeights.get(listIndex)).doubleValue();
		}		
		private void setPrevWeight(int listIndex,double newWeight) {
			prevWeights.set(listIndex,new Double(newWeight));
		}		
		private void addPrevWeight(double newWeight) {
			prevWeights.add(new Double(newWeight));
		}	
		private int getPrevNodeIndex(int listIndex) {
			return ((Integer)prevNodeIndices.get(listIndex)).intValue();
		}		
		private void setPrevNodeIndex(int listIndex, int newNodeIndex) {
			prevNodeIndices.set(listIndex, new Integer(newNodeIndex));
		}		
		private void addPrevNodeIndex(int newNodeIndex) {
			prevNodeIndices.add(new Integer(newNodeIndex));
		}
				
		
		Node getNode(int index) {
			return nodeList[index];
		}
		
		char getId() {
			return id;
		}
		
		double getScore() {
			return score;
		}
		
		void setScore(double score) {
			this.score = score;
		}
		
		int getIndex() {
			return index;
		}
		
		int getTraceBackIndex() {
			return traceBackIndex;
		}
		
		List getPrevNodeIndices() {
			return prevNodeIndices;
		}
		
		List getSuccessorNodeIndices() {
			return prevNodeIndices;
		}
		
		/**
		* computes the previous score using the specified algorithm
		* returns UNDEF if the node has no predecessors;
		*/
		
		double getPrevScore() {			
			double maxWeight=UNDEF;
			Node maxNode = null;
			for (int count=0;count<prevWeights.size();count++) {
				double curWeight = getPrevWeight(count);
				
				Node prevNode = getNode(getPrevNodeIndex(count));
				double prevScore = prevNode.score;
				if (curWeight>maxWeight) {					
					maxWeight = curWeight; // we want the heaviest weighted edge leading to us
					maxNode = prevNode;
				}
				else if (curWeight==maxWeight) { // if we get a tie, break it if the previous node score is biggest
					if (prevScore>maxNode.score) {
						maxNode = prevNode;
					}
				}
			}
			if (maxWeight!=UNDEF) { // this condition will only happen when we have no predecessors	
				traceBackIndex = maxNode.index;
				//Logger.debugln("score to add "+maxWeight);
				//Logger.debugln("Returning score of "+(maxNode.score+maxWeight));
				return maxNode.score+maxWeight;
			}			
			return maxWeight;
		}
		
		/**
		* adding a node to the graph
		*/
		
		void addPrevNode(int nodeIndex) {
			int counter=0;
			boolean found=false;
			// search the previous nodes
			while ((!found)&&(counter<prevNodeIndices.size())) {
				if (this.getPrevNodeIndex(counter)==nodeIndex) {  // if we see a duplicate, just increase the edge weight
					double edge = this.getPrevWeight(counter);
					edge+=WEIGHT_INC;
					this.setPrevWeight(counter,edge);					
					found=true;
				}
				counter++;
			}
			if (!found) {				
				// we we didn't find a duplicate, just add this new node
				//Logger.debugln("Added prev position: "+nodeIndex);
			    this.addPrevNodeIndex(nodeIndex);
				Node node = getNode(nodeIndex);
				this.addPrevWeight(WEIGHT_INC);				
			}			
		}
		
		/**
		* remove a node from the graph
		*/
		
		void removePrev(int prevNodeIndex) {			
			for (int count=0;count<prevNodeIndices.size();count++) {
				int curNodeIndex = getPrevNodeIndex(count);
				double curWeight = getPrevWeight(count);
				if (curNodeIndex == prevNodeIndex) {  					
					double updatedWeight = curWeight - WEIGHT_INC;
					//Logger.debugln("Updated weight is "+updatedWeight);
					if (updatedWeight>0) { // if the edge still has weight keep it and update the weight
						//Logger.debugln("Removing weight from "+curNodeIndex+" and updating to "+updatedWeight);
						setPrevWeight(count,updatedWeight);	
					}
					else {
						//Logger.debugln("Removing edge and weight");
						prevWeights.remove(count);
						prevNodeIndices.remove(count);
					}
				}
			}
		}
	
		// constructor
		
		Node(char id, int index, Node[] nodeList){
			this.nodeList = nodeList;
			if (nodeList==null) {
					Logger.debugln("It's null at position and id !!"+index+" "+id);
				}
			this.id = id;
			this.index = index;
		}
	}
}
