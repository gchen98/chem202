package edu.ucla.chem202.gchen.algorithm;

import edu.ucla.chem202.gchen.io.Logger;

/**
* Encapsulates finished alignment. Should be able to marshall out into various formats
*/

public class Alignment {
	
	private StringBuffer seq1positions,seq2positions;
	public static String GAP_CHAR="_";
	private String seq1id;
	private String seq2id;
	private final int CLUSTAL_COL = 50;
	private final int SEQ_ID_COL = 36;
	
	public Alignment(String seq1id, String seq2id){
		seq1positions=new StringBuffer();
		seq2positions=new StringBuffer();
		this.seq1id = seq1id;
		this.seq2id = seq2id;
	}
	
	public void addPosition(String res1,String res2) {		
		seq1positions.insert(0,res1);
		seq2positions.insert(0,res2);
	}	

	public String getSeq1() {
		return seq1positions.toString();
	}
	
	public String getSeq2() {
		return seq2positions.toString();
	}
	
	public void toClustal() {
		int seqLength = seq1positions.length();
		String seq1 = getSeq1();
		String seq2 = getSeq2();
		int lineNum = 0;
		while (lineNum!=-1) {
			int startPos = lineNum * CLUSTAL_COL;
			int length = CLUSTAL_COL;
			if (startPos+length>seqLength) {
				length = seqLength-startPos;
				lineNum=-1;
			}
			else {
				lineNum++;
			}
			//Logger.debugln("start "+startPos+" length "+length);
			Logger.print(seq1id);
			int seq1pad = SEQ_ID_COL-seq1id.length();
			for (int count=0;count<seq1pad;count++) {
				Logger.print(" ");
			}
			Logger.print(seq1.substring(startPos,startPos+length));
			Logger.println("");
			//Logger.debugln("start "+startPos+" length "+length);
			Logger.print(seq2id);
			int seq2pad = SEQ_ID_COL-seq2id.length();
			for (int count=0;count<seq2pad;count++) {
				Logger.print(" ");
			}
			Logger.print(seq2.substring(startPos,startPos+length));
			Logger.println("");
			Logger.println("");
				
		}		
	}
}