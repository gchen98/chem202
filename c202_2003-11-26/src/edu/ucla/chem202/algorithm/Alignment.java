package edu.ucla.chem202.algorithm;

import edu.ucla.chem202.io.Logger;
import edu.ucla.chem202.DataFormatException;
import java.io.*;
import java.util.Iterator;
import java.util.List;
import java.util.ArrayList;

/**
* Encapsulates finished alignment. Should be able to marshall out into various formats
*/

public class Alignment{
	
	public static final int UNDEF = -1;
	public static final char GAP_CHAR = '-';	
	public static final String CLUSTAL_HEADER = "CLUSTAL W (1.74) multiple sequence alignment";
	public static final int MAX_SEQNAME_LENGTH = 36;
	public static final int MAX_SEQCLUSTAL_LENGTH = 50;

	
	protected int seqLength;
	protected List sequences = new ArrayList();	
	
	public Alignment(List sequences)
	throws DataFormatException {
		Iterator it = sequences.iterator();
		while (it.hasNext()) {
			addSequence((Sequence)it.next());
		}
	}

	protected char getCharAt(int seqNum, int charPosition) {
		Sequence seq = (Sequence)sequences.get(seqNum);
		String rawStr = seq.getRawSeq();
		return rawStr.charAt(charPosition);
	}
	
	public void addSequence(Sequence seq)
	throws DataFormatException{
		sequences.add(seq);
		if ((sequences.size()>1) && (seqLength!=seq.getLength())) {
			throw new DataFormatException("The newly added sequence does not match in size to the previous one");
		}
		seqLength = seq.getLength();
	}

	public void clustalOut(OutputStream out)
	throws IOException{
		BufferedWriter writer = new BufferedWriter(new OutputStreamWriter(out));
		writer.write(CLUSTAL_HEADER);
		writer.newLine();
		writer.newLine();
		int newStart = 0;
		int newEnd = (MAX_SEQCLUSTAL_LENGTH>seqLength) ? seqLength: MAX_SEQCLUSTAL_LENGTH;		
		do {
			for (int seqNum=0;seqNum<sequences.size();seqNum++) {
				Sequence curSeq = (Sequence)sequences.get(seqNum);
				writer.write(curSeq.getName());
				for (int count=0;count<MAX_SEQNAME_LENGTH-curSeq.getName().length();count++) {
					writer.write(" ");	// pad with spaces
				}
				String rawSeq = curSeq.getRawSeq();
				writer.write(rawSeq.substring(newStart,newEnd));				
				writer.newLine();
			}
			writer.newLine();
			writer.newLine();
			newStart+=MAX_SEQCLUSTAL_LENGTH;
			newEnd=(newStart+MAX_SEQCLUSTAL_LENGTH)>seqLength?seqLength:MAX_SEQCLUSTAL_LENGTH;			
	} while ((newEnd-newStart>0) );
		writer.flush();
	}
}
