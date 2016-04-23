/*
 * Author: gchen
 * Created: Sunday, October 12, 2003 3:06:59 PM
 * Modified: Sunday, October 12, 2003 3:06:59 PM
 */

/**
* Handles parsing responsibilies from inputstreams
*/


package edu.ucla.chem202.io;

import edu.ucla.chem202.*;
import edu.ucla.chem202.algorithm.*;
	
import java.io.*;
import java.util.*;


public class Parser
{
	/**
	* This is a Clustal MSA
	* @return: List of sequences
	*/
	
	public List parseClustal(InputStream in)
	throws IOException{	
		List sequenceNames = new ArrayList();		
		Map sequenceMap = new HashMap();
		BufferedReader bufRead=new BufferedReader(new InputStreamReader(in));		
		String seqName = null, seqDesc = null, newLine = null;
		StringBuffer seqBuffer = new StringBuffer();		
		while ((newLine=bufRead.readLine())!=null) {
			if ((newLine.length()>0) && (!newLine.startsWith(" "))) {
				int seqNameEndIndex = newLine.indexOf(" ");
				if (seqNameEndIndex>=-1) {					
					seqName = newLine.substring(0,seqNameEndIndex);
					if (!sequenceNames.contains(seqName)) {
						sequenceNames.add(seqName);
					}
					StringBuffer buf = (StringBuffer)sequenceMap.get(seqName);
					if (buf==null) buf = new StringBuffer();					
					String sequenceStr = newLine.substring(seqNameEndIndex+1,newLine.length());
					buf.append(sequenceStr.replaceAll("\\s",""));
					sequenceMap.put(seqName,buf);
				}				
			}			
		}
		// Capture the last sequence
		Iterator it = sequenceNames.iterator();
		List sequences = new ArrayList();
		while (it.hasNext()) {
			seqName = it.next().toString();
			String seqFull = ((StringBuffer)sequenceMap.get(seqName)).toString();
			sequences.add(new Sequence(seqName,null,seqFull));
		}		
		return sequences;
	}
	
	
	/**
	* This must be a FASTA formatted file
	* Returns a list of sequences
	*/
	
	public List parseFasta(InputStream in)
	throws IOException{
		List sequences = new ArrayList();
		BufferedReader bufRead=new BufferedReader(new InputStreamReader(in));		
		String seqName = null, seqDesc = null, newLine = null;
		StringBuffer seqBuffer = new StringBuffer();		
		while ((newLine=bufRead.readLine())!=null) {
			if (newLine.startsWith(">")) {				
				if (seqName!=null) { // Create the previous sequence
					Sequence sequence = new Sequence(seqName, seqDesc, seqBuffer.toString());
					sequences.add(sequence);
				}				
				int delim = newLine.indexOf(" "); // Prepare the new sequence
				if (delim>1) {
					seqName = newLine.substring(1,delim);
					seqDesc = newLine.substring(delim+1,newLine.length());				
				}
				else {
					seqName = newLine.substring(1,newLine.length());					
				}
				seqBuffer = new StringBuffer();
			}
			else {
				seqBuffer.append(newLine.toUpperCase());
			}
		}
		// Capture the last sequence
		if (seqName!=null) {
			Sequence sequence = new Sequence(seqName, seqDesc, seqBuffer.toString());
			sequences.add(sequence);
		}
		return sequences;
	}
	
	/**
	* parse matrix format.  must be symmetrical
	*/
	
	public ScoringMatrix parseMatrix(InputStream in)
	throws IOException,DataFormatException{		
		BufferedReader bufRead=new BufferedReader(new InputStreamReader(in));		
		String newLine = null;
		newLine=bufRead.readLine().toUpperCase().trim();
		String[] header = newLine.split("\\s+");
		char[] headerChars = new char[header.length];
		for (int count=0;count<header.length;count++) {
			headerChars[count] = header[count].charAt(0);
		}
		ScoringMatrix scoringMatrix = new ScoringMatrix(headerChars);
		while ((newLine=bufRead.readLine())!=null) {			
			newLine = newLine.toUpperCase().trim();
			if (newLine.length()>0) {
				scoringMatrix.parseLine(newLine.split("\\s+"));
			}
		}
		return scoringMatrix;	
	}

}
