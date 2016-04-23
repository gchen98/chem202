/*
 * Author: gchen
 * Created: Sunday, October 12, 2003 3:06:59 PM
 * Modified: Sunday, October 12, 2003 3:06:59 PM
 */

/**
* Encapsulate a sequence and its attributes
*/

package edu.ucla.chem202.algorithm;

import edu.ucla.chem202.algorithm.Alignment;
	
public class Sequence{
	
	private String name;
	private String desc;
	private String seq;
	
	public Sequence(String name, String desc, String seq) {
		this.name = name;
		this.desc = desc;
		seq.replace('_',Alignment.GAP_CHAR); // normalize errors in datainput
		this.seq = seq.toUpperCase();
	}	
	
	public String getName(){
		return name;
	}
	
	public String getDesc() {
		return desc;
	}	
	
	public String getRawSeq() {
		return seq;
	}
	
	public int getLength() {
		return seq.length();
	}
}
