/**
 * 
 */
package dr.evolution.tree;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.List;

import dr.math.MathUtils;
import dr.math.distributions.Distribution;
import dr.math.distributions.NormalDistribution;

/**
 * @author Sebastian Hoehna
 *
 */
public class Clade implements Comparable<Clade> {

    private final BitSet bits;
    private double height;
    private List<Double> heights;
    private final int size;
    
    public Clade(final BitSet bits, final double height) {
        this.bits = bits;
        this.height = height;
        size = bits.cardinality();
        heights = new ArrayList<Double>();
    }

    public BitSet getBits() {
        return bits;
    }

    public double getHeight() {
        return height;
    }
    
    public double getHeight(int index) {
        return heights.get(index);
    }
    
    public int getSampleCount(){
    	return heights.size();
    }

    public int getSize() {
        return size;
    }

    public int compareTo(Clade clade) {
		int setBitIndexI = -1;
		int setBitIndexJ = -1;
		BitSet otherBits = clade.getBits();
		do {
			setBitIndexI = bits.nextSetBit(setBitIndexI + 1);
			setBitIndexJ = otherBits.nextSetBit(setBitIndexJ + 1);
		} while (setBitIndexI == setBitIndexJ && setBitIndexI != -1);

		return (setBitIndexI < setBitIndexJ ? -1
				: (setBitIndexI > setBitIndexJ ? 1 : 0));
	}
    
    public boolean equals(Object o){
    	if (o instanceof Clade){
    		return equals((Clade)o);
    	}
    	return false;
    }
    
    public boolean equals(Clade c){
    	if (size != c.getSize()){
    		return false;
    	}
    	return compareTo(c) == 0;
    }
    
    public String toString(){
    	return bits.toString();
    }
    
    public int hashCode(){
    	return bits.hashCode();
    }
    
    public void addHeight(double height){
    	heights.add(height);
    	setHeight(height);
    }
    
    public void setHeight(double height){
    	this.height = height;
    }
    
}
