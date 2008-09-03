/**
 * 
 */
package dr.evolution.tree;

import java.util.ArrayList;
import java.util.BitSet;
import java.util.List;

/**
 * @author Sebastian Hoehna
 *
 */
public class Clade implements Comparable<Clade> {

    private final BitSet bits;
    private final double height;
    private List<Double> heights;
    private final int size;
    private double lognormalMean;
    private double lognormalStDev;
    private double lognormalVariance;
    private boolean dirty = true;
    
    private final int PRECISION = 100;
    
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
		int setBitIndexI = 0;
		int setBitIndexJ = 0;
		do {
			setBitIndexI = bits.nextSetBit(setBitIndexI + 1);
			setBitIndexJ = clade.getBits().nextSetBit(setBitIndexJ + 1);
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
    	return compareTo(c) == 0;
    }
    
    public String toString(){
    	return bits.toString();
    }
    
    public int hashCode(){
    	return bits.hashCode();
    }
    
    public void addHeight(double height){
    	dirty = true;
    	double logHeight = Math.log(height);
    	heights.add(logHeight);
    }
    
    public double getLognormalMean(){
    	if (dirty){
    		lognormalMean = estimateLognormalMean();
    		lognormalStDev = Math.sqrt(getLognormalVariance());
    		
    		dirty = false;
    	}
    	
    	return lognormalMean;
    }
    
    private double estimateLognormalMean(){
    	double max = getLogMax();
    	double min = getLogMin();
    	
    	double slotSize = (max - min) / PRECISION;
    	int [] frequencies = new int[PRECISION];
    	
    	for (double height : heights){
    		double tmp = (height - min) / slotSize;
    		int index = (int) tmp;
    		index = (index == PRECISION ? PRECISION-1 : index);
    		
    		if (index - 2 >= 0){
        		frequencies[index-2]++; 
        		frequencies[index-1]++;    			
    		}
    		else if (index - 1 >= 0){
        		frequencies[index-1]++;    			
    		}
    		
    		if (index + 2 < PRECISION){
        		frequencies[index+2]++; 
        		frequencies[index+1]++;  
        		frequencies[index]++;   			
    		}
    		else if (index + 1 < PRECISION){ 
        		frequencies[index+1]++;  
        		frequencies[index]++; 		
    		}
    		else if (index < PRECISION){ 
        		frequencies[index]++;  			
    		}
    	}
    	
    	double [] smoothedFrequencies = new double[PRECISION];
    	smoothedFrequencies[0] = frequencies[0]/3.0;
    	smoothedFrequencies[1] = frequencies[1]/4.0;
    	for (int i=2; i<PRECISION-2; i++){
    		smoothedFrequencies[i] = frequencies[i]/5.0;
    	}
    	smoothedFrequencies[PRECISION-2] = frequencies[PRECISION-2]/4.0;
    	smoothedFrequencies[PRECISION-1] = frequencies[PRECISION-1]/3.0;
    	
    	// find slot with the max entries
    	int maxIndex = 0;
    	for (int i=0; i<PRECISION; i++){
    		if (smoothedFrequencies[i] > smoothedFrequencies[maxIndex]){
    			maxIndex = i;
    		}
    	}
    	
    	double mean = min + (slotSize * (maxIndex + 0.5));
    	
    	return mean;
    }
    
    public double getLognormalStDev(){
    	if (dirty){
    		lognormalStDev = Math.sqrt(getLognormalVariance());
    	}
    	
    	return lognormalStDev;
    }
    
    private double getLognormalVariance(){
    	if (dirty){
    		double sum = 0.0;
    		double count = 0.0;
        	for (double height : heights){
        		if (height <= lognormalMean){
        			sum += (lognormalMean - height) * (lognormalMean - height);
        			count++;
        		}
        	}    		
        	lognormalVariance = sum / count;
        	lognormalStDev = Math.sqrt(lognormalVariance);
        	
        	dirty = false;
    	}
    	
    	return lognormalVariance;
    }
    
    public void printLogHeights(){
    	for (double height : heights){
    		System.out.println(height);
    	}
    }
    
    public void printHeights(){
    	for (double height : heights){
    		System.out.println(Math.exp(height));
    	}
    }
    
    public void printLogHistogram(){
    	double max = getLogMax();
    	double min = getLogMin();
    	
    	double slotSize = (max - min) / 100.0;
    	int [] frequencies = new int[100];
    	
    	for (double height : heights){
    		double tmp = (height - min) / slotSize;
    		int index = (int) tmp;
    		index = (index == 100 ? 99 : index);
    		frequencies[index]++;
    	}
    	
    	for (int i=0; i<100; i++){
    		int freq = frequencies[i];
    		System.out.println(slotSize * i + "\t" + freq);
    	}
    }
    
    public void printHistogram(){
    	double max = Math.exp(getLogMax());
    	double min = Math.exp(getLogMin());
    	
    	double slotSize = (max - min) / 100.0;
    	int [] frequencies = new int[100];
    	
    	for (double height : heights){
    		double tmp = (Math.exp(height) - min) / slotSize;
    		int index = (int) tmp;
    		index = (index == 100 ? 99 : index);
    		frequencies[index]++;
    	}
    	
    	for (int i=0; i<100; i++){
    		int freq = frequencies[i];
    		System.out.println(slotSize * i + "\t" + freq);
    	}
    }
    
    public double getLogMax(){
    	double max = Double.NEGATIVE_INFINITY;
    	
    	for (double height : heights){
    		if (max < height){
    			max = height;
    		}
    	}
    	
    	return max;
    }
    
    public double getLogMin(){
    	double min = Double.MIN_VALUE;
    	
    	for (double height : heights){
    		if (min > height){
    			min = height;
    		}
    	}
    	
    	return min;
    }
}
