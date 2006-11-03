// Trace.java

package dr.app.coalgen;



/**
 * A simple class that stores a trace for a single statistic
 *
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @version $Id: Trace.java,v 1.2 2004/10/01 22:40:02 alexei Exp $
 */
public class Trace {

	public static final int INITIAL_SIZE = 1000;
	public static final int INCREMENT_SIZE = 1000;
	
	public Trace(String name) {

		this.name = name;
	}
	
	public Trace(String name, double[] values) {

		this.values = new double[values.length];
		System.arraycopy(values, 0, this.values, 0, values.length);
	}
	
	/**
	 * add a value
	 */
	public void add(double value) {
	
		if (valueCount == values.length) {
			double[] newValues = new double[valueCount + INCREMENT_SIZE];
			System.arraycopy(values, 0, newValues, 0, values.length);
			values = newValues;
		}
		
		values[valueCount] = value;
		valueCount++;
	}
	
	/**
	 * add all the values in an array of doubles
	 */
	public void add(double[] values) {
		for (int i = 0; i < values.length; i++) {
			add(values[i]);
		}
	}
	
	public int getCount() { return valueCount; }
	public double getValue(int index) { return values[index]; }
	public void getValues(int start, double[] destination) { 
		System.arraycopy(values, start, destination, 0, valueCount - start);
	}
	public void getValues(int start, double[] destination, int offset) { 
		System.arraycopy(values, start, destination, offset, valueCount - start);
	}
	      
	public String getName() { return name; }
	
	//************************************************************************
	// private methods
	//************************************************************************
	
	private double[] values = new double[INITIAL_SIZE];
	private int valueCount = 0;
	private String name;
}