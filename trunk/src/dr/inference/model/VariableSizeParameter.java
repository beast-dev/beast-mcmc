package dr.inference.model;

public class VariableSizeParameter extends Parameter.Abstract {


	public VariableSizeParameter(int dimension) {
		this(dimension, 0.0);
	}

	public VariableSizeParameter(double initialValue) {
		values = new double[1];
		values[0] = initialValue;
		this.bounds = null;
	}

	public VariableSizeParameter(int dimension, double initialValue) {
		values = new double[dimension];
		for (int i = 0; i < dimension; i++) {
			values[i] = initialValue;
		}
		this.bounds = null;
	}

	public VariableSizeParameter(double[] values) {
		this.values = new double[values.length];
		for (int i = 0; i < values.length; i++) {
			this.values[i] = values[i];
		}
	}

	public void addBounds(Bounds boundary) {
		if (bounds == null) {
			bounds = new IntersectionBounds(getDimension());
		}
		bounds.addBounds(boundary);

		// can't change dimension after bounds are added!
		// TODO what to do about bounds???
		hasBeenStored = true;
	}

	//********************************************************************
	// GETTERS
	//********************************************************************

	public int getDimension() {
		return values.length;
	}

	public double getParameterValue(int i) {
		return values[i];
	}

	/**
	 * Defensively returns copy of parameter array.
	 *
	 * @return
	 */
	public final double[] getParameterValues() {

		double[] copyOfValues = new double[values.length];
		System.arraycopy(values, 0, copyOfValues, 0, values.length);
		return copyOfValues;
	}

	public double[] inspectParametersValues() {
		throw new RuntimeException("Not yet implemented.");
	}

	public Bounds getBounds() {
		if (bounds == null) {
//				bounds = new IntersectionBounds(getDimension());
			throw new NullPointerException(getParameterName() + " parameter: Bounds not set");
		}
		return bounds;
	}

	public String getParameterName() {
		return getId();
	}

	//********************************************************************
	// SETTERS
	//********************************************************************

	/**
	 * Can only be called before store is called. If it results in new
	 * dimensions, then the value of the first dimension is copied into the new dimensions.
	 */
	public void setDimension(int dim) {
		//if (!hasBeenStored) {
		double[] newValues = new double[dim];
		for (int i = 0; i < values.length; i++) {
			newValues[i] = values[i];
		}
		for (int i = values.length; i < newValues.length; i++) {
			newValues[i] = values[0];
		}
		values = newValues;
		//} else throw new RuntimeException("Can't change dimension after store has been called!");
	}

	/**
	 * Adds an extra dimension to the end of values
	 *
	 * @param value value to save at end of new array
	 */
	public void addDimension(double value) {
		final int n = values.length;
		double[] newValues = new double[n + 1];
		for (int i = 0; i < n; i++)
			newValues[i] = values[i];
		newValues[n] = value;
		values = newValues;
	}


	/**
	 * Removes a single dimension from value array
	 *
	 * @param index Index of dimension to lose
	 */
	public void removeDimension(int index) {
		final int n = values.length;
		double[] newValues = new double[n - 1];
		for (int i = 0; i < index; i++)
			newValues[i] = values[i];
		for (int i = index; i < n; i++)
			newValues[i - 1] = values[i];
		values = newValues;
	}

	public void setParameterValue(int i, double val) {
		values[i] = val;

		fireParameterChangedEvent(i);
	}

	/**
	 * Sets the value of the parameter without firing a changed event.
	 *
	 * @param i
	 * @param val
	 */
	public void setParameterValueQuietly(int i, double val) {
		values[i] = val;
	}

	protected final void storeValues() {
		//hasBeenStored = true;
		if (storedValues == null || storedValues.length != values.length) {
			storedValues = new double[values.length];
		}
		System.arraycopy(values, 0, storedValues, 0, values.length);
	}

	protected final void restoreValues() {

		//swap the arrays
		double[] temp = storedValues;
		storedValues = values;
		values = temp;

		//if (storedValues != null) {
		//	System.arraycopy(storedValues, 0, values, 0, values.length);
		//} else throw new RuntimeException("restore called before store!");
	}

	/**
	 * Nothing to do
	 */
	protected final void acceptValues() {
	}

	protected final void adoptValues(Parameter source) {

		if (getDimension() != source.getDimension()) {
			throw new RuntimeException("The two parameters don't have the same number of dimensions");
		}

		for (int i = 0, n = getDimension(); i < n; i++) {
			values[i] = source.getParameterValue(i);
		}
	}

	private double[] values;

	private double[] storedValues;

	private boolean hasBeenStored = false;
	private IntersectionBounds bounds = null;
}

