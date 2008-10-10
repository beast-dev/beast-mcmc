package dr.evomodel.substmodel;

import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.impl.DenseDoubleMatrix2D;
import cern.colt.matrix.linalg.Algebra;
import cern.colt.matrix.linalg.EigenvalueDecomposition;
import dr.evolution.datatype.*;
import dr.evoxml.DataTypeUtils;
import dr.inference.loggers.LogColumn;
import dr.inference.loggers.Loggable;
import dr.inference.loggers.MatrixEntryColumn;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Likelihood;
import dr.math.matrixAlgebra.Vector;
import dr.xml.*;

import java.util.logging.Logger;

/**
 * <b>A general irreversible class for any
 * data type; allows complex eigenstructures.
 * Eigendecomposition is done in colt.  Finite-time
 * transition probabilities are caches in native memory</b>
 *
 * @author Marc Suchard
 */

public class NativeSubstitutionModel extends AbstractSubstitutionModel
    implements  Loggable, Likelihood {

	public static final String NATIVE_SUBSTITUTION_MODEL = "nativeSubstitutionModel";
	public static final String RATES = "rates";
	public static final String ROOT_FREQUENCIES = "rootFrequencies";
	public static final String INDICATOR = "rateIndicator";

//	private static final double minProb = Property.DEFAULT.tolerance();

	public NativeSubstitutionModel(String name, DataType dataType,
	                                FrequencyModel rootFreqModel, Parameter parameter) {

		super(name, dataType, rootFreqModel);
		this.infinitesimalRates = parameter;

		rateCount = stateCount * (stateCount - 1);

		stateCountSquared = stateCount * stateCount;

		if (rateCount != infinitesimalRates.getDimension()) {
			throw new RuntimeException("Dimension of '" + infinitesimalRates.getId() + "' ("
					+ infinitesimalRates.getDimension() + ") must equal " + rateCount);
		}

		stationaryDistribution = new double[stateCount];
		storedStationaryDistribution = new double[stateCount];

		addParameter(infinitesimalRates);

        illConditionedProbabilities = new double[stateCount*stateCount];

		Logger.getLogger("dr.evomodel.substmodel").info("Trying a native substitution model. Best of luck to you!") ;

    }

    public Model getModel() { return this; }

    public double getLogLikelihood() { // this can be used as a 0/1 prior on complex matrices
	if( isComplex ) 
	    return Double.NEGATIVE_INFINITY;
	else
	    return 0.0;
    }

    public void makeDirty() { }

	protected void handleModelChangedEvent(Model model, Object object, int index) {
		if (model == freqModel)
			return; // freqModel only affects the likelihood calculation at the tree root
		super.handleModelChangedEvent(model, object, index);
	}

	protected void restoreState() {


		double[] tmp3 = storedStationaryDistribution;
		storedStationaryDistribution = stationaryDistribution;
		stationaryDistribution = tmp3;

		normalization = storedNormalization;

		updateMatrix = storedUpdateMatrix;

		nativeRestoreState();
		
		//nativeStoreState();
		
		if (DEBUG_PRINT)
			System.err.println("Restore SM");
		
//		if( firstCall ) {
//		    firstCall = false;
//		    migrateThread();
//		}

	}

    protected void migrateThread() {
		updateMatrix = true;
		eigenInitialised = false;
		setupMatrix();
	}


    protected boolean firstCall = true;

	protected void storeState() {

		storedUpdateMatrix = updateMatrix;

		System.arraycopy(stationaryDistribution, 0, storedStationaryDistribution, 0, stateCount);
		storedNormalization = normalization;

		if (DEBUG_PRINT)
			System.err.println("Storing SM");
		
		if( firstCall ) {
			firstCall = false;
			migrateThread();
		}		
	
		nativeStoreState();


	}

	public void getTransitionProbabilities(double distance, long ptrMatrix) {

//		synchronized (this) {
			if (updateMatrix) {
				setupMatrix();
			}
//		}
//			if ( isComplex )
//			    nativeGetTransitionProbabilities(A
		nativeGetTransitionProbabilities(ptrEigenDecomposition, distance, ptrMatrix, ptrExtraCache, stateCount);

	}

	public void getTransitionProbabilities(double distance, double[] matrix) {
			throw new RuntimeException("Should only get here in debug mode");

	}



	public double[] getStationaryDistribution() {
		return stationaryDistribution;
	}

	protected void computeStationaryDistribution() {

		int i;
		final int end = stateCount - 1;
		for (i = 0; i < end; i++) {
			if (Math.abs(Eval[i]) < 1E-12)
				break;
		}

		double total = 0.0;

		for (int k = 0; k < stateCount; k++) {
			double value = Evec[k][i];
			total += value;
			stationaryDistribution[k] = value;
		}

		for (int k = 0; k < stateCount; k++)
			stationaryDistribution[k] /= total;

	}


	protected double[] getRates() {
		return infinitesimalRates.getParameterValues();
	}


	protected void setupMatrix() {

		if (!eigenInitialised) {

			amat = new double[stateCount][stateCount];
			q = new double[stateCount][stateCount];

			nativeInitialiseEigen();
			eigenInitialised = true;
			updateMatrix = true;

		}

		int i, j, k = 0;

		double[] rates = getRates();

		// Set the instantaneous rate matrix
		for (i = 0; i < stateCount; i++) {
			for (j = 0; j < stateCount; j++) {
				if (i != j)
					amat[i][j] = rates[k++];
			}
		}

		makeValid(amat, stateCount);

		// compute eigenvalues and eigenvectors
		EigenvalueDecomposition eigenDecomp = new EigenvalueDecomposition(new DenseDoubleMatrix2D(amat));

		DoubleMatrix2D eigenV = eigenDecomp.getV();
		DoubleMatrix1D eigenVReal = eigenDecomp.getRealEigenvalues();
		DoubleMatrix1D eigenVImag = eigenDecomp.getImagEigenvalues();

        DoubleMatrix2D eigenVInv;

        try {
            eigenVInv = alegbra.inverse(eigenV);
        } catch (IllegalArgumentException e) {
            wellConditioned = false;
            return;
//            throw e;
        }

        // fill AbstractSubstitutionModel parameters

		Ievc = eigenVInv.toArray();
		Evec = eigenV.toArray();
		Eval = eigenVReal.toArray();
		EvalImag = eigenVImag.toArray();

		checkComplexSolutions();

		// compute normalization and rescale eigenvalues

		computeStationaryDistribution();

		double subst = 0.0;

		for (i = 0; i < stateCount; i++)
			subst += -amat[i][i] * stationaryDistribution[i];

		normalization = subst;

		for (i = 0; i < stateCount; i++) {
			Eval[i] /= subst;
			EvalImag[i] /= subst;
		}

		nativeSetup(ptrEigenDecomposition,Ievc,Evec,Eval,EvalImag,stateCount);

		//int length = 2 * (stateCount*stateCount + stateCount);
		//double[] tmp = new double[length];
		//getNativeMemoryArray(ptrCache,0,tmp,0,length);
		//System.err.println("Cache = "+new Vector(tmp));
		//System.exit(-1);

        wellConditioned = true;
        updateMatrix = false;
	}



	protected void checkComplexSolutions() {
		boolean complex = false;
		for (int i = 0; i < stateCount && !complex; i++) {
			if (EvalImag[i] != 0)
				complex = true;
		}
		isComplex = complex;
		//	if( isComplex && usingGPU ) {
		//  System.err.println("No complex GPU support!!!");
		    //System.exit(-1);
		//}
	}

	public boolean getIsComplex() {
		return isComplex;
	}

	protected void frequenciesChanged() {
	}

	protected void ratesChanged() {
	}

	protected void setupRelativeRates() {
	}

	protected Parameter infinitesimalRates;

	public LogColumn[] getColumns() {

		LogColumn[] columnList = new LogColumn[stateCount * stateCount];
		int index = 0;
		for (int i = 0; i < stateCount; i++) {
			for (int j = 0; j < stateCount; j++)
				columnList[index++] = new MatrixEntryColumn(getId(), i, j, amat);
		}
		return columnList;
	}


	public static void main(String[] arg) {

	    //	Parameter rates = new Parameter.Default(new double[]{5.0, 1.0, 1.0, 0.1, 5.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0});
		//			Parameter rates = new Parameter.Default(new double[] {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0});
//		Parameter rates = new Parameter.Default(new double[] {1.0, 1.0});
		
//		double[] r = new double[4032];
		double[] r = new double[3660];
		for(int i=0; i<r.length; i++)
			r[i] = 1.0;
		Parameter rates = new Parameter.Default(r);

		NativeSubstitutionModel substModel = new NativeSubstitutionModel("test",
//				TwoStates.INSTANCE,
//				Nucleotides.INSTANCE,
				Codons.UNIVERSAL,
				null,
				rates);

		ComplexSubstitutionModel testModel = new ComplexSubstitutionModel("test",
										  Codons.UNIVERSAL, null, rates);

		int states = substModel.getDataType().getStateCount();
		
		System.err.println("States = "+states);

		double[] finiteTimeProbs = new double[substModel.getDataType().getStateCount() * substModel.getDataType().getStateCount()];
		double time = 1.0;
//		substModel.getTransitionProbabilities(time, finiteTimeProbs);

		long ptrMatrix = substModel.allocateNativeMemoryArray(substModel.getDataType().getStateCount()*substModel.getDataType().getStateCount());


				int len = 2 * (states*states + states);
				//		double[] mem = new double[len];
				//		substModel.getNativeMemoryArray(substModel.ptrCache,0,mem,0,len);
				//		System.out.println("Cache = "+new Vector(mem));
		
		if(substModel.getIsComplex()) {
		    System.err.println("GPU current does not support complex stuff.");
		    System.exit(-1);

		}

		substModel.getTransitionProbabilities(time, ptrMatrix);

		System.out.println("Results:");
//		System.out.println(new Vector(finiteTimeProbs));
//		finiteTimeProbs = new double[substModel.getDataType().getStateCount()*substModel.getDataType().getStateCount()];

		substModel.getNativeMemoryArray(ptrMatrix,0,finiteTimeProbs,0,substModel.getDataType().getStateCount()*substModel.getDataType().getStateCount());
		double[] print = new double[10];
		
		int count = 0;
		for(int i=0; i<finiteTimeProbs.length; i++) {
			if(Double.isNaN(finiteTimeProbs[i]))
				System.err.println("BAD NUMBER!!!");
			if ( finiteTimeProbs[i] == 1E-10)
				count++;
		}
		
		System.out.println("Zeroes = "+count);
		//System.out.println(new Vector(finiteTimeProbs));
		//System.exit(0);
		for(int i=0; i<10; i++)
//			print[i] = finiteTimeProbs[finiteTimeProbs.length-i-1];
			print[i] = finiteTimeProbs[states-i-1];
		System.out.println(new Vector(print));
		for(int i=0; i<10; i++)
//			print[i] = finiteTimeProbs[finiteTimeProbs.length-i-1];
			print[i] = finiteTimeProbs[2*states+i];
		System.out.println(new Vector(print));
		for(int i=0; i<10; i++)
			print[i] = finiteTimeProbs[finiteTimeProbs.length-i-1];
		System.out.println(new Vector(print));
		
		System.out.println();
		testModel.getTransitionProbabilities(time,finiteTimeProbs);
		for(int i=0; i<10; i++)
//			print[i] = finiteTimeProbs[finiteTimeProbs.length-i-1];
			print[i] = finiteTimeProbs[states-i-1];	
		System.out.println(new Vector(print));
		for(int i=0; i<10; i++)
//			print[i] = finiteTimeProbs[finiteTimeProbs.length-i-1];
			print[i] = finiteTimeProbs[2*states+i];	
		System.out.println(new Vector(print));
		for(int i=0; i<10; i++)
			print[i] = finiteTimeProbs[finiteTimeProbs.length-i-1];
		System.out.println(new Vector(print));

	}



	public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

		public String getParserName() {
			return NATIVE_SUBSTITUTION_MODEL;
		}

		public Object parseXMLObject(XMLObject xo) throws XMLParseException {

//			DataType dataType = null;
//
//			if (xo.hasAttribute(DataType.DATA_TYPE)) {
//				String dataTypeStr = xo.getStringAttribute(DataType.DATA_TYPE);
//				if (dataTypeStr.equals(Nucleotides.DESCRIPTION)) {
//					dataType = Nucleotides.INSTANCE;
//				} else if (dataTypeStr.equals(AminoAcids.DESCRIPTION)) {
//					dataType = AminoAcids.INSTANCE;
//				} else if (dataTypeStr.equals(Codons.DESCRIPTION)) {
//					dataType = Codons.UNIVERSAL;
//				} else if (dataTypeStr.equals(TwoStates.DESCRIPTION)) {
//					dataType = TwoStates.INSTANCE;
//				}
//			}
//
//			if (dataType == null) dataType = (DataType) xo.getChild(DataType.class);
			
	         DataType dataType = DataTypeUtils.getDataType(xo);

			XMLObject cxo = (XMLObject) xo.getChild(RATES);

			Parameter ratesParameter = (Parameter) cxo.getChild(Parameter.class);

			int rateCount = (dataType.getStateCount() - 1) * dataType.getStateCount();

			if (ratesParameter.getDimension() != rateCount) {
				throw new XMLParseException("Rates parameter in " + getParserName() + " element should have " + (rateCount) + " dimensions.  However parameter dimension is " + ratesParameter.getDimension());
			}


			cxo = (XMLObject) xo.getChild(ROOT_FREQUENCIES);
			FrequencyModel rootFreq = (FrequencyModel) cxo.getChild(FrequencyModel.class);

			if (dataType != rootFreq.getDataType()) {
				throw new XMLParseException("Data type of " + getParserName() + " element does not match that of its rootFrequencyModel.");
			}

			Parameter indicators = null;

			if (xo.hasChildNamed(INDICATOR)) {
				indicators = (Parameter) ((XMLObject) xo.getChild(INDICATOR)).getChild(Parameter.class);
				if (ratesParameter.getDimension() != indicators.getDimension())
					throw new XMLParseException("Rate parameter dimension must match indicator parameter dimension");
			}

			if (indicators == null)
				return new NativeSubstitutionModel(xo.getId(), dataType, rootFreq, ratesParameter);
			else
				return new NativeSubstitutionModel(xo.getId(), dataType, rootFreq, ratesParameter);

		}

		//************************************************************************
		// AbstractXMLObjectParser implementation
		//************************************************************************

		public String getParserDescription() {
			return "A general reversible model of sequence substitution for any data type with stochastic variable selection.";
		}

		public Class getReturnType() {
			return SubstitutionModel.class;
		}

		public XMLSyntaxRule[] getSyntaxRules() {
			return rules;
		}

		private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
				 new XORRule(
	                        new StringAttributeRule(DataType.DATA_TYPE, "The type of sequence data",
	                                DataType.getRegisteredDataTypeNames(), false),
	                        new ElementRule(DataType.class)
	                ),
				new ElementRule(ROOT_FREQUENCIES, FrequencyModel.class),
				new ElementRule(RATES,
						new XMLSyntaxRule[]{
								new ElementRule(Parameter.class)}
				),
				new ElementRule(INDICATOR,
						new XMLSyntaxRule[]{
								new ElementRule(Parameter.class)
						}),
		};

	};

	private boolean isComplex = false;
	private double[] stationaryDistribution = null;
	private double[] storedStationaryDistribution;
	private Double normalization;
	private Double storedNormalization;

	protected double[] EvalImag;

    private boolean wellConditioned = true;
    private double[] illConditionedProbabilities;

    private static final Algebra alegbra = new Algebra();

	protected void nativeStoreState() {
		copyNativeMemoryArray(ptrEigenDecomposition,ptrStoredEigenDecomposition, 2*(stateCount + stateCountSquared));
	}

	protected void nativeRestoreState() {
		long tmpPtr = ptrEigenDecomposition;
		ptrEigenDecomposition = ptrStoredEigenDecomposition;
		ptrStoredEigenDecomposition = tmpPtr;
	}

	protected boolean DEBUG_PRINT = false;
	int count = 0;
	
	
	
	
	protected void nativeInitialiseEigen() {

		if (DEBUG_PRINT)
			System.err.println("Allocating substitution matrix memory");
		ptrEigenDecomposition = allocateNativeMemoryArray(2*stateCountSquared + 2*stateCount);
		ptrStoredEigenDecomposition = allocateNativeMemoryArray(2*stateCountSquared + 2*stateCount);
		if( stateCount > 32)
			ptrExtraCache = allocateNativeMemoryArray(stateCountSquared);
		if (DEBUG_PRINT) {
			System.err.println("Done allocating matrix memory");
			System.err.println("Cache = "+ptrEigenDecomposition+" to "+(ptrEigenDecomposition +(2*stateCountSquared + 2*stateCount)*getNativeRealSize()));	
			if (++count == 3)
				System.exit(-1);
		}
	}


	protected native void nativeSetup(long ptr, double[][] Ievc, double[][] Evec, double[] Eval,
	                                  double[] EvalImag, int stateCount);

	protected native long getContext();
	
	protected native void migrate();
	
	protected native long allocateNativeMemoryArray(int length);
	protected native void copyNativeMemoryArray(long from, long to, int length);

	protected native int getNativeRealSize();

	protected native void getNativeMemoryArray(long from, int fromOffset, double[] to, int toOffset, int length);
	
	protected native void nativeGetTransitionProbabilities(long ptrEigendecomposition, double distance,
	                                                           long ptrMatrix, long extraCache, int stateCount);

    protected native void nativeGetTransitionProbabilitiesComplex(long ptrCache, double distance, long ptrMatrix, int stateCount);

    private static boolean usingGPU = false;

	static {

		try {
			System.loadLibrary("GPUMemoryLikelihoodCore");
			usingGPU = true;
		} catch (UnsatisfiedLinkError e) {
		    try {
		    
			System.loadLibrary("NativeMemoryLikelihoodCore");
		    } catch (UnsatisfiedLinkError e2) {
		    }
		}
	}

	int stateCountSquared;

	long ptrEigenDecomposition;
	long ptrExtraCache;
//  long ptrIevc;
//	long ptrEvec;
//	long ptrEval;
//	long ptrEvalImag;

	long ptrStoredEigenDecomposition;
//	long ptrStoredIevc;
//	long ptrStoredEvec;
//	long ptrStoredEval;
//	long ptrStoredEvalImag;

}