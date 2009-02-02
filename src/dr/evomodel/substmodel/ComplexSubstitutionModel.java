package dr.evomodel.substmodel;

import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.impl.DenseDoubleMatrix2D;
import cern.colt.matrix.linalg.Algebra;
import cern.colt.matrix.linalg.EigenvalueDecomposition;
import cern.colt.matrix.linalg.Property;
import dr.evolution.datatype.DataType;
import dr.evolution.datatype.Nucleotides;
import dr.evoxml.DataTypeUtils;
import dr.inference.loggers.LogColumn;
import dr.inference.loggers.Loggable;
import dr.inference.loggers.MatrixEntryColumn;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.math.matrixAlgebra.Matrix;
import dr.math.matrixAlgebra.Vector;
import dr.xml.*;

/**
 * <b>A general irreversible class for any
 * data type; allows complex eigenstructures.</b>
 *
 * @author Marc Suchard
 */

public class ComplexSubstitutionModel extends AbstractSubstitutionModel implements Loggable {

    public static final String COMPLEX_SUBSTITUTION_MODEL = "complexSubstitutionModel";
    public static final String RATES = "rates";
    public static final String ROOT_FREQUENCIES = "rootFrequencies";
    public static final String INDICATOR = "rateIndicator";


    public ComplexSubstitutionModel(String name, DataType dataType,
                                    FrequencyModel rootFreqModel, Parameter parameter) {

        super(name, dataType, rootFreqModel);
        this.infinitesimalRates = parameter;

        rateCount = stateCount * (stateCount - 1);

        if (rateCount != infinitesimalRates.getDimension()) {
            throw new RuntimeException("Dimension of '" + infinitesimalRates.getId() + "' ("
                    + infinitesimalRates.getDimension() + ") must equal " + rateCount);
        }

        stationaryDistribution = new double[stateCount];
        storedStationaryDistribution = new double[stateCount];

        addParameter(infinitesimalRates);

        illConditionedProbabilities = new double[stateCount * stateCount];
        for (int i = 0; i < stateCount * stateCount; i++)
            illConditionedProbabilities[i] = 1.0 / stateCount;

    }

    protected void handleModelChangedEvent(Model model, Object object, int index) {
        if (model == freqModel)
            return; // freqModel only affects the likelihood calculation at the tree root
        super.handleModelChangedEvent(model, object, index);
    }

    protected void restoreState() {

        // To restore all this stuff just swap the pointers...
//		DoubleMatrix2D tmp = eigenD;
//		eigenD = storedEigenD;
//		storedEigenD = tmp;
//
//		tmp = eigenV;
//		eigenV = storedEigenV;
//		storedEigenV = tmp;

        double[] tmp3 = storedEvalImag;
        storedEvalImag = EvalImag;
        EvalImag = tmp3;

        tmp3 = storedStationaryDistribution;
        storedStationaryDistribution = stationaryDistribution;
        stationaryDistribution = tmp3;

        normalization = storedNormalization;

        // Inherited
        updateMatrix = storedUpdateMatrix;
        wellConditioned = storedWellConditioned;

        double[] tmp1 = storedEval;
        storedEval = Eval;
        Eval = tmp1;

        double[][] tmp2 = storedIevc;
        storedIevc = Ievc;
        Ievc = tmp2;

        tmp2 = storedEvec;
        storedEvec = Evec;
        Evec = tmp2;

    }

    protected void storeState() {

        storedUpdateMatrix = updateMatrix;
        storedWellConditioned = wellConditioned;

        //         storedEigenDecomp = eigenDecomp.


//		if( storedEigenD == null )
//				storedEigenD = new DenseDoubleMatrix2D(eigenD.rows(),eigenD.columns());
//		if( storedEigenV == null )
//				storedEigenV = new DenseDoubleMatrix2D(eigenV.rows(),eigenV.columns());
//		storedEigenD.assign(eigenD);
//		storedEigenV.assign(eigenV);
        // todo assign values instead of allocating new object

//		storedEigenD = eigenD.copy();
//		storedEigenV = eigenV.copy();
//		storedEigenVReal = eigenVReal.copy();
//		storedEigenVImag = eigenVImag.copy();
//

        System.arraycopy(stationaryDistribution, 0, storedStationaryDistribution, 0, stateCount);
        System.arraycopy(EvalImag, 0, storedEvalImag, 0, stateCount);
        storedNormalization = normalization;

        // Inherited
        System.arraycopy(Eval, 0, storedEval, 0, stateCount);
        for (int i = 0; i < stateCount; i++) {
            System.arraycopy(Ievc[i], 0, storedIevc[i], 0, stateCount);
            System.arraycopy(Evec[i], 0, storedEvec[i], 0, stateCount);
        }

    }

    public void getTransitionProbabilities(double distance, double[] matrix) {

        double temp;

        int i, j, k;

        synchronized (this) {
            if (updateMatrix) {
                setupMatrix();
            }
        }

        if (!wellConditioned) {
            System.arraycopy(illConditionedProbabilities, 0, matrix, 0, stateCount * stateCount);
            return;
        }

// Eigenvalues and eigenvectors of a real matrix A.
//
// If A is symmetric, then A = V*D*V' where the eigenvalue matrix D is diagonal
// and the eigenvector matrix V is orthogonal. I.e. A = V D V^t and V V^t equals
// the identity matrix.
//
// If A is not symmetric, then the eigenvalue matrix D is block diagonal with
// the real eigenvalues in 1-by-1 blocks and any complex eigenvalues,
// lambda + i*mu, in 2-by-2 blocks, [lambda, mu; -mu, lambda]. The columns
// of V represent the eigenvectors in the sense that A*V = V*D. The matrix
// V may be badly conditioned, or even singular, so the validity of the
// equation A = V D V^{-1} depends on the conditioning of V.

        double[][] iexp = popiexp();

        for (i = 0; i < stateCount; i++) {

            if (EvalImag[i] == 0) {
                // 1x1 block
                temp = Math.exp(distance * Eval[i]);
                for (j = 0; j < stateCount; j++) {
                    iexp[i][j] = Ievc[i][j] * temp;
                }
            } else {
                // 2x2 conjugate block
                // If A is 2x2 with complex conjugate pair eigenvalues a +/- bi, then
                // exp(At) = exp(at)*( cos(bt)I + \frac{sin(bt)}{b}(A - aI)).
                int i2 = i + 1;
                double b = EvalImag[i];
                double expat = Math.exp(distance * Eval[i]);
                double expatcosbt = expat * Math.cos(distance * b);
                double expatsinbt = expat * Math.sin(distance * b);

                for (j = 0; j < stateCount; j++) {
//                    try {
                        iexp[i][j] = expatcosbt * Ievc[i][j] + expatsinbt * Ievc[i2][j];
                        iexp[i2][j] = expatcosbt * Ievc[i2][j] - expatsinbt * Ievc[i][j];
//                    } catch (Exception e) {
//                        System.err.println("Exception: " + e.getMessage());
//                        DoubleMatrix1D eigenVReal = eigenDecomp.getRealEigenvalues();
//                        DoubleMatrix1D eigenVImag = eigenDecomp.getImagEigenvalues();
//                        System.err.println("Real = " + eigenVReal);
//                        System.err.println(new Vector(Eval));
//                        System.err.println("Imag = " + eigenVImag);
//                        System.err.println(new Vector(EvalImag));
//                        System.exit(-1);
//                    }
                }
                i++; // processed two conjugate rows
            }
        }

        int u = 0;
        for (i = 0; i < stateCount; i++) {
            for (j = 0; j < stateCount; j++) {
                temp = 0.0;
                for (k = 0; k < stateCount; k++) {
                    temp += Evec[i][k] * iexp[k][j];
                }
                if (temp < 0.0)
                    matrix[u] = minProb;
                else
                    matrix[u] = temp;
                u++;
            }
        }
        pushiexp(iexp);
    }

    public double[] getStationaryDistribution() {
        return stationaryDistribution;
    }

    protected void computeStationaryDistribution() {
        stationaryDistribution = freqModel.getFrequencies();       
    }


    protected double[] getRates() {
        return infinitesimalRates.getParameterValues();
    }

    protected void setupMatrix() {

        if (!eigenInitialised) {
            initialiseEigen();
            storedEvalImag = new double[stateCount];
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


        //

        DoubleMatrix2D eigenVInv = null;

        updateMatrix = false;

        try {
            eigenVInv = alegbra.inverse(eigenV);
            wellConditioned = true;
        } catch (IllegalArgumentException e) {
            wellConditioned = false;
//            updateMatrix = false;
//            count++;
//            if (count > 10000) {
//                System.err.println("Too many ill-conditioned matrices");
//                 System.exit(-1);
//            }

//            Eval[0] = 42;
//            System.err.println("bad");
            return;
//            throw e;
        }

//              System.err.println("good");
        // fill AbstractSubstitutionModel parameters

        Ievc = eigenVInv.toArray();
        Evec = eigenV.toArray();
        Eval = eigenVReal.toArray();
        EvalImag = eigenVImag.toArray();

        checkComplexSolutions();

        // Check for valid decomposition
        for (i = 0; i < stateCount; i++) {
            if (Eval[i] == Double.NaN || EvalImag[i] == Double.NaN) {
                System.err.println("Got here");
                wellConditioned = false;
                return;
            }
        }

        // compute normalization and rescale eigenvalues

        computeStationaryDistribution();

        double subst = 0.0;

        for (i = 0; i < stateCount; i++)
            subst += -amat[i][i] * stationaryDistribution[i];

        normalization = subst;

//        if (subst == Double.NaN || subst == Double.POSITIVE_INFINITY || subst == Double.NEGATIVE_INFINITY) {
        if (subst == 0 ) {
            System.err.println("subst = "+subst);
            System.exit(-1);
        }

        for (i = 0; i < stateCount; i++) {
            Eval[i] /= subst;
            EvalImag[i] /= subst;
        }

//		printDebugSetupMatrix();
//        wellConditioned = true;
//        updateMatrix = false;
    }

    private void printDebugSetupMatrix() {
        System.out.println("Normalized infinitesimal rate matrix:");
        System.out.println(new Matrix(amat));
        System.out.println(new Matrix(amat).toStringOctave());
        System.out.println("Normalization = " + normalization);
        System.out.println("Values in setupMatrix():");
//		System.out.println(eigenV);
//		System.out.println(eigenVInv);
//		System.out.println(eigenVReal);
    }

    protected void checkComplexSolutions() {
        boolean complex = false;
        for (int i = 0; i < stateCount && !complex; i++) {
            if (EvalImag[i] != 0)
                complex = true;
        }
        isComplex = complex;
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

        Parameter rates = new Parameter.Default(new double[]{5.0, 1.0, 1.0, 0.1, 5.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0});
//		Parameter rates = new Parameter.Default(new double[] {5.0, 1.0, 1.0, 1.0, 5.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0});
//		Parameter rates = new Parameter.Default(new double[] {1.0, 1.0});

        ComplexSubstitutionModel substModel = new ComplexSubstitutionModel("test",
//				TwoStates.INSTANCE,
                Nucleotides.INSTANCE,
                null,
                rates);

        double[] finiteTimeProbs = new double[substModel.getDataType().getStateCount() * substModel.getDataType().getStateCount()];
        double time = 1.0;
        substModel.getTransitionProbabilities(time, finiteTimeProbs);

        System.out.println("Results:");
        System.out.println(new Vector(finiteTimeProbs));

//		System.out.println("COLT value:");
//		 This should work, matches 'octave' results
//		DoubleMatrix2D result = alegbra.mult(substModel.eigenV, alegbra.mult(blockDiagonalExponential(1.0, substModel.eigenD), substModel.eigenVInv));
//
//		System.out.println(result);

    }

    private static DoubleMatrix2D blockDiagonalExponential(double distance, DoubleMatrix2D mat) {
        for (int i = 0; i < mat.rows(); i++) {
            if ((i + 1) < mat.rows() && mat.getQuick(i, i + 1) != 0) {
                double a = mat.getQuick(i, i);
                double b = mat.getQuick(i, i + 1);
                double expat = Math.exp(distance * a);
                double cosbt = Math.cos(distance * b);
                double sinbt = Math.sin(distance * b);
                mat.setQuick(i, i, expat * cosbt);
                mat.setQuick(i + 1, i + 1, expat * cosbt);
                mat.setQuick(i, i + 1, expat * sinbt);
                mat.setQuick(i + 1, i, -expat * sinbt);
                i++; // processed two entries in loop
            } else
                mat.setQuick(i, i, Math.exp(distance * mat.getQuick(i, i))); // 1x1 block
        }
        return mat;
    }

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return COMPLEX_SUBSTITUTION_MODEL;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            DataType dataType = DataTypeUtils.getDataType(xo);

            if (dataType == null) dataType = (DataType) xo.getChild(DataType.class);

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
                return new ComplexSubstitutionModel(xo.getId(), dataType, rootFreq, ratesParameter);
            else
                return new SVSComplexSubstitutionModel(xo.getId(), dataType, rootFreq, ratesParameter, indicators);

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

//	private DoubleMatrix2D eigenD;
    //	private DoubleMatrix2D storedEigenD;
//	private DoubleMatrix2D eigenV;
//	private DoubleMatrix2D eigenVInv;
//	private DoubleMatrix2D storedEigenV;
    //	private DoubleMatrix2D storedEigenVInv;
//	private DoubleMatrix1D eigenVReal;
//	private DoubleMatrix1D eigenVImag;
//	private DoubleMatrix1D storedEigenVReal;
//	private DoubleMatrix1D storedEigenVImag;

    private boolean isComplex = false;
    private double[] stationaryDistribution = null;
    private double[] storedStationaryDistribution;
    private Double normalization;
    private Double storedNormalization;

    protected double[] EvalImag;
    protected double[] storedEvalImag;

//	private boolean normalizationAmat = false;

    private boolean wellConditioned = true;
    private boolean storedWellConditioned;
    private double[] illConditionedProbabilities;

    private static final double minProb = Property.DEFAULT.tolerance();
//    private static final double minProb = 1E-20;
    //    private static final double minProb = Property.ZERO.tolerance();
    private static final Algebra alegbra = new Algebra(minProb);
    EigenvalueDecomposition eigenDecomp;
    EigenvalueDecomposition storedEigenDecomp;

    private int count = 0;
}