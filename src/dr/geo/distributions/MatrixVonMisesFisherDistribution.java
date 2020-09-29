package dr.geo.distributions;

import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treedatalikelihood.continuous.IntegratedFactorAnalysisLikelihood;
import dr.inference.operators.factorAnalysis.FactorAnalysisOperatorAdaptor;
import dr.math.MathUtils;
import dr.math.ModifiedBesselFirstKind;
import dr.math.distributions.MultivariateDistribution;
import dr.math.distributions.RandomGenerator;
import dr.math.matrixAlgebra.Vector;
import dr.xml.*;
import org.ejml.data.DenseMatrix64F;
import org.ejml.factory.DecompositionFactory;
import org.ejml.interfaces.decomposition.EigenDecomposition;
import org.ejml.interfaces.decomposition.SingularValueDecomposition;
import org.ejml.ops.CommonOps;

import static dr.math.MathUtils.nextDouble;

public class MatrixVonMisesFisherDistribution implements RandomGenerator, MultivariateDistribution, Reportable {

    private final FactorAnalysisOperatorAdaptor adaptor; //TODO: way fewer buffers
    private final DenseMatrix64F C;
    private final int nRows;
    private final int nColumns;
    private final DenseMatrix64F mkBuffer1;
    private final DenseMatrix64F mkBuffer2;
    private final DenseMatrix64F kkBuffer1;
    private final DenseMatrix64F kkBuffer2;
    private final DenseMatrix64F kkBuffer3;
    private final DenseMatrix64F mmBuffer;
    private final DenseMatrix64F mBuffer1;
    private final DenseMatrix64F mBuffer2;


    private final DenseMatrix64F D;
    private final DenseMatrix64F V;
    private final DenseMatrix64F F;
    private final DenseMatrix64F Y;
    private final DenseMatrix64F H;
    private static final int MAX_REJECTS = 100;

    public MatrixVonMisesFisherDistribution(FactorAnalysisOperatorAdaptor adaptor) { //TODO:remove adaptor?

        this.adaptor = adaptor;
        this.nRows = adaptor.getNumberOfTraits();
        this.nColumns = adaptor.getNumberOfFactors();
        this.C = new DenseMatrix64F(nRows, nColumns);

        this.mkBuffer1 = new DenseMatrix64F(nRows, nColumns);
        this.mkBuffer2 = new DenseMatrix64F(nRows, nColumns);
        this.kkBuffer1 = new DenseMatrix64F(nColumns, nColumns);
        this.kkBuffer2 = new DenseMatrix64F(nColumns, nColumns);
        this.kkBuffer3 = new DenseMatrix64F(nColumns, nColumns);
        this.mmBuffer = new DenseMatrix64F(nRows, nRows);
        this.mBuffer1 = new DenseMatrix64F(nRows, 1);
        this.mBuffer2 = new DenseMatrix64F(nRows, 1);

        this.V = new DenseMatrix64F(adaptor.getNumberOfTraits(), adaptor.getNumberOfFactors());
        this.D = new DenseMatrix64F(adaptor.getNumberOfFactors(), adaptor.getNumberOfFactors());
        this.F = new DenseMatrix64F(adaptor.getNumberOfTaxa(), adaptor.getNumberOfFactors());
        this.Y = new DenseMatrix64F(adaptor.getNumberOfTaxa(), adaptor.getNumberOfTraits());
        this.H = new DenseMatrix64F(nRows, nColumns);

    }


    @Override
    public double[] nextRandom() {
        updateC();
        return nextRandomNoUpdate();
    }

    public double[] nextRandomNoUpdate() {
        SingularValueDecomposition svd = DecompositionFactory.svd(C.numRows, C.numCols, true, true, true);
        svd.decompose(C);
        double[] singularValues = svd.getSingularValues();
        svd.getU(H, false);
        for (int i = 0; i < nRows; i++) {
            for (int j = 0; j < nColumns; j++) {
                H.set(i, j, H.get(i, j) * singularValues[j]);
            }
        }

        int rejects = 0;
        boolean valid = false;


        double ratio = 1.0;

        while (!valid && rejects < MAX_REJECTS) {
            double u = MathUtils.nextDouble();

            DenseMatrix64F x1 = transferColumns(H, 0, 1);
            double[] draw1 = nextVectorVonMisesFisher(x1.data);
            x1.setData(draw1);
            transferColumns(x1, 0, mkBuffer1, 0, 1); //mkBuffer1 stores matrix draw

            ratio = 1.0;
            double besselOrder = (C.numRows - C.numCols - 1) / 2;

            for (int i = 1; i < C.numCols; i++) {
                DenseMatrix64F previousDraws = new DenseMatrix64F(C.numRows, i);
                transferColumns(mkBuffer1, 0, previousDraws, 0, i);
                SingularValueDecomposition nullSvd = DecompositionFactory.svd(C.numRows, i, true, false, false);
                nullSvd.decompose(previousDraws);
                nullSvd.getU(mmBuffer, false);
                DenseMatrix64F nullMatrix = transferColumns(mmBuffer, i, C.numRows - i);
                transferColumns(H, i, mBuffer1, 0, 1);
                double normHr = computeNorm(mBuffer1.data);

                DenseMatrix64F z = new DenseMatrix64F(C.numRows - i, 1);
                CommonOps.multTransA(nullMatrix, mBuffer1, z);
                double normNtHr = computeNorm(z.data);
                double[] zBuffer = nextVectorVonMisesFisher(z.data);
                z.setData(zBuffer);
                CommonOps.mult(nullMatrix, z, mBuffer1);
                transferColumns(mBuffer1, 0, mkBuffer1, i, 1);

                double iHr = ModifiedBesselFirstKind.bessi(normHr, besselOrder);
                double iNtHr = ModifiedBesselFirstKind.bessi(normNtHr, besselOrder);
                ratio *= (iHr / iNtHr) * Math.pow(normHr / normNtHr, besselOrder);
            }

            if (u < ratio) {
                valid = true;
            }

            System.out.println(rejects);
            rejects++;

        }

        if (!valid) {
            System.err.println("Didn't work.");
            return null;
        }

        svd.getV(kkBuffer1, true);

        DenseMatrix64F X = new DenseMatrix64F(C.numRows, C.numCols);
        CommonOps.multTransB(mkBuffer1, kkBuffer1, X);
        return X.data;
    }

    public void transferColumns(DenseMatrix64F src, int srcStart, DenseMatrix64F dest, int destStart, int nColumns) {
        for (int i = 0; i < src.numRows; i++) {
            for (int j = 0; j < nColumns; j++) {
                dest.set(i, j + destStart, src.get(i, j + srcStart));
            }
        }
    }

    public DenseMatrix64F transferColumns(DenseMatrix64F src, int srcStart, int nColumns) {
        DenseMatrix64F mat = new DenseMatrix64F(src.numRows, nColumns);
        transferColumns(src, srcStart, mat, 0, nColumns);
        return mat;
    }


    private double[] nextVectorVonMisesFisher(double[] c) {
        double[] u = c.clone();
        double norm = makeUnit(u);
        double[] unitDraw = nextVectorVonMisesFisherUnitMode(c.length, norm);
        SingularValueDecomposition svd = DecompositionFactory.svd(c.length, 1, true, false, false);
        DenseMatrix64F mBuffer = new DenseMatrix64F(c.length, 1);
        DenseMatrix64F mmBuffer = new DenseMatrix64F(c.length, c.length);
        DenseMatrix64F mBuffer2 = new DenseMatrix64F(c.length, 1);

        mBuffer.setData(u);
        svd.decompose(mBuffer); //TODO: probably better way to construct orthogonal matrix including u
        svd.getU(mmBuffer, false);

        for (int i = 0; i < u.length; i++) {
            double f = mmBuffer.get(i, 0);
            double l = mmBuffer.get(i, u.length - 1);
            mmBuffer.set(i, 0, l);
            mmBuffer.set(i, u.length - 1, f);

        }

        mBuffer.setData(unitDraw);
        CommonOps.mult(mmBuffer, mBuffer, mBuffer2);

        return mBuffer2.getData();
    }

    //Wood, Andrew TA. "Simulation of the von Mises Fisher distribution." Communications in statistics-simulation and computation 23.1 (1994): 157-164.
    //draws from von Mises-Fisher with mode (0, ..., 0, 1)^t and concentration `norm`
    private double[] nextVectorVonMisesFisherUnitMode(int m, double norm) {

        int rejects = 0;

        while (rejects < MAX_REJECTS) {

            double mMinusOne = m - 1;
            double b = -2 * norm + Math.sqrt(4 * norm * norm + mMinusOne * mMinusOne);
            b /= mMinusOne;

            double x0 = (1 - b) / (1 + b);
            double c = norm * x0 + mMinusOne * Math.log(1 - x0 * x0);

            double z = MathUtils.nextBeta(mMinusOne / 2, mMinusOne / 2);
            double u = MathUtils.nextDouble();

            double w = (1 - (1 + b) * z) / (1 - (1 - b) * z);

            if (norm * w + mMinusOne * Math.log(1 - x0 * w) - c > Math.log(u)) {
                double[] v = nextUniformVector(m - 1);
                double[] draw = new double[m];
                double norm2 = Math.sqrt(1 - w * w);
                for (int i = 0; i < v.length; i++) {
                    draw[i] = norm2 * v[i];
                }
                draw[m - 1] = w;
                return draw;
            }

            rejects++;
        }

        return null;
    }

    private double[] nextUniformVector(int dim) {
        double[] draw = new double[dim];
        for (int i = 0; i < dim; i++) {
            draw[i] = MathUtils.nextGaussian();
        }
        makeUnit(draw);
        return draw;
    }

    private double makeUnit(double[] x) {
        double norm = computeNorm(x);
        double invNorm = 1.0 / norm;
        for (int i = 0; i < x.length; i++) {
            x[i] *= invNorm;
        }
        return norm;
    }

    private double computeNorm(double[] x) {
        double sumSqares = 0;
        for (int i = 0; i < x.length; i++) {
            sumSqares += x[i] * x[i];
        }
        return Math.sqrt(sumSqares);
    }

    private double[] slowNextRandom() {
        updateC();
        int rejects = 0;

        while (rejects < MAX_REJECTS) {
            DenseMatrix64F uniformDraw = nextUniform();

            SingularValueDecomposition svd = new DecompositionFactory().svd(C.numRows, C.numCols, false, false, true);
            svd.decompose(C);
            double[] singularValues = svd.getSingularValues();
            System.out.println("Rejects: " + rejects);
            CommonOps.multTransA(C, uniformDraw, kkBuffer1); //TODO: just need the trace, super inefficient
            double trace = 0;
            for (int i = 0; i < C.numCols; i++) {
                trace += kkBuffer1.get(i, i) - singularValues[i];
            }

            System.out.println("ExpTrace: " + Math.exp(trace));

            if (nextDouble() < Math.exp(trace)) {
                return uniformDraw.getData();
            }
            rejects++;
            System.out.println("");
        }

        throw new RuntimeException("Rejection sampler failed.");//TODO: handle better
    }

    private DenseMatrix64F nextUniform() {
        double[] X = new double[C.getNumElements()];
        for (int i = 0; i < X.length; i++) {
            X[i] = MathUtils.nextGaussian();
        }
        mkBuffer1.setData(X);
        CommonOps.multTransA(mkBuffer1, mkBuffer1, kkBuffer1);
        EigenDecomposition<DenseMatrix64F> eig = DecompositionFactory.eig(nColumns, true, true);
        eig.decompose(kkBuffer1);
        for (int i = 0; i < nColumns; i++) {
            DenseMatrix64F vector = eig.getEigenVector(i);
            for (int j = 0; j < nColumns; j++) {
                kkBuffer2.set(j, i, vector.get(j, 0));
            }
            double value = eig.getEigenvalue(i).getReal(); //TODO: should be real
            kkBuffer3.set(i, i, 1.0 / Math.sqrt(value));
        }

        CommonOps.mult(kkBuffer2, kkBuffer3, kkBuffer1); // kkBuffer1 = VD
        CommonOps.multTransB(kkBuffer1, kkBuffer2, kkBuffer3); // kkBuffer3 = VDVt
        CommonOps.mult(mkBuffer1, kkBuffer3, mkBuffer2);
        return mkBuffer2;
    }


    @Override
    public double logPdf(Object x) {
        throw new RuntimeException("Not yet implemented.");
    }

    public void setC(double[] values) {
        C.setData(values);
    }

    public void updateC() {
        splitLoadings();
        fillFactors();
        fillTraits();
        double maxPrecision = getMaximumPrecision();


        CommonOps.multTransA(Y, F, C);
        for (int i = 0; i < adaptor.getNumberOfFactors(); i++) {
            double scaledNorm = D.get(i, i) * maxPrecision;
            for (int j = 0; j < adaptor.getNumberOfTraits(); j++) {
                C.set(j, i, C.get(j, i) * scaledNorm);
            }
        }

    }

    private void splitLoadings() {
        int offset1 = 0;
        int offset2 = 0;
        for (int i = 0; i < adaptor.getNumberOfFactors(); i++) {
            double sumSquares = 0;
            for (int j = 0; j < adaptor.getNumberOfTraits(); j++) {
                sumSquares += adaptor.getLoadingsValue(offset1);
                offset1++;
            }
            double norm = Math.sqrt(sumSquares);
            D.set(i, i, norm);
            double invNorm = 1.0 / norm;
            for (int j = 0; j < adaptor.getNumberOfTraits(); j++) {
                V.set(offset2, adaptor.getLoadingsValue(offset2) * invNorm);
                offset2++;
            }
        }
    }

    private void fillFactors() {
        adaptor.drawFactors();
        for (int i = 0; i < adaptor.getNumberOfTaxa(); i++) {
            for (int j = 0; j < adaptor.getNumberOfFactors(); j++) {
                F.set(i, j, adaptor.getFactorValue(j, i));
            }
        }
    }

    private void fillTraits() {
        //TODO: fill in missing traits
        for (int i = 0; i < adaptor.getNumberOfTaxa(); i++) {
            for (int j = 0; j < adaptor.getNumberOfTraits(); j++) {
                Y.set(i, j, adaptor.getDataValue(j, i));
            }
        }
    }

    private double getMaximumPrecision() {
        double maxPrec = 0;
        for (int i = 0; i < adaptor.getNumberOfTraits(); i++) {
            if (adaptor.getColumnPrecision(i) > maxPrec) {
                maxPrec = adaptor.getColumnPrecision(i);
            }
        }
        return maxPrec;
    }

    @Override
    public double logPdf(double[] x) {
        updateC();
        // TODO: normalizing constant (actually need it for changing other parameters via MH)
        mkBuffer1.setData(x);
        CommonOps.multTransA(C, mkBuffer1, kkBuffer1); //TODO: inefficient (only need trace)
        double trace = 0;
        for (int i = 0; i < kkBuffer1.numCols; i++) {
            trace += kkBuffer1.get(i, i);
        }
        return trace;
    }

    @Override
    public double[][] getScaleMatrix() {
        throw new RuntimeException("Not yet implemented.");
    }

    @Override
    public double[] getMean() {
        throw new RuntimeException("Not yet implemented.");
    }

    @Override
    public String getType() {
        return "MatrixVonMises-Fisher";
    }

    @Override
    public String getReport() {
        int repeats = 100;
        int dim = C.numRows;

        double[] c = nextUniformVector(dim);

        double norm = 10000;

        for (int i = 0; i < dim; i++) {
            c[i] *= norm;
        }

        double[] mean = new double[dim];
        double[] var = new double[dim];
        for (int i = 0; i < repeats; i++) {
            double[] draw = nextVectorVonMisesFisher(c);

            for (int j = 0; j < dim; j++) {
                mean[j] += draw[j];
                var[j] += draw[j] * draw[j];
            }
            double drawNorm = makeUnit(draw);
            if (!MathUtils.isClose(drawNorm, 1.0, 1e-8)) {
                System.err.println("Norm: " + drawNorm);
            }
        }

        for (int i = 0; i < dim; i++) {
            mean[i] /= repeats;
            var[i] /= repeats;
            var[i] -= mean[i] * mean[i];
        }

        StringBuilder sb = new StringBuilder("matrix von Mises-Fisher distribution:\n");
        makeUnit(c);
        sb.append("original: " + new Vector(c) + "\n");
        sb.append("mean: " + new Vector(mean));
        sb.append("\n");
        sb.append("variance: " + new Vector(var));
        sb.append("\n\n");

        DenseMatrix64F newC = nextUniform();
        for (int i = 0; i < newC.data.length; i++) {
            C.data[i] = newC.data[i] * norm;
        }

        double[] matDraw = nextRandomNoUpdate();
        sb.append(new Vector(matDraw));
        sb.append("\n");
        sb.append(new Vector(newC.data));

        return sb.toString();


    }

    public static final AbstractXMLObjectParser PARSER = new AbstractXMLObjectParser() {
        private static final String MATRIX_VON_MISES_FISHER_DISTRIBUTION = "matrixVonMisesFisherDistribution";


        @Override
        public Object parseXMLObject(XMLObject xo) throws XMLParseException {
            TreeDataLikelihood treeLikelihood = (TreeDataLikelihood) xo.getChild(TreeDataLikelihood.class);
            IntegratedFactorAnalysisLikelihood factorLikelihood =
                    (IntegratedFactorAnalysisLikelihood) xo.getChild(IntegratedFactorAnalysisLikelihood.class);
            FactorAnalysisOperatorAdaptor.IntegratedFactors adaptor =
                    new FactorAnalysisOperatorAdaptor.IntegratedFactors(factorLikelihood, treeLikelihood);
            return new MatrixVonMisesFisherDistribution(adaptor);
        }

        @Override
        public XMLSyntaxRule[] getSyntaxRules() {
            return new XMLSyntaxRule[]{
                    new ElementRule(TreeDataLikelihood.class),
                    new ElementRule(IntegratedFactorAnalysisLikelihood.class)
            };
        }

        @Override
        public String getParserDescription() {
            return null;
        }


        @Override
        public Class getReturnType() {
            return MatrixVonMisesFisherDistribution.class;
        }

        @Override
        public String getParserName() {
            return MATRIX_VON_MISES_FISHER_DISTRIBUTION;
        }
    };
}
