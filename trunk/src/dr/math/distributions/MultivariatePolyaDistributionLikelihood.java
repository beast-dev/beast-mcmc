package dr.math.distributions;

import dr.inference.loggers.LogColumn;
import dr.inference.loggers.NumberColumn;
import dr.inference.model.*;
import dr.math.GammaFunction;
import dr.xml.*;

/**
 * Package: MultivariatePolyaDistributionLikelihood
 * Description:
 * this class provides a model for over-dispersed multinomial counts. The model follows Dirichlet-Multinomial distribution with
 * multinomial parameters integrated out analytically. This model is also known as Multivariate Polya distribution.
 * Standard parameterization involves k intensities a_i's. This implementation reparametrizes those as frequencies (k-1 df)
 * and dispersion parameters, where a = sum_i=1^k a_i is dispersion and f_i = a_i/a
 * <p/>
 * Created by
 *
 * @author Alexander V. Alekseyenko (alexander.alekseyenko@gmail.com)
 *         Date: Sep 23, 2010
 *         Time: 12:25:14 PM
 */
public class MultivariatePolyaDistributionLikelihood extends AbstractModel implements Likelihood {

    protected Parameter frequencies;
    protected Parameter dispersion;
    protected MatrixParameter data;
    protected double fixedNorm;
    protected double variableNorm;
    protected double storedFixedNorm;
    protected double storedVariableNorm;
    protected double logLikelihood;
    protected double storedLogLikelihood;
    protected boolean isLogLikelihoodKnown;
    protected boolean isFixedNormKnown;
    protected boolean isVariableNormKnown;
    protected double rowSums[];


    public MultivariatePolyaDistributionLikelihood(String modelID, MatrixParameter data, Parameter frequencies, Parameter dispersion) {
        super(modelID);
        this.frequencies = frequencies;
        this.dispersion = dispersion;
        this.data = data;
        isFixedNormKnown = false;
        isVariableNormKnown = false;
        addVariable(frequencies);
        addVariable(dispersion);
        addVariable(data);
        if (frequencies.getDimension() != data.getColumnDimension()) {
            System.err.println("Dimensions of the frequncy vector and number of columns do not match!");
        }
    }

    public MultivariatePolyaDistributionLikelihood(String modelID) {
        super(modelID);
    }

    public double calculateLogLikelihood() {
        // R code for this function:
        // //assuming X[,1] is row totals, lfactX is log factorial of X //
        // logLikes = sapply(1:n, function(i)  lfactX[subset,1][i] - sum(lfactX[subset,][i,2:p1]) +
        // lgamma(sum(alpha)) - lgamma(X[subset,1][i] + sum(alpha)) + sum(lgamma(X[subset,2:p1][i,] + alpha)) - sum(lgamma(alpha)));
        // sum(logLikes)

        if (!isFixedNormKnown) {
            computeFixedNorm();
        }
        if (!isVariableNormKnown) {
            computeVariableNorm();
        }
        double logP = fixedNorm + variableNorm;

        double freqs[] = frequencies.getParameterValues();
        double ad = dispersion.getParameterValue(0);

        for (int i = 0; i < data.getRowDimension(); ++i) {
            for (int j = 0; j < data.getColumnDimension(); ++j) {
                logP += GammaFunction.lnGamma(data.getParameterValue(i, j) + ad * freqs[j]);
            }
            logP -= GammaFunction.lnGamma(rowSums[i] + dispersion.getParameterValue(0));
        }
        return logP;
    }

    protected void computeFixedNorm() {
        rowSums = new double[data.getRowDimension()];
        for (int i = 0; i < data.getRowDimension(); ++i) {
            rowSums[i] = 0;
            for (int j = 0; j < data.getColumnDimension(); ++j) {
                rowSums[i] += data.getParameterValue(i, j);
            }
        }

        fixedNorm = 0;
        for (int i = 0; i < data.getRowDimension(); ++i) {
            for (int j = 0; j < data.getColumnDimension(); ++j) {
                fixedNorm -= GammaFunction.lnGamma(data.getParameterValue(i, j) + 1);
            }
            fixedNorm += GammaFunction.lnGamma(rowSums[i] + 1);
        }

        isFixedNormKnown = true;
    }

    protected void computeVariableNorm() {
        variableNorm = GammaFunction.lnGamma(dispersion.getParameterValue(0));

        double ad = dispersion.getParameterValue(0);
        for (int i = 0; i < frequencies.getDimension(); ++i) {
            variableNorm -= GammaFunction.lnGamma(frequencies.getParameterValue(i) * ad);
        }
        variableNorm *= data.getRowDimension();
    }

    protected void handleModelChangedEvent(Model model, Object object, int index) {
    }

    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        if (variable.getVariableName().equals(frequencies.getVariableName()) || variable.getVariableName().equals(dispersion.getVariableName())) {
            isVariableNormKnown = false;
        } else if (variable.getVariableName().equals(data.getVariableName())) {
            isFixedNormKnown = false;
        }
    }

    protected void storeState() {
        storedVariableNorm = variableNorm;
        storedFixedNorm = fixedNorm;
        storedLogLikelihood = logLikelihood;
    }

    protected void restoreState() {
        variableNorm = storedVariableNorm;
        fixedNorm = storedFixedNorm;
        logLikelihood = storedLogLikelihood;
    }

    protected void acceptState() {
    }

    public Model getModel() {
        return this;
    }

    public double getLogLikelihood() {
        if (!isLogLikelihoodKnown) {
            logLikelihood = calculateLogLikelihood();
        }
        return logLikelihood;
    }

    public void makeDirty() {
        isLogLikelihoodKnown = false;
        isVariableNormKnown = false;
        isFixedNormKnown = false;
    }

    public String prettyName() {
        return "Multivariate Polya Distribution Likelihood";
    }

    public void setUsed() {
    }


    public LogColumn[] getColumns() {
        return new LogColumn[]{
                new NumberColumn(this.getId()) {
                    public double getDoubleValue() {
                        return getLogLikelihood();
                    }
                }
        };
    }

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return MVPLIKE;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {
            MatrixParameter data;
            Parameter dispersion;
            Parameter frequencies;

            if (xo.hasChildNamed(DATA)) {
                data = (MatrixParameter) xo.getChild(DATA).getChild(MatrixParameter.class);
            } else {
                throw new XMLParseException("Missing data element!");
            }

            if (xo.hasChildNamed(DISPERSION)) {
                dispersion = (Parameter) xo.getChild(DISPERSION).getChild(Parameter.class);
            } else {
                throw new XMLParseException("Missing dispersion element!");
            }

            if (xo.hasChildNamed(FREQ)) {
                frequencies = (Parameter) xo.getChild(FREQ).getChild(Parameter.class);
            } else {
                throw new XMLParseException("Missing frequencies element!");
            }

            if (dispersion.getDimension() != 1) {
                throw new XMLParseException("Dispersion parameter must be of dimmension exactly 1!");
            }

            if (frequencies.getDimension() != data.getColumnDimension()) {
                throw new XMLParseException("The number of data columns must match the dimension" +
                        " of frequencies parameter (" +
                        data.getColumnDimension() + " != " + frequencies.getDimension() + "!");
            }

            return new MultivariatePolyaDistributionLikelihood(MVPLIKE, data, frequencies, dispersion);
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return "A matrix parameter constructed from its component parameters.";
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private final XMLSyntaxRule[] rules = {
                new ElementRule(DATA, new XMLSyntaxRule[]{new ElementRule(MatrixParameter.class)}, false),
                new ElementRule(DISPERSION, new XMLSyntaxRule[]{new ElementRule(Parameter.class)}, false),
                new ElementRule(FREQ, new XMLSyntaxRule[]{new ElementRule(Parameter.class)}, false),
        };

        public Class getReturnType() {
            return MatrixParameter.class;
        }
    };
    public static final String MVPLIKE = "mvPolyaLikelihood";
    public static final String DATA = "data";
    public static final String DISPERSION = "dispersion";
    public static final String FREQ = "frequencies";
}
