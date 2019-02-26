package dr.util;

import dr.inference.model.Parameter;
import dr.inference.loggers.LogColumn;
import dr.inference.model.Statistic;
import dr.inference.model.Variable;
import dr.inference.model.VariableListener;
import dr.xml.ElementRule;
import dr.xml.XMLObject;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;

public interface CrossValidationProvider {

    Parameter getTrueParameter();

    Parameter getInferredParameter();

    int[] getRelevantDimensions();

    String getName(int dim);





    public class CrossValidator extends Statistic.Abstract implements VariableListener {
        final CrossValidationProvider provider;
        final double[] squaredErrors;
        final int[] relevantDims;
        final Parameter truthParameter;
        final Parameter inferredParameter;
        final int dimStat;
        boolean statKnown = false;

        CrossValidator(CrossValidationProvider provider) {
            this.provider = provider;
            this.relevantDims = provider.getRelevantDimensions();

            this.dimStat = relevantDims.length;
            this.squaredErrors = new double[dimStat];
            this.truthParameter = provider.getTrueParameter();
            this.inferredParameter = provider.getInferredParameter();

            inferredParameter.addParameterListener(this);


        }

        private void updateSquaredErrors() {
            for (int i = 0; i < dimStat; i++) {
                double truth = truthParameter.getParameterValue(relevantDims[i]);
                double inferred = inferredParameter.getParameterValue(relevantDims[i]);
                double error = truth - inferred;
                squaredErrors[i] = error * error;
            }
        }


        @Override
        public String getDimensionName(int dim) {
            return provider.getName(dim);
        }

        @Override
        public int getDimension() {
            return dimStat;
        }


        @Override
        public double getStatisticValue(int dim) {

            if (!statKnown) {
                updateSquaredErrors();
            }

            return squaredErrors[dim];
        }


        @Override
        public void variableChangedEvent(Variable variable, int index, Variable.ChangeType type) {
            statKnown = false;
        }
    }
}
