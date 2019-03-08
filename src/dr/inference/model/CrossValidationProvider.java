package dr.inference.model;

public interface CrossValidationProvider {

    Parameter getTrueParameter();

    Parameter getInferredParameter();

    int[] getRelevantDimensions();

    String getName(int dim);

    String getNameSum(int dim);


    class CrossValidator extends Statistic.Abstract {
        protected final CrossValidationProvider provider;
        private final double[] squaredErrors;
        private final int[] relevantDims;
        private Parameter truthParameter;
        private Parameter inferredParameter;
        private final int dimStat;
//        boolean statKnown = false;

        CrossValidator(CrossValidationProvider provider) {
            this.provider = provider;
            this.relevantDims = provider.getRelevantDimensions();

            this.dimStat = relevantDims.length;
            this.squaredErrors = new double[dimStat];
//            this.truthParameter = provider.getTrueParameter();
//            this.inferredParameter = provider.getInferredParameter();


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

            //TODO: add variable listeners as needed
            if (dim == 0) {
                this.inferredParameter = provider.getInferredParameter();
                this.truthParameter = provider.getTrueParameter();
                updateSquaredErrors();
            }

            return squaredErrors[dim];
        }


    }

    class CrossValidatorSum extends CrossValidator {

        CrossValidatorSum(CrossValidationProvider provider) {
            super(provider);
        }

        @Override
        public String getDimensionName(int dim) {
            return provider.getNameSum(dim);
        }

        @Override
        public int getDimension() {
            return 1;
        }


        @Override
        public double getStatisticValue(int dim) {
            double sum = 0;
            for (int i = 0; i < super.getDimension(); i++) {
                sum += super.getStatisticValue(i);
            }

            return sum;
        }


    }
}
