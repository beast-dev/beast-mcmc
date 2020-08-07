package dr.inference.model;

import dr.xml.Reportable;

public interface CrossValidationProvider {

    double[] getTrueValues();

    double[] getInferredValues();

    int[] getRelevantDimensions();

    String getName(int dim);

    String getNameSum(int dim);


    class CrossValidator extends Statistic.Abstract implements Reportable {
        protected final CrossValidationProvider provider;
        private final double[] statValues;
        private final int[] relevantDims;
        private double[] truthValues;
        private double[] inferredValues;
        //        private Parameter truthParameter;
//        private Parameter inferredParameter;
        private final int dimStat;
        //        boolean statKnown = false;
        private final ValidationType validationType;

        public CrossValidator(CrossValidationProvider provider, ValidationType validationType) {
            this.provider = provider;
            this.relevantDims = provider.getRelevantDimensions();

            this.dimStat = relevantDims.length;
            this.statValues = new double[dimStat];
            this.validationType = validationType;
//            this.truthParameter = provider.getTrueParameter();
//            this.inferredParameter = provider.getInferredParameter();


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
                validationType.updateValues(this);
            }

            return statValues[dim];
        }

        @Override
        public String getReport() {
            int reps = 1000000;
            double[] vals = new double[dimStat];
            for (int i = 0; i < reps; i++) {
                for (int j = 0; j < dimStat; j++) {
                    vals[j] += getStatisticValue(j);
                }
            }

            StringBuilder sb = new StringBuilder();
            sb.append("Cross Validation Report:\n\n");
            for(int j = 0; j < dimStat; j++) {
                sb.append("\t" + getDimensionName(j) + ": " + vals[j] / reps + "\n");
            }
            sb.append("\n");

            return sb.toString();
        }
    }

    class CrossValidatorSum extends CrossValidator {

        public CrossValidatorSum(CrossValidationProvider provider, ValidationType validationType) {
            super(provider, validationType);
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


    public enum ValidationType {
        SQUARED_ERROR("squaredError") {
            @Override
            void updateValues(CrossValidator crossValidator) {
                BIAS.updateValues(crossValidator);
                for (int i = 0; i < crossValidator.dimStat; i++) {
                    double error = crossValidator.statValues[i];
                    crossValidator.statValues[i] = error * error;
                }
            }
        },

        BIAS("bias") {
            @Override
            void updateValues(CrossValidator crossValidator) {
                crossValidator.truthValues = crossValidator.provider.getTrueValues();
                crossValidator.inferredValues = crossValidator.provider.getInferredValues();


                for (int i = 0; i < crossValidator.dimStat; i++) {
                    double truth = crossValidator.truthValues[crossValidator.relevantDims[i]];
                    double inferred = crossValidator.inferredValues[crossValidator.relevantDims[i]];
                    double error = truth - inferred;
                    crossValidator.statValues[i] = error;
                }
            }
        },

        VALUE("value") {
            @Override
            void updateValues(CrossValidator crossValidator) {
                crossValidator.inferredValues = crossValidator.provider.getInferredValues();
                for (int i = 0; i < crossValidator.dimStat; i++) {
                    double inferred = crossValidator.inferredValues[crossValidator.relevantDims[i]];
                    crossValidator.statValues[i] = inferred;
                }
            }
        };

        private final String name;

        ValidationType(String name) {
            this.name = name;
        }

        abstract void updateValues(CrossValidator crossValidator);

        public String getName() {
            return this.name;
        }
    }
}
