package dr.evomodel.antigenic;

import dr.inference.model.Parameter;

public interface AntigenicGradientWrtParameter {

    boolean requiresLocationGradient();

    boolean requiresObservationGradient();

    int getSize();

    void getGradient(double[] gradient, int offset,
                     double[] locationGradient,
                     double[] observationGradient);

    Parameter getParameter();

    class VirusLocations extends Locations {

        VirusLocations(int viruses, int sera, int mdsDim, Parameter parameter, NewAntigenicLikelihood.Layout layout) {
            super(viruses, sera, mdsDim, parameter, layout);
        }

        @Override
        int getLocationOffset() {
            return layout.getVirusLocationOffset();
        }

        @Override
        public int getSize() {
            return viruses * mdsDim;
        }
    }

    class SerumLocations extends Locations {

        SerumLocations(int viruses, int sera, int mdsDim, Parameter parameter, NewAntigenicLikelihood.Layout layout) {
            super(viruses, sera, mdsDim, parameter, layout);
        }

        @Override
        int getLocationOffset() {
            return layout.getSerumLocationOffset();
        }

        @Override
        public int getSize() {
            return sera * mdsDim;
        }
    }

    abstract class Locations extends Base {

        final NewAntigenicLikelihood.Layout layout;

        Locations(int viruses, int sera, int mdsDim,
                  Parameter parameter,
                  NewAntigenicLikelihood.Layout layout) {
            super(viruses, sera, mdsDim, parameter);
            this.layout = layout;
        }

        abstract int getLocationOffset();

        abstract public int getSize();

        @Override
        public boolean requiresLocationGradient() {
            return true;
        }

        @Override
        public boolean requiresObservationGradient() {
            return false;
        }

        @Override
        public void getGradient(double[] gradient, int offset,
                                double[] locationGradient,
                                double[] observationGradient) {
            System.arraycopy(locationGradient, getLocationOffset(), gradient, offset, getSize());
        }
    }

    class Drift extends Locations {

        private final Parameter virusTime;
        private final Parameter serumTime;

        Drift(int viruses, int sera, int mdsDim,
              Parameter locationDrift, Parameter virusTime, Parameter serumTime,
              NewAntigenicLikelihood.Layout layout) {
            super(viruses, sera, mdsDim, locationDrift, layout);
            this.virusTime = virusTime;
            this.serumTime = serumTime;
        }

        @Override
        int getLocationOffset() {
            throw new RuntimeException("Should not be called");
        }

        @Override
        public int getSize() { return 1; }

        @Override
        public void getGradient(double[] gradient, int offset,
                                double[] locationGradient, double[] observationGradient) {

            double derivative = 0;

            int virusOffset = layout.getVirusLocationOffset();
            for (int i = 0; i < viruses; ++i) {
                derivative += locationGradient[virusOffset + i * mdsDim] * virusTime.getParameterValue(i);
            }

            int serumOffset = layout.getSerumLocationOffset();
            for (int i = 0; i < sera; ++i) {
                derivative += locationGradient[serumOffset + i * mdsDim] * serumTime.getParameterValue(i);
            }

            gradient[offset] = derivative;
        }
    }

    abstract class Base implements AntigenicGradientWrtParameter {

        Base(int viruses, int sera, int mdsDim, Parameter parameter) {
            this.viruses = viruses;
            this.sera = sera;
            this.mdsDim = mdsDim;
            this.parameter = parameter;
        }

        @Override
        public Parameter getParameter() {
            return parameter;
        }

        final int viruses;
        final int sera;
        final int mdsDim;
        final Parameter parameter;
    }
}
