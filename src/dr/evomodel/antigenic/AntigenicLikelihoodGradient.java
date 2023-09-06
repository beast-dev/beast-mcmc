package dr.evomodel.antigenic;

import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.model.*;
import dr.inference.multidimensionalscaling.MultiDimensionalScalingCore;
import dr.xml.Reportable;

import java.util.List;

/**
 * @author Marc A. Suchard
 * @author Philippe Lemey
 */

public class AntigenicLikelihoodGradient extends AbstractModel implements GradientWrtParameterProvider, Reportable {

    private final NewAntigenicLikelihood likelihood;
    private final MultiDimensionalScalingCore mdsCore;
    private final List<AntigenicGradientWrtParameter> wrtList;

    private final int numViruses;
    private final int numSera;
    private final int mdsDim;

    private final Parameter parameter;

    private boolean locationGradientKnown;
    private boolean savedLocationGradientKnown;

    private double[] locationGradient;
    private double[] savedLocationGradient;

    private boolean observationGradientKnown;
    private boolean savedObservationGradientKnown;

    private double[] observationGradient;
    private double[] savedObservationGradient;

    public AntigenicLikelihoodGradient(NewAntigenicLikelihood likelihood, List<AntigenicGradientWrtParameter> wrtList) {
        super("AntigenicLikelihoodGradient");
        this.likelihood = likelihood;
        this.mdsCore = likelihood.getCore();
        this.wrtList = wrtList;

        this.numViruses = likelihood.getNumberOfViruses();
        this.numSera = likelihood.getNumberOfSera();
        this.mdsDim = likelihood.getMdsDimension();

        addModel(likelihood);

        if (wrtList.size() == 1) {
            this.parameter = wrtList.get(0).getParameter();
        } else {
            CompoundParameter cp = new CompoundParameter("AntigenicLikelihoodGradient");
            for (AntigenicGradientWrtParameter wrt : wrtList) {
                cp.addParameter(wrt.getParameter());
            }
            this.parameter = cp;
        }

        locationGradientKnown = false;
    }

    @Override
    public Parameter getParameter() {
        return parameter;
    }

    @Override
    public int getDimension() {
        return parameter.getDimension();
    }

    @Override
    public double[] getGradientLogDensity() {

        if (!locationGradientKnown && requiresLocationGradient()) {
            getLocationGradients();
            locationGradientKnown = true;
        }

        if (!observationGradientKnown && requiresObservationGradient()) {
            getObservationGradients();
            observationGradientKnown = true;
        }

        double[] gradient = new double[getGradientSize()];

        int offset = 0;
        for (AntigenicGradientWrtParameter wrt : wrtList) {
            wrt.getGradient(gradient, offset, locationGradient, observationGradient);
            offset += wrt.getSize();
        }

        return gradient;
    }

    @Override
    public String getReport() {
        return GradientWrtParameterProvider.getReportAndCheckForError(this,
                Double.NEGATIVE_INFINITY,
                Double.POSITIVE_INFINITY, tolerance);
    }

    private final static double tolerance = 1E-3;

    private boolean requiresLocationGradient() {
        for (AntigenicGradientWrtParameter wrt : wrtList) {
            if (wrt.requiresLocationGradient()) {
                return true;
            }
        }
        return false;
    }

    private boolean requiresObservationGradient() {
        for (AntigenicGradientWrtParameter wrt : wrtList) {
            if (wrt.requiresObservationGradient()) {
                return true;
            }
        }
        return false;
    }

    private int getGradientSize() {
        int size = 0;
        for (AntigenicGradientWrtParameter wrt : wrtList) {
            size += wrt.getSize();
        }
        return size;
    }

    private void getLocationGradients() {
        if (locationGradient == null) {
            locationGradient = new double[(numViruses + numSera) * mdsDim];
        }

        mdsCore.getLocationGradient(locationGradient);
    }

    private void getObservationGradients() {
        if (observationGradient == null) {
            observationGradient = new double[numViruses * numSera];
        }

        mdsCore.getObservationGradient(observationGradient);
    }

    @Override
    protected void handleModelChangedEvent(Model model, Object object, int index) {
        if (model == likelihood) {
            locationGradientKnown = false;
            observationGradientKnown = false;
        } else {
            throw new IllegalArgumentException("Unknown model");
        }
    }

    @Override
    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        // Do nothing
    }

    @Override
    protected void storeState() {
        savedLocationGradientKnown = locationGradientKnown;
        savedObservationGradientKnown = observationGradientKnown;

        if (locationGradientKnown) {
            if (savedLocationGradient == null) {
                savedLocationGradient = new double[locationGradient.length];
            }
            System.arraycopy(locationGradient, 0, savedLocationGradient, 0, locationGradient.length);
        }

        if (observationGradientKnown) {
            if (savedObservationGradient == null) {
                savedObservationGradient = new double[observationGradient.length];
            }
            System.arraycopy(observationGradient, 0, savedObservationGradient, 0, observationGradient.length);
        }
    }

    @Override
    protected void restoreState() {
        locationGradientKnown = savedLocationGradientKnown;
        observationGradientKnown = savedObservationGradientKnown;

        if (locationGradientKnown) {
            double[] swap = locationGradient;
            locationGradient = savedLocationGradient;
            savedLocationGradient = swap;
        }

        if (observationGradientKnown) {
            double[] swap = observationGradient;
            observationGradient = savedObservationGradient;
            savedObservationGradient = swap;
        }
    }

    @Override
    protected void acceptState() {
        // Do nothing
    }

    @Override
    public Likelihood getLikelihood() {
        return likelihood;
    }
}
