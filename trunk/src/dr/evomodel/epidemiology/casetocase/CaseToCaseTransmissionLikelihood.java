package dr.evomodel.epidemiology.casetocase;

import dr.inference.model.*;

/**
 * A likelihood function for transmission between identified epidemiological cases
 *
 * Timescale must be in days. Python scripts to write XML for it and analyse the posterior set of networks exist;
 * contact MH. @todo make timescale not just in days
 *
 * @author Matthew Hall
 * @version $Id: $
 */

public class CaseToCaseTransmissionLikelihood extends AbstractModelLikelihood {

    private AbstractOutbreak outbreak;
    private CaseToCaseTreeLikelihood treeLikelihood;
    private AbstractKernelFunction spatialKernel;
    private Parameter transmissionRate;

    public CaseToCaseTransmissionLikelihood(AbstractOutbreak outbreak, // CaseToCaseTreeLikelihood treeLikelihood,
                                            AbstractKernelFunction spatialKernal, Parameter transmissionRate,
                                            Parameter infectionTimes, String name){
        super(name);
        this.outbreak = outbreak;
        this.treeLikelihood = treeLikelihood;
        this.spatialKernel = spatialKernal;
        this.transmissionRate = transmissionRate;
    }


    @Override
    protected void handleModelChangedEvent(Model model, Object object, int index) {
        //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    protected void storeState() {
        //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    protected void restoreState() {
        //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    protected void acceptState() {
        //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public Model getModel() {
        return null;  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public double getLogLikelihood() {
        return 0;  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public void makeDirty() {
        //To change body of implemented methods use File | Settings | File Templates.
    }
}
