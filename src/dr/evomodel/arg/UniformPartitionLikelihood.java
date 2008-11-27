package dr.evomodel.arg;

import dr.evomodel.arg.operators.ARGAddRemoveEventOperator;
import dr.inference.loggers.LogColumn;
import dr.inference.model.*;
import dr.xml.*;

import java.util.logging.Logger;

public class UniformPartitionLikelihood extends AbstractModel implements Likelihood {

    public static final String UNIFORM_PARTITION_LIKELIHOOD = "uniformPartitionLikelihood";
    private double logPartitionNumber; //Transformed initially for computational reasons
    private ARGModel arg;

    public UniformPartitionLikelihood(ARGModel arg) {
        super("");
        this.arg = arg;

        addModel(arg);

        if (arg.isRecombinationPartitionType()) {
            logPartitionNumber = -Math.log(arg.getNumberOfPartitions() - 1);
        } else {
            logPartitionNumber = -(arg.getNumberOfPartitions() - 1) * ARGAddRemoveEventOperator.LOG_TWO;
        }
    }


    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserDescription() {
            return "Provides a uniform prior for partitions";
        }

        public Class getReturnType() {
            return UniformPartitionLikelihood.class;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return new XMLSyntaxRule[]{
                    new ElementRule(ARGModel.class),
            };
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {
            Logger.getLogger("dr.evomodel").info("Creating " + UNIFORM_PARTITION_LIKELIHOOD);

            return new UniformPartitionLikelihood((ARGModel) xo.getChild(ARGModel.class));
        }

        public String getParserName() {
            return UNIFORM_PARTITION_LIKELIHOOD;
        }

    };

    public double getLogLikelihood() {
        return logPartitionNumber * arg.getReassortmentNodeCount();
    }

    public Model getModel() {
        // TODO Auto-generated method stub
        return null;
    }

    public void makeDirty() {
        // TODO Auto-generated method stub

    }

    public LogColumn[] getColumns() {
        // TODO Auto-generated method stub
        return null;
    }

    public String getId() {
        // TODO Auto-generated method stub
        return null;
    }

    public void setId(String id) {
        // TODO Auto-generated method stub

    }

    @Override
    protected void acceptState() {
        // TODO Auto-generated method stub

    }

    @Override
    protected void handleModelChangedEvent(Model model, Object object, int index) {
        // TODO Auto-generated method stub

    }

    @Override
    protected void handleParameterChangedEvent(Parameter parameter, int index, ParameterChangeType type) {
        // TODO Auto-generated method stub

    }

    @Override
    protected void restoreState() {
        // TODO Auto-generated method stub

    }

    @Override
    protected void storeState() {
        // TODO Auto-generated method stub

    }


}
