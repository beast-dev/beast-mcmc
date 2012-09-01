package dr.evomodel.operators;

import dr.evomodel.speciation.AlloppSpeciesBindings;
import dr.evomodel.speciation.AlloppSpeciesNetworkModel;
import dr.evomodelxml.operators.AlloppHybPopSizesScaleParser;
import dr.inference.operators.OperatorFailedException;
import dr.inference.operators.SimpleMCMCOperator;
import dr.math.MathUtils;

/**
 * Created with IntelliJ IDEA.
 * User: Graham
 * Date: 03/08/12
 */
public class AlloppHybPopSizesScale  extends SimpleMCMCOperator {

    private final AlloppSpeciesNetworkModel apspnet;
    private final AlloppSpeciesBindings apsp;
    private final double scalingFactor;

    public AlloppHybPopSizesScale(AlloppSpeciesNetworkModel apspnet, AlloppSpeciesBindings apsp, double scalingFactor, double weight) {
        this.apspnet = apspnet;
        this.apsp = apsp;
        this.scalingFactor = scalingFactor;
        setWeight(weight);
    }


    public String getPerformanceSuggestion() {
        return "None";
    }

    @Override
    public String getOperatorName() {
        return AlloppHybPopSizesScaleParser.HYB_POP_SIZES_SCALE + "(" + apspnet.getId() +
                "," + apsp.getId() + ")";
    }

    @Override
    public double doOperation() throws OperatorFailedException {
        apspnet.beginNetworkEdit();
        double b = (1.0-scalingFactor) * (1.0-scalingFactor) / scalingFactor;
        double c = scalingFactor / (1.0-scalingFactor);
        double y = MathUtils.nextDouble();
        double s = b * (y+c) * (y+c);
        int i = MathUtils.nextInt(apspnet.getNumberOfTetraTrees());
        apspnet.setOneHybPopValue(i, s * apspnet.getOneHybPopValue(i));
        apspnet.endNetworkEdit();
        return 0.0;  // this way of scaling, with proposal proportional to x^-(1/2) has hastings ratio 1
    }

}
