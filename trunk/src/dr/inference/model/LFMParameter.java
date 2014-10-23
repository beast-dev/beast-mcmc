package dr.inference.model;

/**
 * Created by max on 10/22/14.
 */
//Designed to return a data matrix post computation if asked. Designed for latent liabilities
public class LFMParameter extends CompoundParameter {
    LatentFactorModel LFM;

    public LFMParameter(LatentFactorModel LFM) {
        super(null);
        this.LFM = LFM;
    }

    @Override
    public int getDimension(){
        return LFM.getData().getDimension();
    }

    public Parameter getParameter(int PID) {
        return LFM.getDataWithLatent(PID);
    }
};

