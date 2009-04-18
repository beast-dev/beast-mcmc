package dr.inference.model;

/**
 * @author joseph
 *         Date: 16/04/2009
 */
public abstract class AbstractModelLikelihood extends AbstractModel implements Likelihood {
    /**
     * @param name Model Name
     */
    public AbstractModelLikelihood(String name) {
        super(name);
    }

    public String prettyName() {
         return Likelihood.Abstract.getPrettyName(this);
     }
}
