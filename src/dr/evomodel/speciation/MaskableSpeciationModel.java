package dr.evomodel.speciation;

/**
 * Created by IntelliJ IDEA.
 * User: adru001
 * Date: Nov 2, 2010
 * Time: 10:56:03 AM
 * To change this template use File | Settings | File Templates.
 */
public abstract class MaskableSpeciationModel extends SpeciationModel {

    public MaskableSpeciationModel(String name, Type units) {
        super(name, units);
    }

    // a model specific implementation that allows this speciation model
    // to be partially masked by another -- useful in model averaging applications
    public abstract void mask(SpeciationModel mask);

    public abstract void unmask();
}
