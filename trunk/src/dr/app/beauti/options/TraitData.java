package dr.app.beauti.options;

/**
 *
 */
public class TraitData extends PartitionData {
    private TraitOptions traitOptions;

    private String traitName = TraitOptions.Traits.TRAIT_SPECIES.toString();
    private TraitOptions.TraitType traitType = TraitOptions.TraitType.DISCRETE;


    public TraitData(String traitName, String fileName, TraitOptions.TraitType traitType) {
        super(traitName, fileName, null);
        this.traitType = traitType;

        createTraitOptions();
    }   

    private void createTraitOptions(){
        if (traitType == TraitOptions.TraitType.DISCRETE) {
            traitOptions = new DiscreteTraitOptions(this);
        } else {
            traitOptions = null; //TODO integer and continuous
        }
    }

    /////////////////////////////////////////////////////////////////////////
    public String getTraitName() {
        return traitName;
    }

    public void setTraitName(String traitName) {
        this.traitName = traitName;
    }

    public TraitOptions.TraitType getTraitType() {
        return traitType;
    }

    public void setTraitType(TraitOptions.TraitType traitType) {
        this.traitType = traitType;
        createTraitOptions();
    }

    public TraitOptions getTraitOptions() {
        return traitOptions;
    }
    
}
