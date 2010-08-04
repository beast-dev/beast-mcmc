package dr.app.beauti.options;

import dr.evolution.util.Taxa;
import dr.evolution.util.Taxon;

import java.util.ArrayList;
import java.util.List;

/**
 *
 */
public class TraitData extends PartitionData {
    public static final String TRAIT_SPECIES = "species";

    public static enum TraitType {
        DISCRETE,
        INTEGER,
        CONTINUOUS
    }

//    private String traitName = TraitOptions.Traits.TRAIT_SPECIES.toString();
    private TraitType traitType = TraitType.DISCRETE;


    public TraitData(BeautiOptions options, String traitName, String fileName, TraitType traitType) {
        super(options, traitName, fileName, null);
        this.traitType = traitType;

//        createTraitOptions();
    }   

//    private void createTraitOptions(){
//        if (traitType == TraitOptions.TraitType.DISCRETE) {
//            traitOptions = new DiscreteTraitOptions(this);
//        } else {
//            traitOptions = null; //TODO integer and continuous
//        }
//    }

    /////////////////////////////////////////////////////////////////////////

    public TraitType getTraitType() {
        return traitType;
    }

    public void setTraitType(TraitType traitType) {
        this.traitType = traitType;
    }

//    public TraitOptions getTraitOptions() {
//        return traitOptions;
//    }

    public int getSiteCount() {
        return 0;
    }

    public int getTaxaCount() {
        return options.taxonList.getTaxonCount();
    }

    public String getDataType() {
        return getTraitType().toString();
    }


    public boolean isSpecifiedTraitAnalysis(String traitName) {
        return  getName().equalsIgnoreCase(traitName);
    }

    public List<String> getStatesListOfTrait(Taxa taxonList) {
        List<String> states = new ArrayList<String>();
        String attr;

        if (taxonList != null) {
            for (int i = 0; i < taxonList.getTaxonCount(); i++) {
                Taxon taxon = taxonList.getTaxon(i);
                attr = (String) taxon.getAttribute(getName());

                if (attr == null) return null;

                if (!states.contains(attr)) {
                    states.add(attr);
                }
            }
            return states;
        } else {
            return null;
        }
    }

    public static String getPhylogeographicDescription() {
        return "Discrete phylogeographic inference in BEAST (PLoS Comput Biol. 2009 Sep;5(9):e1000520)";
    }

    
    public static List<String> getStatesListOfTrait(Taxa taxonList, String traitName) {
        List<String> states = new ArrayList<String>();
        String attr;

        if (taxonList != null) {
            for (int i = 0; i < taxonList.getTaxonCount(); i++) {
                Taxon taxon = taxonList.getTaxon(i);
                attr = (String) taxon.getAttribute(traitName);

                if (attr == null) return null;

                if (!states.contains(attr)) {
                    states.add(attr);
                }
            }
            return states;
        } else {
            return null;
        }
    }
    
}
