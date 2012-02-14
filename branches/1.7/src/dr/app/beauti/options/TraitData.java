package dr.app.beauti.options;

import dr.evolution.datatype.ContinuousDataType;
import dr.evolution.datatype.DataType;
import dr.evolution.datatype.GeneralDataType;
import dr.evolution.util.Taxa;
import dr.evolution.util.Taxon;

import java.io.Serializable;
import java.util.Set;
import java.util.TreeSet;

/**
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @author Walter Xie
 */
public class TraitData implements Serializable {
    public static final String TRAIT_SPECIES = "species";

    public static enum TraitType {
        DISCRETE,
        INTEGER,
        CONTINUOUS;

        public String toString() {
            return name().toLowerCase();
        }
    }

    private TraitType traitType = TraitType.DISCRETE;

    private final String fileName;
    private String name;

    protected final BeautiOptions options;

    public TraitData(BeautiOptions options, String name, String fileName, TraitType traitType) {
        this.options = options;
        this.name = name;
        this.fileName = fileName;
        this.traitType = traitType;
    }

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

    public Taxon getTaxon(int i) {
        return options.taxonList.getTaxon(i);
    }

    public boolean hasValue(int i) {
        if (options.taxonList == null || options.taxonList.getTaxon(i) == null
                || options.taxonList.getTaxon(i).getAttribute(getName()) == null) return false;
        return options.taxonList.getTaxon(i).getAttribute(getName()).toString().trim().length() > 0;
    }

    public DataType
    getDataType() {
        switch (traitType) {
            case DISCRETE:
                return GeneralDataType.INSTANCE;
            case INTEGER:
                return GeneralDataType.INSTANCE;
            case CONTINUOUS:
                return ContinuousDataType.INSTANCE;
        }
        throw new IllegalArgumentException("Unknown trait type");
    }

    public String getFileName() {
        return fileName;
    }

    public String getName() {
        return name;
    }

    public void setName(final String name) {
        this.name = name;
    }

    public Set<String> getStatesOfTrait(Taxa taxonList) {
        return getStatesListOfTrait(taxonList, getName());
    }


    public static Set<String> getStatesListOfTrait(Taxa taxonList, String traitName) {
        Set<String> states = new TreeSet<String>();

        if (taxonList == null) {
            throw new IllegalArgumentException("taxon list is null");
        }

        for (int i = 0; i < taxonList.getTaxonCount(); i++) {
            Taxon taxon = taxonList.getTaxon(i);
            String attr = (String) taxon.getAttribute(traitName);

            // ? is used to denote missing data so is not a state...
            if (attr != null && !attr.equals("?")) {
                states.add(attr);
            }
        }


        return states;
    }

    public static int getEmptyStateIndex(Taxa taxonList, String traitName) {
        if (taxonList == null) {
            throw new IllegalArgumentException("taxon list is null");
        }

        for (int i = 0; i < taxonList.getTaxonCount(); i++) {
            Taxon taxon = taxonList.getTaxon(i);
            String attr = (String) taxon.getAttribute(traitName);

            // ? is used to denote missing data so is not a state...
            if (attr == null || attr.equals("?")) {
                return i;
            }
        }

        return -1;
    }

    public String toString() {
        return name;
    }
}
