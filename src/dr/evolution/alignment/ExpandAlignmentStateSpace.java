package dr.evolution.alignment;

import dr.evolution.datatype.*;
import dr.evolution.util.TaxonList;

import java.util.Arrays;

public class ExpandAlignmentStateSpace extends WrappedAlignment {

    /**
     * Constructor.
     */
    public ExpandAlignmentStateSpace(DataType dataType) {
        this(dataType, null, null);
    }

    /**
     * Constructor.
     */
    public ExpandAlignmentStateSpace(DataType dataType, PatternList trait) {
        this(dataType, null, trait);
    }

    /**
     * Constructor.
     */
    public ExpandAlignmentStateSpace(DataType dataType, Alignment alignment) {
        this(dataType, alignment, null);
    }

    /**
     * Constructor.
     */
    public ExpandAlignmentStateSpace(DataType dataType, Alignment alignment, PatternList trait) {
        super(alignment);
        if ( alignment.getTaxonCount() != trait.getTaxonCount() ) {throw new RuntimeException("Alignment and trait must have same number of taxa.");}
        if ( !checkTaxa(alignment, trait) ) {throw new RuntimeException("Alignment and trait must have same taxa.");}
        setDataType(dataType);
        setExpansionIndices(trait);
        setAlignment(alignment);
    }

    /**
     * Sets the expansionIndices of this alignment.
     */
    public void setExpansionIndices(PatternList trait) {

        int[] traitValues = new int[trait.getTaxonCount()];

        for (int i = 0; i < traitValues.length; i++) {
            traitValues[i] = trait.getPatternState(i,0);
        }

        this.numDependentStates = trait.getStateCount();
        this.expansionIndices = traitValues;
    }

    /**
     * Sets the contained.
     */
    public void setAlignment(Alignment alignment) {
        if (dataType == null) {
            dataType = alignment.getDataType();
        }

        this.alignment = alignment;

        int originalType = alignment.getDataType().getType();

        if ( ((StateDependentDataType)dataType).getDependentClassCount() != numDependentStates ) {
            throw new RuntimeException("Number of dependent states found (" + numDependentStates + ") incompatible with data type " + dataType);
        }

        if (originalType == DataType.NUCLEOTIDES) {

            if ( !(dataType instanceof StateDependentNucleotides) ) {
                throw new RuntimeException("Incompatible alignment DataType for ExpandedStateSpaceAlignment");
            }

        } else {
            // TODO generalize to other datatypes

            throw new RuntimeException("Incompatible alignment DataType for ExpandedStateSpaceAlignment");

        }//END: original type check
    }

    /**
     * Sets the dataType of this alignment. This can be different from
     * the dataTypes of the contained alignment - they will be translated
     * as required.
     */
    public void setDataType(DataType dataType) {
        this.dataType = dataType;
    }

    /**
     * @return the DataType of this siteList
     */
    public DataType getDataType() {
        return dataType;
    }

    /**
     * @return number of sites
     */
    public int getSiteCount() {
        if (alignment == null) throw new RuntimeException("ExpandAlignmentStateSpace has no alignment");

        return alignment.getSiteCount();
    }

    /**
     * @return the sequence state at (taxon, site)
     */
    public int getState(int taxonIndex, int siteIndex) {
        if (alignment == null) throw new RuntimeException("ExpandAlignmentStateSpace has no alignment");

        int newType = dataType.getType();
        int originalType = alignment.getDataType().getType();

        int state = 0;

        if (originalType == DataType.NUCLEOTIDES) {
            int nucState = alignment.getState(taxonIndex, siteIndex);
            state = ((StateDependentNucleotides)dataType).getState(nucState,expansionIndices[taxonIndex]);
        }

        return state;
    }

    private boolean checkTaxa(Alignment alignment, PatternList trait) {
        boolean allMatch = true;

        String[] traitTaxa = new String[alignment.getTaxonCount()];
        for ( int i = 0; i < alignment.getTaxonCount(); i++ ) {
            traitTaxa[i] = trait.getTaxon(i).getId();
        }

        for ( int i = 0; i < alignment.getTaxonCount(); i++ ) {
            String taxon = alignment.getTaxon(i).getId();
            // Best guess, alignments in same order so our match is here
            if ( !(taxon == traitTaxa[i]) ) {
                // Fall back to search
                if ( !Arrays.asList(traitTaxa).contains(taxon) ) {
                    allMatch = false;
                    break;
                }
            }
        }

        return allMatch;
    }

    private DataType dataType = null;
    private int numDependentStates;
    private int[] expansionIndices = null;
}
