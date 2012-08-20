package dr.app.seqgen;

import dr.evolution.tree.Tree;
import dr.evolution.tree.NodeRef;
import dr.evolution.datatype.Microsatellite;
import dr.evolution.alignment.Alignment;
import dr.evolution.alignment.Patterns;
import dr.evolution.sequence.Sequence;
import dr.evolution.util.Taxa;
import dr.evomodel.sitemodel.SiteModel;
import dr.evomodel.sitemodel.GammaSiteModel;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.substmodel.MicrosatelliteModel;
import dr.math.MathUtils;

/**
 * @author Chieh-Hsi Wu
 *
 * Simulates a pattern of microsatellites given a tree and microsatellite model
 *
 */
public class MicrosatelliteSimulator extends SequenceSimulator{
    private Taxa taxa;
    private Microsatellite dataType;

    public MicrosatelliteSimulator(
            Microsatellite dataType,
            Taxa taxa,
            Tree tree,
            MicrosatelliteModel msatModel,
            BranchRateModel branchRateModel){
        this(dataType, taxa, tree, new GammaSiteModel(msatModel), branchRateModel);

    }

    public MicrosatelliteSimulator(
            Microsatellite dataType,
            Taxa taxa,
            Tree tree,
            SiteModel siteModel,
            BranchRateModel branchRateModel) {

    	super(tree, siteModel, branchRateModel, 1);
        this.dataType = dataType;
        this.taxa = taxa;
    }

    /**
     * Convert integer representation of microsatellite length to string.
     */
	Sequence intArray2Sequence(int [] seq, NodeRef node) {
    	String sSeq = ""+seq[0];
		return new Sequence(m_tree.getNodeTaxon(node), sSeq);
    } // intArray2Sequence

    /**
     * Convert an alignment to a pattern
     */
    public Patterns simulateMsatPattern(){
        Alignment align = simulate();

        int[] pattern = new int[align.getTaxonCount()];
        for(int i = 0; i < pattern.length; i++){
            String taxonName = align.getSequence(i).getTaxon().getId();
            int index = taxa.getTaxonIndex(taxonName);
            pattern[index] = Integer.parseInt(align.getSequence(i).getSequenceString());
        }
        Patterns patterns = new Patterns(dataType,taxa);
        patterns.addPattern(pattern);
        for(int i = 0; i < pattern.length;i++){
            System.out.print(pattern[i]+",");
        }
        System.out.println();
        
        return patterns;
    }

}
