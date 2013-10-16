package dr.app.seqgen;

import dr.evolution.util.Taxa;
import dr.evolution.datatype.Microsatellite;
import dr.evolution.tree.Tree;
import dr.evolution.tree.NodeRef;
import dr.evolution.sequence.Sequence;
import dr.evomodel.sitemodel.SiteModel;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.math.MathUtils;
import dr.inference.model.Parameter;

import java.util.ArrayList;

/**
 * @author Chieh-Hsi Wu
 */
public class GeoDiffusionSimulator extends SequenceSimulator{
    public static final int LATITUDE_INDEX = 0;
    public static final int LONGITUDE_INDEX = 1;
    private Taxa taxa;
    private Microsatellite dataType;
    private double maxLat;
    private double minLat;
    private double maxLong;
    private double minLong;

    public GeoDiffusionSimulator(
            Microsatellite dataType,
            Taxa taxa,
            Tree tree,
            SiteModel siteModel,
            BranchRateModel branchRateModel,
            double maxLat,
            double minLat,
            double maxLong,
            double minLong) {

    	super(tree, siteModel, branchRateModel, 1);
        this.dataType = dataType;
        this.taxa = taxa;
        this.maxLat = maxLat;
        this.minLat = minLat;
        this.maxLong = maxLong;
        this.minLong = minLong;

    }

    /**
     * Convert integer representation of microsatellite length to string.
     */
	Sequence intArray2Sequence(int [] seq, NodeRef node) {
    	String sSeq = ""+seq[0];
		return new Sequence(m_tree.getNodeTaxon(node), sSeq);
    } // intArray2Sequence

    public double[][] simulateLocations() {
    	NodeRef root =  m_tree.getRoot();

    	//assume uniform
    	double[][] latLongs = new double[m_tree.getNodeCount()][2];
        double rootLat = MathUtils.nextDouble()*(maxLat-minLat)+minLat;
        double rootLong = MathUtils.nextDouble()*(maxLong-minLong)+minLong;
        int rootNum = root.getNumber();
        latLongs[rootNum] [LATITUDE_INDEX] = rootLat;
        latLongs[rootNum] [LONGITUDE_INDEX] = rootLong;
    	traverse(root, latLongs[rootNum], latLongs);



    	return latLongs;
    }

    void traverse(NodeRef node, double [] parentSequence, double[][] latLongs) {
		for (int iChild = 0; iChild < m_tree.getChildCount(node); iChild++) {
			NodeRef child = m_tree.getChild(node, iChild);

            //find the branch length
            final double branchRate = m_branchRateModel.getBranchRate(m_tree, child);
            final double branchLength = branchRate * (m_tree.getNodeHeight(node) - m_tree.getNodeHeight(child));
            if (branchLength < 0.0) {
                        throw new RuntimeException("Negative branch length: " + branchLength);
            }

            double childLat = MathUtils.nextGaussian()*Math.sqrt(branchLength)+parentSequence[LATITUDE_INDEX];
            double childLong = MathUtils.nextGaussian()*Math.sqrt(branchLength)+parentSequence[LONGITUDE_INDEX];
            int childNum = child.getNumber();

        	latLongs[childNum][LATITUDE_INDEX] = childLat;
            latLongs[childNum][LONGITUDE_INDEX] = childLong;

			traverse(m_tree.getChild(node, iChild), latLongs[childNum], latLongs);
		}
	}

    /**
     * Convert an alignment to a pattern
     */
    public ArrayList simulateGeoAttr(){
        double[][] locations = simulateLocations();

        ArrayList<Parameter> locationList = new ArrayList<Parameter>();
        for(int i = 0; i < m_tree.getExternalNodeCount(); i++){
            NodeRef node = m_tree.getNode(i);
            String taxaName = m_tree.getTaxon(node.getNumber()).getId();
            Parameter location = new Parameter.Default(locations[i]);
            System.out.println("taxon: "+taxaName+", lat: "+locations[i][0]+", long: "+locations[i][1]);
            locationList.add(location);
        }
        
        return locationList;
    }

}
