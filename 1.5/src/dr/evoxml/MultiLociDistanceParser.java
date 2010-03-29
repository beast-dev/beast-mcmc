package dr.evoxml;

import dr.xml.*;
import dr.evolution.alignment.Patterns;
import dr.evolution.distance.DistanceMatrix;
import java.util.ArrayList;

/**
 * @author Chieh-Hsi Wu
 *
 * Calculates overall mean pairwise distances across multiple loci
 * 
 */
public class MultiLociDistanceParser extends AbstractXMLObjectParser {

    public static final String MULTI_LOCI_DISTANCE = "multiLociDistance";
    public String getParserName() {
        return MULTI_LOCI_DISTANCE; 
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        int childNum = xo.getChildCount();
        ArrayList<DistanceMatrix> distMatList = new ArrayList<DistanceMatrix>();
        for(int i = 0; i < childNum; i++){
            distMatList.add((DistanceMatrix)xo.getChild(i));
        }
        double[] meanDists = new double[childNum];
        for(int i = 0; i < childNum; i++){
            meanDists[i] = distMatList.get(i).getMeanDistance();
        }

        double sum = 0.0;
        for(int i = 0; i < childNum; i++){
            sum = sum+meanDists[i];
        }
        double mean = sum/(double)childNum;

        printMeans(distMatList,mean);

        return mean;

    }

    public void printMeans(ArrayList<DistanceMatrix> distMat, double mean){
        System.out.println("Individual mean distances:");
        for(int i = 0; i < distMat.size(); i++){
            System.out.print(distMat.get(i).getId()+": ");
            System.out.println(distMat.get(i).getMeanDistance());
        }
        System.out.println("overallMean: "+mean);
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
        new ElementRule(DistanceMatrix.class, 1, Integer.MAX_VALUE)
    };

    public String getParserDescription() {
        return "Constructs a distance matrix from a pattern list or alignment";
    }

    public Class getReturnType() { return Double.class; }
}
