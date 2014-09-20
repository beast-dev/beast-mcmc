package dr.geo;

import dr.evolution.util.Taxa;
import dr.evolution.util.Taxon;
import dr.geo.math.SphericalPolarCoordinates;
import dr.inference.model.Statistic;
import dr.evoxml.TaxaParser;
import dr.evoxml.TaxonParser;
import dr.evoxml.DateParser;
import dr.xml.XMLParseException;
import dr.xml.XMLParser;
import dr.xml.AttributeParser;
import org.xml.sax.SAXException;

import javax.xml.parsers.ParserConfigurationException;
import java.io.IOException;
import java.io.File;
import java.io.FileReader;

/**
 * @author Alexei Drummond
 */
public class GreatCircleDistances {

    double[][] distances;

    public GreatCircleDistances(Taxa taxa, String attributeName) {

        distances = new double[taxa.getTaxonCount()][taxa.getTaxonCount()];

        for (int i = 0; i < taxa.getTaxonCount(); i++) {

            Taxon taxon = taxa.getTaxon(i);

            String attr = (String)taxon.getAttribute(attributeName);
            String[] loc = attr.split(" ");
            double latitude = Double.parseDouble(loc[0]);
            double longitude = Double.parseDouble(loc[1]);


            SphericalPolarCoordinates coord = new SphericalPolarCoordinates(latitude, longitude);
            //System.out.println(coord);

            for (int j = i+1; j < taxa.getTaxonCount(); j++) {
                Taxon taxon2 = taxa.getTaxon(j);
                attr = (String)taxon2.getAttribute(attributeName);
                String[] loc2 = attr.split(" ");
                latitude = Double.parseDouble(loc2[0]);
                longitude = Double.parseDouble(loc2[1]);

                SphericalPolarCoordinates coord2 = new SphericalPolarCoordinates(latitude, longitude);

                distances[i][j] = distances[j][i] = coord.distance(coord2);
            }
        }
    }

    class DistancesStatistic extends Statistic.Abstract {

        double[] dists = new double[distances.length*distances.length];

        public DistancesStatistic(boolean logTransformed) {

            int k = 0;
            for (double[] distance : distances) {
                for (int j = 0; j < distances.length; j++) {
                    dists[k] = distance[j];
                    if (logTransformed) dists[k] = Math.log(dists[k]);
                    k += 1;
                }
            }
        }

        public int getDimension() {
            return dists.length;
        }

        public double getStatisticValue(int dim) {
            return dists[dim];
        }
    }

    public Statistic getDistanceStatistic(boolean logTransformed) {
        return new DistancesStatistic(logTransformed);
    }

    public static void main(String[] args) throws ParserConfigurationException, IOException, SAXException, XMLParseException {

        XMLParser parser = new XMLParser(true, true);
        parser.addXMLObjectParser(new TaxonParser());
        parser.addXMLObjectParser(new TaxaParser());
        parser.addXMLObjectParser(new AttributeParser());
        parser.addXMLObjectParser(new DateParser());

        parser.parse(new FileReader(new File(args[0])), true);

        Taxa taxa = (Taxa)parser.getRoot().getChild(0);

        System.out.println("Found " + taxa.getTaxonCount() + " taxa");

        GreatCircleDistances distances = new GreatCircleDistances(taxa, "location");

        Statistic statistic = distances.getDistanceStatistic(true);

        for (int i = 0; i < statistic.getDimension(); i++) {
            System.out.println(statistic.getStatisticValue(i));
        }

    }
}
