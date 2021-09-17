package dr.inferencexml.operators.hmc;

import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.model.Parameter;
import dr.inference.operators.AdaptationMode;
import dr.inference.operators.hmc.GeodesicHamiltonianMonteCarloOperator;
import dr.inference.operators.hmc.HamiltonianMonteCarloOperator;
import dr.inference.operators.hmc.MassPreconditionScheduler;
import dr.inference.operators.hmc.MassPreconditioner;
import dr.util.Transform;
import dr.xml.XMLObject;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;

import java.util.ArrayList;

/**
 * @author Gabriel Hassler
 * @author Marc A. Suchard
 */

public class GeodesicHamiltonianMonteCarloOperatorParser extends HamiltonianMonteCarloOperatorParser {
    public final static String OPERATOR_NAME = "geodesicHamiltonianMonteCarloOperator";
    public final static String ORTHOGONALITY_STRUCTURE = "orthogonalityStructure";
    public final static String ROWS = "rows";

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        GeodesicHamiltonianMonteCarloOperator hmc = (GeodesicHamiltonianMonteCarloOperator) super.parseXMLObject(xo);
        if (xo.hasChildNamed(ORTHOGONALITY_STRUCTURE)) {
            XMLObject cxo = xo.getChild(ORTHOGONALITY_STRUCTURE);
            ArrayList<ArrayList<Integer>> orthogonalityStructure = new ArrayList<>();
            for (int i = 0; i < cxo.getChildCount(); i++) {
                XMLObject group = (XMLObject) cxo.getChild(i);
                int[] rows = group.getIntegerArrayAttribute(ROWS);
                ArrayList<Integer> rowList = new ArrayList<>();

                for (int j = 0; j < rows.length; j++) {
                    rowList.add(rows[j] - 1);
                }

                orthogonalityStructure.add(rowList);
            }

            hmc.setOrthogonalityStructure(orthogonalityStructure);
        }

        return hmc;
    }

    @Override
    protected HamiltonianMonteCarloOperator factory(AdaptationMode adaptationMode, double weight, GradientWrtParameterProvider derivative,
                                                    Parameter parameter, Transform transform, Parameter mask,
                                                    HamiltonianMonteCarloOperator.Options runtimeOptions,
                                                    MassPreconditioner preconditioner, MassPreconditionScheduler.Type schedulerType) {
        return new GeodesicHamiltonianMonteCarloOperator(adaptationMode, weight, derivative,
                parameter, transform, mask,
                runtimeOptions, preconditioner);
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    } //TODO: add orthogonality structure rules


    @Override
    public String getParserDescription() {
        return "Returns a geodesic Hamiltonian Monte Carlo transition kernel";
    }

    @Override
    public Class getReturnType() {
        return GeodesicHamiltonianMonteCarloOperator.class;
    }

    @Override
    public String getParserName() {
        return OPERATOR_NAME;
    }
}
