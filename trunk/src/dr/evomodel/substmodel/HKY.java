/*
 * HKY.java
 *
 * Copyright (C) 2002-2006 Alexei Drummond and Andrew Rambaut
 *
 * This file is part of BEAST.
 * See the NOTICE file distributed with this work for additional
 * information regarding copyright ownership and licensing.
 *
 * BEAST is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 *  BEAST is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with BEAST; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA  02110-1301  USA
 */

package dr.evomodel.substmodel;

import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Statistic;
import dr.xml.*;

import java.util.logging.Logger;

/**
 * Hasegawa-Kishino-Yano model of nucleotide evolution
 *
 * @version $Id: HKY.java,v 1.42 2005/09/23 13:17:59 rambaut Exp $
 *
 * @author Alexei Drummond
 * @author Andrew Rambaut
 */
public class HKY extends AbstractNucleotideModel
{
    public static final String HKY_MODEL = "hkyModel";
    public static final String KAPPA = "kappa";
    public static final String FREQUENCIES = "frequencies";

    /** tsTv */
    private double tsTv;

    private Parameter kappaParameter = null;

    private boolean updateIntermediates = true;

    /** Used for precalculations */
    private double 				beta, A_R, A_Y;
    private double 				tab1A, tab2A, tab3A;
    private double 				tab1C, tab2C, tab3C;
    private double 				tab1G, tab2G, tab3G;
    private double 				tab1T, tab2T, tab3T;

    /**
     * Constructor
     */
    public HKY(Parameter kappaParameter, FrequencyModel freqModel) {

        super(HKY_MODEL, freqModel);
        this.kappaParameter = kappaParameter;
        addParameter(kappaParameter);
        kappaParameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, 1));
        updateIntermediates = true;

        addStatistic(tsTvStatistic);
    }

    /**
     * set kappa
     */
    public void setKappa(double kappa)
    {
        kappaParameter.setParameterValue(0, kappa);
        updateMatrix = true;
    }

    /**
     * @return kappa
     */
    public final double getKappa() { return kappaParameter.getParameterValue(0); }

    /**
     * set ts/tv
     */
    public void setTsTv(double tsTv) {
        this.tsTv = tsTv;
          calculateFreqRY();
        setKappa((tsTv*freqR*freqY)/(freqA*freqG+freqC*freqT));
    }

    /**
     * @return tsTv
     */
    public double getTsTv()
    {
        calculateFreqRY();
        tsTv = (getKappa() * (freqA*freqG + freqC*freqT))/(freqR*freqY);

        return tsTv;
    }

    protected void frequenciesChanged() {
        // frequencyModel changed
        updateIntermediates = true;
    }

    private void calculateIntermediates() {

        calculateFreqRY();

        tab1A = freqA*((1/freqR)-1);
        tab2A = (freqR-freqA)/freqR;
        tab3A = freqA/freqR;
        tab1C = freqC*((1/freqY)-1);
        tab2C = (freqY-freqC)/freqY;
        tab3C = freqC/freqY;
        tab1G = freqG*((1/freqR)-1);
        tab2G = (freqR-freqG)/freqR;
        tab3G = freqG/freqR;
        tab1T = freqT*((1/freqY)-1);
        tab2T = (freqY-freqT)/freqY;
        tab3T = freqT/freqY;

        updateMatrix = true;
        updateIntermediates = false;
    }

    /**
     * get the complete transition probability matrix for the given distance
     *
     * @param distance the expected number of substitutions
     * @param matrix an array to store the matrix
     */
    public void getTransitionProbabilities(double distance, double[] matrix)
    {
        if (updateIntermediates) {
            calculateIntermediates();
        }

        if (updateMatrix) {
            setupMatrix();
        }

        final double xx = beta * distance;
        final double bbR = Math.exp(xx*A_R);
        final double bbY = Math.exp(xx*A_Y);

        final double aa = Math.exp(xx);
        final double oneminusa = (1 - aa);

        matrix[0] =	freqA+(tab1A*aa)+(tab2A*bbR);

        matrix[1] =	freqC* oneminusa;
        matrix[2] =	freqG+(tab1G*aa)-(tab3G*bbR);
        matrix[3] =	freqT* oneminusa;

        matrix[4] =	freqA* oneminusa;
        matrix[5] =	freqC+(tab1C*aa)+(tab2C*bbY);
        matrix[6] =	freqG* oneminusa;
        matrix[7] =	freqT+(tab1T*aa)-(tab3T*bbY);

        matrix[8] =	freqA+(tab1A*aa)-(tab3A*bbR);
        matrix[9] =	matrix[1];
        matrix[10] =freqG+(tab1G*aa)+(tab2G*bbR);
        matrix[11] =matrix[3];

        matrix[12] =matrix[4];
        matrix[13] =freqC+(tab1C*aa)-(tab3C*bbY);
        matrix[14] =matrix[6];
        matrix[15] =freqT+(tab1T*aa)+(tab2T*bbY);
    }

    /**
     * setup substitution matrix
     */
    protected void setupMatrix()
    {
        double kappa = getKappa();
        beta = -1.0 / (2.0*(freqR*freqY + kappa*(freqA*freqG +
                                                    freqC*freqT)));

        A_R = 1.0 +freqR*(kappa-1);
        A_Y = 1.0 +freqY*(kappa-1);

        updateMatrix = false;
    }

     protected void setupRelativeRates() { }

    // *****************************************************************
    // Interface Model
    // *****************************************************************

    /**
     * Restore the stored state
     */
    public void restoreState() {
        super.restoreState();
        updateIntermediates = true;
    }

    /**
     * Restore the stored state
     */
    public void adoptState(Model source) {
        super.adoptState(source);
        updateIntermediates = true;
    }

    // **************************************************************
    // XHTMLable IMPLEMENTATION
    // **************************************************************

    public String toXHTML() {
        StringBuffer buffer = new StringBuffer();

        buffer.append("<em>HKY Model</em> Ts/Tv = ");
        buffer.append(getTsTv());
        buffer.append(" (kappa = ");
        buffer.append(getKappa());
        buffer.append(")");

        return buffer.toString();
    }

    /**
     * Parses an element from an DOM document into a DemographicModel. Recognises
     * ConstantPopulation and ExponentialGrowth.
     */
    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() { return HKY_MODEL; }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            Parameter kappaParam = (Parameter)xo.getSocketChild(KAPPA);
            FrequencyModel freqModel = (FrequencyModel)xo.getSocketChild(FREQUENCIES);

            Logger.getLogger("dr.evomodel").info("Creating HKY substitution model. Initial kappa = " + kappaParam.getParameterValue(0));

            return new HKY(kappaParam, freqModel);
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return "This element represents an instance of the HKY85 (Hasegawa, Kishino & Yano, 1985) model of nucleotide evolution.";
        }

        public Class getReturnType() { return HKY.class; }

        public XMLSyntaxRule[] getSyntaxRules() { return rules; }

        private XMLSyntaxRule[] rules = new XMLSyntaxRule[] {
            new ElementRule(FREQUENCIES,
                new XMLSyntaxRule[] { new ElementRule(FrequencyModel.class) }),
            new ElementRule(KAPPA,
                new XMLSyntaxRule[] { new ElementRule(Parameter.class) })
        };

    };

    //
    // Private stuff
    //

    private Statistic tsTvStatistic = new Statistic.Abstract() {

        public String getStatisticName() {
            return "tsTv";
        }

        public int getDimension() { return 1; }

        public double getStatisticValue(int dim) {
            return getTsTv();
        }

    };
}
