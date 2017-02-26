/*
 * LewisMkSubstitutionModelParser.java
 *
 * Copyright (c) 2002-2016 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.evomodelxml.substmodel;

import dr.evomodel.substmodel.FrequencyModel;
import dr.evomodel.substmodel.GeneralSubstitutionModel;
import dr.evolution.datatype.DataType;
import dr.inference.model.Parameter;
import dr.xml.*;

import java.util.Stack;

/**
 * Package: LewisMkSubstitutionModelParser
 * Description:
 * <p/>
 * <p/>
 * Created by
 *
 * @author Alexander V. Alekseyenko (alexander.alekseyenko@gmail.com)
 *         Date: Oct 13, 2009
 *         Time: 4:01:38 PM
 */
public class LewisMkSubstitutionModelParser extends AbstractXMLObjectParser {
    public static final String LEWIS_MK_MODEL = "lewisMk";
    public static final String TOTAL_ORDER = "totalOrder";
    public static final String FREQUENCIES = "frequencies";
    public static final String ORDER = "order";
    public static final String STATE = "state";
    public static final String ADJACENT = "adjacentTo";

    //public static XMLObjectParser PARSER=new LewisMkSubstitutionModelParser();

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        XMLObject cxo = xo.getChild(FREQUENCIES);
        FrequencyModel freqModel = (FrequencyModel) cxo.getChild(FrequencyModel.class);
        DataType dataType = freqModel.getDataType();
        int k = dataType.getStateCount();
        System.err.println("Number of states " + k);
        Parameter ratesParameter;
        if (xo.hasAttribute(TOTAL_ORDER) && xo.getBooleanAttribute(TOTAL_ORDER)) { //TOTAL ORDERING OF THE STATES BASED ON DATATYPE
            ratesParameter = new Parameter.Default(k * (k - 1) / 2, 0);
            int j = k - 1;
            for (int i = 0; i < (k - 1) * k / 2; i = i + j + 1) {
                ratesParameter.setParameterValue(i, 1);
                j -= 1;
            }
        } else if (xo.hasChildNamed(ORDER)) { // USER-SPECIFIED ORDERING OF THE STATES
            ratesParameter = new Parameter.Default(k * (k - 1) / 2, 0);
            for (int i = 0; i < xo.getChildCount(); ++i) {
                if (xo.getChildName(i).equals(ORDER)) {
                    cxo = (XMLObject) xo.getChild(i);
                    if (cxo.getName().equals(ORDER)) {
                        int from = dataType.getState(cxo.getStringAttribute(STATE).charAt(0));
                        int to = dataType.getState(cxo.getStringAttribute(ADJACENT).charAt(0));
                        if (from > to) {//SWAP: from should have the smaller state number
                            to += from;
                            from = to - from;
                            to -= from;
                        }
                        int ratesIndex = (from * (2 * k - 3) - from * from) / 2 + to - 1;
                        ratesParameter.setParameterValue(ratesIndex, 1);
                    }
                }
            }
        } else {
            ratesParameter = new Parameter.Default(k * (k - 1) / 2, 1);
        }
        System.err.println(ratesParameter.toString());
        System.err.println("Infinitesimal matrix:");
        for (int i = 0; i < k; ++i) {
            for (int j = 0; j < k; ++j) {
                int from, to;
                if (i < j) {
                    from = i;
                    to = j;
                } else {
                    from = j;
                    to = i;
                }


                int ratesIndex = (from * (2 * k - 3) - from * from) / 2 + to - 1;    //This is right now!!! Thanks, Marc!

                if (i != j)
                    System.err.print(Double.toString(ratesParameter.getValue(ratesIndex)) + "\t(" + ratesIndex + ")\t");
                else System.err.print("-\t\t");

            }
            System.err.println("");//newline
        }
        System.err.println("");

        if (!checkConnected(ratesParameter.getValues(), k)) {
            throw (new XMLParseException("The state transitions form a disconnected graph! This model is not suited for this case."));
        }

        return new GeneralSubstitutionModel(LEWIS_MK_MODEL, dataType, freqModel, ratesParameter, -1);
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private boolean checkConnected(Double rates[], int states) {
        boolean[] visited = new boolean[states];
        Stack<Integer> open = new Stack<Integer>();

        open.push(0);

        for (int i = 1; i < states; ++i) {
            visited[i] = false;
        }
        visited[0] = true;

        while (!open.empty()) {
            int current = open.pop();
            for (int j = 0; j < states; ++j) {
                int rateIndex;
                if (current < j) rateIndex = (current * (2 * states - 3) - current * current) / 2 + j - 1;
                else rateIndex = (j * (2 * states - 3) - j * j) / 2 + current - 1;
                if (current != j && !visited[j] && rates[rateIndex] != 0.0) {
                    visited[j] = true;
                    open.push(j);
                }
            }
        }
        for (int i = 0; i < states; ++i) {
            if (!visited[i]) return false;
        }
        return true;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(FREQUENCIES, FrequencyModel.class),
            AttributeRule.newBooleanRule(TOTAL_ORDER, true),
            new ElementRule(ORDER,
                    new XMLSyntaxRule[]{AttributeRule.newStringRule(STATE, false),
                            AttributeRule.newStringRule(ADJACENT, false)}, 0, Integer.MAX_VALUE)
    };

    public String getParserDescription() {
        return "A parser for Lewis's Mk model";
    }

    public Class getReturnType() {
        return GeneralSubstitutionModel.class;
    }

    public String getParserName() {
        return LEWIS_MK_MODEL;
    }
}
