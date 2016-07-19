/*
 * MsatBMAParser.java
 *
 * Copyright (c) 2002-2015 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.oldevomodelxml.substmodel;

import dr.xml.*;
import dr.evolution.datatype.Microsatellite;
import dr.inference.model.Parameter;
import dr.oldevomodel.substmodel.MsatBMA;

import java.util.HashMap;
import java.util.ArrayList;

/**
 * @author Chieh-Hsi Wu
 *
 * Parser of MsatAveragingSubsetModel
 */
public class MsatBMAParser extends AbstractXMLObjectParser{

    public static final String MODELS = "models";
    public static final String MODEL = "model";
    public static final String BINARY = "binary";
    public static final String CODE = "code";
    public static final String LOGIT = "logit";
    public static final String RATE_PROPS = "rateProps";
    public static final String RATE_PROP = "rateProp";
    public static final String RATE_QUADS = "rateQuads";
    public static final String RATE_QUAD = "rateQuad";
    public static final String BIAS_CONSTS = "biasConsts";
    public static final String BIAS_CONST = "biasConst";
    public static final String BIAS_LINS = "biasLins";
    public static final String BIAS_LIN = "biasLin";
    public static final String GEOS = "geos";
    public static final String GEO = "geo";
    public static final String PHASE_PROBS = "phaseProbs";
    public static final String PHASE_PROB = "phaseProb";
    public static final String IN_MODELS = "inModels";
    public static final String MODEL_CHOOSE = "modelChoose";
    public static final String MODEL_INDICATOR = "modelIndicator";
    public static final String MSAT_BMA = "msatBMA";
    public static final int PROP_RATES_MAX_COUNT = 6;
    public static final int BIAS_CONST_MAX_COUNT = 8;
    public static final int BIAS_LIN_MAX_COUNT = 4;
    public static final int GEO_MAX_COUNT = 6;



    public String getParserName(){
        return MSAT_BMA;
    }


    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        //get microsatellite data type
        Microsatellite dataType = (Microsatellite)xo.getChild(Microsatellite.class);

        //whether mutational bias is in logit space
        boolean logit = xo.getAttribute(LOGIT,true);

        XMLObject modelsXO = xo.getChild(MODELS);
        int modelCount = modelsXO.getChildCount();

        HashMap<Integer, Integer> modelBitIndMap = new HashMap<Integer, Integer>(modelCount);
        for(int i = 0; i < modelCount; i++){

            XMLObject modelXO = (XMLObject)modelsXO.getChild(i);
            String bitVec = modelXO.getStringAttribute(BINARY);
            int bitVecVal = Integer.parseInt(bitVec,2);
            int modelCode = modelXO.getIntegerAttribute(CODE);
            modelBitIndMap.put(bitVecVal,modelCode);

        }

        Parameter[][] paramModelMap = new Parameter[6][modelCount];



        XMLObject propRatesXO = xo.getChild(RATE_PROPS);
        ArrayList<Parameter> rateProps =
                processParameters(
                        propRatesXO,
                        paramModelMap,
                        MsatBMA.PROP_INDEX
                );

        XMLObject quadRatesXO = xo.getChild(RATE_QUADS);
        ArrayList<Parameter> rateQuads =
                processParameters(
                        quadRatesXO,
                        paramModelMap,
                        MsatBMA.QUAD_INDEX
                );

        XMLObject biasConstsXO = xo.getChild(BIAS_CONSTS);
        ArrayList<Parameter> biasConsts =
                processParameters(
                        biasConstsXO,
                        paramModelMap,
                        MsatBMA.BIAS_CONST_INDEX
                );

        XMLObject biasLinsXO = xo.getChild(BIAS_LINS);
        ArrayList<Parameter> biasLins =
                processParameters(
                        biasLinsXO,
                        paramModelMap,
                        MsatBMA.BIAS_LIN_INDEX
                );

        XMLObject geosXO = xo.getChild(GEOS);
        ArrayList<Parameter> geos =
                processParameters(
                        geosXO,
                        paramModelMap,
                        MsatBMA.GEO_INDEX
                );

        XMLObject phaseProbXO = xo.getChild(PHASE_PROBS);
        ArrayList<Parameter> phaseProbs =
                processParameters(
                        phaseProbXO,
                        paramModelMap,
                        MsatBMA.PHASE_PROB_INDEX
                );

        Parameter modelChoose = (Parameter) xo.getElementFirstChild(MODEL_CHOOSE);
        Parameter modelIndicator = (Parameter) xo.getElementFirstChild(MODEL_INDICATOR);

        printParameters(paramModelMap);
        return new MsatBMA(
                dataType,
                logit,
                rateProps,
                rateQuads,
                biasConsts,
                biasLins,
                geos,
                phaseProbs,
                paramModelMap,
                modelChoose,
                modelIndicator,
                modelBitIndMap
        );
    }

    public ArrayList<Parameter> processParameters(
            XMLObject paramsXO,
            Parameter[][] paramModelMap,
            int paramIndex)throws XMLParseException{

        ArrayList<Parameter> paramList = new ArrayList<Parameter>();
        int paramsCount = paramsXO.getChildCount();

        for(int i = 0; i < paramsCount; i++){

            XMLObject paramXO = (XMLObject) paramsXO.getChild(i);
            int[] inModels = paramXO.getIntegerArrayAttribute(IN_MODELS);
            Parameter param = (Parameter)paramXO.getChild(Parameter.class);

            for(int j = 0; j < inModels.length; j++){
                if(paramModelMap[paramIndex][inModels[j]] == null){
                    paramModelMap[paramIndex][inModels[j]] = param;
                }else{
                    throw new RuntimeException("Different objects cannot be assigned to the same parameter in a model");
                }
            }

            paramList.add(param);

        }
        return paramList;
    }

    public void printParameters(Parameter[][] paramModelMap){
        for(int i = 0; i < paramModelMap.length; i++){
            for(int j = 0; j < paramModelMap[i].length; j++){
                System.out.print(paramModelMap[i][j]+" ");
            }
            System.out.println();
        }
    }




    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            new ElementRule(Microsatellite.class),
            AttributeRule.newBooleanRule(LOGIT,true),
            new ElementRule(
                    MODELS,
                    new XMLSyntaxRule[]{
                            new ElementRule(
                                    MODEL,
                                    new XMLSyntaxRule[]{
                                            AttributeRule.newStringRule(BINARY),
                                            AttributeRule.newStringRule(CODE)
                                    },
                                    1,
                                    27
                            )
                    }
            ),
            new ElementRule(
                    RATE_PROPS,
                    new XMLSyntaxRule[]{
                        new ElementRule(
                                RATE_PROP,
                                new XMLSyntaxRule[]{
                                        AttributeRule.newIntegerArrayRule(IN_MODELS, false),
                                        new ElementRule(Parameter.class)
                                },
                                1,
                                18
                        )
                    }
            ),
            new ElementRule(
                    RATE_QUADS,
                    new XMLSyntaxRule[]{
                        new ElementRule(
                                RATE_QUAD,
                                new XMLSyntaxRule[]{
                                        AttributeRule.newIntegerArrayRule(IN_MODELS, false),
                                        new ElementRule(Parameter.class)
                                },
                                1,
                                9
                        )
                    }
            ),
            new ElementRule(
                    BIAS_CONSTS,
                    new XMLSyntaxRule[]{
                        new ElementRule(
                                BIAS_CONST,
                                new XMLSyntaxRule[]{
                                        AttributeRule.newIntegerArrayRule(IN_MODELS, false),
                                        new ElementRule(Parameter.class)
                                },
                                1,
                                18
                        )
                    }
            ),
            new ElementRule(
                    BIAS_LINS,
                    new XMLSyntaxRule[]{
                        new ElementRule(
                                BIAS_LIN,
                                new XMLSyntaxRule[]{
                                        AttributeRule.newIntegerArrayRule(IN_MODELS, false),
                                        new ElementRule(Parameter.class)
                                },
                                1,
                                9
                        )
                    }
            ),
            new ElementRule(
                    GEOS,
                    new XMLSyntaxRule[]{
                        new ElementRule(
                                GEO,
                                new XMLSyntaxRule[]{
                                        AttributeRule.newIntegerArrayRule(IN_MODELS, false),
                                        new ElementRule(Parameter.class)
                                },
                                1,
                                18
                        )
                    }
            ),
            new ElementRule(
                    PHASE_PROBS,
                    new XMLSyntaxRule[]{
                        new ElementRule(
                                PHASE_PROB,
                                new XMLSyntaxRule[]{
                                        AttributeRule.newIntegerArrayRule(IN_MODELS, false),
                                        new ElementRule(Parameter.class)
                                },
                                1,
                                9
                        )
                    }
            ),
            new ElementRule(MODEL_CHOOSE, new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
            new ElementRule(MODEL_INDICATOR, new XMLSyntaxRule[]{new ElementRule(Parameter.class)})

    };

    public String getParserDescription() {
        return "This element represents an instance of the Microsatellite Averaging Model of microsatellite evolution.";
    }

    public Class getReturnType(){
        return MsatBMA.class;
    }
}
