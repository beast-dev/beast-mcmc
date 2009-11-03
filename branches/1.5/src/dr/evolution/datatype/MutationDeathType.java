package dr.evolution.datatype;

import dr.evoxml.DataTypeUtils;
import dr.xml.*;

import java.util.logging.Logger;

/**
 * Package: MutationDeathType
 * Description:
 * <p/>
 * <p/>
 * Created by
 * Alexander V. Alekseyenko (alexander.alekseyenko@gmail.com)
 * Date: Mar 5, 2008
 * Time: 1:09:40 PM
 */
public class MutationDeathType extends DataType {
    protected static String DESCRIPTION = "MutationDeathType";
    protected static int UNKNOWN_STATE = 128;
    protected int[][] codes;
    protected char[] stateCodes;
    public int DEATHSTATE;

    public MutationDeathType(char deathCode, char extantCode) { // Constructor for pure death type
        super();
        initialize_internals();
        codes[extantCode] = new int[]{0};
        stateCodes[0] = extantCode;

        codes[deathCode] = new int[]{1};
        stateCodes[1] = deathCode;
        stateCount = 2;
        DEATHSTATE = 1;
        ambiguousStateCount = 0;
    }

    public MutationDeathType(DataType x, char deathCode) { // constructor for extention type
        super();
        int i;
        char stateCode;

        initialize_internals();
        for (i = 0; i < x.getStateCount(); ++i) {  /* Copy unique codes */
            stateCode = x.getCode(i).charAt(0);
            this.codes[stateCode] = new int[]{i};
            stateCodes[i] = stateCode;
        }
        this.codes[deathCode] = new int[]{i}; /* Append the state space with the death state */
        stateCodes[i] = deathCode;
        DEATHSTATE = i;
        stateCount = i + 1;

        for (i = 0; i < 128; ++i) {
            int state = x.getState((char) i);
            if (state > 0 && state < 128 && x.isAmbiguousState(state) && i != deathCode) {
                if (!x.isUnknownState(state)) {
                    int[] states = x.getStates(state);
                    this.codes[i] = new int[states.length];
                    System.arraycopy(states, 0, this.codes[i], 0, states.length);
                }
            }

        }
        ambiguousStateCount = x.getAmbiguousStateCount() + 1;
    }

    public void addAmbiguity(char ambiguityCode, String s) {
        if (s.length() == 0) {
            this.codes[ambiguityCode] = new int[stateCount];
            for (int i = 0; i < stateCount; ++i) {
                this.codes[ambiguityCode][i] = i;
            }
        } else {
            this.codes[ambiguityCode] = new int[s.length()];
            for (int i = 0; i < s.length(); ++i) {
                this.codes[ambiguityCode][i] = getState(s.charAt(i));
            }
        }
        ambiguousStateCount += 1;
    }

    private void initialize_internals() {
        this.codes = new int[128][]; /* stores states (w/ ambiguities) corresponding to codes */
        this.stateCodes = new char[128]; /* Stores characters corresponding to unique state codes*/
    }

    /**
     * Get state corresponding to a character
     *
     * @param c character
     * @return state
     */
    public int getState(char c) {
        if (codes[c] != null && codes[c].length == 1) {
            return codes[c][0];
        } else {
            return c;
        }
    }

    /**
     * Get state corresponding to an unknown
     *
     * @return state
     */
    public int getUnknownState() {
        return UNKNOWN_STATE;
    }


    /**
     * Get character corresponding to a given state
     *
     * @param state state
     *              <p/>
     *              return corresponding character
     */
    public char getChar(int state) {
        if (state < stateCount)
            return stateCodes[state];
        return super.getChar(state);
    }

    /**
     * returns an array containing the non-ambiguous states that this state represents.
     */
    public int[] getStates(int state) {
        if (state < stateCount)
            return codes[stateCodes[state]];
        else
            return codes[state];
    }

    /**
     * returns anarray of indicatiors for non-ambiguous states that this state represents.
     */
    public boolean[] getStateSet(int state) {
        boolean[] stateSet = new boolean[stateCount];
        int states[];
        int i;
        for (i = 0; i < stateCount; ++i)
            stateSet[i] = false;
        states = getStates(state);

        for (i = 0; states != null && i < states.length; ++i) {
            stateSet[states[i]] = true;
        }
        return stateSet;
    }

    /**
     * @return true if this character is an ambiguous state
     */
    public boolean isAmbiguousChar(char c) {
        return codes[c] != null && codes[c].length > 1;
    }

    /**
     * @return true if this character is a gap
     */
    public boolean isUnknownChar(char c) {
        return codes[c] == null;
    }

    /**
     * returns true if this state is an ambiguous state.
     */
    public boolean isAmbiguousState(int state) {
        return state >= stateCount;
    }

    /**
     * @return true if this state is an unknown state
     */
    public boolean isUnknownState(int state) {
        return state >= stateCount && codes[state] == null;
    }

    public String getDescription() {
        return DESCRIPTION;
    }

    public int getType() {
        return 314;
    }

    public static final String MODEL_NAME = "extendedDataType";
    public static final String STATE = "deathState";
    public static final String CODE = "code";
    public static final String STATES = "states";
    public static final String AMBIGUITY = "ambiguity";
    public static final String EXTANT = "extantState";

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return MODEL_NAME;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {
            Logger.getLogger("dr.evolution").info("\nCreating an extended data type.");
            Logger.getLogger("dr.evolution").info("\tIf you publish results using this model, please reference Alekseyenko, Lee and Suchard (in submision).\n");

            DataType dataType = DataTypeUtils.getDataType(xo);

            char extantChar = '\0';
            XMLObject cxo = xo.getChild(EXTANT);
            if (cxo != null) {
                extantChar = cxo.getStringAttribute(CODE).charAt(0);
            }


            cxo = xo.getChild(STATE);
            char stateChar;
            if (cxo != null) {
                stateChar = cxo.getStringAttribute(CODE).charAt(0);
            } else {
                stateChar = dataType.getChar(dataType.getGapState());
            }

            Logger.getLogger("dr.evolution").info("\tNon-existent code: " + stateChar);

            if (dataType == null && extantChar == '\0')
                throw new XMLParseException("In " + xo.getName() + " you must either provide a data type or a code for extant state");

            MutationDeathType mdt;
            if (dataType != null) {
                Logger.getLogger("dr.evolution").info("\tBase type: " + dataType.toString());
                mdt = new MutationDeathType(dataType, stateChar);
            } else {
                Logger.getLogger("dr.evolution").info("\tExtant code: " + extantChar);
                mdt = new MutationDeathType(stateChar, extantChar);
            }

            char ambiguityChar = '\0';
            String states;
            Object dxo;
            for (int i = 0; i < xo.getChildCount(); ++i) {
                dxo = xo.getChild(i);
                if (dxo instanceof XMLObject) {
                    cxo = (XMLObject) dxo;
                    if (cxo.getName().equals(AMBIGUITY)) {
                        ambiguityChar = cxo.getStringAttribute(CODE).charAt(0);
                        if (cxo.hasAttribute(STATES)) {
                            states = cxo.getStringAttribute(STATES);
                        } else {
                            states = "";
                        }
                        mdt.addAmbiguity(ambiguityChar, states);
                        Logger.getLogger("dr.evolution").info("\tAmbiguity code: " + ambiguityChar);
                    }
                }
            }
            return mdt;
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return "This element represents an instance of the MutationDeathType which extends a base datatype with an additional \"death\" state.";
        }

        public Class getReturnType() {
            return MutationDeathType.class;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private final XMLSyntaxRule[] rules = {
                new StringAttributeRule(DataType.DATA_TYPE, "Base datatype name", DataType.getRegisteredDataTypeNames(), true),
                new ElementRule(DataType.class, true),
                new ElementRule(STATE, new XMLSyntaxRule[]{AttributeRule.newStringRule(CODE)}, true),
                new ElementRule(EXTANT, new XMLSyntaxRule[]{AttributeRule.newStringRule(CODE)}, true),
                new ContentRule("<ambiguity code=\"Z\" states=\"XY\"/>")
        };

    };
}
