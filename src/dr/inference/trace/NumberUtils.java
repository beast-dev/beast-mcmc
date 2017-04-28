package dr.inference.trace;

/**
 * Detect number type from a string.
 * <code>isNumber</code> is modified from org.apache.commons.lang3.math.NumberUtils
 * to handle NaN and Inf.
 *
 * @author Walter Xie
 */
public class NumberUtils {

    /**
     * <p>Checks whether the <code>String</code> does not contain decimal point.</p>
     * <p>
     * <p><code>Null</code> and empty String will return <code>false</code>.</p>
     * <p><code>Double.NaN</code>, <code>Double.NEGATIVE_INFINITY</code>,
     * and <code>Double.POSITIVE_INFINITY</code> will return <code>true</code>.</p>
     *
     * @param str the <code>String</code> to check
     * @return <code>true</code> if str contains only Unicode numeric
     */
    public static boolean hasDecimalPoint(String str) {
        if (str == null || str.length() == 0)
            return false;
        if (str.equalsIgnoreCase(Double.toString(Double.NaN)) ||
                str.equalsIgnoreCase(Double.toString(Double.NEGATIVE_INFINITY)) ||
                str.equalsIgnoreCase(Double.toString(Double.POSITIVE_INFINITY)))
            return true;
        return str.contains(".");
    }

    /**
     * <p>Checks whether the String a valid Java number.</p>
     * <p>
     * <p>Valid numbers include hexadecimal marked with the <code>0x</code>
     * qualifier, scientific notation and numbers marked with a type
     * qualifier (e.g. 123L).</p>
     * <p>
     * <p><code>Null</code> and empty String will return <code>false</code>.</p>
     * <p><code>Double.NaN</code>, <code>Double.NEGATIVE_INFINITY</code>,
     * and <code>Double.POSITIVE_INFINITY</code> will return <code>true</code>.</p>
     *
     * @param str the <code>String</code> to check
     * @return <code>true</code> if the string is a correctly formatted number
     */
    public static boolean isNumber(String str) {
        if (str == null || str.length() == 0) {
            return false;
        }
// modified
        if (str.equalsIgnoreCase(Double.toString(Double.NaN)) ||
                str.equalsIgnoreCase(Double.toString(Double.NEGATIVE_INFINITY)) ||
                str.equalsIgnoreCase(Double.toString(Double.POSITIVE_INFINITY)))
            return true;

        char[] chars = str.toCharArray();
        int sz = chars.length;
        boolean hasExp = false;
        boolean hasDecPoint = false;
        boolean allowSigns = false;
        boolean foundDigit = false;
// deal with any possible sign up front
        int start = (chars[0] == '-') ? 1 : 0;
        if (sz > start + 1 && chars[start] == '0' && chars[start + 1] == 'x') {
            int i = start + 2;
            if (i == sz) {
                return false; // str == "0x"
            }
// checking hex (it can't be anything else)
            for (; i < chars.length; i++) {
                if ((chars[i] < '0' || chars[i] > '9')
                        && (chars[i] < 'a' || chars[i] > 'f')
                        && (chars[i] < 'A' || chars[i] > 'F')) {
                    return false;
                }
            }
            return true;
        }
        sz--; // don't want to loop to the last char, check it afterwords
// for type qualifiers
        int i = start;
// loop to the next to last char or to the last char if we need another digit to
// make a valid number (e.g. chars[0..5] = "1234E")
        while (i < sz || (i < sz + 1 && allowSigns && !foundDigit)) {
            if (chars[i] >= '0' && chars[i] <= '9') {
                foundDigit = true;
                allowSigns = false;
            } else if (chars[i] == '.') {
                if (hasDecPoint || hasExp) {
// two decimal points or dec in exponent
                    return false;
                }
                hasDecPoint = true;
            } else if (chars[i] == 'e' || chars[i] == 'E') {
// we've already taken care of hex.
                if (hasExp) {
// two E's
                    return false;
                }
                if (!foundDigit) {
                    return false;
                }
                hasExp = true;
                allowSigns = true;
            } else if (chars[i] == '+' || chars[i] == '-') {
                if (!allowSigns) {
                    return false;
                }
                allowSigns = false;
                foundDigit = false; // we need a digit after the E
            } else {
                return false;
            }
            i++;
        }
        if (i < chars.length) {
            if (chars[i] >= '0' && chars[i] <= '9') {
// no type qualifier, OK
                return true;
            }
            if (chars[i] == 'e' || chars[i] == 'E') {
// can't have an E at the last byte
                return false;
            }
            if (chars[i] == '.') {
                if (hasDecPoint || hasExp) {
// two decimal points or dec in exponent
                    return false;
                }
// single trailing decimal point after non-exponent is ok
                return foundDigit;
            }
            if (!allowSigns
                    && (chars[i] == 'd'
                    || chars[i] == 'D'
                    || chars[i] == 'f'
                    || chars[i] == 'F')) {
                return foundDigit;
            }
            if (chars[i] == 'l'
                    || chars[i] == 'L') {
// not allowing L with an exponent or decimal point
                return foundDigit && !hasExp && !hasDecPoint;
            }
// last character is illegal
            return false;
        }
// allowSigns is true iff the val ends in 'E'
// found digit it to make sure weird stuff like '.' and '1E-' doesn't pass
        return !allowSigns && foundDigit;
    }

}
