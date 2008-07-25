package dr.app.beauti.options;

/**
 * @author Alexei Drummond
 */
public enum OperatorType {

    SCALE("scale"),
    RANDOM_WALK("randomWalk"),
    INTEGER_RANDOM_WALK("integerRandomWalk"),
    UP_DOWN("upDown"),
    SCALE_ALL("scaleAll"),
    CENTERED_SCALE("centeredScale"),
    DELTA_EXCHANGE("deltaExchange"),
    INTEGER_DELTA_EXCHANGE("integerDeltaExchange"),
    SWAP("swap"),
    BITFLIP("bitFlip"),
    TREE_BIT_MOVE("treeBitMove"),
    SAMPLE_NONACTIVE("sampleNoneActiveOperator"),
    SCALE_WITH_INDICATORS("scaleWithIndicators"),
    UNIFORM("uniform"),
    INTEGER_UNIFORM("integerUniform"),
    SUBTREE_SLIDE("subtreeSlide"),
    NARROW_EXCHANGE("narrowExchange"),
    WIDE_EXCHANGE("wideExchange"),
    WILSON_BALDING("wilsonBalding");

    OperatorType(String name) {
        this.name = name;
    }

    public String getName() {
        return name;
    }

    private String name;
}
