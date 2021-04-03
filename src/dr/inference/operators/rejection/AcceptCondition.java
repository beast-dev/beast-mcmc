package dr.inference.operators.rejection;

public interface AcceptCondition {

    boolean satisfiesCondition(double[] values);

    enum SimpleAcceptCondition implements AcceptCondition {
        DescendingAbsoluteValue("descendingAbsoluteValue") {
            @Override
            public boolean satisfiesCondition(double[] values) {
                for (int i = 1; i < values.length; i++) {
                    if (Math.abs(values[i - 1]) < Math.abs(values[i])) {
                        return false;
                    }
                }
                return true;
            }
        },

//        DescendingAbsoluteValueSpaced("descendingAbsoluteValueSpaced") {
//            @Override
//            public boolean satisfiesCondition(double[] values) {
//                for (int i = 1; i < values.length; i++) {
//                    if (0.9 * Math.abs(values[i - 1]) < Math.abs(values[i])) {
//                        return false;
//                    }
//                }
//                return true;
//            }
//        },

        AlternatingSigns("descendingAlternatingSigns") {
            @Override
            public boolean satisfiesCondition(double[] values) {
                for (int i = 1; i < values.length; i++) {
                    Boolean signa = (values[i] > 0);
                    Boolean signb = (values[i - 1] > 0);
                    if (Math.abs(values[i - 1]) < Math.abs(values[i]) || signa == signb) {
                        return false;
                    }
                }
                return true;
            }
        };

        private final String name;

        SimpleAcceptCondition(String name) {
            this.name = name;
        }

        public String getName() {
            return name;
        }

        @Override
        public abstract boolean satisfiesCondition(double[] values);
    }

}


