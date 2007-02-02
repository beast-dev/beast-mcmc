package dr.app.beauti;

/**
 * Created by IntelliJ IDEA.
 * User: msuchard
 * Date: Jan 31, 2007
 * Time: 3:56:10 PM
 * To change this template use File | Settings | File Templates.
 */
public class BeautiARG {


    static public void main(String[] args) {


        if (args.length < 3) {
            System.err.println("Usage: beautiARG <template_file> <output_file> <alignment_file1> <alignment_file2> ...");
            return;
        }

        //String inputFileName = args[0];
        String templateFileName = args[0];
        String outputFileName = args[1];
        String[] inputFileNames = new String[args.length - 2];
        for (int i = 0; i < args.length - 2; i++)
            inputFileNames[i] = args[i + 2];

        CommandLineBeauti beauti = new CommandLineBeauti(inputFileNames, templateFileName, outputFileName);
    }


}
