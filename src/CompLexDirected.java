import lexvoltsdirected.*;

import java.io.File;
import java.io.FileWriter;
import java.util.Scanner;


public class CompLexDirected {

    public static void main(String args[]) throws java.io.IOException {

        if (args.length < 3) {
            System.err.print("usage: CompLexDirected graph conditions outputFile\n" +
            "  graph is the name of the file containing the input directed graph\n" +
                    "  conditions is a list of the boundary conditions (fixed voltages)\n" +
                    "  outputFile is the name of the file to which output should be written\n" +
                    "\n" +
                    "Graph format: the first line should be numVertices numEdges\n" +
                    "  each line after that should list one edge, as destination source length\n" +
                    "  The lowest numbered vertex should be 0.\n" +
                    "\n" +
                    "Conditions format: the first line should give the number of boundary nodes\n" +
                    "  this should be followed by lines in the format vertexNumber value\n" +
                    "\n" +
                    "To generate examples of these files in Matlab, see the routine lexDirectedFromFiles.m\n" +
                    "For more info, see this URL\n");
            return;
        }

        String graphName = args[0];
        String conditionsName = args[1];
        String outfile = args[2];

        // Read the graph

        Scanner sc = new Scanner(new File(graphName));
        int nv = sc.nextInt();
        int ne = sc.nextInt();


        // these are what we read from the file.
        // we will sort them into appropriate order later
        int[] ai = new int[ne];
        int[] aj = new int[ne];
        double[] av = new double[ne];


        for (int i = 0; i < ne; i++) {
            int dest = sc.nextInt();
            int src = sc.nextInt();
            double weight = sc.nextDouble();
            ai[i] = dest;
            aj[i] = src;
            av[i] = weight;
        }


        // Read the boundary conditions

        sc = new Scanner(new File(conditionsName));
        int bdrySize = sc.nextInt();
        int[] where = new int[bdrySize];
        double[] what = new double[bdrySize];
        for (int i = 0; i < bdrySize; i++) {
            where[i] = sc.nextInt();
            what[i] = sc.nextDouble();
        }


        RandLexVoltDirected rlv = new RandLexVoltDirected(ai,aj,av, System.currentTimeMillis());
        rlv.setTerminals(where,what);
        rlv.computeVoltages();

        FileWriter fw = new FileWriter(outfile);
        for (int i = 0; i < rlv.nv; i++)
            fw.write(rlv.vertices[i].voltage + "\n");
        fw.close();




    }

}