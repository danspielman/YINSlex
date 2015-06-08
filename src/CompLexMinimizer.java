import lexvolts.EdgeListGraph;
import lexvolts.RandLexVolt;

import java.io.File;
import java.io.FileWriter;
import java.util.Scanner;


public class CompLexMinimizer {

    public static void main(String args[]) throws java.io.IOException {

        if (args.length < 3) {
            System.err.print("usage: CompLexMinimizer graph conditions outputFile\n" +
            "  graph is the name of the file containing the input graph\n" +
                    "  conditions is a list of the boundary conditions (fixed voltages)\n" +
                    "  outputFile is the name of the file to which output should be written\n" +
                    "\n" +
                    "Graph format: the first line should be numVertices numEdges\n" +
                    "  each line after that should list one edge, as source destination length\n" +
                    "  The graph is treated as undirected, and each edge should be listed just once\n" +
                    "  The lowest numbered vertex should be 0.\n" +
                    "\n" +
                    "Conditions format: the first line should give the number of boundary nodes\n" +
                    "  this should be followed by lines in the format vertexNumber value\n" +
                    "\n" +
                    "To generate examples of these files in Matlab, see the route lexFromFiles.m\n" +
                    "For more info, see this URL\n");
            return;
        }

        String graphName = args[0];
        String conditionsName = args[1];
        String outfile = args[2];

        EdgeListGraph g = new EdgeListGraph(graphName);

        Scanner sc = new Scanner(new File(conditionsName));
        int bdrySize = sc.nextInt();
        int[] where = new int[bdrySize];
        double[] what = new double[bdrySize];
        for (int i = 0; i < bdrySize; i++) {
            where[i] = sc.nextInt();
            what[i] = sc.nextDouble();
        }


        RandLexVolt rlv = new RandLexVolt(g);
        rlv.setTerminals(where,what);
        rlv.computeVoltages5();

        FileWriter fw = new FileWriter(outfile);
        for (int i = 0; i < g.nv; i++)
            fw.write(rlv.voltages[i] + "\n");
        fw.close();




    }

}