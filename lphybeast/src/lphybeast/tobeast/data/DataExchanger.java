package lphybeast.tobeast.data;

import beast.evolution.alignment.Alignment;
import lphy.core.LPhyParser;
import lphy.evolution.DataFrame;
import lphy.evolution.tree.TimeTree;
import lphy.evolution.tree.TimeTreeNode;
import lphy.graphicalModel.Value;

import java.io.PrintStream;
import java.nio.file.Path;
import java.util.List;
import java.util.Map;

/**
 * @author Walter Xie
 */
@Deprecated
public class DataExchanger {

    protected NexusParser nexusParser;

    // if null, then 1 partition
    protected Map<String, String> partitionMap; // LPhy <=> Nex

    // parse partition(s) from CMD --nex
    // create Path
    public DataExchanger(Map<String, String> partitionMap) {
        Path[] nexfile = getNexFiles(partitionMap);

//        nexusParser = new NexusParser(nexfile);
    }

    private Path[] getNexFiles(Map<String, String> partitionMap) {
        return new Path[0];
    }

    // Single file, if partitionMap null, then 1 partition
    public DataExchanger(Path nexfile, Map<String, String> partitionMap) {
        this.partitionMap = partitionMap;
        nexusParser = new NexusParser(nexfile);
    }

    public Alignment getBEASTAlignment(String partitionName) {
        if (partitionMap != null) {
            //TODO
        }
        return nexusParser.getAlignment();
    }
    @Deprecated
    public Alignment getAlignment() {
        return nexusParser.getAlignment();
    }
//    public Alignment getAlignment(String algID) {
//        return nexusParser.getAlignment(algID); //TODO
//    }

    //****** before LPhyParser ******//

//    public void preloadArgs() {
//        // run ArgI.putArgument before LPhy
//        for (Map.Entry<String, String> entry : varMap.entrySet()) {
//            String lphyVar = entry.getKey();
//            String nexusVar = entry.getValue();
//            String nexusValue = nexusParser.getVal(nexusVar);
//
//            ArgI.putArgument(lphyVar, Integer.parseInt(nexusValue));
//        }
//    }


    // replace dataframe to nexus
    public void updateDF(LPhyParser parser) {

        for (Map.Entry<String, Value<?>> entry : parser.getDictionary().entrySet()) {

            String lphyVar = entry.getKey();
            DataFrame df = (DataFrame) entry.getValue().value();



            if (partitionMap.containsKey(lphyVar)) {
                String partName = partitionMap.get(lphyVar);

//                if (df.get.hasParts()) {
//
//                }



                Alignment beastAlign = getBEASTAlignment(partName);
                int ntax = beastAlign.getTaxonCount();
                int nchar = beastAlign.getSiteCount();


            }

        }



//        Matcher matcher = pattern.matcher(line);
//        while (matcher.find()) {
//            // Get the group matched
//            String regx = matcher.group();
//            String replacement = args.get(regx);
//            line = line.replaceAll(regx, replacement);
//
//            System.out.println("Find argument " + matcher.group() + " in LPhy and replace to " + replacement);
//        }
//        return line;

    }



    //****** after LPhyParser ******//

    // replace taxa names in TimeTree TODO how to map to alignment
    public void replaceTaxaNamesByOrder(TimeTree timeTree) {
        List<String> beastTaxaNm = nexusParser.getAlignment().getTaxaNames();
        //TODO allow taxa diff ?
        if (beastTaxaNm.size() != timeTree.getTaxaNames().length)
            throw new IllegalArgumentException("The given taxa have to match the taxa in the LPhy model !\n"
                    + beastTaxaNm.size() + " != " + timeTree.getTaxaNames().length);

        int i= 0;
        for (TimeTreeNode node : timeTree.getNodes()) {
            if (node.isLeaf()) {
                node.setId(beastTaxaNm.get(i));
                i++;
            }
        }
    }

    public void printPartMap(PrintStream out){
        if (partitionMap != null) {
            out.println("Map LPhy script partitions to Nexus data : ");
            for (Map.Entry<String, String> entry : partitionMap.entrySet()) {
                out.println(entry.getKey() + " => " + entry.getValue());
            }
            out.println();
        }
    }

    public boolean containsDF(String line) {
        return line.trim().toLowerCase().contains("dataframe(");
    }

// TODO map taxa names


    //*** in dev ***//

    //    public Alignment getAlignment() {
//        Sequence human = new Sequence("0human", "AAAACCCCGGGGTTTT");
//        Sequence chimp = new Sequence("1chimp", "ACGTACGTACGTACGT");
//        Sequence gorilla = new Sequence("2gorilla", "ATGTACGTACGTACTT");
//
//        Alignment data = new Alignment();
//        data.initByName("sequence", human, "sequence", chimp, "sequence", gorilla,
//                "dataType", "nucleotide"
//        );
//        return data;
//    }

//    public void replaceAlignment(final Set<BEASTInterface> elements, final SortedMap<String, Taxon> allTaxa) {
//        Alignment alg = null;
//        for (BEASTInterface bo : elements) {
//            if (bo instanceof Alignment) {
//                Alignment newAlg = getAlignment();
//                List<Sequence> newSeqs = newAlg.sequenceInput.get();
//                alg = (Alignment) bo;
//
//                // validate SiteCount
//                if (newAlg.getSiteCount() != alg.getSiteCount()) {
//                    System.err.println("Warning: the number of sites ("+ alg.getSiteCount() +
//                            ") defined in LPhy != " + " sites (" + newAlg.getSiteCount() +
//                            ") in the given alignment !");
//                }
//
//                alg.sequenceInput.get().clear();
//                // this will add list, not replace
//                alg.sequenceInput.setValue(newSeqs, alg);
//                alg.initAndValidate();
//
//            }
//        }
//
//        assert alg != null;
//
//        replaceAllTaxaBy(allTaxa, alg); // TODO cannot access taxonset in Tree
//
//    }
//
//    private void replaceAllTaxaBy(final SortedMap<String, Taxon> allTaxa, final Alignment alg) {
//        assert allTaxa.size() == alg.getTaxonCount();
//
//        List<String> tn = alg.getTaxaNames();
//        allTaxa.clear();//TODO
//        for (String taxonID : tn) {
//            if (!allTaxa.containsKey(taxonID)) {
//                allTaxa.put(taxonID, new Taxon(taxonID));
//            }
//        }
//
//        for (Map.Entry<String, Taxon> pair : allTaxa.entrySet()) {
////             pair.getKey()  pair.getValue();
//        }
//
//    }
}
