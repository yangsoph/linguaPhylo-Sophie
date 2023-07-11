package lphy.core.simulator;

import lphy.core.io.FileConfig;
import lphy.core.logger.ValueLoggerListener;
import lphy.core.model.Value;

import java.io.File;
import java.io.IOException;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class Simulator {

    SimulatorListener simulatorListener;

    Sampler sampler;

    public Simulator() {
        simulatorListener = new ValueLoggerListener();
    }

    public Map<Integer, List<Value>> simulateAndSaveResults(File lphyFile, int numReplicates, Long seed) throws IOException {
        simulatorListener.start(List.of(numReplicates, lphyFile));
        return simulate(lphyFile, numReplicates, seed);
    }

    public Map<Integer, List<Value>> simulateAndSaveResults(FileConfig fileConfig) throws IOException {
        simulatorListener.start(List.of(fileConfig));

        File lphyFile = fileConfig.lphyFile;
        int numReplicates = fileConfig.numReplicates;
        Long seed = fileConfig.seed;

        return simulate(lphyFile, numReplicates, seed);
    }

    private Map<Integer, List<Value>> simulate(File lphyFile, int numReplicates, Long seed) throws IOException {
        // TODO duplicate to maps in ValueLoggerListener
        Map<Integer, List<Value>> simResMap = new HashMap<>();

        // create Sampler given a lphy script file
        sampler = Sampler.createSampler(lphyFile);

        long start = System.currentTimeMillis();

        for (int i = SimulatorListener.REPLICATES_START_INDEX; i < numReplicates; i++) {
            List<Value> values = sampler.sample(seed);

            simulatorListener.replicate(i, values);

            simResMap.put(i, values);
        }
        simulatorListener.complete();

        long end = System.currentTimeMillis();
        System.out.println("Sampled " + lphyFile + " at " + numReplicates + (numReplicates >1?" times":" time") +
                " which takes " + (end - start) + " ms.");

        return simResMap;
    }


}