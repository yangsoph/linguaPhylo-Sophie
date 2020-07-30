package lphy2beast.tobeast.generators;

import beast.core.BEASTInterface;
import beast.core.parameter.RealParameter;
import beast.math.distributions.LogNormalDistributionModel;
import lphy.core.distributions.LogNormal;
import lphy.core.distributions.Normal;
import lphy2beast.BEASTContext;
import lphy2beast.GeneratorToBEAST;

public class NormalToBEAST implements GeneratorToBEAST<Normal> {
    @Override
    public BEASTInterface generatorToBEAST(Normal generator, BEASTInterface value, BEASTContext context) {
        beast.math.distributions.Normal normal = new beast.math.distributions.Normal();
        normal.setInputValue("mean", context.getBEASTObject(generator.getMean()));
        normal.setInputValue("sigma", context.getBEASTObject(generator.getSd()));
        normal.initAndValidate();

        return BEASTContext.createPrior(normal, (RealParameter)value);
    }

    @Override
    public Class<Normal> getGeneratorClass() {
        return Normal.class;
    }
}
