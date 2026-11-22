#include <fstream>
#include <sstream>
#include <map>

#include "../../../src/IndelSimulator.h"
#include "../../../src/SubstitutionSimulator.h"
#include "../../../libs/sfc/sfc64.h"


static void simulate(bool toDisk, const std::string& outPath, size_t& msaLength, size_t& numLeaves) {
    tree tree_("../../trees/normalbranches_nLeaves10000.treefile");
    SimulationContext<SFC64> simContext(&tree_, 1);
    simContext.setCacheBranchProbs(true);
    simContext.setSaveAll();
    simContext.setOutputLeavesOnly(true);

    size_t numBranches = tree_.getNodesNum() - 1;
    // Zipf distribution with 'a' parameter of 1.7 (truncated at 50 for a short test)
    std::vector<double> zipf(50); double sum = 0;
    for (int i = 0; i < 50; i++) { zipf[i] = std::pow(i + 1, -1.7); sum += zipf[i]; }
    for (auto& z : zipf) z /= sum;
    DiscreteDistribution d1(zipf);
    std::vector<DiscreteDistribution*> insertionDists(numBranches, &d1);
    std::vector<DiscreteDistribution*> deletionDists(numBranches, &d1);
    std::vector<double> insertionRates(numBranches, 0.03);
    std::vector<double> deletionRates(numBranches, 0.09);

    SimulationProtocol protocol(numBranches);
    simContext.setProtocol(&protocol);
    protocol.setInsertionLengthDistributions(insertionDists);
    protocol.setDeletionLengthDistributions(deletionDists);
    protocol.setInsertionRates(insertionRates);
    protocol.setDeletionRates(deletionRates);
    protocol.setMaxInsertionLength(20);
    protocol.setMinSequenceSize(10);
    protocol.setSiteRateModel(SiteRateModel::SIMPLE);
    protocol.setSequenceSize(30000);

    IndelSimulator<SFC64> indelSim(simContext, &protocol);
    auto eventMap = indelSim.generateSimulation();
    std::cout << "simulated indel events\n";

    MSA<SFC64> msa(eventMap, simContext);
    msaLength = msa.getMSAlength();
    numLeaves = msa.getNumberOfSequences();
    std::cout << "output MSA: " << numLeaves << " rows, length " << msaLength
              << " (simulation width " << msa.getSimulationLength() << ")\n";

    modelFactory mFac;
    mFac.setSiteRateModel({0.01, 0.5, 2.0, 4.0}, {0.25, 0.25, 0.25, 0.25});
    mFac.buildModel(alphabetCode::AMINOACID, modelCode::WAG);

    SubstitutionSimulator<SFC64, 20> substitutionSim(mFac, simContext);
    substitutionSim.setAlignedSequenceMap(msa);

    // NOTE: the substitution simulator needs the simulation width, not the (trimmed) output length
    if (toDisk) {
        substitutionSim.simulateAndWriteSubstitutions(msa.getSimulationLength(), outPath, "", msa.getRootPositionsInMsa());
    } else {
        auto fullContainer = substitutionSim.simulateSubstitutions(msa.getSimulationLength(), "", msa.getRootPositionsInMsa());
        msa.fillSubstitutions(fullContainer);
        // msa.writeFullMsa(outPath.c_str());
    }
}

int main() {
    size_t len1, rows1, len2, rows2;
    // simulate(true,  "leafonly_disk.fasta", len1, rows1);
    simulate(false, "leafonly_mem.fasta",  len2, rows2);

    return 0;
}
