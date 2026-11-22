#ifndef _MSA
#define _MSA

#include <vector>
#include <string>
#include <cmath>
#include <cstdlib>
#include <unordered_map>
#include <iostream>
#include <fstream>
#include <filesystem>

#include "../libs/Phylolib/includes/tree.h"

#include "SimulationContext.h"
#include "SimulationProtocol.h"
#include "Sequence.h"
#include "Event.h"



typedef std::vector<std::vector<int>> SparseMSA;

template<typename RngType = std::mt19937_64>
class MSA
{
public:
	// using iteratorType = typename std::list<typename SuperSequence<RngType, BlockTree>::columnContainer>::iterator;
	using iteratorType = typename plf::list<typename SuperSequence<RngType, BlockTree>::columnContainer>::iterator;


    MSA<RngType>(EventMap &eventmap, SimulationContext<RngType> &simContext): 
        _simContext(simContext),  _numberOfSequences(simContext.getNumberOfNodesToSave()), _msaLength(0) {

        tree::nodeP rootNode = simContext.getRootNode();
        const std::vector<bool>& nodesToSave = simContext.getNodesToSave();
        // Retrieve root sequence size from dummy insertion event in event map
        size_t sequenceSize = eventmap.at(rootNode->id())[0].length;

        // DISPATCH based on protocol setting
        auto* protocol = simContext.getProtocol();
        if (protocol->getSiteRateModel() == SiteRateModel::SIMPLE) {
            buildMsa<BlockTree>(eventmap, rootNode, nodesToSave, sequenceSize);
        } else {
            buildMsa<BlockTreeWithRates>(eventmap, rootNode, nodesToSave, sequenceSize);
        }

    }


    template<typename BlockTreeType>
    void buildMsa(EventMap &eventmap, tree::nodeP rootNode,
                  const std::vector<bool>& nodesToSave, size_t sequenceSize) {
        // Create SuperSequence with the appropriate BlockTree type
        SuperSequence<RngType, BlockTreeType> superSequence(sequenceSize, _simContext);
        
        const bool trackLeaves = _simContext.getOutputLeavesOnly();
        Sequence rootSequence(superSequence, nodesToSave[rootNode->id()], trackLeaves && rootNode->isLeaf());
        rootSequence.initSequence();

        std::vector<CompressedSequence> finalSequences;
        finalSequences.reserve(_numberOfSequences);

        buildMsaRecursively(finalSequences, eventmap, *rootNode, 
                           superSequence, rootSequence, nodesToSave, trackLeaves);
        eventmap.clear();
        
        prepareLazyRows(finalSequences, superSequence);
    }


    template<typename BlockTreeType>
    void buildMsaRecursively(std::vector<CompressedSequence> &finalSequences,
                             EventMap &eventmap, const tree::TreeNode &parrentNode,
                             SuperSequence<RngType, BlockTreeType> &superSequence,
                             const Sequence<RngType, BlockTreeType>& parentSequence, 
                             const std::vector<bool>& nodesToSave, bool trackLeaves) {
        if ((nodesToSave)[parrentNode.id()]) finalSequences.emplace_back(parentSequence.compress());
        if (parrentNode.isLeaf()) return;

        for (size_t i = 0; i < parrentNode.getNumberOfSons(); i++) {
            tree::TreeNode* childNode = parrentNode.getSon(i);
            Sequence currentSequence(superSequence, nodesToSave[childNode->id()], trackLeaves && childNode->isLeaf());

            auto events = eventmap.at(childNode->id());
            currentSequence.generateSequence(events, &parentSequence);
            buildMsaRecursively(finalSequences, eventmap, *childNode, superSequence, currentSequence, nodesToSave, trackLeaves);
        }
        
    }

    // PROTOTYPE: lazy rows. Keep compressed runs + flat position->column array; build a sparse row on demand.
    template<typename BlockTreeType>
    void prepareLazyRows(std::vector<CompressedSequence> &sequences, SuperSequence<RngType, BlockTreeType> &superSeq) {
        _numberOfSequences = superSeq.getNumSequences();
        _msaLength = superSeq.getMsaSequenceLength();
        _rootPositionsInMsa = superSeq.setAbsolutePositions();
        _perSiterateCategories = superSeq.getMsaRateCategories();
        _absByPos = superSeq.getAbsByPos();
        if (_simContext.getOutputLeavesOnly()) {
            _leavesOnlyOutput = true;
            _outByPos = superSeq.getOutByPos();
            _outLength = superSeq.getOutputSequenceLength();
            _outputRows = _simContext.getLeafRows();
        }
        _lazyRuns = std::move(sequences);
        _lazy = true;
    }

    // Build the sparse (run-length) row for `row` into `out` (reused buffer).
    // trimmed=false: simulation coordinates (full MSA width). trimmed=true: output coordinates (columns that are gap in all leaves removed; leaves-only output mode).
    void fillRow(size_t row, std::vector<int>& out, bool trimmed = false) const {
        if (!_lazy) { out = _alignedSequenceMap->at(row); return; }
        out.clear();
        const CompressedSequence& cs = _lazyRuns[row];
        trimmed = trimmed && _leavesOnlyOutput;
        const std::vector<uint32_t>& colByPos = trimmed ? _outByPos : _absByPos;
        const int msaLen = trimmed ? (int)_outLength : (int)_msaLength;
        if (cs.uncompressedSize == 0) { out.push_back(-msaLen); return; }
        const auto& runs = cs.runs;
        int totalSize = 0, cumulatedDifference = 1;
        int lastPosition = (int)colByPos[runs[0].first];
        if (lastPosition > 0) { out.push_back(-lastPosition); totalSize += lastPosition; }
        for (size_t r = 0; r < runs.size(); ++r) {
            const size_t start = runs[r].first, len = runs[r].second;
            for (size_t k = (r == 0 ? 1 : 0); k < len; ++k) {
                int currentPosition = (int)colByPos[start + k];
                int positionDifference = currentPosition - lastPosition - 1;
                if (positionDifference == 0) cumulatedDifference++;
                if (positionDifference > 0) {
                    out.push_back(cumulatedDifference);
                    out.push_back(-positionDifference);
                    totalSize += (cumulatedDifference + positionDifference);
                    cumulatedDifference = 1;
                }
                if (totalSize > msaLen) errorMsg::reportError("sequence lengths mismatch in fillRow");
                lastPosition = currentPosition;
            }
        }
        if (cumulatedDifference > 0 && (totalSize != msaLen)) { out.push_back(cumulatedDifference); totalSize += cumulatedDifference; }
        if (totalSize < msaLen) out.push_back(-(msaLen - totalSize));
    }

    template<typename BlockTreeType>
    void fillMSA(std::vector<CompressedSequence> &sequences, SuperSequence<RngType, BlockTreeType> &superSeq) {
		_numberOfSequences = superSeq.getNumSequences();
		_msaLength = superSeq.getMsaSequenceLength();
		_rootPositionsInMsa = superSeq.setAbsolutePositions();
        _perSiterateCategories = superSeq.getMsaRateCategories();
		// create AlignedSequenceMap shared pointer, can be shared with SubstitutionSimulator
        _alignedSequenceMap = std::make_shared<SparseMSA>();
        auto & _alignedSequence = *_alignedSequenceMap;
        _alignedSequence.resize(_numberOfSequences, std::vector<int>());


        int totalSize = 0;
		int currentPosition = 0;
		int lastPosition = 0;
		int positionDifference = 0;
		int cumulatedDifference = 1;
        
        const std::vector<uint32_t>& absByPos = superSeq.getAbsByPos();
        for(size_t i = 0; i < sequences.size(); ++i) {
            const auto& runs = sequences[i].runs;
            auto& currentSequence = _alignedSequence[i];
            // if the sequence is only made up of gaps:
            if (sequences[i].uncompressedSize == 0) {
                currentSequence.push_back(-_msaLength);
                continue;
            }
            lastPosition = absByPos[runs[0].first];
            if (lastPosition > 0) {
                currentSequence.push_back(-lastPosition);
                totalSize += lastPosition;
            }

            for (size_t r = 0; r < runs.size(); ++r) {
                const size_t start = runs[r].first, len = runs[r].second;
                for (size_t k = (r == 0 ? 1 : 0); k < len; ++k) {
                    currentPosition = absByPos[start + k];
                    positionDifference = currentPosition - lastPosition - 1;

                    if (positionDifference == 0) cumulatedDifference++;
                    if (positionDifference > 0) {
                        currentSequence.push_back(cumulatedDifference);
                        if (i == 0) _seq0NonGapSites += cumulatedDifference;
                        currentSequence.push_back(-(positionDifference));
                        totalSize += (cumulatedDifference + positionDifference);
                        cumulatedDifference = 1;
                    }
                    if (totalSize > _msaLength) errorMsg::reportError("sequence lengths mismatch in fillMSA");

                    lastPosition = currentPosition;
                }
            }
			if (cumulatedDifference > 0 && (totalSize != _msaLength)) {
                currentSequence.push_back(cumulatedDifference);
                if (i == 0) _seq0NonGapSites += cumulatedDifference; 
                totalSize += cumulatedDifference;
            }
            if (totalSize < _msaLength) currentSequence.push_back(-(_msaLength - totalSize));
			cumulatedDifference = 1;
            lastPosition = 0;
            totalSize = 0;
			// rowInMSA++;
        }
    };


    void fillSubstitutions(std::shared_ptr<const SparseSequenceContainer> _seqContainer) {
        _substitutions = _seqContainer;
    }

    void setCharLen(size_t charLen) {
        assert(charLen == 1 || charLen == 3);
        _charLen = charLen;
    }

	MSA<RngType>(size_t msaLength, SimulationContext<RngType> &simContext):
        _simContext(simContext), _msaLength(msaLength), 
        _numberOfSequences(simContext.getNumberOfNodesToSave()) {
            // create dummy aligned sequence map without gaps.
            _alignedSequenceMap = std::make_shared<SparseMSA>();
            auto & _alignedSequence = *_alignedSequenceMap;
            _alignedSequence.resize(_numberOfSequences, std::vector<int>());
            for (size_t i = 0; i < _numberOfSequences; i++) {
                _alignedSequence[i].push_back(msaLength);
            }
            if (simContext.getOutputLeavesOnly()) {
                _leavesOnlyOutput = true;
                _outputRows = simContext.getLeafRows();
                _outLength = msaLength; // no indels -> nothing to trim
            }
        };

	MSA<RngType>(const MSA<RngType> &msa){
        _simContext = msa._simContext;
		_numberOfSequences = msa._numberOfSequences;
		_msaLength = msa._msaLength;
        _alignedSequenceMap = msa._alignedSequenceMap;
	};

	// Output length (leaves-only output mode drops columns that are gap in all leaves)
	int getMSAlength() const {return _leavesOnlyOutput ? _outLength : _msaLength;}
	// Length in simulation coordinates: this is what the substitution simulator needs
	int getSimulationLength() const {return _msaLength;}
	int getNumberOfSequences() const {return _leavesOnlyOutput ? _outputRows.size() : _numberOfSequences;} 

	void printMSAInfo() {
		std::cout << _numberOfSequences << "x" << _numberOfSequences << "\n";
		std::cout << _msaLength << "\n";
	}





    std::string generateMsaRowString(size_t row) {
        std::string msaString;
        msaString.reserve((getMSAlength() * _charLen) + 64);

        msaString.append(">");
        const size_t internalRow = _leavesOnlyOutput ? _outputRows[row] : row; // `row` indexes output rows
        msaString.append(_simContext.getNodeToSaveNames()[internalRow]);
        msaString.append("\n");

        std::string currentSeq;

        if (_substitutions == nullptr) {
            currentSeq = "";
        } else {
            currentSeq = (*_substitutions)[row];
        }

        int passedSeq = 0;
        std::vector<int> alignedSeqRow; fillRow(internalRow, alignedSeqRow, true);
        for (size_t col = 0; col < alignedSeqRow.size(); col++) {
            int strSize = alignedSeqRow[col];
            if (strSize < 0) {
                strSize = -strSize;
                msaString.append(strSize* _charLen, '-');
            } else {
                if (_substitutions == nullptr) {
                    msaString.append(strSize * _charLen, 'A');
                } else {
                    msaString.append(currentSeq.substr(passedSeq, strSize * _charLen));
                }
                passedSeq += (strSize * _charLen);
            }
        }
        msaString.append("\n");

        return msaString;
    }

    // print the Msa to console line by line, to avoid creating a large string in memory that represents the whole MSA.
    void printFullMsa() {
        for (size_t row = 0; row < (size_t)getNumberOfSequences(); row++) {
            std::cout << generateMsaRowString(row);
        }
    }

    // overload "<<" operator to print MSA to output stream
    friend std::ostream& operator<<(std::ostream& os, MSA<RngType>& msa) {
        for (size_t row = 0; row < (size_t)msa.getNumberOfSequences(); row++) {
            os << msa.generateMsaRowString(row);
        }
        return os;
    }


    void writeFullMsa(const char * filePath) {
        ofstream msafile (filePath);
        if (msafile.is_open()) {
            for (size_t row = 0; row < (size_t)getNumberOfSequences(); row++) {
                msafile << generateMsaRowString(row);
            }
            msafile.close();
        }
        else cout << "Unable to open file";
    }

    std::shared_ptr<const SparseMSA> getSparseMSA() {
        if (_alignedSequenceMap == nullptr) {   // lazy mode: materialize all output rows on request
            const size_t n = getNumberOfSequences();
            auto m = std::make_shared<SparseMSA>(); m->resize(n);
            for (size_t r = 0; r < n; ++r) fillRow(_leavesOnlyOutput ? _outputRows[r] : r, (*m)[r], true);
            _alignedSequenceMap = m;
        } else if (!_lazy && _leavesOnlyOutput) {   // no-indel MSA: keep leaf rows only
            auto m = std::make_shared<SparseMSA>();
            for (size_t r : _outputRows) m->push_back(_alignedSequenceMap->at(r));
            return m;
        }
        return _alignedSequenceMap;
    }

    std::shared_ptr<std::vector<uint8_t>> getPerSiteRateCategories() const {
        return _perSiterateCategories;
    }

    const std::vector<size_t>& getRootPositionsInMsa() const { return _rootPositionsInMsa; }

	
	~MSA() {
	}

private:
    SimulationContext<RngType>& _simContext;
	size_t _numberOfSequences; // NUMBER OF SEQUENCES IN THE FINAL MSA
    size_t _msaLength; // Length of the MSA
    std::shared_ptr<const SparseSequenceContainer> _substitutions;

	std::shared_ptr<SparseMSA> _alignedSequenceMap = nullptr;
    std::shared_ptr<std::vector<uint8_t>> _perSiterateCategories = nullptr;
    // IDs of sequences to save are accessed via _simContext->getNodesToSaveIndices()
    // Names of sequences to save are accessed via _simContext->getNodeToSaveNames()
    std::vector<size_t> _rootPositionsInMsa;

    size_t _seq0NonGapSites = 0;
    // leaves-only output mode: internal row ids of the leaves, position-id -> output column, output length
    bool _leavesOnlyOutput = false;
    std::vector<size_t> _outputRows;
    std::vector<uint32_t> _outByPos;
    size_t _outLength = 0;
    std::vector<CompressedSequence> _lazyRuns;
    std::vector<uint32_t> _absByPos;
    bool _lazy = false;
    size_t _charLen = 1; // default to 1, will be set to 3 if alphabet size > 20
};
#endif
