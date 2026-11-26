#pragma once
#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>
#include <random>
#include <utility>
#include <vector>

// Alias-method sampler over N outcomes (N is a power of two).
// Each bin is one 64-bit word: [ threshold (high bits) | alias (low bits) ].
// One random number picks a bin (its low bits) and decides between
// that bin and its alias (its high bits).
template<uint8_t N>
class DiscreteNDistribution
{
    static_assert(N >= 2 && (N & (N - 1)) == 0, "N must be a power of two");

private:
    static constexpr int countTrailingZeros(unsigned value) {
        return (value & 1) ? 0 : 1 + countTrailingZeros(value >> 1);
    }

    static constexpr int kIndexBits = countTrailingZeros(N);   // bits needed for a bin number
    static constexpr uint64_t kIndexMask = N - 1;              // selects those bits
    static constexpr double kThresholdScale = static_cast<double>(1ULL << (64 - kIndexBits));

    // During construction, probability is counted in integer "mass": one bin holds kBinCapacity.
    // All N bins together hold 2^62, which fits safely in a signed 64-bit integer.
    static constexpr int64_t kBinCapacity = 1LL << (62 - kIndexBits);
    // Shift that turns a mass into the high bits of a table word.
    static constexpr int kMassToThresholdShift = kIndexBits + 2;

    std::array<uint64_t, N> table_;

    // bin that keeps `ownMass` of its capacity and gives the rest to `alias`
    static uint64_t makeEntry(int64_t ownMass, int alias) {
        return (static_cast<uint64_t>(ownMass) << kMassToThresholdShift) | static_cast<uint64_t>(alias);
    }

    // full bin: maximum threshold, alias = itself
    static uint64_t makeFullEntry(int bin) {
        return (~0ULL << kIndexBits) | static_cast<uint64_t>(bin);
    }

public:
    DiscreteNDistribution(const std::array<double, N> &probabilities, double normalizingFactor = 1.0) {
        const double probabilityToMass =
            static_cast<double>(N) * normalizingFactor * static_cast<double>(kBinCapacity);

        int64_t mass[N];         // mass still to be placed, per outcome
        uint8_t smallBins[N];    // stack: outcomes with less than one bin of mass
        uint8_t largeBins[N];    // stack: outcomes with one bin of mass or more
        int numSmall = 0, numLarge = 0;

        for (int bin = 0; bin < N; bin++) {
            // max() guards against tiny negative probabilities
            mass[bin] = static_cast<int64_t>(std::max(probabilities[bin], 0.0) * probabilityToMass);

            // write to both stacks, advance only the right one (avoids a branch)
            const bool isSmall = mass[bin] < kBinCapacity;
            smallBins[numSmall] = static_cast<uint8_t>(bin);
            largeBins[numLarge] = static_cast<uint8_t>(bin);
            numSmall += isSmall;
            numLarge += !isSmall;
        }

        // Shortcut: with a single large outcome, every small bin must alias to it.
        if (numLarge == 1) {
            const int largeBin = largeBins[0];
            for (int i = 0; i < numSmall; ++i) {
                const int smallBin = smallBins[i];
                table_[smallBin] = makeEntry(mass[smallBin], largeBin);
            }
            table_[largeBin] = makeFullEntry(largeBin);
            return;
        }

        // General case: fill each small bin's gap with mass from a large outcome.
        if (numSmall > 0 && numLarge > 0) {
            int largeBin = largeBins[--numLarge];
            int64_t largeMass = mass[largeBin];
            while (true) {
                const int smallBin = smallBins[--numSmall];
                table_[smallBin] = makeEntry(mass[smallBin], largeBin);
                largeMass -= kBinCapacity - mass[smallBin];

                if (largeMass < kBinCapacity) {
                    // the large outcome became small: take the next large one
                    mass[largeBin] = largeMass;
                    smallBins[numSmall++] = static_cast<uint8_t>(largeBin);
                    if (numLarge == 0) break;
                    largeBin = largeBins[--numLarge];
                    largeMass = mass[largeBin];
                } else if (numSmall == 0) {
                    largeBins[numLarge++] = static_cast<uint8_t>(largeBin);
                    break;
                }
            }
        }

        // Whatever is left fills its own bin (up to rounding).
        while (numLarge > 0) {
            const int bin = largeBins[--numLarge];
            table_[bin] = makeFullEntry(bin);
        }
        while (numSmall > 0) {
            const int bin = smallBins[--numSmall];
            table_[bin] = makeFullEntry(bin);
        }
    }

    // returns an outcome in 1..N
    template<typename RngType = std::mt19937_64>
    int drawSample(RngType &rng) const {
        const uint64_t random = rng();
        const uint64_t entry = table_[random & kIndexMask];

        // Index bits are set to 1 for the comparison only, so a tie picks the alias.
        // keepMask is all ones (keep the bin) or all zeros (take the alias).
        const uint64_t keepMask = -static_cast<uint64_t>((random | kIndexMask) < entry);
        const uint64_t chosen = (random & keepMask) | (entry & ~keepMask);

        return static_cast<int>(chosen & kIndexMask) + 1;
    }

    // chance that bin i keeps its own outcome
    double getProb(int i) const {
        return static_cast<double>(table_[i] >> kIndexBits) / kThresholdScale;
    }

    std::vector<std::pair<double, int>> getTable() const {
        std::vector<std::pair<double, int>> out;
        for (int i = 0; i < N; i++) {
            out.emplace_back(getProb(i), static_cast<int>(table_[i] & kIndexMask));
        }
        return out;
    }

    void printTable() const {
        for (auto &entry : getTable()) std::cout << entry.first << " ";
        std::cout << "\n";
        for (auto &entry : getTable()) std::cout << entry.second << " ";
        std::cout << "\n";
    }
};