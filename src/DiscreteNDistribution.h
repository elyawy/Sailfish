#pragma once
#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>
#include <random>
#include <utility>
#include <vector>

template<uint8_t N>
class DiscreteNDistribution
{
    static_assert(N >= 2 && (N & (N - 1)) == 0, "N must be a power of two");

private:
    static constexpr int ctz(unsigned v) { return (v & 1) ? 0 : 1 + ctz(v >> 1); }
    static constexpr int kBits = ctz(N);            // low bits: bin / alias
    static constexpr uint64_t kMask = N - 1;
    static constexpr double kScale = static_cast<double>(1ULL << (64 - kBits));

    // each word: [ threshold | alias ]
    std::array<uint64_t, N> table_;

    static uint64_t pack(double prob, int alias) {
        // max() guards against tiny negative rounding leftovers
        uint64_t threshold = static_cast<uint64_t>(std::max(prob, 0.0) * kScale);
        return (threshold << kBits) | static_cast<uint64_t>(alias);
    }

    // full bin: maximum threshold, alias = itself
    static uint64_t packFull(int i) {
        return (~0ULL << kBits) | static_cast<uint64_t>(i);
    }

public:
    DiscreteNDistribution(const std::array<double, N> &probabilities, double normalizingFactor=1.0) {
        std::pair<int, double> small_buf[N];
        std::pair<int, double> large_buf[N];
        int small_top = 0, large_top = 0;

        for(int i = 0; i < N; i++) {
            double scaled_prob = N * probabilities[i] * normalizingFactor;
            if (scaled_prob < 1.0) {
                small_buf[small_top++] = {i, scaled_prob};
            } else {
                large_buf[large_top++] = {i, scaled_prob};
            }
        }

        while (small_top > 0 && large_top > 0) {
            auto s = small_buf[--small_top];
            auto l = large_buf[--large_top];

            table_[s.first] = pack(s.second, l.first);

            l.second = (l.second + s.second) - 1.0;
            if (l.second < 1.0) {
                small_buf[small_top++] = l;
            } else {
                large_buf[large_top++] = l;
            }
        }

        while (large_top > 0) {
            int i = large_buf[--large_top].first;
            table_[i] = packFull(i);
        }
        while (small_top > 0) {
            int i = small_buf[--small_top].first;
            table_[i] = packFull(i);
        }
    }

    template<typename RngType = std::mt19937_64>
    int drawSample(RngType &rng) const {
        uint64_t v = rng();
        uint64_t p = table_[v & kMask];

        // bin bits set to 1 for the comparison only, so a tie picks the alias
        uint64_t keep = -static_cast<uint64_t>((v | kMask) < p);  // all 1s or all 0s
        uint64_t picked = (v & keep) | (p & ~keep);

        return static_cast<int>(picked & kMask) + 1;
    }

    double getProb(int i) const {
        return static_cast<double>(table_[i] >> kBits) / kScale;
    }

    std::vector<std::pair<double, int>> getTable() const {
        std::vector<std::pair<double, int>> out;
        for (int i = 0; i < N; i++) {
            out.emplace_back(getProb(i), static_cast<int>(table_[i] & kMask));
        }
        return out;
    }

    void printTable() const {
        for (auto &e : getTable()) std::cout << e.first << " ";
        std::cout << "\n";
        for (auto &e : getTable()) std::cout << e.second << " ";
        std::cout << "\n";
    }
};