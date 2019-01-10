/**********************************************************************
 * Copyright (c) 2018 Pieter Wuille, Greg Maxwell, Gleb Naumenko      *
 * Distributed under the MIT software license, see the accompanying   *
 * file LICENSE or http://www.opensource.org/licenses/mit-license.php.*
 **********************************************************************/

#include "../include/minisketch.h"
#include <string.h>
#include <memory>
#include <vector>
#include <algorithm>
#include <random>
#include <iostream>
#include <thread>
#include "util.h"

uint64_t Combination(uint64_t n, uint64_t k) {
    if (n - k < k) k = n - k;
    uint64_t ret = 1;
    for (uint64_t i = 1; i <= k; ++i) {
        ret = (ret * n) / i;
        --n;
    }
    return ret;
}

void TestAll(int bits, int impl, int count, uint32_t threadid, uint32_t threads, std::vector<uint64_t>& ret) {
    minisketch* state = minisketch_create(bits, impl, count);
    if (!state) return;

    // Iterate over all (bits)-bit sketches with (count) syndromes.
    for (uint64_t x = threadid; (x >> (bits * count)) == 0; x += threads) {
        // Construct the serialization and load it.
        unsigned char ser[8];
        ser[0] = x;
        ser[1] = x >> 8;
        ser[2] = x >> 16;
        ser[3] = x >> 24;
        ser[4] = x >> 32;
        ser[5] = x >> 40;
        ser[6] = x >> 48;
        ser[7] = x >> 56;

        minisketch_deserialize(state, ser);

        // Compute all the solutions.
        uint64_t roots[64];
        int num_roots = minisketch_decode(state, 64, roots);

        // If there are solutions:
        if (num_roots >= 0) {
            // Asking for one root less should fail.
            CHECK(num_roots < 1 || minisketch_decode(state, num_roots - 1, roots) == -1);
            // Reconstruct the sketch from the solutions.
            minisketch* state2 = minisketch_create(bits, 0, count);
            for (int i = 0; i < num_roots; ++i) {
                minisketch_add_uint64(state2, roots[i]);
            }
            // Serialize it.
            unsigned char nser[8] = {0};
            minisketch_serialize(state2, nser);
            // Compare it to the original.
            CHECK(memcmp(ser, nser, 8) == 0);
            // Count it.
            if (num_roots +1 >= (int)ret.size()) ret.resize(num_roots + 2);
            ret[num_roots + 1]++;
            minisketch_destroy(state2);
        } else {
            if (ret.size() == 0) ret.resize(1);
            ret[0]++;
        }
    }
    minisketch_destroy(state);
}

std::vector<uint64_t> TestAll(int bits, int impl, int count, uint32_t threads) {
    std::vector<std::vector<uint64_t>> outputs;
    std::vector<std::thread> thread_list;
    thread_list.reserve(threads);
    outputs.resize(threads);
    for (uint32_t i = 0; i < threads; ++i) {
        thread_list.emplace_back([=,&outputs](){ TestAll(bits, impl, count, i, threads, outputs[i]); });
    }
    std::vector<uint64_t> ret;
    for (uint32_t i = 0; i < threads; ++i) {
        thread_list[i].join();
        if (ret.size() < outputs[i].size()) ret.resize(outputs[i].size());
        for (size_t j = 0; j < outputs[i].size(); ++j) {
            ret[j] += outputs[i][j];
        }
    }
    if (ret.size()) {
        for (int i = 1; i <= count + 1; ++i) {
            CHECK(ret[i] == Combination((uint64_t(1) << bits) - 1, i - 1));
        }
    }
    return ret;
}

void TestRand(int bits, int impl, int count, int iter) {
    std::vector<uint64_t> elems(count);
    std::vector<uint64_t> roots(count + 1);
    std::random_device rnd;
    std::uniform_int_distribution<uint64_t> dist(1, bits == 64 ? -1 : ((uint64_t(1) << bits) - 1));

    for (int i = 0; i < iter; ++i) {
        bool overfill = iter & 1; // Test some cases with overfull sketches that may not decode.
        minisketch* state = minisketch_create(bits, impl, count);
        if (!state) return;
        minisketch* basestate = minisketch_create(bits, 0, count);
        for (int j = 0; j < count + 3*overfill; ++j) {
            uint64_t r = dist(rnd);
            if (!overfill) elems[j] = r;
            minisketch_add_uint64(state, r);
            minisketch_add_uint64(basestate, r);
        }
        roots.assign(count + 1, 0);
        std::vector<unsigned char> data, basedata;
        basedata.resize(((count + 1) * bits + 7) / 8);
        data.resize(((count + 1) * bits + 7) / 8);
        minisketch_serialize(basestate, basedata.data());
        minisketch_serialize(state, data.data());
        CHECK(data == basedata);
        minisketch_deserialize(state, basedata.data());
        int num_roots = minisketch_decode(state, count + 1, roots.data());
        CHECK(overfill || num_roots >= 0);
        CHECK(num_roots < 1 || minisketch_decode(state, num_roots - 1, roots.data()) == -1); // Decoding with a too-low maximum should fail.
        if (!overfill) {
            std::sort(roots.begin(), roots.begin() + num_roots);
//            fprintf(stderr, "Solut: ");
//            for (int j = 0; j < num_roots; ++j) {
//                fprintf(stderr, "%016llx ", (unsigned long long)roots[j]);
//            }
//            fprintf(stderr, "\n");
            std::sort(elems.begin(), elems.end());
            int expected = elems.size();
            for (size_t pos = 0; pos < elems.size(); ++pos) {
                if (pos + 1 < elems.size() && elems[pos] == elems[pos + 1]) {
                    expected -= 2;
                    elems[pos] = 0;
                    elems[pos + 1] = 0;
                    ++pos;
                }
            }
            CHECK(num_roots == expected);
            std::sort(elems.begin(), elems.end());
//            fprintf(stderr, "Elems: ");
//            for (int j = 0; j < expected; ++j) {
//                fprintf(stderr, "%016llx ", (unsigned long long)elems[j + elems.size() - expected]);
//            }
//            fprintf(stderr, "\n");
            CHECK(std::equal(roots.begin(), roots.begin() + num_roots, elems.end() - expected));
        }
        minisketch_destroy(state);
        minisketch_destroy(basestate);
    }
}

int main(void) {
    for (int j = 2; j <= 64; j += 1) {
        fprintf(stderr, "%i random tests with %i bits:\n", 500 / j, j);
        TestRand(j, 0, 150, 500 / j);
        TestRand(j, 1, 150, 500 / j);
        TestRand(j, 2, 150, 500 / j);
        fprintf(stderr, "%i random tests with %i bits: done\n", 500 / j, j);
    }

    int counts[65] = {0};
    // Initialize capacities to 1 because a 0 capacity is not allowed.
    for (int bits = 0; bits < 65; ++bits) {
        counts[bits] = 1;
    }
    for (int weight = 0; weight <= 40; weight += 1) {
        for (int bits = 2; bits <= 32; ++bits) {
            int count = counts[bits];
            while (count < (1 << bits) && count * bits <= weight) {
                auto ret = TestAll(bits, 0, count, 4);
                auto ret2 = TestAll(bits, 1, count, 4);
                auto ret3 = TestAll(bits, 2, count, 4);
                CHECK(ret2.empty() || ret == ret2);
                CHECK(ret3.empty() || ret == ret3);
                fprintf(stderr, "bits=%i count=%i below_bound=[", bits, count);
                for (int i = 0; i <= count; ++i) {
                    if (i) fprintf(stderr,  ",");
                    fprintf(stderr, "%llu", (unsigned long long)ret[i + 1]);
                }
                fprintf(stderr, "] above_bound=[");
                for (int i = count + 1; i + 1 < (int)ret.size(); ++i) {
                    if (i > count + 1) fprintf(stderr,  ",");
                    fprintf(stderr, "%llu/%llu", (unsigned long long)ret[i + 1], (unsigned long long)Combination((uint64_t(1) << bits) - 1, i));
                }
                fprintf(stderr, "] nodecode=[%g]\n", (double)ret[0] * pow(0.5, bits * count));
                ++count;
            }
            counts[bits] = count;
        }
    }

    return 0;
}
