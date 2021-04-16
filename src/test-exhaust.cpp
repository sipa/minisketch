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

std::vector<uint64_t> TestAll(int bits, int impl, int capacity) {
    bool supported = minisketch_implementation_supported(bits, impl);
    minisketch* state = minisketch_create(bits, impl, capacity);
    CHECK(supported == (state != nullptr));
    std::vector<uint64_t> ret;
    if (!state) return ret;

    // Iterate over all (bits)-bit sketches with (capacity) syndromes.
    for (uint64_t x = 0; (x >> (bits * capacity)) == 0; ++x) {
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
            minisketch* state2 = minisketch_create(bits, 0, capacity);
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

    for (int i = 1; i - 1 <= capacity && (i - 1) >> bits; ++i) {
        CHECK(ret[i] == Combination((uint64_t(1) << bits) - 1, i - 1));
    }
    return ret;
}

void TestRand(int bits, int impl, int capacity, int iter) {
    std::vector<uint64_t> elems(capacity);
    std::vector<uint64_t> roots(capacity + 1);
    std::random_device rnd;
    std::uniform_int_distribution<uint64_t> dist(1, bits == 64 ? -1 : ((uint64_t(1) << bits) - 1));

    for (int i = 0; i < iter; ++i) {
        bool overfill = iter & 1; // Test some cases with overfull sketches that may not decode.
        minisketch* state = minisketch_create(bits, impl, capacity);
        if (!state) return;
        minisketch* basestate = minisketch_create(bits, 0, capacity);
        for (int j = 0; j < capacity + 3*overfill; ++j) {
            uint64_t r = dist(rnd);
            if (!overfill) elems[j] = r;
            minisketch_add_uint64(state, r);
            minisketch_add_uint64(basestate, r);
        }
        roots.assign(capacity + 1, 0);
        std::vector<unsigned char> data, basedata;
        basedata.resize(((capacity + 1) * bits + 7) / 8);
        data.resize(((capacity + 1) * bits + 7) / 8);
        minisketch_serialize(basestate, basedata.data());
        minisketch_serialize(state, data.data());
        CHECK(data == basedata);
        minisketch_deserialize(state, basedata.data());
        int num_roots = minisketch_decode(state, capacity + 1, roots.data());
        CHECK(overfill || num_roots >= 0);
        CHECK(num_roots < 1 || minisketch_decode(state, num_roots - 1, roots.data()) == -1); // Decoding with a too-low maximum should fail.
        if (!overfill) {
            std::sort(roots.begin(), roots.begin() + num_roots);
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
            CHECK(std::equal(roots.begin(), roots.begin() + num_roots, elems.end() - expected));
        }
        minisketch_destroy(state);
        minisketch_destroy(basestate);
    }
}

void TestComputeFunctions() {
    for (uint32_t bits = 0; bits <= 256; ++bits) {
        for (uint32_t fpbits = 0; fpbits <= 512; ++fpbits) {
            std::vector<size_t> table_max_elements(1025);
            for (size_t capacity = 0; capacity <= 1024; ++capacity) {
                table_max_elements[capacity] = minisketch_compute_max_elements(bits, capacity, fpbits);
                // Exception for bits==0
                if (bits == 0) CHECK(table_max_elements[capacity] == 0);
                // A sketch with capacity N cannot guarantee decoding more than N elements.
                CHECK(table_max_elements[capacity] <= capacity);
                // When asking for N bits of false positive protection, either no solution exists, or no more than ceil(N / bits) excess capacity should be needed.
                if (bits > 0) CHECK(table_max_elements[capacity] == 0 || capacity - table_max_elements[capacity] <= (fpbits + bits - 1) / bits);
                // Increasing capacity by one, if there is a solution, should always increment the max_elements by at least one as well.
                if (capacity > 0) CHECK(table_max_elements[capacity] == 0 || table_max_elements[capacity] > table_max_elements[capacity - 1]);
            }

            std::vector<size_t> table_capacity(513);
            for (size_t max_elements = 0; max_elements <= 512; ++max_elements) {
                table_capacity[max_elements] = minisketch_compute_capacity(bits, max_elements, fpbits);
                // Exception for bits==0
                if (bits == 0) CHECK(table_capacity[max_elements] == 0);
                // To be able to decode N elements, capacity needs to be at least N.
                if (bits > 0) CHECK(table_capacity[max_elements] >= max_elements);
                // A sketch of N bits in total cannot have more than N bits of false positive protection;
                if (bits > 0) CHECK(bits * table_capacity[max_elements] >= fpbits);
                // When asking for N bits of false positive protection, no more than ceil(N / bits) excess capacity should be needed.
                if (bits > 0) CHECK(table_capacity[max_elements] - max_elements <= (fpbits + bits - 1) / bits);
                // Increasing max_elements by one can only increment the capacity by 0 or 1.
                if (max_elements > 0 && fpbits < 256) CHECK(table_capacity[max_elements] == table_capacity[max_elements - 1] || table_capacity[max_elements] == table_capacity[max_elements - 1] + 1);
                // Check round-tripping max_elements->capacity->max_elements (only a lower bound)
                CHECK(table_capacity[max_elements] <= 1024);
                CHECK(table_max_elements[table_capacity[max_elements]] == 0 || table_max_elements[table_capacity[max_elements]] >= max_elements);
            }

            for (size_t capacity = 0; capacity <= 512; ++capacity) {
                // Check round-tripping capacity->max_elements->capacity (exact, if it exists)
                CHECK(table_max_elements[capacity] <= 512);
                CHECK(table_max_elements[capacity] == 0 || table_capacity[table_max_elements[capacity]] == capacity);
            }
        }
    }
}

int main(void) {
    TestComputeFunctions();

    for (int j = 2; j <= 64; j += 1) {
        TestRand(j, 0, 150, 500 / j);
        TestRand(j, 1, 150, 500 / j);
        TestRand(j, 2, 150, 500 / j);
    }

    for (int weight = 2; weight <= 40; weight += 1) {
        for (int bits = 2; bits <= 32 && bits <= weight; ++bits) {
            int capacity = weight / bits;
            if (capacity * bits != weight) continue;
            auto ret = TestAll(bits, 0, capacity);
            auto ret2 = TestAll(bits, 1, capacity);
            auto ret3 = TestAll(bits, 2, capacity);
            CHECK(ret2.empty() || ret == ret2);
            CHECK(ret3.empty() || ret == ret3);
        }
    }

    return 0;
}
