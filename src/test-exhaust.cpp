/**********************************************************************
 * Copyright (c) 2018,2021 Pieter Wuille, Greg Maxwell, Gleb Naumenko *
 * Distributed under the MIT software license, see the accompanying   *
 * file LICENSE or http://www.opensource.org/licenses/mit-license.php.*
 **********************************************************************/

#include "../include/minisketch.h"
#include <vector>
#include <algorithm>
#include <stdexcept>
#include <random>
#include <stdio.h>
#include "util.h"

namespace {

uint64_t Combination(uint64_t n, uint64_t k) {
    if (n - k < k) k = n - k;
    uint64_t ret = 1;
    for (uint64_t i = 1; i <= k; ++i) {
        ret = (ret * n) / i;
        --n;
    }
    return ret;
}

/** Create a vector with Minisketch objects, one for each implementation. */
std::vector<Minisketch> CreateSketches(uint32_t bits, size_t capacity) {
    if (!Minisketch::BitsSupported(bits)) return {};
    std::vector<Minisketch> ret;
    for (uint32_t impl = 0; impl <= Minisketch::MaxImplementation(); ++impl) {
        if (Minisketch::ImplementationSupported(bits, impl)) {
            ret.push_back(Minisketch(bits, impl, capacity));
            CHECK((bool)ret.back());
        } else {
            CHECK(impl != 0); // implementation 0 must always work
        }
    }
    return ret;
}

/** Test properties by exhaustively decoding all 2**(bits*capacity) sketches
 *  with specified capacity and bits. */
void TestExhaustive(uint32_t bits, size_t capacity) {
    auto sketches = CreateSketches(bits, capacity);
    CHECK(!sketches.empty());
    auto sketches_rebuild = CreateSketches(bits, capacity);

    std::vector<unsigned char> serialized;
    std::vector<unsigned char> serialized_empty;
    std::vector<uint64_t> counts; //!< counts[i] = number of results with i elements
    std::vector<uint64_t> elements_0; //!< Result vector for elements for impl=0
    std::vector<uint64_t> elements_other; //!< Result vector for elements for other impls
    std::vector<uint64_t> elements_too_small; //!< Result vector that's too small

    counts.resize(capacity + 1);
    serialized.resize(sketches[0].GetSerializedSize());
    serialized_empty.resize(sketches[0].GetSerializedSize());

    // Iterate over all (bits)-bit sketches with (capacity) syndromes.
    for (uint64_t x = 0; (x >> (bits * capacity)) == 0; ++x) {
        // Construct the serialization.
        for (size_t i = 0; i < serialized.size(); ++i) {
            serialized[i] = (x >> (i * 8)) & 0xFF;
        }

        // Compute all the solutions
        sketches[0].Deserialize(serialized);
        elements_0.resize(64);
        bool decodable_0 = sketches[0].Decode(elements_0);
        std::sort(elements_0.begin(), elements_0.end());

        // Verify that decoding with other implementations agrees.
        for (size_t impl = 1; impl < sketches.size(); ++impl) {
            sketches[impl].Deserialize(serialized);
            elements_other.resize(64);
            bool decodable_other = sketches[impl].Decode(elements_other);
            CHECK(decodable_other == decodable_0);
            std::sort(elements_other.begin(), elements_other.end());
            CHECK(elements_other == elements_0);
        }

        // If there are solutions:
        if (decodable_0) {
            if (!elements_0.empty()) {
                // Decoding with limit one less than the number of elements should fail.
                elements_too_small.resize(elements_0.size() - 1);
                for (size_t impl = 0; impl < sketches.size(); ++impl) {
                    CHECK(!sketches[impl].Decode(elements_too_small));
                }
            }

            // Reconstruct the sketch from the solutions.
            for (size_t impl = 0; impl < sketches.size(); ++impl) {
                // Clear the sketch.
                sketches_rebuild[impl].Deserialize(serialized_empty);
                // Load all decoded elements into it.
                for (uint64_t elem : elements_0) {
                    CHECK(elem != 0);
                    CHECK(elem >> bits == 0);
                    sketches_rebuild[impl].Add(elem);
                }
                // Reserialize the result
                auto serialized_rebuild = sketches_rebuild[impl].Serialize();
                // Compare
                CHECK(serialized == serialized_rebuild);
                // Count it
                if (elements_0.size() <= capacity) ++counts[elements_0.size()];
            }
        }
    }

    // Verify that the number of decodable sketches with given elements is expected.
    for (uint64_t i = 0; i <= capacity && i >> bits; ++i) {
        CHECK(counts[i] == Combination((uint64_t{1} << bits) - 1, i));
    }
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

} // namespace

int main(int argc, char** argv) {
    uint64_t test_complexity = 4;
    if (argc > 1) {
        size_t len = 0;
        std::string arg{argv[1]};
        try {
            test_complexity = 0;
            long long complexity = std::stoll(arg, &len);
            if (complexity >= 1 && len == arg.size() && ((uint64_t)complexity <= std::numeric_limits<uint64_t>::max() >> 7)) {
                test_complexity = complexity;
            }
        } catch (const std::logic_error&) {}
        if (test_complexity == 0) {
            fprintf(stderr, "Invalid complexity specified: `%s'\n", arg.c_str());
            return 1;
        }
    }

#ifdef MINISKETCH_VERIFY
    const char* mode = " in verify mode";
#else
    const char* mode = "";
#endif
    printf("Running libminisketch tests%s with complexity=%llu\n", mode, (unsigned long long)test_complexity);

    TestComputeFunctions();

    for (unsigned j = 2; j <= 64; j += 1) {
        TestRand(j, 0, 150, (test_complexity << 7) / j);
        TestRand(j, 1, 150, (test_complexity << 7) / j);
        TestRand(j, 2, 150, (test_complexity << 7) / j);
    }

    for (int weight = 2; weight <= 40; ++weight) {
        for (int bits = 2; bits <= 32 && bits <= weight; ++bits) {
            int capacity = weight / bits;
            if (capacity * bits != weight) continue;
            TestExhaustive(bits, capacity);
        }
        if (weight >= 16 && test_complexity >> (weight - 16) == 0) break;
    }

    printf("All tests successful.\n");
    return 0;
}
