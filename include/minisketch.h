#ifndef _MINISKETCH_H_
#define _MINISKETCH_H_ 1

#ifdef __cplusplus
extern "C" {
#endif

#include <stdint.h>
#include <stdlib.h>
#include <unistd.h>

/** Opaque type for decoded sketches. */
typedef struct minisketch minisketch;

/** Determine whether support for elements of `bits` bits was compiled in. */
int minisketch_bits_supported(uint32_t bits);

/** Determine the maximum number of implementations available.
 *
 * Multiple implementations may be available for a given element size, with
 * different performance characteristics on different hardware.
 *
 * Each implementation is identified by a number from 0 to the output of this
 * function call, inclusive. Note that not every combination of implementation
 * and element size may exist (see further).
*/
uint32_t minisketch_implementation_max(void);

/** Construct a sketch for a given element size, implementation and capacity.
 *
 * If the combination of `bits` and `implementation` is unavailable, or if
 * `capacity` is 0, NULL is returned.
 * If the result is not NULL, it must be destroyed using minisketch_destroy.
 */
minisketch* minisketch_create(uint32_t bits, uint32_t implementation, size_t capacity);

/** Get the element size of a sketch in bits. */
uint32_t minisketch_bits(const minisketch* sketch);

/** Get the capacity of a sketch. */
size_t minisketch_capacity(const minisketch* sketch);

/** Get the implementation of a sketch. */
uint32_t minisketch_implementation(const minisketch* sketch);

/** Set the seed for randomizing algorithm choices to a fixed value.
 *
 * By default, sketches are initialized with a random seed. This is important
 * to avoid scenarios where an attacker could force worst-case behavior.
 *
 * This function initializes the seed to a user-provided value (any 64-bit
 * integer is acceptable, regardless of field size).
 *
 * When seed is -1, a fixed internal value with predictable behavior is
 * used. It is only intended for testing.
 */
void minisketch_set_seed(minisketch* sketch, uint64_t seed);

/** Clone a sketch.
 *
 * The result must be destroyed using minisketch_destroy.
 */
minisketch* minisketch_clone(const minisketch* sketch);

/** Destroy a sketch.
 *
 * The pointer that was passed in may not be used anymore afterwards.
 */
void minisketch_destroy(minisketch* sketch);

/** Compute the size in bytes for serializing a given sketch. */
size_t minisketch_serialized_size(const minisketch* sketch);

/** Serialize a sketch to bytes. */
void minisketch_serialize(const minisketch* sketch, unsigned char* output);

/** Deserialize a sketch from bytes. */
void minisketch_deserialize(minisketch* sketch, const unsigned char* input);

/** Add an element to a sketch.
 * 
 * If the element to be added is too large for the sketch, the most significant
 * bits of the element are dropped. More precisely, if the element size of
 * `sketch` is b bits, then this function adds the unsigned integer represented
 * by the b least significant bits of `element` to `sketch`.
 * 
 * If the element to be added is 0 (after potentially dropping the most significant
 * bits), then this function is a no-op. Sketches cannot contain an element with
 * the value 0.
 */
void minisketch_add_uint64(minisketch* sketch, uint64_t element);

/** Merge the elements of another sketch into this sketch.
 *
 * After merging, `sketch` will contain every element that existed in one but not
 * both of the input sketches. It can be seen as an exclusive or operation on
 * the set elements.  If the capacity of `other_sketch` is lower than `sketch`'s,
 * merging reduces the capacity of `sketch` to that of `other_sketch`.
 *
 * This function returns the capacity of `sketch` after merging has been performed
 * (where this capacity is at least 1), or 0 to indicate that merging has failed because
 * the two input sketches differ in their element size or implementation. If 0 is
 * returned, `sketch` (and its capacity) have not been modified.
 *
 * It is also possible to perform this operation directly on the serializations
 * of two sketches with the same element size and capacity by performing a bitwise XOR
 * of the serializations.
 */
size_t minisketch_merge(minisketch* sketch, const minisketch* other_sketch);

/** Decode a sketch.
 *
 * `output` is a pointer to an array of `max_element` uint64_t's, which will be
 * filled with the elements in this sketch.
 *
 * The return value is the number of decoded elements, or -1 if decoding failed.
 */
ssize_t minisketch_decode(const minisketch* sketch, size_t max_elements, uint64_t* output);

/** Compute the capacity needed to achieve a certain rate of false positives.
 *
 * A sketch with capacity c and no more than c elements can always be decoded
 * correctly. However, if it has more than c elements, or contains just random
 * bytes, it is possible that it will still decode, but the result will be
 * nonsense. This can be counteracted by increasing the capacity slightly.
 *
 * Given a field size bits, an intended number of elements that can be decoded
 * max_elements, and a false positive probability of 1 in 2**fpbits, this
 * function computes the necessary capacity. It is only guaranteed to be
 * accurate up to fpbits=256.
 */
size_t minisketch_compute_capacity(uint32_t bits, size_t max_elements, uint32_t fpbits);

/** Compute what max_elements can be decoded for a certain rate of false positives.
 *
 * This is the inverse operation of minisketch_compute_capacity. It determines,
 * given a field size bits, a capacity of a sketch, and an acceptable false
 * positive probability of 1 in 2**fpbits, what the maximum allowed
 * max_elements value is. If no value of max_elements would give the desired
 * false positive probability, 0 is returned.
 *
 * Note that this is not an exact inverse of minisketch_compute_capacity. For
 * example, with bits=32, fpbits=16, and max_elements=8,
 * minisketch_compute_capacity will return 9, as capacity 8 would only have a
 * false positive chance of 1 in 2^15.3. Increasing the capacity to 9 however
 * decreases the fp chance to 1 in 2^47.3, enough for max_elements=9 (with fp
 * chance of 1 in 2^18.5). Therefore, minisketch_compute_max_elements with
 * capacity=9 will return 9.
 */
size_t minisketch_compute_max_elements(uint32_t bits, size_t capacity, uint32_t fpbits);

#ifdef __cplusplus
}
#endif

#endif
