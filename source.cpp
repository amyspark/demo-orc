#include <stdint.h>
#include <stddef.h>
#include <immintrin.h>

#define ORC_RESTRICT

typedef uint8_t guint8;
typedef uint8_t orc_uint8;
typedef int8_t orc_int8;
typedef uint16_t orc_uint16;
typedef int16_t orc_int16;
typedef uint32_t orc_uint32;
typedef int32_t orc_int32;
typedef int64_t orc_int64;

typedef union { orc_int16 i; orc_int8 x2[2]; } orc_union16;

typedef union
{
  int32_t i;
  float f;
  int16_t x2[2];
  int8_t x4[4];
} orc_union32;

typedef union
{
  int64_t i;
  double f;
  int16_t x2[4];
  int8_t x4[8];
} orc_union64;

#define ORC_UINT64_C(x) UINT64_C(x)
#define ORC_CLAMP(x,a,b) ((x)<(a) ? (a) : ((x)>(b) ? (b) : (x)))
#define ORC_CLAMP_UL(x) ORC_CLAMP(x,0,UINT32_MAX)
#define ORC_MIN(a,b) ((a)<(b) ? (a) : (b))
#define ORC_ABS(a) ((a)<0 ? -(a) : (a))
#define ORC_SB_MAX 127
#define ORC_SB_MIN (-1-ORC_SB_MAX)

#define ORC_STATIC_OPCODE_N_SRC 4
#define ORC_STATIC_OPCODE_N_DEST 2

typedef struct _OrcOpcodeExecutor OrcOpcodeExecutor;
typedef void (*OrcOpcodeEmulateNFunc)(OrcOpcodeExecutor *ex, int index, int n);

/**
 * OrcOpcodeExecutor:
 *
 * The OrcOpcodeExecutor structure has no public members
 */
struct _OrcOpcodeExecutor {
  /*< private >*/
  int src_values[ORC_STATIC_OPCODE_N_SRC];
  int dest_values[ORC_STATIC_OPCODE_N_DEST];

  OrcOpcodeEmulateNFunc emulateN;

  void *src_ptrs[ORC_STATIC_OPCODE_N_SRC];
  void *dest_ptrs[ORC_STATIC_OPCODE_N_DEST];
  int shift;
};

void
emulate_convusswb (OrcOpcodeExecutor *ex, int offset, int n)
{
  int i;
  orc_int8 * ORC_RESTRICT ptr0;
  const orc_union16 * ORC_RESTRICT ptr4;
  orc_union16 var32;
  orc_int8 var33;

  ptr0 = (orc_int8 *)ex->dest_ptrs[0];
  ptr4 = (orc_union16 *)ex->src_ptrs[0];


  for (i = 0; i < n; i++) {
    /* 0: loadw */
    var32 = ptr4[i];
    /* 1: convusswb */
    var33 = ORC_MIN((orc_uint16)var32.i,ORC_SB_MAX);
    /* 2: storeb */
    ptr0[i] = var33;
  }

}

void
emulate_convusswb_2 (const int8_t *src, uint16_t *dst, int n)
{
    for (int i = 0; i < n; i++) {
        *dst = ORC_MIN((orc_uint16)src[i],ORC_SB_MAX);
    }
}

void
emulate_convusswb_3 (const int8_t *src, uint16_t *dst, int n)
{
    for (int i = 0; i < n; i+=8) {
        __m128i x = _mm_loadu_si128((const __m128i*)src);
        const __m128i j = _mm_set1_epi8(127u);
        const __m128i k = _mm_setzero_si128();
        x = _mm_min_epu8(x, j); // shave off everything negative
        _mm_storeu_si128((__m128i*)dst, _mm_unpacklo_epi8(x, k));
        _mm_storeu_si128((__m128i*)dst+8, _mm_unpackhi_epi8(x, k));
    }
}

void
emulate_convusswb_4 (const int8_t *src, uint16_t *dst, int n)
{
    for (int i = 0; i < n; i+=16) {
        __m256i x = _mm256_loadu_si256((const __m256i*)src);
        const __m256i j = _mm256_set1_epi8(127u);
        const __m256i k = _mm256_setzero_si256();
        x = _mm256_min_epu8(x, j); // shave off everything negative
        // now we have a conundrum; we need to exchange
        // the high lane of the first result 
        // with the low lane of the second
        const __m256i low = _mm256_unpacklo_epi8(x, k);
        const __m256i high = _mm256_unpackhi_epi8(x, k);
        _mm256_storeu_si256((__m256i*)dst, _mm256_permute2f128_si256(low, high, 0b00100000));
        _mm256_storeu_si256((__m256i*)dst+16, _mm256_permute2f128_si256(low, high, 0b00110001));
    }
}

void
emulate_splitql (OrcOpcodeExecutor *ex, int offset, int n)
{
  int i;
  orc_union32 * ORC_RESTRICT ptr0;
  orc_union32 * ORC_RESTRICT ptr1;
  const orc_union64 * ORC_RESTRICT ptr4;
  orc_union64 var32;
  orc_union32 var33;
  orc_union32 var34;

  ptr0 = (orc_union32 *)ex->dest_ptrs[0];
  ptr1 = (orc_union32 *)ex->dest_ptrs[1];
  ptr4 = (orc_union64 *)ex->src_ptrs[0];


  for (i = 0; i < n; i++) {
    /* 0: loadq */
    var32 = ptr4[i];
    /* 1: splitql */
    {
       orc_union64 _src;
       _src.i = var32.i;
       var33.i = _src.x2[1];
       var34.i = _src.x2[0];
    }
    /* 2: storel */
    ptr0[i] = var33;
    /* 3: storel */
    ptr1[i] = var34;
  }

}

void
emulate_ldreslinl (OrcOpcodeExecutor *ex, int offset, int n)
{
  int i;
  orc_union32 * ORC_RESTRICT ptr0;
  const orc_union32 * ORC_RESTRICT ptr4;
  orc_union32 var32;

  ptr0 = (orc_union32 *)ex->dest_ptrs[0];
  ptr4 = (orc_union32 *)ex->src_ptrs[0];


  for (i = 0; i < n; i++) {
    /* 0: ldreslinl */
    {
    auto tmp = ((orc_union64 *)(ex->src_ptrs[1]))->i + (offset + i) * ((orc_union64 *)(ex->src_ptrs[2]))->i;
    orc_union32 a = ptr4[tmp>>16];
    orc_union32 b = ptr4[(tmp>>16)+1];
    // dest <- a * (256 - tmp >> 8)  + b * (tmp >> 8)
    // dest <- a * 256 - a * (tmp >> 8) + b * (tmp >> 8)
    // dest <- a * 256 + (b - a) * (tmp >> 8)
    auto x1 = ((orc_uint8)a.x4[0] * (256-((tmp>>8)&0xff)) + (orc_uint8)b.x4[0] * ((tmp>>8)&0xff))>>8;
    var32.x4[0] = x1;
    var32.x4[1] = ((orc_uint8)a.x4[1] * (256-((tmp>>8)&0xff)) + (orc_uint8)b.x4[1] * ((tmp>>8)&0xff))>>8;
    var32.x4[2] = ((orc_uint8)a.x4[2] * (256-((tmp>>8)&0xff)) + (orc_uint8)b.x4[2] * ((tmp>>8)&0xff))>>8;
    var32.x4[3] = ((orc_uint8)a.x4[3] * (256-((tmp>>8)&0xff)) + (orc_uint8)b.x4[3] * ((tmp>>8)&0xff))>>8;

    auto ab = _mm_loadl_epi64((const __m128i*)(ptr4 + (tmp >> 16)));
    // sign extend channel to 16 bits
    const auto zero = _mm_setzero_si128();
    ab = _mm_unpacklo_epi8(ab, zero);
    const auto b_only = _mm_shuffle_epi32(ab, 0b11101110);
    // codegen'd as accessing the upper nibble of the low word
    auto tmpv = _mm_set1_epi16(tmp);
    tmpv = _mm_srli_epi16(tmpv, 8);
    // end of cursed codegen
    const auto factor = _mm_mullo_epi16 (_mm_sub_epi16(b_only, ab), tmpv);
    const auto result = _mm_packs_epi16(_mm_srai_epi16(factor, 8), zero);

    var32.i = _mm_cvtsi128_si32 (result);
    }
    /* 1: storel */
    ptr0[i] = var32;
  }

}
