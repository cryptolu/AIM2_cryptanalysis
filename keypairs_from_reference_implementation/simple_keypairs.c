#include <stdio.h>
#include <stdint.h>
#include "params.h"
#include "aim2.h"
#include "rng.h"

static void print_hex(const char *label, const uint8_t *buf, size_t len)
{
  size_t i;
  printf("%s", label);
  for (i = 0; i < len; i++)
  {
    printf("%02x", buf[i]);
  }
  printf("\n");
}

int main(void)
{
  uint8_t sk[AIM2_NUM_BYTES_FIELD];
  uint8_t iv[AIM2_IV_SIZE];
  uint8_t pk_ct[AIM2_NUM_BYTES_FIELD];

  randombytes(sk, AIM2_NUM_BYTES_FIELD);
  print_hex("sk = ", sk, sizeof(sk));

  for (int i = 0; i < 100; i++)
  {
    randombytes(iv, AIM2_IV_SIZE);

    #ifndef PARAMS
      // NIST
      aim2(sk, iv, pk_ct);
    #else
      // New API
      aim2(pk_ct, sk, iv);
    #endif

    print_hex("iv = ", iv, sizeof(iv));
    print_hex("pk = ", pk_ct, sizeof(pk_ct));
    printf("\n");
  }

  return 0;
}

// For NIST package
// $ cd aimer128f && gcc simple_keypairs.c aim2.c field.c hash.c shake/*.c rng.c aes.c && ./a.out >keypairs128.txt
// $ cd aimer192f && gcc simple_keypairs.c aim2.c field.c hash.c shake/*.c rng.c aes.c && ./a.out >keypairs192.txt
// $ cd aimer256f && gcc simple_keypairs.c aim2.c field.c hash.c shake/*.c rng.c aes.c && ./a.out >keypairs256.txt

// For  https://github.com/samsungsds-opensource/AIMer/blob/main/AIMer-package.zip
// $ gcc -DPARAMS=128f simple_keypairs.c aim2.c field128.c hash.c -Icommon common/*.c && ./a.out >keypairs128.txt
// $ gcc -DPARAMS=192f simple_keypairs.c aim2.c field192.c hash.c -Icommon common/*.c && ./a.out >keypairs192.txt
// $ gcc -DPARAMS=256f simple_keypairs.c aim2.c field256.c hash.c -Icommon common/*.c && ./a.out >keypairs256.txt