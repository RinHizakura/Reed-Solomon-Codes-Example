#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include "rs.h"

int main()
{
    static struct rs_control *rs_decoder;
    rs_decoder = init_rs(10);

    uint16_t par[10];
    uint8_t data8[] = {0x40, 0xd2, 0x75, 0x47, 0x76, 0x17, 0x32, 0x06,
                       0x27, 0x26, 0x96, 0xc6, 0xc6, 0x96, 0x70, 0xec};
    /* Initialize the parity buffer */
    memset(par, 0, sizeof(par));
    encode_rs8(rs_decoder, data8, 16, par);

    for (int i = 0; i < 16; i++)
        printf("%x ", data8[i]);
    for (int i = 0; i < 10; i++)
        printf("%x ", par[i]);
    printf("\n");

    data8[0] = 0x50;
    data8[1] = 0x50;
    data8[2] = 0x50;

    int numerr = decode_rs8(rs_decoder, data8, 16, par);
    printf("numerr: %d\n", numerr);

    for (int i = 0; i < 16; i++)
        printf("%x ", data8[i]);
    for (int i = 0; i < 10; i++)
        printf("%x ", par[i]);
    printf("\n");
}
