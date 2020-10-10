#ifndef _RS_H
#define _RS_H

#include <stdint.h>

struct rs_codec {
    int mm;
    int nn;
    uint16_t *alpha_to;
    uint16_t *index_of;
    uint16_t *genpoly;
    int nroots;
    int fcr;
    int prim;
    int iprim;
    int gfpoly;
};

struct rs_control {
    struct rs_codec *codec;
    uint16_t buffers[];
};


static inline int rs_modnn(struct rs_codec *rs, int x)
{
    while (x >= rs->nn) {
        x -= rs->nn;
        x = (x >> rs->mm) + (x & rs->nn);
    }
    return x;
}
struct rs_control *init_rs(int nroots);
int encode_rs8(struct rs_control *rsc, uint8_t *data, int len, uint16_t *par);
int decode_rs8(struct rs_control *rsc, uint8_t *data, int len, uint16_t *par);

#endif
