#include "rs.h"
#include <assert.h>
#include <stdlib.h>
#include <string.h>

#ifndef min
#define min(x, y) ((x) < (y) ? (x) : (y))
#endif

enum {
    RS_DECODE_LAMBDA,
    RS_DECODE_SYN,
    RS_DECODE_B,
    RS_DECODE_T,
    RS_DECODE_OMEGA,
    RS_DECODE_ROOT,
    RS_DECODE_REG,
    RS_DECODE_LOC,
    RS_DECODE_NUM_BUFFERS
};
static struct rs_codec *codec_init(int symsize,
                                   int gfpoly,
                                   int fcr,
                                   int prim,
                                   int nroots)
{
    int i, j, sr, root, iprim;
    struct rs_codec *rs;

    rs = malloc(sizeof(*rs));
    if (!rs)
        return NULL;


    rs->mm = symsize;
    rs->nn = (1 << symsize) - 1;
    rs->fcr = fcr;
    rs->prim = prim;
    rs->nroots = nroots;
    rs->gfpoly = gfpoly;

    /* Allocate the arrays */
    rs->alpha_to = malloc((rs->nn + 1) * sizeof(uint16_t));
    if (rs->alpha_to == NULL)
        goto err;

    rs->index_of = malloc((rs->nn + 1) * sizeof(uint16_t));
    if (rs->index_of == NULL)
        goto err;

    rs->genpoly = malloc((rs->nroots + 1) * sizeof(uint16_t));
    if (rs->genpoly == NULL)
        goto err;

    /* Generate Galois field lookup tables */
    rs->index_of[0] = rs->nn; /* log(zero) = -inf */
    rs->alpha_to[rs->nn] = 0; /* alpha**-inf = 0 */
    if (gfpoly) {
        sr = 1;
        for (i = 0; i < rs->nn; i++) {
            rs->index_of[sr] = i;
            rs->alpha_to[i] = sr;
            sr <<= 1;
            if (sr & (1 << symsize))
                sr ^= gfpoly;
            sr &= rs->nn;
        }
    } else {
        goto err;
    }
    /* If it's not primitive, exit */
    if (sr != rs->alpha_to[0])
        goto err;


    /* Find prim-th root of 1, used in decoding */
    for (iprim = 1; (iprim % prim) != 0; iprim += rs->nn)
        ;
    /* prim-th root of 1, index form */
    rs->iprim = iprim / prim;



    /* Form RS code generator polynomial from its roots */
    rs->genpoly[0] = 1;
    for (i = 0, root = fcr * prim; i < nroots; i++, root += prim) {
        rs->genpoly[i + 1] = 1;

        /* Multiply rs->genpoly[] by  @**(root + x) */
        for (j = i; j > 0; j--) {
            if (rs->genpoly[j] != 0) {
                rs->genpoly[j] = rs->genpoly[j - 1] ^
                                 rs->alpha_to[rs_modnn(
                                     rs, rs->index_of[rs->genpoly[j]] + root)];
            } else
                rs->genpoly[j] = rs->genpoly[j - 1];
        }
        /* rs->genpoly[0] can never be zero */
        rs->genpoly[0] =
            rs->alpha_to[rs_modnn(rs, rs->index_of[rs->genpoly[0]] + root)];
    }

    /* convert rs->genpoly[] to index form for quicker encoding */
    for (i = 0; i <= nroots; i++)
        rs->genpoly[i] = rs->index_of[rs->genpoly[i]];

    return rs;

err:
    free(rs->genpoly);
    free(rs->index_of);
    free(rs->alpha_to);
    free(rs);
    return NULL;
}

struct rs_control *init_rs(int nroots)
{
    int symsize = 8;
    int gfpoly = 0x11d;
    int fcr = 0;
    int prim = 1;


    struct rs_control *rs;
    unsigned int bsize;

    /* Sanity checks */
    if (symsize < 1)
        return NULL;
    if (fcr < 0 || fcr >= (1 << symsize))
        return NULL;
    if (prim <= 0 || prim >= (1 << symsize))
        return NULL;
    if (nroots < 0 || nroots >= (1 << symsize))
        return NULL;

    /*
     * The decoder needs buffers in each control struct instance to
     * avoid variable size or large fixed size allocations on
     * stack. Size the buffers to arrays of [nroots + 1].
     */
    bsize = sizeof(uint16_t) * RS_DECODE_NUM_BUFFERS * (nroots + 1);
    rs = malloc(sizeof(*rs) + bsize);
    if (!rs)
        return NULL;

    /* Create a new one */
    rs->codec = codec_init(symsize, gfpoly, fcr, prim, nroots);
    if (!rs->codec) {
        free(rs);
        rs = NULL;
    }

    return rs;
}

int encode_rs8(struct rs_control *rsc, uint8_t *data, int len, uint16_t *par)
{
    struct rs_codec *rs = rsc->codec;
    int i, j, pad;
    int nn = rs->nn;
    int nroots = rs->nroots;
    uint16_t *alpha_to = rs->alpha_to;
    uint16_t *index_of = rs->index_of;
    uint16_t *genpoly = rs->genpoly;
    uint16_t fb;
    uint16_t msk = (uint16_t) rs->nn;

    /* Check length parameter for validity */
    pad = nn - nroots - len;
    if (pad < 0 || pad >= nn)
        return -1;

    for (i = 0; i < len; i++) {
        fb = index_of[(((uint16_t) data[i]) & msk) ^ par[0]];
        /* feedback term is non-zero */
        if (fb != nn) {
            for (j = 1; j < nroots; j++) {
                par[j] ^= alpha_to[rs_modnn(rs, fb + genpoly[nroots - j])];
            }
        }
        /* Shift */
        memmove(&par[0], &par[1], sizeof(uint16_t) * (nroots - 1));
        if (fb != nn) {
            par[nroots - 1] = alpha_to[rs_modnn(rs, fb + genpoly[0])];
        } else {
            par[nroots - 1] = 0;
        }
    }

    return 0;
}


int decode_rs8(struct rs_control *rsc, uint8_t *data, int len, uint16_t *par)
{
    struct rs_codec *rs = rsc->codec;
    int deg_lambda, el, deg_omega;
    int i, j, r, k, pad;
    int nn = rs->nn;
    int nroots = rs->nroots;
    int fcr = rs->fcr;
    int prim = rs->prim;
    int iprim = rs->iprim;
    uint16_t *alpha_to = rs->alpha_to;
    uint16_t *index_of = rs->index_of;
    uint16_t q, tmp, num1, num2, den, discr_r, syn_error;
    int count = 0;
    int num_corrected;
    uint16_t msk = (uint16_t) rs->nn;

    /*
     * The decoder buffers are in the rs control struct. They are
     * arrays sized [nroots + 1]
     */
    uint16_t *lambda = rsc->buffers + RS_DECODE_LAMBDA * (nroots + 1);
    uint16_t *syn = rsc->buffers + RS_DECODE_SYN * (nroots + 1);
    uint16_t *b = rsc->buffers + RS_DECODE_B * (nroots + 1);
    uint16_t *t = rsc->buffers + RS_DECODE_T * (nroots + 1);
    uint16_t *omega = rsc->buffers + RS_DECODE_OMEGA * (nroots + 1);
    uint16_t *root = rsc->buffers + RS_DECODE_ROOT * (nroots + 1);
    uint16_t *reg = rsc->buffers + RS_DECODE_REG * (nroots + 1);
    uint16_t *loc = rsc->buffers + RS_DECODE_LOC * (nroots + 1);

    /* Check length parameter for validity */
    pad = nn - nroots - len;
    assert(!(pad < 0 || pad >= nn - nroots));

    /* form the syndromes; i.e., evaluate data(x) at roots of
     * g(x) */
    for (i = 0; i < nroots; i++)
        syn[i] = ((uint16_t) data[0]) & msk;

    for (j = 1; j < len; j++) {
        for (i = 0; i < nroots; i++) {
            if (syn[i] == 0) {
                syn[i] = ((uint16_t) data[j]) & msk;
            } else {
                syn[i] =
                    (((uint16_t) data[j]) & msk) ^
                    alpha_to[rs_modnn(rs, index_of[syn[i]] + (fcr + i) * prim)];
            }
        }
    }

    for (j = 0; j < nroots; j++) {
        for (i = 0; i < nroots; i++) {
            if (syn[i] == 0) {
                syn[i] = ((uint16_t) par[j]) & msk;
            } else {
                syn[i] =
                    (((uint16_t) par[j]) & msk) ^
                    alpha_to[rs_modnn(rs, index_of[syn[i]] + (fcr + i) * prim)];
            }
        }
    }

    uint16_t *s = syn;

    /* Convert syndromes to index form, checking for nonzero condition */
    syn_error = 0;
    for (i = 0; i < nroots; i++) {
        syn_error |= s[i];
        s[i] = index_of[s[i]];
    }

    if (!syn_error) {
        /* if syndrome is zero, data[] is a codeword and there are no
         * errors to correct. So return data[] unmodified
         */
        return 0;
    }

    memset(&lambda[1], 0, nroots * sizeof(lambda[0]));
    lambda[0] = 1;

    for (i = 0; i < nroots + 1; i++)
        b[i] = index_of[lambda[i]];

    /*
     * Begin Berlekamp-Massey algorithm to determine error+erasure
     * locator polynomial
     */
    r = 0;
    el = 0;
    while (++r <= nroots) { /* r is the step number */
        /* Compute discrepancy at the r-th step in poly-form */
        discr_r = 0;
        for (i = 0; i < r; i++) {
            if ((lambda[i] != 0) && (s[r - i - 1] != nn)) {
                discr_r ^=
                    alpha_to[rs_modnn(rs, index_of[lambda[i]] + s[r - i - 1])];
            }
        }
        discr_r = index_of[discr_r]; /* Index form */
        if (discr_r == nn) {
            /* 2 lines below: B(x) <-- x*B(x) */
            memmove(&b[1], b, nroots * sizeof(b[0]));
            b[0] = nn;
        } else {
            /* 7 lines below: T(x) <-- lambda(x)-discr_r*x*b(x) */
            t[0] = lambda[0];
            for (i = 0; i < nroots; i++) {
                if (b[i] != nn) {
                    t[i + 1] =
                        lambda[i + 1] ^ alpha_to[rs_modnn(rs, discr_r + b[i])];
                } else
                    t[i + 1] = lambda[i + 1];
            }
            if (2 * el <= r - 1) {
                el = r - el;
                /*
                 * 2 lines below: B(x) <-- inv(discr_r) *
                 * lambda(x)
                 */
                for (i = 0; i <= nroots; i++) {
                    b[i] =
                        (lambda[i] == 0)
                            ? nn
                            : rs_modnn(rs, index_of[lambda[i]] - discr_r + nn);
                }
            } else {
                /* 2 lines below: B(x) <-- x*B(x) */
                memmove(&b[1], b, nroots * sizeof(b[0]));
                b[0] = nn;
            }
            memcpy(lambda, t, (nroots + 1) * sizeof(t[0]));
        }
    }

    /* Convert lambda to index form and compute deg(lambda(x)) */
    deg_lambda = 0;
    for (i = 0; i < nroots + 1; i++) {
        lambda[i] = index_of[lambda[i]];
        if (lambda[i] != nn)
            deg_lambda = i;
    }

    if (deg_lambda == 0) {
        /*
         * deg(lambda) is zero even though the syndrome is non-zero
         * => uncorrectable error detected
         */
        return -4;
    }

    /* Find roots of error+erasure locator polynomial by Chien search */
    memcpy(&reg[1], &lambda[1], nroots * sizeof(reg[0]));
    count = 0; /* Number of roots of lambda(x) */
    for (i = 1, k = iprim - 1; i <= nn; i++, k = rs_modnn(rs, k + iprim)) {
        q = 1; /* lambda[0] is always 0 */
        for (j = deg_lambda; j > 0; j--) {
            if (reg[j] != nn) {
                reg[j] = rs_modnn(rs, reg[j] + j);
                q ^= alpha_to[reg[j]];
            }
        }
        if (q != 0)
            continue; /* Not a root */

        if (k < pad) {
            /* Impossible error location. Uncorrectable error. */
            return -3;
        }

        /* store root (index-form) and error location number */
        root[count] = i;
        loc[count] = k;
        /* If we've already found max possible roots,
         * abort the search to save time
         */
        if (++count == deg_lambda)
            break;
    }
    if (deg_lambda != count) {
        /*
         * deg(lambda) unequal to number of roots => uncorrectable
         * error detected
         */
        return -2;
    }
    /*
     * Compute err+eras evaluator poly omega(x) = s(x)*lambda(x) (modulo
     * x**nroots). in index form. Also find deg(omega).
     */
    deg_omega = deg_lambda - 1;
    for (i = 0; i <= deg_omega; i++) {
        tmp = 0;
        for (j = i; j >= 0; j--) {
            if ((s[i - j] != nn) && (lambda[j] != nn))
                tmp ^= alpha_to[rs_modnn(rs, s[i - j] + lambda[j])];
        }
        omega[i] = index_of[tmp];
    }

    /*
     * Compute error values in poly-form. num1 = omega(inv(X(l))), num2 =
     * inv(X(l))**(fcr-1) and den = lambda_pr(inv(X(l))) all in poly-form
     * Note: we reuse the buffer for b to store the correction pattern
     */
    num_corrected = 0;
    for (j = count - 1; j >= 0; j--) {
        num1 = 0;
        for (i = deg_omega; i >= 0; i--) {
            if (omega[i] != nn)
                num1 ^= alpha_to[rs_modnn(rs, omega[i] + i * root[j])];
        }

        if (num1 == 0) {
            /* Nothing to correct at this position */
            b[j] = 0;
            continue;
        }

        num2 = alpha_to[rs_modnn(rs, root[j] * (fcr - 1) + nn)];
        den = 0;

        /* lambda[i+1] for i even is the formal derivative
         * lambda_pr of lambda[i] */
        for (i = min(deg_lambda, nroots - 1) & ~1; i >= 0; i -= 2) {
            if (lambda[i + 1] != nn) {
                den ^= alpha_to[rs_modnn(rs, lambda[i + 1] + i * root[j])];
            }
        }

        b[j] = alpha_to[rs_modnn(
            rs, index_of[num1] + index_of[num2] + nn - index_of[den])];
        num_corrected++;
    }

    /*
     * We compute the syndrome of the 'error' and check that it matches
     * the syndrome of the received word
     */
    for (i = 0; i < nroots; i++) {
        tmp = 0;
        for (j = 0; j < count; j++) {
            if (b[j] == 0)
                continue;

            k = (fcr + i) * prim * (nn - loc[j] - 1);
            tmp ^= alpha_to[rs_modnn(rs, index_of[b[j]] + k)];
        }

        if (tmp != alpha_to[s[i]])
            return -1;
    }

    /* Apply error to data and parity */
    for (i = 0; i < count; i++) {
        if (loc[i] < (nn - nroots))
            data[loc[i] - pad] ^= b[i];
        else
            par[loc[i] - pad - len] ^= b[i];
    }

    return num_corrected;
}
