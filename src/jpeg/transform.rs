use simd::{i16x8, i32x4, u8x16};

// The forward dct's output coefficients are scaled by 8
// The inverse dct's output samples are clamped to the range [0, 255]
/*
idct and fdct are Rust translations of jfdctint.c and jidctint.c from the
Independent JPEG Group's libjpeg version 9a
obtained from http://www.ijg.org/files/jpegsr9a.zip
They come with the following conditions of ditstribution and use:

	In plain English:

	1. We don't promise that this software works.  (But if you find any bugs,
		please let us know!)
	2. You can use this software for whatever you want.  You don't have to pay us.
	3. You may not pretend that you wrote this software.  If you use it in a
	   program, you must acknowledge somewhere in your documentation that
	   you've used the IJG code.

	In legalese:

	The authors make NO WARRANTY or representation, either express or implied,
	with respect to this software, its quality, accuracy, merchantability, or
	fitness for a particular purpose.  This software is provided "AS IS", and you,
	its user, assume the entire risk as to its quality and accuracy.

	This software is copyright (C) 1991-2014, Thomas G. Lane, Guido Vollbeding.
	All Rights Reserved except as specified below.

	Permission is hereby granted to use, copy, modify, and distribute this
	software (or portions thereof) for any purpose, without fee, subject to these
	conditions:
	(1) If any part of the source code for this software is distributed, then this
	README file must be included, with this copyright and no-warranty notice
	unaltered; and any additions, deletions, or changes to the original files
	must be clearly indicated in accompanying documentation.
	(2) If only executable code is distributed, then the accompanying
	documentation must state that "this software is based in part on the work of
	the Independent JPEG Group".
	(3) Permission for use of this software is granted only if the user accepts
	full responsibility for any undesirable consequences; the authors accept
	NO LIABILITY for damages of any kind.

	These conditions apply to any software derived from or based on the IJG code,
	not just to the unmodified library.  If you use our work, you ought to
	acknowledge us.

	Permission is NOT granted for the use of any IJG author's name or company name
	in advertising or publicity relating to this software or products derived from
	it.  This software may be referred to only as "the Independent JPEG Group's
	software".

	We specifically permit and encourage the use of this software as the basis of
	commercial products, provided that all warranty or liability claims are
	assumed by the product vendor.
*/

const CONST_BITS: i32 = 13;
const PASS1_BITS: i32 = 2;

const FIX_0_298631336: i32 = 2446;
const FIX_0_390180644: i32 = 3196;
const FIX_0_541196100: i32 = 4433;
const FIX_0_765366865: i32 = 6270;
const FIX_0_899976223: i32 = 7373;
const FIX_1_175875602: i32 = 9633;
const FIX_1_501321110: i32 = 12299;
const FIX_1_847759065: i32 = 15137;
const FIX_1_961570560: i32 = 16069;
const FIX_2_053119869: i32 = 16819;
const FIX_2_562915447: i32 = 20995;
const FIX_3_072711026: i32 = 25172;

pub fn fdct(samples: &[u8], coeffs: &mut [i32]) {
    // Pass 1: process rows.
    // Results are scaled by sqrt(8) compared to a true DCT
    // furthermore we scale the results by 2**PASS1_BITS
    for y in (0usize..8) {
        let y0 = y * 8;

        // Even part
        let t0 = samples[y0 + 0] as i32 + samples[y0 + 7] as i32;
        let t1 = samples[y0 + 1] as i32 + samples[y0 + 6] as i32;
        let t2 = samples[y0 + 2] as i32 + samples[y0 + 5] as i32;
        let t3 = samples[y0 + 3] as i32 + samples[y0 + 4] as i32;

        let t10 = t0 + t3;
        let t12 = t0 - t3;
        let t11 = t1 + t2;
        let t13 = t1 - t2;

        let t0 = samples[y0 + 0] as i32 - samples[y0 + 7] as i32;
        let t1 = samples[y0 + 1] as i32 - samples[y0 + 6] as i32;
        let t2 = samples[y0 + 2] as i32 - samples[y0 + 5] as i32;
        let t3 = samples[y0 + 3] as i32 - samples[y0 + 4] as i32;

        // Apply unsigned -> signed conversion
        coeffs[y0 + 0] = (t10 + t11 - 8 * 128) << PASS1_BITS as usize;
        coeffs[y0 + 4] = (t10 - t11) << PASS1_BITS as usize;

        let mut z1 = (t12 + t13) * FIX_0_541196100;
        // Add fudge factor here for final descale
        z1 += 1 << (CONST_BITS - PASS1_BITS - 1) as usize;

        coeffs[y0 + 2] = (z1 + t12 * FIX_0_765366865) >> (CONST_BITS - PASS1_BITS) as usize;
        coeffs[y0 + 6] = (z1 - t13 * FIX_1_847759065) >> (CONST_BITS - PASS1_BITS) as usize;

        // Odd part
        let t12 = t0 + t2;
        let t13 = t1 + t3;

        let mut z1 = (t12 + t13) * FIX_1_175875602;
        // Add fudge factor here for final descale
        z1 += 1 << (CONST_BITS - PASS1_BITS - 1) as usize;

        let mut t12 = t12 * (-FIX_0_390180644);
        let mut t13 = t13 * (-FIX_1_961570560);
        t12 += z1;
        t13 += z1;

        let z1 = (t0 + t3) * (-FIX_0_899976223);
        let mut t0 = t0 * FIX_1_501321110;
        let mut t3 = t3 * FIX_0_298631336;
        t0 += z1 + t12;
        t3 += z1 + t13;

        let z1 = (t1 + t2) * (-FIX_2_562915447);
        let mut t1 = t1 * FIX_3_072711026;
        let mut t2 = t2 * FIX_2_053119869;
        t1 += z1 + t13;
        t2 += z1 + t12;

        coeffs[y0 + 1] = t0 >> (CONST_BITS - PASS1_BITS) as usize;
        coeffs[y0 + 3] = t1 >> (CONST_BITS - PASS1_BITS) as usize;
        coeffs[y0 + 5] = t2 >> (CONST_BITS - PASS1_BITS) as usize;
        coeffs[y0 + 7] = t3 >> (CONST_BITS - PASS1_BITS) as usize;
    }

    // Pass 2: process columns
    // We remove the PASS1_BITS scaling but leave the results scaled up an
    // overall factor of 8
    for x in (0usize..8).rev() {
        // Even part
        let t0 = coeffs[x + 8 * 0] + coeffs[x + 8 * 7];
        let t1 = coeffs[x + 8 * 1] + coeffs[x + 8 * 6];
        let t2 = coeffs[x + 8 * 2] + coeffs[x + 8 * 5];
        let t3 = coeffs[x + 8 * 3] + coeffs[x + 8 * 4];

        // Add fudge factor here for final descale
        let t10 = t0 + t3 + (1 << (PASS1_BITS - 1) as usize);
        let t12 = t0 - t3;
        let t11 = t1 + t2;
        let t13 = t1 - t2;

        let t0 = coeffs[x + 8 * 0] - coeffs[x + 8 * 7];
        let t1 = coeffs[x + 8 * 1] - coeffs[x + 8 * 6];
        let t2 = coeffs[x + 8 * 2] - coeffs[x + 8 * 5];
        let t3 = coeffs[x + 8 * 3] - coeffs[x + 8 * 4];

        coeffs[x + 8 * 0] = (t10 + t11) >> PASS1_BITS as usize;
        coeffs[x + 8 * 4] = (t10 - t11) >> PASS1_BITS as usize;

        let mut z1 = (t12 + t13) * FIX_0_541196100;
        // Add fudge factor here for final descale
        z1 += 1 << (CONST_BITS + PASS1_BITS - 1) as usize;

        coeffs[x + 8 * 2] = (z1 + t12 * FIX_0_765366865) >> (CONST_BITS + PASS1_BITS) as usize;
        coeffs[x + 8 * 6] = (z1 - t13 * FIX_1_847759065) >> (CONST_BITS + PASS1_BITS) as usize;

        // Odd part
        let t12 = t0 + t2;
        let t13 = t1 + t3;

        let mut z1 = (t12 + t13) * FIX_1_175875602;
        // Add fudge factor here for final descale
        z1 += 1 << (CONST_BITS - PASS1_BITS - 1) as usize;

        let mut t12 = t12 * (-FIX_0_390180644);
        let mut t13 = t13 * (-FIX_1_961570560);
        t12 += z1;
        t13 += z1;

        let z1 = (t0 + t3) * (-FIX_0_899976223);
        let mut t0 = t0 * FIX_1_501321110;
        let mut t3 = t3 * FIX_0_298631336;
        t0 += z1 + t12;
        t3 += z1 + t13;

        let z1 = (t1 + t2) * (-FIX_2_562915447);
        let mut t1 = t1 * FIX_3_072711026;
        let mut t2 = t2 * FIX_2_053119869;
        t1 += z1 + t13;
        t2 += z1 + t12;

        coeffs[x + 8 * 1] = t0 >> (CONST_BITS + PASS1_BITS) as usize;
        coeffs[x + 8 * 3] = t1 >> (CONST_BITS + PASS1_BITS) as usize;
        coeffs[x + 8 * 5] = t2 >> (CONST_BITS + PASS1_BITS) as usize;
        coeffs[x + 8 * 7] = t3 >> (CONST_BITS + PASS1_BITS) as usize;
    }
}

#[inline(always)]
fn interleave(x: i16x8, y: i16x8) -> (i16x8, i16x8) {
    (i16x8::new(x.extract(0), y.extract(0), x.extract(1), y.extract(1),
                x.extract(2), y.extract(2), x.extract(3), y.extract(3)),
     i16x8::new(x.extract(4), y.extract(4), x.extract(5), y.extract(5),
                x.extract(6), y.extract(6), x.extract(7), y.extract(7)))
}

fn as_i32(x: i16x8) -> (i32x4, i32x4) {
    let (a, b): (i32x4, i32x4) = unsafe {::std::mem::transmute(interleave(x, i16x8::splat(0)))};
    (a << 16 >> 16, b << 16 >> 16)
}
pub fn transpose(matrix: &[i16x8; 8]) -> (i16x8, i16x8, i16x8, i16x8,
                                          i16x8, i16x8, i16x8, i16x8)
{
    let (a, b) = interleave(matrix[0], matrix[4]);
    let (c, d) = interleave(matrix[1], matrix[5]);
    let (e, f) = interleave(matrix[2], matrix[6]);
    let (g, h) = interleave(matrix[3], matrix[7]);

    let (a1, b1) = interleave(a, e);
    let (c1, d1) = interleave(b, f);
    let (e1, f1) = interleave(c, g);
    let (g1, h1) = interleave(d, h);

    let (a2, b2) = interleave(a1, e1);
    let (c2, d2) = interleave(b1, f1);
    let (e2, f2) = interleave(c1, g1);
    let (g2, h2) = interleave(d1, h1);

    (a2, b2, c2, d2, e2, f2, g2, h2)
}

#[test]
fn transpose_() {
    let (a, b, c, d, e, f, g, h) = transpose(&[i16x8::splat(1),
                                               i16x8::splat(2),
                                               i16x8::splat(3),
                                               i16x8::splat(4),
                                               i16x8::splat(5),
                                               i16x8::splat(6),
                                               i16x8::splat(7),
                                               i16x8::splat(8)]);

    let ret = i16x8::new(1, 2, 3, 4, 5, 6, 7, 8);
    assert!(a.eq(ret).all());
    assert!(b.eq(ret).all());
    assert!(c.eq(ret).all());
    assert!(d.eq(ret).all());
    assert!(e.eq(ret).all());
    assert!(f.eq(ret).all());
    assert!(g.eq(ret).all());
    assert!(h.eq(ret).all());
}

fn transpose_bytes(x: u8x16, y: u8x16, z: u8x16, w: u8x16) -> (u8x16, u8x16, u8x16, u8x16) {
    #[inline(always)]
    fn interleave_8(x: u8x16, y: u8x16) -> (u8x16, u8x16) {
        (u8x16::new(x.extract(0), y.extract(0), x.extract(1), y.extract(1),
                    x.extract(2), y.extract(2), x.extract(3), y.extract(3),
                    x.extract(4), y.extract(4), x.extract(5), y.extract(5),
                    x.extract(6), y.extract(6), x.extract(7), y.extract(7)),
         u8x16::new(x.extract(8 + 0), y.extract(8 + 0), x.extract(8 + 1), y.extract(8 + 1),
                    x.extract(8 + 2), y.extract(8 + 2), x.extract(8 + 3), y.extract(8 + 3),
                    x.extract(8 + 4), y.extract(8 + 4), x.extract(8 + 5), y.extract(8 + 5),
                    x.extract(8 + 6), y.extract(8 + 6), x.extract(8 + 7), y.extract(8 + 7)))

    }

    let (a, b) = interleave_8(x, z);
    let (c, d) = interleave_8(y, w);
    let (a1, b1) = interleave_8(a, c);
    let (c1, d1) = interleave_8(b, d);
    let (a2, b2) = interleave_8(a1, c1);
    let (c2, d2) = interleave_8(b1, d1);
    (a2, b2, c2, d2)
}

#[test]
fn transpose_bytes_() {
    let (a, b, c, d) = transpose_bytes(u8x16::new(1, 1, 1, 1, 1, 1, 1, 1,
                                                  2, 2, 2, 2, 2, 2, 2, 2),
                                       u8x16::new(3, 3, 3, 3, 3, 3, 3, 3,
                                                  4, 4, 4, 4, 4, 4, 4, 4),
                                       u8x16::new(5, 5, 5, 5, 5, 5, 5, 5,
                                                  6, 6, 6, 6, 6, 6, 6, 6),
                                       u8x16::new(7, 7, 7, 7, 7, 7, 7, 7,
                                                  8, 8, 8, 8, 8, 8, 8, 8));
    let e = u8x16::new(1, 2, 3, 4, 5, 6, 7, 8,
                       1, 2, 3, 4, 5, 6, 7, 8);
    assert!(a.eq(e).all());
    assert!(b.eq(e).all());
    assert!(c.eq(e).all());
    assert!(d.eq(e).all());
}

#[inline(never)]
pub fn idct(coeffs: &[i16x8; 8], samples: &mut [u8]) {
    use simd::x86::sse2::*;

    let mut tmp = [i16x8::splat(0); 64 / 8];

    //let (c0, c1, c2, c3, c4, c5, c6, c7) = transpose(coeffs);
    let (c0, c1, c2, c3, c4, c5, c6, c7) = (coeffs[0], coeffs[1], coeffs[2], coeffs[3],
                                            coeffs[4], coeffs[5], coeffs[6], coeffs[7]);
    let zero = i16x8::splat(0);

    #[inline(always)]
    fn cf(x: i32, y: i32, z: i32) -> (i16x8, i16x8) {
        let (x, y, z) = (x as i16, y as i16, z as i16);
        let xy = x + y;
        let xz = x + z;
        (i16x8::new(xy, x, xy, x, xy, x, xy, x),
         i16x8::new(x, xz, x, xz, x, xz, x, xz))
    }
    let (f05_07_05, f05_05_m18) = cf(FIX_0_541196100, FIX_0_765366865, -FIX_1_847759065);
    let (f11_m19_11, f11_11_m03) = cf(FIX_1_175875602, -FIX_1_961570560, -FIX_0_390180644);
    let (fm08_02_m08, fm08_m08_15) = cf(-FIX_0_899976223, FIX_0_298631336, FIX_1_501321110);
    let (fm25_20_m25, fm25_m25_30) = cf(-FIX_2_562915447, FIX_2_053119869, FIX_3_072711026);

    if (c0 | c1 | c2 | c3 | c4 | c5 | c6 | c7).eq(zero).all() {
        tmp = [c0 << PASS1_BITS; 64 / 8];
    } else {
        // Even part: reverse the even part of the forward DCT
        let z2 = c2;
        let z3 = c6;

        //println!("z2: {:?} z3: {:?}", z2, z3);
        /*
        let z1 = (z2 as i32 + z3 as i32) * FIX_0_541196100;
        let t2 = z1 + z2 as i32 * FIX_0_765366865;
        let t3 = z1 - z3 as i32 * FIX_1_847759065;
        */
        //let t2 = z2 as i32 * (FIX_0_765366865 + FIX_0_541196100) + z3 as i32 * FIX_0_541196100;
        //let t3 = z2 as i32 * FIX_0_541196100 + z3 as i32 * (FIX_0_541196100 - FIX_1_847759065);
        let (z2z3l, z2z3h) = interleave(z2, z3);
        let t2l = z2z3l.madd(f05_07_05);
        let t2h = z2z3h.madd(f05_07_05);
        let t3l = z2z3l.madd(f05_05_m18);
        let t3h = z2z3h.madd(f05_05_m18);
        //println!("\tt2: {:?} {:?} t3: {:?} {:?}", t2l, t2h, t3l, t3h);

        /*
        let mut z2 = coeffs[x + 8 * 0] as i32;
        let mut z3 = coeffs[x + 8 * 4] as i32;
        z2 <<= CONST_BITS as usize;
        */
        let (mut z2l, mut z2h) = as_i32(c0);
        z2l = z2l << CONST_BITS;
        z2h = z2h << CONST_BITS;
        let (mut z3l, mut z3h) = as_i32(c4);
        z3l = z3l << CONST_BITS;
        z3h = z3h << CONST_BITS;

        z2l = z2l + i32x4::splat(1 << (CONST_BITS - PASS1_BITS - 1));
        z2h = z2h + i32x4::splat(1 << (CONST_BITS - PASS1_BITS - 1));

        //let t0 = z2 as i32 + z3 as i32;
        //let t1 = z2 as i32 - z3 as i32;
        let t0l = z2l + z3l;
        let t0h = z2h + z3h;
        let t1l = z2l - z3l;
        let t1h = z2h - z3h;

        //let t10 = t0 + t2;
        //let t13 = t0 - t2;
        //let t11 = t1 + t3;
        //let t12 = t1 - t3;
        let t10l = t0l + t2l;
        let t10h = t0h + t2h;
        let t13l = t0l - t2l;
        let t13h = t0h - t2h;
        let t11l = t1l + t3l;
        let t11h = t1h + t3h;
        let t12l = t1l - t3l;
        let t12h = t1h - t3h;

        let t0 = c7;
        let t1 = c5;
        let t2 = c3;
        let t3 = c1;

        let z2 = t0 + t2;
        let z3 = t1 + t3;

        /*
        let z1 = (z2 as i32 + z3 as i32) * FIX_1_175875602;
        let mut z2 = z2 as i32 * (-FIX_1_961570560);
        let mut z3 = z3 as i32 * (-FIX_0_390180644);
        z2 += z1;
        z3 += z1;
         */
        //let z2_ = z2 as i32 * (FIX_1_175875602 - FIX_1_961570560) + z3 as i32 * FIX_1_175875602;
        //let z3 = z2 as i32 * FIX_1_175875602 + z3 as i32 * (FIX_1_175875602 - FIX_0_390180644);
        //let z2 = z2_;
        let (z2z3l, z2z3h) = interleave(z2, z3);
        let z2l = z2z3l.madd(f11_m19_11);
        let z2h = z2z3h.madd(f11_m19_11);
        let z3l = z2z3l.madd(f11_11_m03);
        let z3h = z2z3h.madd(f11_11_m03);

        /*
        let z1 = (t0 as i32 + t3 as i32) * (-FIX_0_899976223);
        let mut t0 = t0 as i32 * FIX_0_298631336;
        let mut t3 = t3 as i32 * FIX_1_501321110;
        t0 += z1 + z2;
        t3 += z1 + z3;
        */
        //let t0_ = t0 as i32 * (-FIX_0_899976223 + FIX_0_298631336) + t3 as i32 * (-FIX_0_899976223) + z2 as i32;
        //let t3 = t0 as i32 * (-FIX_0_899976223) + t3 as i32 * (-FIX_0_899976223 + FIX_1_501321110) + z3 as i32;
        //let t0 = t0_;
        let (t0t3l, t0t3h) = interleave(t0, t3);
        let t0l = t0t3l.madd(fm08_02_m08) + z2l;
        let t0h = t0t3h.madd(fm08_02_m08) + z2h;
        let t3l = t0t3l.madd(fm08_m08_15) + z3l;
        let t3h = t0t3h.madd(fm08_m08_15) + z3h;
        /*
        let z1 = (t1 as i32 + t2 as i32) * (-FIX_2_562915447);
        let mut t1 = t1 as i32 * FIX_2_053119869;
        let mut t2 = t2 as i32 * FIX_3_072711026;
        t1 += z1 + z3 as i32;
        t2 += z1 + z2 as i32;
         */
        //let t1_ = t1 as i32 * (-FIX_2_562915447 + FIX_2_053119869) + t2 as i32 * (-FIX_2_562915447) + z3 as i32;
        //let t2 = t1 as i32 * (-FIX_2_562915447) + t2 as i32 * (-FIX_2_562915447 + FIX_3_072711026) + z2 as i32;
        //let t1 = t1_;
        let (t1t2l, t1t2h) = interleave(t1, t2);
        let t1l = t1t2l.madd(fm25_20_m25) + z3l;
        let t1h = t1t2h.madd(fm25_20_m25) + z3h;
        let t2l = t1t2l.madd(fm25_m25_30) + z2l;
        let t2h = t1t2h.madd(fm25_m25_30) + z2h;

        fn combine(xl: i32x4, xh: i32x4,
                   yl: i32x4, yh: i32x4) -> (i16x8, i16x8) {
            const SHIFT: i32 = CONST_BITS - PASS1_BITS;

            (((xl + yl) >> SHIFT).packs((xh + yh) >> SHIFT),
             ((xl - yl) >> SHIFT).packs((xh - yh) >> SHIFT))

        }

        //println!("\t{:?} {:?} {:?} {:?} {:?} {:?} {:?} {:?}", t10l, t11l, t12l, t13l, t0l, t1l, t2l, t3l);
        //println!("\t{:?} {:?} {:?} {:?} {:?} {:?} {:?} {:?}", t10h, t11h, t12h, t13h, t0h, t1h, t2h, t3h);
        let (a, b) = combine(t10l, t10h, t3l, t3h);
        tmp[0] = a;
        tmp[7] = b;
        let (a, b) = combine(t11l, t11h, t2l, t2h);
        tmp[1] = a;
        tmp[6] = b;
        let (a, b) = combine(t12l, t12h, t1l, t1h);
        tmp[2] = a;
        tmp[5] = b;
        let (a, b) = combine(t13l, t13h, t0l, t0h);
        tmp[3] = a;
        tmp[4] = b;
    }

    //println!("{:?}", &tmp);

    let (y0, y1, y2, y3, y4, y5, y6, y7) = transpose(&tmp);

    let z2 = y2;
    let z3 = y6;

    /*
    let z1 = (z2 as i32 + z3 as i32) * FIX_0_541196100;
    let t2 = z1 + z2 as i32 * FIX_0_765366865;
    let t3 = z1 - z3 as i32 * FIX_1_847759065;
     */
    //let t2 = z2 as i32 * (FIX_0_541196100 + FIX_0_765366865) + z3 as i32 * FIX_0_541196100;
    //let t3 = z2 as i32 * (FIX_0_541196100) + z3 as i32 * (FIX_0_541196100 - FIX_1_847759065);
    let (z2z3l, z2z3h) = interleave(z2, z3);
    let t2l = z2z3l.madd(f05_07_05);
    let t2h = z2z3h.madd(f05_07_05);
    let t3l = z2z3l.madd(f05_05_m18);
    let t3h = z2z3h.madd(f05_05_m18);

    let z2 = y0 + i16x8::splat(1 << (PASS1_BITS + 2));
    let z3 = y4;


    //let t0 = (z2 + z3) << CONST_BITS as usize;
    //let t1 = (z2 - z3) << CONST_BITS as usize;
    let (mut t0l, mut t0h) = as_i32(z2 + z3);
    t0l = t0l << CONST_BITS;
    t0h = t0h << CONST_BITS;
    let (mut t1l, mut t1h) = as_i32(z2 - z3);
    t1l = t1l << CONST_BITS;
    t1h = t1h << CONST_BITS;

    let t10l = t0l + t2l;
    let t10h = t0h + t2h;
    let t13l = t0l - t2l;
    let t13h = t0h - t2h;
    let t11l = t1l + t3l;
    let t11h = t1h + t3h;
    let t12l = t1l - t3l;
    let t12h = t1h - t3h;

    let t0 = y7;
    let t1 = y5;
    let t2 = y3;
    let t3 = y1;

    let z2 = t0 + t2;
    let z3 = t1 + t3;

    /*
    let z1 = (z2 as i32 + z3 as i32) * FIX_1_175875602;
    let mut z2 = z2 as i32 * (-FIX_1_961570560);
    let mut z3 = z3 as i32 * (-FIX_0_390180644);
    z2 += z1;
    z3 += z1;
     */
    //let z2_ = z2 as i32 * (FIX_1_175875602 - FIX_1_961570560) + z3 as i32 * FIX_1_175875602;
    //let z3 = z2 as i32 * FIX_1_175875602 + z3 as i32 * (FIX_1_175875602 - FIX_0_390180644);
    //let z2 = z2_;
    let (z2z3l, z2z3h) = interleave(z2, z3);
    let z2l = z2z3l.madd(f11_m19_11);
    let z2h = z2z3h.madd(f11_m19_11);
    let z3l = z2z3l.madd(f11_11_m03);
    let z3h = z2z3h.madd(f11_11_m03);

    /*
    let z1 = (t0 as i32 + t3 as i32) * (-FIX_0_899976223);
    let mut t0 = t0 as i32 * FIX_0_298631336;
    let mut t3 = t3 as i32 * FIX_1_501321110;
    t0 += z1 + z2 as i32;
    t3 += z1 + z3 as i32;
     */
    //let t0_ = t0 as i32 * (-FIX_0_899976223 + FIX_0_298631336) + t3 as i32 * (-FIX_0_899976223) + z2;
    //let t3 = t0 as i32 * (-FIX_0_899976223) + t3 as i32 * (-FIX_0_899976223 + FIX_1_501321110) + z3;
    //let t0 = t0_;
    let (t0t3l, t0t3h) = interleave(t0, t3);
    let t0l = t0t3l.madd(fm08_02_m08) + z2l;
    let t0h = t0t3h.madd(fm08_02_m08) + z2h;
    let t3l = t0t3l.madd(fm08_m08_15) + z3l;
    let t3h = t0t3h.madd(fm08_m08_15) + z3h;

    /*
    let z1 = (t1 as i32 + t2 as i32) * (-FIX_2_562915447);
    let mut t1 = t1 as i32 * FIX_2_053119869;
    let mut t2 = t2 as i32 * FIX_3_072711026;
    t1 += z1 + z3 as i32;
    t2 += z1 + z2 as i32;
     */
    //let t1_ = t1 as i32 * (-FIX_2_562915447 + FIX_2_053119869) + t2 as i32 * (-FIX_2_562915447) + z3;
    //let t2 = t1 as i32 * (-FIX_2_562915447) + t2 as i32 * (-FIX_2_562915447 + FIX_3_072711026) + z2;
    //let t1 = t1_;
    let (t1t2l, t1t2h) = interleave(t1, t2);
    let t1l = t1t2l.madd(fm25_20_m25) + z3l;
    let t1h = t1t2h.madd(fm25_20_m25) + z3h;
    let t2l = t1t2l.madd(fm25_m25_30) + z2l;
    let t2h = t1t2h.madd(fm25_m25_30) + z2h;

    fn sp(a: i32x4, b: i32x4) -> i16x8 {
        const SHIFT: i32 = CONST_BITS + PASS1_BITS + 3;
        (a >> SHIFT).packs(b >> SHIFT)
    }
    fn ps(a: i16x8, b: i16x8) -> u8x16 {
        let bytes: u8x16 = unsafe {::std::mem::transmute(a.packs(b))};
        bytes + u8x16::splat(128)
    }

    let a0 = sp(t10l + t3l, t10h + t3h);
    let a7 = sp(t10l - t3l, t10h - t3h);
    let a1 = sp(t11l + t2l, t11h + t2h);
    let a6 = sp(t11l - t2l, t11h - t2h);
    let a2 = sp(t12l + t1l, t12h + t1h);
    let a5 = sp(t12l - t1l, t12h - t1h);
    let a3 = sp(t13l + t0l, t13h + t0h);
    let a4 = sp(t13l - t0l, t13h - t0h);

    let x = ps(a0, a1);
    let y = ps(a2, a3);
    let z = ps(a4, a5);
    let w = ps(a6, a7);

    let (a, b, c, d) = transpose_bytes(x, y, z, w);

    a.store(samples, 0 * 16);
    b.store(samples, 1 * 16);
    c.store(samples, 2 * 16);
    d.store(samples, 3 * 16);
}
