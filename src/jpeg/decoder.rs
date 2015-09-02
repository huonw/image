use std::cmp;
use std::io::Read;
use std::default::Default;
use std::collections::HashMap;
use byteorder::{ReadBytesExt, BigEndian};
use num::range_step;

use simd::{i16x8, u8x16};

use color;
use super::transform;

use super::entropy:: {
    HuffTable,
    HuffDecoder,
    derive_tables,
};

use image;
use image::ImageResult;
use image::ImageDecoder;
use math::utils::clamp;

/// The permutation of dct coefficients.
pub static UNZIGZAG: [u8; 64] = [
    0,  1,  8, 16,  9,  2,  3, 10,
    17, 24, 32, 25, 18, 11,  4,  5,
    12, 19, 26, 33, 40, 48, 41, 34,
    27, 20, 13,  6,  7, 14, 21, 28,
    35, 42, 49, 56, 57, 50, 43, 36,
    29, 22, 15, 23, 30, 37, 44, 51,
    58, 59, 52, 45, 38, 31, 39, 46,
    53, 60, 61, 54, 47, 55, 62, 63,
];

/// A representation of a JPEG component
#[derive(Copy, Clone)]
pub struct Component {
    /// The Component's identifier
    pub id: u8,

    /// Horizontal sampling factor
    pub h: u8,

    /// Vertical sampling factor
    pub v: u8,

    /// The quantization table selector
    pub tq: u8,

    /// Index to the Huffman DC Table
    pub dc_table: u8,

    /// Index to the AC Huffman Table
    pub ac_table: u8,

    /// The dc prediction of the component
    pub dc_pred: i32
}

// Markers
// Baseline DCT
const SOF0: u8 = 0xC0;
// Progressive DCT
const SOF2: u8 = 0xC2;
// Huffman Tables
const DHT: u8 = 0xC4;
// Restart Interval start and End (standalone)
const RST0: u8 = 0xD0;
const RST7: u8 = 0xD7;
// Start of Image (standalone)
const SOI: u8 = 0xD8;
// End of image (standalone)
const EOI: u8 = 0xD9;
// Start of Scan
const SOS: u8 = 0xDA;
// Quantization Tables
const DQT: u8 = 0xDB;
// Number of lines
const DNL: u8 = 0xDC;
// Restart Interval
const DRI: u8 = 0xDD;
// Application segments start and end
const APP0: u8 = 0xE0;
const APPF: u8 = 0xEF;
// Comment
const COM: u8 = 0xFE;
// Reserved
const TEM: u8 = 0x01;

#[derive(PartialEq)]
enum JPEGState {
    Start,
    HaveSOI,
    HaveFirstFrame,
    HaveFirstScan,
    #[allow(dead_code)]
    End
}

/// The representation of a JPEG decoder
///
/// Does not support decoding progressive JPEG images
pub struct JPEGDecoder<R> {
    r: R,

    qtables: [u8; 64 * 4],
    dctables: [HuffTable; 2],
    actables: [HuffTable; 2],

    h: HuffDecoder,

    height: u16,
    width: u16,

    num_components: u8,
    scan_components: Vec<u8>,
    components: HashMap<usize, Component>,     // TODO: replace by `VecMap`

    mcu_row: Vec<u8>,
    mcu: Vec<u8>,
    hmax: u8,
    vmax: u8,

    interval: u16,
    mcucount: u32,
    expected_rst: u8,

    row_count: u8,
    decoded_rows: u32,
    padded_width: usize,
    state: JPEGState,
}

impl<R: Read>JPEGDecoder<R> {
    /// Create a new decoder that decodes from the stream ```r```
    pub fn new(r: R) -> JPEGDecoder<R> {
        let h: HuffTable  = Default::default();

        JPEGDecoder {
            r: r,

            qtables: [0u8; 64 * 4],
            dctables: [h.clone(), h.clone()],
            actables: [h.clone(), h.clone()],

            h: HuffDecoder::new(),

            height: 0,
            width: 0,

            num_components: 0,
            scan_components: Vec::new(),
            components: HashMap::new(),

            mcu_row: Vec::new(),
            mcu: Vec::new(),
            hmax: 0,
            vmax: 0,

            interval: 0,
            mcucount: 0,
            expected_rst: RST0,

            row_count: 0,
            decoded_rows: 0,
            state: JPEGState::Start,
            padded_width: 0
        }
    }

    fn decode_mcu_row(&mut self) -> ImageResult<()> {
        let bytesperpixel = self.num_components as usize;

        for x0 in range_step(0, self.padded_width * bytesperpixel, bytesperpixel * 8 * self.hmax as usize) {

            let _ = try!(self.decode_mcu());

            upsample_mcu (
                &mut self.mcu_row,
                x0,
                self.padded_width,
                bytesperpixel,
                &self.mcu,
                self.hmax,
                self.vmax
            );
        }

        Ok(())
    }

    fn decode_mcu(&mut self) -> ImageResult<()> {
        let mut i = 0;
        let tmp = self.scan_components.clone();

        for id in tmp.iter() {
            let mut c = self.components.get(&(*id as usize)).unwrap().clone();

            for _ in (0..c.h * c.v) {
                let pred  = try!(self.decode_block(i, c.dc_table, c.dc_pred, c.ac_table, c.tq));
                c.dc_pred = pred;
                i += 1;
            }

            self.components.insert(*id as usize, c);
        }

        self.mcucount += 1;
        self.read_restart()
    }

    fn decode_block(&mut self, i: usize, dc: u8, pred: i32, ac: u8, q: u8) -> ImageResult<i32> {
        let zz   = &mut self.mcu[i * 64..i * 64 + 64];
        let mut tmp = [i16x8::splat(0); 64 / 8];

        let dc_;
        {
            let tmp = unsafe {&mut *(tmp.as_mut_ptr() as *mut [i16; 64])};

            let dctable = &self.dctables[dc as usize];
            let actable = &self.actables[ac as usize];
            let qtable  = &self.qtables[64 * q as usize..64 * q as usize + 64];

            let t     = try!(self.h.decode_symbol(&mut self.r, dctable));

            let diff  = if t > 0 {
                try!(self.h.receive(&mut self.r, t))
            } else {
                0
            };

            // Section F.2.1.3.1
            let diff = extend(diff, t);
            let dc = diff + pred;
            tmp[0] = dc as i16;
            dc_ = dc;

            let mut k = 1usize;
            while k < 64 {
                let rs = try!(self.h.decode_symbol(&mut self.r, actable));

                let ssss = rs & 0x0F;
                let rrrr = rs >> 4;

                if ssss == 0 {
                    if rrrr != 15 {
                        break
                    }

                    k += 16;
                } else {
                    k += rrrr as usize;

                    // Figure F.14
                    let t = try!(self.h.receive(&mut self.r, ssss));

                    tmp[UNZIGZAG[k] as usize] = extend(t, ssss) as i16 * qtable[k] as i16;
                    k += 1;
                }
            }
        }

        transform::idct(&tmp, zz);

        Ok(dc_)
    }

    fn read_metadata(&mut self) -> ImageResult<()> {
        while self.state != JPEGState::HaveFirstScan {
            let byte = try!(self.r.read_u8());

            if byte != 0xFF {
                continue;
            }

            let marker = try!(self.r.read_u8());

            match marker {
                SOI => self.state = JPEGState::HaveSOI,
                DHT => try!(self.read_huffman_tables()),
                DQT => try!(self.read_quantization_tables()),
                SOF0 => {
                    let _ = try!(self.read_frame_header());
                    self.state = JPEGState::HaveFirstFrame;
                }
                SOS => {
                    let _ = try!(self.read_scan_header());
                    self.state = JPEGState::HaveFirstScan;
                }
                DRI => try!(self.read_restart_interval()),
                APP0 ... APPF | COM => {
                    let length = try!(self.r.read_u16::<BigEndian>());
                    let mut buf = Vec::with_capacity((length - 2) as usize);
                    try!(self.r.by_ref().take((length - 2) as u64).read_to_end(&mut buf));
                }
                TEM  => continue,
                SOF2 => return Err(image::ImageError::UnsupportedError("Marker SOF2 ist not supported.".to_string())),
                DNL  => return Err(image::ImageError::UnsupportedError("Marker DNL ist not supported.".to_string())),
                marker => return Err(image::ImageError::FormatError(format!("Unkown marker {} encountered.", marker))),
            }
        }

        Ok(())
    }

    fn read_frame_header(&mut self) -> ImageResult<()> {
        let _frame_length = try!(self.r.read_u16::<BigEndian>());
        let sample_precision = try!(self.r.read_u8());

        if sample_precision != 8 {
            return Err(image::ImageError::UnsupportedError(format!(
                "A sample precision of {} is not supported",
                sample_precision
            )))
        }

        self.height         = try!(self.r.read_u16::<BigEndian>());
        self.width          = try!(self.r.read_u16::<BigEndian>());
        self.num_components = try!(self.r.read_u8());

        if self.height == 0 || self.width == 0 {
            return Err(image::ImageError::DimensionError)
        }

        if self.num_components != 1 && self.num_components != 3 {
            return Err(image::ImageError::UnsupportedError(format!(
                "Frames with {} components are not supported",
                self.num_components
            )))
        }

        self.padded_width = 8 * ((self.width as usize + 7) / 8);

        let num_components = self.num_components;
        self.read_frame_components(num_components)
    }

    fn read_frame_components(&mut self, n: u8) -> ImageResult<()> {
        let mut blocks_per_mcu = 0;

        for _ in (0..n) {
            let id = try!(self.r.read_u8());
            let hv = try!(self.r.read_u8());
            let tq = try!(self.r.read_u8());

            let c = Component {
                id: id,
                h: hv >> 4,
                v: hv & 0x0F,
                tq: tq,
                dc_table: 0,
                ac_table: 0,
                dc_pred: 0
            };

            blocks_per_mcu += (hv >> 4) * (hv & 0x0F);
            self.components.insert(id as usize, c);
        }

        let (hmax, vmax) = self.components.iter().fold((0, 0), | (h, v), (_, c) | {
            (cmp::max(h, c.h), cmp::max(v, c.v))
        });

        self.hmax = hmax;
        self.vmax = vmax;

        // only 1 component no interleaving
        if n == 1 {
        for (_, c) in self.components.iter_mut() {
                c.h = 1;
                c.v = 1;
            }

            blocks_per_mcu = 1;
            self.hmax = 1;
            self.vmax = 1;
        }

        self.mcu = vec![0; blocks_per_mcu as usize * 64];

        let mcus_per_row = (self.width as f32 / (8 * hmax) as f32).ceil() as usize;
        let mcu_row_len = (hmax as usize * vmax as usize) * self.mcu.len() * mcus_per_row;

        self.mcu_row = vec![0; mcu_row_len];

        Ok(())
    }

    fn read_scan_header(&mut self) -> ImageResult<()> {
        let _scan_length = try!(self.r.read_u16::<BigEndian>());

        let num_scan_components = try!(self.r.read_u8());

        self.scan_components = Vec::new();

        for _ in (0..num_scan_components as usize) {
            let id = try!(self.r.read_u8());
            let tables = try!(self.r.read_u8());

            let c = self.components.get_mut(&(id as usize)).unwrap();

            c.dc_table = tables >> 4;
            c.ac_table = tables & 0x0F;

            self.scan_components.push(id);
        }

        let _spectral_end   = try!(self.r.read_u8());
        let _spectral_start = try!(self.r.read_u8());

        let approx = try!(self.r.read_u8());

        let _approx_high = approx >> 4;
        let _approx_low  = approx & 0x0F;

        Ok(())
    }

    fn read_quantization_tables(&mut self) -> ImageResult<()> {
        let mut table_length = try!(self.r.read_u16::<BigEndian>()) as i32;
        table_length -= 2;

        while table_length > 0 {
            let pqtq = try!(self.r.read_u8());
            let pq = pqtq >> 4;
            let tq = pqtq & 0x0F;

            if pq != 0 || tq > 3 {
                return Err(image::ImageError::FormatError("Quantization table malformed.".to_string()))
            }

            let slice = &mut self.qtables[64 * tq as usize..64 * tq as usize + 64];

            for i in (0usize..64) {
                slice[i] = try!(self.r.read_u8());
            }

            table_length -= 1 + 64;
        }

        Ok(())
    }

    fn read_huffman_tables(&mut self) -> ImageResult<()> {
        let mut table_length = try!(self.r.read_u16::<BigEndian>());
        table_length -= 2;

        while table_length > 0 {
            let tcth = try!(self.r.read_u8());
            let tc = tcth >> 4;
            let th = tcth & 0x0F;

            if tc != 0 && tc != 1 {
                return Err(image::ImageError::UnsupportedError(format!(
                    "Huffman table class {} is not supported", tc
                )))
            }

            let mut bits = Vec::with_capacity(16);
            try!(self.r.by_ref().take(16).read_to_end(&mut bits));
            let len = bits.len();

            let mt = bits.iter().fold(0, | a, b | a + *b);
            let mut huffval = Vec::with_capacity(mt as usize);
            try!(self.r.by_ref().take(mt as u64).read_to_end(&mut huffval));

            if tc == 0 {
                self.dctables[th as usize] = derive_tables(bits, huffval);
            } else {
                self.actables[th as usize] = derive_tables(bits, huffval);
            }

            table_length -= 1 + len as u16 + mt as u16;
        }

        Ok(())
    }


    fn read_restart_interval(&mut self) -> ImageResult<()> {
        let _length = try!(self.r.read_u16::<BigEndian>());
        self.interval = try!(self.r.read_u16::<BigEndian>());

        Ok(())
    }

    fn read_restart(&mut self) -> ImageResult<()> {
        let w = ((self.width + 7) / (self.hmax * 8) as u16) as u32;
        let h = ((self.height + 7) / (self.vmax * 8) as u16) as u32;

        if self.interval != 0  &&
           self.mcucount % (self.interval as u32) == 0 &&
           self.mcucount < w * h {

            let rst = try!(self.find_restart_marker());

            if rst == self.expected_rst {
                self.reset();
                self.expected_rst += 1;

                if self.expected_rst > RST7 {
                    self.expected_rst = RST0;
                }
            } else {
                return Err(image::ImageError::FormatError(format!(
                    "Unexpected restart maker {} found", rst
                )))
            }
        }

        Ok(())
    }

    fn find_restart_marker(&mut self) -> ImageResult<u8> {
        if self.h.marker != 0 {
            let m = self.h.marker;
            self.h.marker = 0;

            return Ok(m);
        }

        let mut b;
        loop {
            b = try!(self.r.read_u8());

            if b == 0xFF {
            b = try!(self.r.read_u8());
                match b {
                    RST0 ... RST7 => break,
                    EOI => return Err(image::ImageError::FormatError("Restart marker not found.".to_string())),
                    _   => continue
                }
            }
        }

        Ok(b)
    }

    fn reset(&mut self) {
        self.h.bits = 0;
        self.h.num_bits = 0;
        self.h.end = false;
        self.h.marker = 0;

        for (_, c) in self.components.iter_mut() {
            c.dc_pred = 0;
        }
    }
}

impl<R: Read> ImageDecoder for JPEGDecoder<R> {
    fn dimensions(&mut self) -> ImageResult<(u32, u32)> {
        if self.state == JPEGState::Start {
            let _ = try!(self.read_metadata());
        }

        Ok((self.width as u32, self.height as u32))
    }

    fn colortype(&mut self) -> ImageResult<color::ColorType> {
        if self.state == JPEGState::Start {
            let _ = try!(self.read_metadata());
        }

        let ctype = if self.num_components == 1 {
            color::ColorType::Gray(8)
        } else {
            color::ColorType::RGB(8)
        };

        Ok(ctype)
    }

    fn row_len(&mut self) -> ImageResult<usize> {
        if self.state == JPEGState::Start {
            let _ = try!(self.read_metadata());
        }

        let len = self.width as usize * self.num_components as usize;

        Ok(len)
    }

    fn read_scanline(&mut self, buf: &mut [u8]) -> ImageResult<u32> {
        if self.state == JPEGState::Start {
            let _ = try!(self.read_metadata());
        }

        if self.row_count == 0 {
            let _ = try!(self.decode_mcu_row());
        }

        let len   = self.padded_width * self.num_components as usize;
        let slice = &self.mcu_row[self.row_count as usize * len..
        self.row_count as usize * len + buf.len()];

        ::copy_memory(slice, buf);

        self.row_count = (self.row_count + 1) % (self.vmax * 8);
        self.decoded_rows += 1;

        Ok(self.decoded_rows)
    }

    fn read_image(&mut self) -> ImageResult<image::DecodingResult> {
        if self.state == JPEGState::Start {
            let _ = try!(self.read_metadata());
        }

        let row = try!(self.row_len());
        let mut buf = vec![0u8; row * self.height as usize];

        for chunk in buf.chunks_mut(row) {
            let _len = try!(self.read_scanline(chunk));
        }

        Ok(image::DecodingResult::U8(buf))
    }
}
#[allow(dead_code)]
fn upsample_mcu_scalar(out: &mut [u8], xoffset: usize, width: usize, bpp: usize, mcu: &[u8], h: u8, v: u8) {
    if mcu.len() == 64 {
        for y in (0usize..8) {
            for x in (0usize..8) {
                out[xoffset + x + (y * width)] = mcu[x + y * 8]
            }
        }
    } else {
        let y_blocks = h * v;

        let y_blocks = &mcu[..y_blocks as usize * 64];
        let cb = &mcu[y_blocks.len()..y_blocks.len() + 64];
        let cr = &mcu[y_blocks.len() + cb.len()..];

        let mut k = 0;

        for by in (0..v as usize) {
            let y0 = by * 8;

            for bx in (0..h as usize) {
                let x0 = xoffset + bx * 8 * bpp;

                for y in (0usize..8) {
                    for x in (0usize..8) {
                        let (a, b, c) = (y_blocks[k * 64 + x + y * 8], cb[x + y * 8], cr[x + y * 8]);
                        let (r, g, b) = ycbcr_to_rgb_scalar(a , b , c );

                        let offset = (y0 + y) * (width * bpp) + x0 + x * bpp;
                        out[offset + 0] = r;
                        out[offset + 1] = g;
                        out[offset + 2] = b;
                    }
                }

                k += 1;
            }
        }
    }
}

#[allow(dead_code)]
fn ycbcr_to_rgb_scalar(y: u8, cb: u8, cr: u8) -> (u8, u8, u8) {
    let y = y as f32;
    let cr = cr as f32;
    let cb = cb as f32;

    let r1 = y + 1.402f32 * (cr - 128f32) ;
    let g1 = y - 0.34414f32 * (cb - 128f32) - 0.71414f32 * (cr - 128f32);
    let b1 = y + 1.772f32 * (cb - 128f32);

    let r = clamp(r1 as i32, 0, 255) as u8;
    let g = clamp(g1 as i32, 0, 255) as u8;
    let b = clamp(b1 as i32, 0, 255) as u8;

    (r, g, b)
}

fn upsample_mcu(out: &mut [u8], xoffset: usize, width: usize, bpp: usize, mcu: &[u8], h: u8, v: u8) {
    if mcu.len() == 64 {
        for y in (0usize..8) {
            for x in (0usize..8) {
                out[xoffset + x + (y * width)] = mcu[x + y * 8]
            }
        }
    } else {
        let y_blocks = h * v;

        let y_blocks = &mcu[..y_blocks as usize * 64];
        let cb = &mcu[y_blocks.len()..y_blocks.len() + 64];
        let cr = &mcu[y_blocks.len() + cb.len()..];

        let mut k = 0;

        for by in (0..v as usize) {
            let y0 = by * 8;

            for bx in (0..h as usize) {
                let x0 = xoffset + bx * 8 * bpp;

                let mut y = 0;
                while y < 8 {
                    let a = u8x16::load(y_blocks, k * 64 + y * 8);
                    let b = u8x16::load(cb, y * 8);
                    let c = u8x16::load(cr, y * 8);
                    let (r, g, b) = ycbcr_to_rgb(a , b , c );

                    let mut rr = [0; 16];
                    r.store(&mut rr, 0);
                    let mut gg = [0; 16];
                    g.store(&mut gg, 0);
                    let mut bb = [0; 16];
                    b.store(&mut bb, 0);
                    let mut offset = (y0 + y) * (width * bpp) + x0;
                    for i in 0..8 {
                        let (r1, g1, b1) = (rr[i], gg[i], bb[i]);
                        let (r2, g2, b2) = (rr[i + 8], gg[i + 8], bb[i + 8]);
                        out[offset + 0] = r1;
                        out[offset + 1] = g1;
                        out[offset + 2] = b1;
                        out[offset + width * bpp + 0] = r2;
                        out[offset + width * bpp + 1] = g2;
                        out[offset + width * bpp + 2] = b2;

                        offset += bpp;
                    }

                    y += 2;
                }

                k += 1;
            }
        }
    }
}

fn ycbcr_to_rgb(y: u8x16, cb: u8x16, cr: u8x16) -> (u8x16, u8x16, u8x16) {
    use simd::{f32x4, i32x4};
    use simd::x86::sse2::*;

    fn to_f32(a: u8x16) -> (f32x4, f32x4, f32x4, f32x4) {
        let a1 = u8x16::new(a.extract(0), 0, 0, 0,
                            a.extract(1), 0, 0, 0,
                            a.extract(2), 0, 0, 0,
                            a.extract(3), 0, 0, 0);
        let a2 = u8x16::new(a.extract(4 + 0), 0, 0, 0,
                            a.extract(4 + 1), 0, 0, 0,
                            a.extract(4 + 2), 0, 0, 0,
                            a.extract(4 + 3), 0, 0, 0);
        let a3 = u8x16::new(a.extract(8 + 0), 0, 0, 0,
                            a.extract(8 + 1), 0, 0, 0,
                            a.extract(8 + 2), 0, 0, 0,
                            a.extract(8 + 3), 0, 0, 0);
        let a4 = u8x16::new(a.extract(12 + 0), 0, 0, 0,
                            a.extract(12 + 1), 0, 0, 0,
                            a.extract(12 + 2), 0, 0, 0,
                            a.extract(12 + 3), 0, 0, 0);

        fn c(a: u8x16) -> f32x4 {
            unsafe {::std::mem::transmute::<_, i32x4>(a)}.to_f32()
        }

        (c(a1), c(a2), c(a3), c(a4))
    }
    let (y1, y2, y3, y4) = to_f32(y);
    let (cr1, cr2, cr3, cr4) = to_f32(cr);
    let (cb1, cb2, cb3, cb4) = to_f32(cb);

    let _128 = f32x4::splat(128.0);
    let r1 = y1 + f32x4::splat(1.402f32) * (cr1 - _128) ;
    let r2 = y2 + f32x4::splat(1.402f32) * (cr2 - _128) ;
    let g1 = y1 - f32x4::splat(0.34414f32) * (cb1 - _128) - f32x4::splat(0.71414f32) * (cr1 - _128);
    let g2 = y2 - f32x4::splat(0.34414f32) * (cb2 - _128) - f32x4::splat(0.71414f32) * (cr2 - _128);
    let b1 = y1 + f32x4::splat(1.772f32) * (cb1 - _128);
    let b2 = y2 + f32x4::splat(1.772f32) * (cb2 - _128);

    let r3 = y3 + f32x4::splat(1.402f32) * (cr3 - _128) ;
    let r4 = y4 + f32x4::splat(1.402f32) * (cr4 - _128) ;
    let g3 = y3 - f32x4::splat(0.34414f32) * (cb3 - _128) - f32x4::splat(0.71414f32) * (cr3 - _128);
    let g4 = y4 - f32x4::splat(0.34414f32) * (cb4 - _128) - f32x4::splat(0.71414f32) * (cr4 - _128);
    let b3 = y3 + f32x4::splat(1.772f32) * (cb3 - _128);
    let b4 = y4 + f32x4::splat(1.772f32) * (cb4 - _128);

    (r1.to_i32().packs(r2.to_i32()).packus(r3.to_i32().packs(r4.to_i32())),
     g1.to_i32().packs(g2.to_i32()).packus(g3.to_i32().packs(g4.to_i32())),
     b1.to_i32().packs(b2.to_i32()).packus(b3.to_i32().packs(b4.to_i32())))

}

// Section F.2.2.1
// Figure F.12
fn extend(v: i32, t: u8) -> i32 {
    // FIXME check if wrapping sub is what we want
    let mut vt: i32 = 0;
    let shift = (t as usize).wrapping_sub(1);
    if shift < 31 {
        vt = 1 << shift;
    }

    if v < vt {
    v + ((-1) << t as usize) + 1
    }
    else {
        v
    }
}
