#![feature(duration_span)]

extern crate image;
use std::time::Duration;
// use std::io::prelude::*; use std::fs::File;
const IMAGE: &'static [u8] = include_bytes!("aurora.jpg");
const REFERENCE: &'static [u8] = include_bytes!("bytes.bin");

fn main() {
    let mut img = None;
    println!("{:?}", Duration::span(|| {
        for _ in 0..10 {
            img = Some(image::load_from_memory_with_format(IMAGE, image::JPEG).unwrap());
        }
    }));

    assert!(img.unwrap().raw_pixels() == REFERENCE);
    //File::create("examples/bytes.bin").unwrap().write(&image::load_from_memory_with_format(IMAGE, image::JPEG).unwrap().raw_pixels()).unwrap();
}
