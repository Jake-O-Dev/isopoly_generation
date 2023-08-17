use std::{time::Instant, fs};

use crate::{
  algebraic_types::{Matrix, generate_iso_polynomials},
  polynomials::{Polynomial, generate_transform_lut}
};

#[allow(non_snake_case)]
#[allow(non_camel_case_types)]
mod algebraic_types;
#[allow(non_snake_case)]
#[allow(non_camel_case_types)]
mod polynomials;

const DEGREE: usize = 5;
const FIELD_ORDER: usize = 3;

const COEFF_BIT_SIZES: [usize; 5] = [1,1,1,2,2];
const COEFF_BIT_SIZE: usize = COEFF_BIT_SIZES[FIELD_ORDER];

// (q^(d+2 choose 2) - 1) / (q - 1)
const DPLUS2_CHOOSE_2: usize = ((DEGREE+2) * (DEGREE+1)) / 2;
const POLYNOMIALS: usize = (FIELD_ORDER.pow(DPLUS2_CHOOSE_2 as u32) - 1) / (FIELD_ORDER - 1);


const CHUNK_SIZE: usize = 1024*64;

const FILE_NAME: &str = "./output.txt";

#[derive(Debug,Clone,Copy,PartialEq)]
struct CustomChunk {
  pub start: usize,
  pub end: usize,
}

#[allow(unreachable_code, unused_variables)]
fn main() {
  let start_time = Instant::now();

  // Matrices
  println!("Generating matrices...");
  let pgl3 = Matrix::generate_pgl3();
  println!("Number of matrices: {}", pgl3.len());
  println!();

  // Lookup Tables
  println!("Generating lookup tables...");
  let normal = Polynomial::generate_default_lut();
  let transform_lut = generate_transform_lut(&pgl3, &normal);
  let lookup_time = Instant::now();
  println!("Generating took: {:?}", (lookup_time-start_time));
  println!();
  
  let results = generate_iso_polynomials(&transform_lut);


  // Counting the polys for verification
  let mut sum: u32 = 0;
  for isopoly in &results {
    let (_, size) = isopoly.deconstruct();
    sum += size;
  }
  let poly_time = Instant::now();
  println!("Total polynomials: {}", sum);
  println!("Generating took: {:?}", (poly_time-lookup_time));
  println!();
  
  
  let a: Vec<String> = results.iter().map(|t| t.to_string(&normal)).collect();
  let begin = "# This file give you represants of all isomorphism classes of homogeneous projective polynomials over P_2 by applying PGL_3 and their class size.\n";
  let begin_2 = "# Homogeneous Degree | Field Order\n";
  let begin_3 = format!("{} | {}\n", DEGREE, FIELD_ORDER);
  let begin_4 = "# Constant_(xpower)(ypower)(zpower) .... | Isomorphism Class size\n";
  let b = a.join("\n");

  let c = begin.to_owned() + &begin_2.to_owned() + &begin_3 + &begin_4.to_owned() + &b;
  fs::write(FILE_NAME, c).expect("Unable to write file");

  println!("Finished printing.");
  println!("Total time: {:?}", (Instant::now()-start_time));
  
}

