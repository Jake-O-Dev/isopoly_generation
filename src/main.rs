use std::{time::Instant, fs};

use crate::{
  algebraic_types::{generate_iso_polynomials, Matrix},
  polynomials::{Polynomial, generate_transform_lut}
};

#[allow(non_snake_case)]
#[allow(non_camel_case_types)]
mod algebraic_types;
#[allow(non_snake_case)]
#[allow(non_camel_case_types)]
mod polynomials;

const DEGREE: usize = 5;
const FIELD_SIZE: usize = 3;

const COEFF_BIT_SIZES: [usize; 5] = [1,1,1,2,2];
const COEFF_BIT_SIZE: usize = COEFF_BIT_SIZES[FIELD_SIZE];

// (q^(d+2 choose 2) - 1) / (q - 1)
const DPLUS2_CHOOSE_2: usize = ((DEGREE+2) * (DEGREE+1)) / 2;
const POLYNOMIALS: usize = (FIELD_SIZE.pow(DPLUS2_CHOOSE_2 as u32) - 1) / (FIELD_SIZE - 1);

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
  println!("Generate matrices");
  let pgl3 = Matrix::generate_pgl3();
  println!("Number of matrices: {}", pgl3.len());
  println!();

  // Lookup Tables
  println!("Generate lookup stuff");
  let normal = Polynomial::generate_default_lut();
 
  /* <> test transform lut <>
  // 
  // let i = 2;
  // let j = 6;
  // let matrix1 = pgl3[i];
  // let term1 = normal[j].create_similar(1);
  // //let transformed = transform_lut[i][j];
  // matrix1.print();
  // println!("{}", term1.str());
  // println!();
  // let transformed_bits = term1.transform_by_matrix(&matrix1, &normal);
  // let trans_poly = Polynomial::new(transformed_bits);
  // println!();
  // println!("{transformed_bits:b}");
  // trans_poly.print(&normal);
  // //println!("{:b} {:b}", transformed, term1.transform_by_matrix(&matrix1, &normal));

  // //Polynomial::new(transformed).print(&normal);
  // return;
  */

  let transform_lut = generate_transform_lut(&pgl3, &normal); // REMOVE THIS LATER



  let lookup_time = Instant::now();
  println!("Generating took: {:?}", (lookup_time-start_time));
  println!();


  //
  // Polynomials
  //
  println!("Generate isomorphic polynomials");
  

  let iso_polys = generate_iso_polynomials(&transform_lut);
  
  println!("Generated {} isomorphic polynomials", iso_polys.len());
  
  // Counting the polys for verification
  let mut sum: u32 = 0;
  for isopoly in &iso_polys {
    let (_, size) = isopoly.deconstruct();
    sum += size;
  }
  let poly_time = Instant::now();
  println!("Total polynomials: {}", sum);
  println!("Generating took: {:?}", (poly_time-lookup_time));
  println!();
  
  
  let a: Vec<String> = iso_polys.iter().map(|t| t.to_string(&normal)).collect();
  let b = a.join("\n");
  fs::write(FILE_NAME, b).expect("Unable to write file");

  
}