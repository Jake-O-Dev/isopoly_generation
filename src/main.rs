use std::{time::Instant, fs, sync::{mpsc, Arc, Mutex, RwLock}, thread};

use algebraic_types::IsoPolynomial;

use crate::{
  algebraic_types::{Matrix, PackedBool, find_isomorphisms, generate_iso_polynomials},
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


const MULTI_THREADING: bool = false;
const PRINTING: bool = true;
const NUM_THREADS: usize = 8;
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
  


  let results = if MULTI_THREADING {
    //
    // Chunks
    //
    println!("Generate chunks, start threads and count smooth polynomials!");
    
    let mut chunks = Vec::new();
    let mut start = 0;
    while start < POLYNOMIALS { // TODO: CHECK THIS RANGE
      chunks.push(CustomChunk {
          start,
          end: std::cmp::min(start + CHUNK_SIZE, POLYNOMIALS), // end is exclusive
      });
      start += CHUNK_SIZE;
    }
    let chunk_length = chunks.len();
    chunks.reverse();
  
    println!("Amount of polynomials to verify: {} | Amount of chunks: {} | Amount of threads: {}", POLYNOMIALS, chunks.len(), NUM_THREADS);
    println!();
                  
  
    //
    // Polynomials
    //
    println!("Allocating memory");
    let mut verified_polynomial = PackedBool::new(usize::pow(3, 21));
    verified_polynomial.set(0, true);
    println!("Allocated memory, do some time checking");
  
  
    // Thread arc stuff
    let (tx, rx) = mpsc::channel();
  
    let arc_transform_lut = Arc::new(transform_lut);
    let arc_verified = Arc::new(RwLock::new(verified_polynomial));
    let arc_chunks = Arc::new(Mutex::new(chunks));
    
  
    
    for _ in 0..NUM_THREADS {
        // Clone the sender to move into each thread
        let a_tx = tx.clone();
  
        // Clone the recomputed results to move into each thread locally
        let local_transform_lut = arc_transform_lut.clone();
        let local_verified = arc_verified.clone();
        let local_chunks = arc_chunks.clone();
  
        // Spawn a new thread
        thread::spawn(move || {
          // let lol = &local_chunks.lock().unwrap().pop();
          loop {
            let (start,end, index);
            {
              let chunk_vec = &mut local_chunks.lock().unwrap();
              let chunk = chunk_vec.pop();
              match chunk {
                Some(t) => { start = t.start; end = t.end; index = chunk_vec.len()}
                None => {return;}
              }
            }
            
            let result =  
            find_isomorphisms(start, end, &local_transform_lut, &local_verified);
            if PRINTING {
              println!("Result size: {} | Chunks left: {index} | Total Chunks: {chunk_length} | Estimated time: {:.2}", result.len(), index as f64 * (Instant::now() - lookup_time).as_secs_f64() / (chunk_length - index) as f64);
            }
            a_tx.send(result).unwrap();      
          }
        });
      }
  
    drop(tx);
  
    // let mut ismorphisms: [usize; MAX_FIELD_EXT] = [0; MAX_FIELD_EXT];
    let mut results = Vec::new();
    let mut checked_polynomials = 0;
    for mut result in rx {
      checked_polynomials += result.iter().fold(0, |acc, iso| acc + iso.size);
      println!("Checked polynomials: {}/{} | Percentage: {:.3}%", checked_polynomials, POLYNOMIALS,  checked_polynomials as f64 / POLYNOMIALS as f64 * 100.0 );
      results.append(&mut result);
    }
    results
  } else {
    generate_iso_polynomials(&transform_lut)
  };
  
  



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

