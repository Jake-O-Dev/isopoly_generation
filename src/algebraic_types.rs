use crate::{POLYNOMIALS, COEFF_BIT_SIZE, FIELD_SIZE};
use crate::polynomials::{Term, Polynomial};


#[derive(Debug, Copy, Clone, PartialEq, Hash, Eq)]
pub struct Matrix {
  pub data: [[u8;3];3]
}

impl Matrix {
  pub fn determinant(&self) -> i32 {
    (self.data[0][0] as i32 * ((self.data[1][1] * self.data[2][2]) as i32 - (self.data[1][2] * self.data[2][1]) as i32) ) as i32 -
    (self.data[0][1] as i32 * ((self.data[1][0] * self.data[2][2]) as i32 - (self.data[1][2] * self.data[2][0]) as i32)) as i32 +
    (self.data[0][2] as i32 * ((self.data[1][0] * self.data[2][1]) as i32 -( self.data[1][1] * self.data[2][0]) as i32)) as i32
  }

  pub fn new(data: [[u8;3];3]) -> Matrix 
  {
    Matrix { data: data }
  }

  #[allow(dead_code)]
  pub fn print(&self) {
    println!("Matrix, det: {}", self.determinant());
    for i in 0..3 {
      for j in 0..3 {
        print!("{} ", self.data[i][j]);
      }
      println!("");
    }
  }


  pub fn generate_pgl3() -> Vec<Matrix> {
    match FIELD_SIZE {
      2 => Matrix::generate_pgl3_f2(),
      3 => Matrix::generate_pgl3_f3(),
      _ => panic!("Field size not supported"),
    }
  }

  pub fn generate_pgl3_f2() -> Vec<Matrix> {
    let mut pgl3_f2: Vec<Matrix> = Vec::new();
    for i in 0..(1<<9) {
      let mut data: [[u8;3];3] = [[0;3];3];
      for j in 0..9 {
        data[j/3][j%3] = ((i >> j) & 1) as u8;
      }
      let matrix = Matrix::new(data);
      if matrix.determinant() % 2 != 0 {
        pgl3_f2.push(matrix);
      }
    }
    pgl3_f2
  }

  pub fn generate_pgl3_f3() -> Vec<Matrix> {
    let mut pgl3_f3: Vec<Matrix> = Vec::new();
    for i in 0..19683 {
      let mut data: [[u8;3];3] = [[0;3];3];
      for j in 0..9 {
        data[j/3][j%3] = Matrix::get_ternary_digit(i, j) as u8;
      }
      let matrix = Matrix::new(data);
      if i == 728 {
        println!("");
      }

      if matrix.determinant() % 3 == 1 || matrix.determinant() % 3 == -2 {
        pgl3_f3.push(matrix);
      }
    }
    pgl3_f3
  }

  fn get_ternary_digit(digit: u64, bit_index: usize) -> u64 {
    const DIGIT_LOOKUP: [u64;9]= [1, 3, 9, 27, 81, 243, 729, 2187, 6561];
    let div = digit / DIGIT_LOOKUP[bit_index];
    div % 3 
  }
}



#[derive(Debug, Copy, Clone, PartialEq)]
pub struct IsoPolynomial {
  pub representative: Polynomial,
  pub size: u32,
}

impl IsoPolynomial {
  pub fn deconstruct(self) -> (Polynomial, u32) {
    (self.representative, self.size)
  }

  pub fn to_string(&self, normal: &Vec<Term>) -> String {
    format!("{} | {}", self.representative.str(normal), self.size)
  }
}


#[derive(Debug,Clone)]
pub struct PackedBool {
  pub data: Vec<u8>,
}

impl PackedBool {
  pub fn new(size: usize) -> PackedBool {
    let mut pack = PackedBool {
      data: Vec::with_capacity(size/8)
    };
    for _ in 0..(size/8) {
        pack.data.push(0);
    }
    pack
  }

  pub fn get(&self, index: usize) -> bool {
    let partial = self.data[index / 8];
    partial >> (index % 8) & 1 == 1
  }
  
  pub fn set(&mut self, index: usize, value: bool) {
    let partial = self.data[index / 8];
    if value {
      self.data[index / 8] = partial | (1 << index % 8)
    } else {      
      self.data[index / 8] = partial & (0xFF ^ (1 << index % 8))
    }
  }
}

fn f3_bijection_inverse(index: u64) -> u64 {
  // const DIGIT_LOOKUP: [u64;16]= [1, 3, 9, 27, 81, 243, 729, 2187, 6561,19683,59049,177147,531441,1594323,4782969,14348907];
  let mut mult = 1;
  let mut res = 0;
  for i in 0..32 {
    res += mult * ((index >> 2*i) % 3); 
    mult *= 3;
  }
  res
}

fn f3_bijection(mut index: u64) -> u64 {
  let mut res = 0;
  for i in 0..32 {
    res += (index % 3) << COEFF_BIT_SIZE*i; 
    index /= 3;
  }
  res
}

pub fn generate_iso_polynomials(transform_lut: &Vec<Vec<u64>>) -> Vec<IsoPolynomial>{
  let mut things = PackedBool::new(POLYNOMIALS + 8);

  things.set(0, true);

  let mut iso_polys = Vec::new();

  for i in 1..POLYNOMIALS/1000 {
    if things.get(i) == false {
      things.set(i, true);
      let poly = Polynomial::new(f3_bijection(i as u64));
      let mut count = 1;
      let mut smallest_poly = poly;
      for i in 0..transform_lut.len() { // loop over matrices
        let perm_poly = poly.transform_by_matrix(&transform_lut[i]);
        if things.get(f3_bijection_inverse(perm_poly.bits) as usize) == false {
          count += 1;
          things.set(f3_bijection_inverse(perm_poly.bits) as usize, true);
          if perm_poly.bits.count_ones() < smallest_poly.bits.count_ones() {
            smallest_poly = perm_poly;
          }
        }
      }
      iso_polys.push(IsoPolynomial { representative: smallest_poly, size: count});
    }
  }
  iso_polys
}
